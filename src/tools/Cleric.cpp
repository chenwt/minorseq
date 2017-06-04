// Copyright (c) 2016-2017, Pacific Biosciences of California, Inc.
//
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted (subject to the limitations in the
// disclaimer below) provided that the following conditions are met:
//
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//
//  * Redistributions in binary form must reproduce the above
//    copyright notice, this list of conditions and the following
//    disclaimer in the documentation and/or other materials provided
//    with the distribution.
//
//  * Neither the name of Pacific Biosciences nor the names of its
//    contributors may be used to endorse or promote products derived
//    from this software without specific prior written permission.
//
// NO EXPRESS OR IMPLIED LICENSES TO ANY PARTY'S PATENT RIGHTS ARE
// GRANTED BY THIS LICENSE. THIS SOFTWARE IS PROVIDED BY PACIFIC
// BIOSCIENCES AND ITS CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
// WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
// OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL PACIFIC BIOSCIENCES OR ITS
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
// USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
// OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
// SUCH DAMAGE.

// Author: Armin Töpfer

// Inspired by work of David Seifert

#include <boost/algorithm/string.hpp>

#include <pbcopper/utility/FileUtils.h>

#include <pbbam/DataSet.h>
#include <pbbam/MD5.h>

#include <pacbio/align/SimdAlignment.h>
#include <pacbio/io/BamParser.h>

#include <pacbio/cleric/Cleric.h>

namespace PacBio {
namespace Cleric {
void Cleric::Align(const std::string& fromReference, const std::string& toReference,
                   std::string* fromReferenceAligned, std::string* toReferenceAligned)
{
    auto align = Align::PariwiseAlignmentFasta(fromReference, toReference);

    *fromReferenceAligned = align.Target;
    *toReferenceAligned = align.Query;
}

void Cleric::Convert(std::string outputFile)
{
    using BAM::Cigar;
    using BAM::CigarOperation;
    using BAM::CigarOperationType;

    // Get data source
    auto query = IO::BamQuery(alignmentPath_);
    std::unique_ptr<BAM::BamWriter> out;

    auto ProcessHeaderAndCreateBamWriter = [this, &outputFile, &out](const BAM::BamRecord& read) {
        BAM::BamHeader h = read.Header().DeepCopy();

        if (h.Sequences().empty())
            throw std::runtime_error("Could not find reference sequence name");

        const auto localFromReferenceName = h.Sequences().begin()->Name();
        if (localFromReferenceName != fromReferenceName_)
            throw std::runtime_error("Internal error. Reference name mismatches");

        const auto RemoveGaps = [](const std::string& input) {
            std::string seq = input;
            seq.erase(std::remove(seq.begin(), seq.end(), '-'), seq.end());
            return seq;
        };
        toReferenceGapless_ = RemoveGaps(toReferenceSequence_);
        fromReferenceGapless_ = RemoveGaps(fromReferenceSequence_);

        for (size_t pos = 0, i = 0; i < fromReferenceSequence_.size(); ++i) {
            if (fromReferenceSequence_.at(i) != '-') {
                if (refToSourcePos_.find(i) == refToSourcePos_.cend()) {
                    refToSourcePos_.emplace(pos, i);
                }
                ++pos;
            }
        }

        for (size_t pos = 0, i = 0; i < toReferenceSequence_.size(); ++i) {
            if (toReferenceSequence_.at(i) != '-') {
                if (sourceToRefPos_.find(i) == sourceToRefPos_.cend()) {
                    sourceToRefPos_.emplace(i, pos);
                }
                ++pos;
            }
        }

        h.ClearSequences();
        auto bamRefSequence =
            BAM::SequenceInfo(toReferenceName_, std::to_string(toReferenceGapless_.size()));
        bamRefSequence.Checksum(BAM::MD5Hash(toReferenceSequence_));
        h.AddSequence(bamRefSequence);

        const bool isXml = Utility::FileExtension(outputFile) == "xml";
        if (isXml) boost::replace_last(outputFile, ".consensusalignmentset.xml", ".bam");

        // Write Dataset
        using BAM::DataSet;
        const std::string metatype = "PacBio.AlignmentFile.AlignmentBamFile";
        DataSet clericSet(DataSet::TypeEnum::ALIGNMENT);
        BAM::ExternalResource resource(metatype, outputFile);

        BAM::FileIndex pbi("PacBio.Index.PacBioIndex", outputFile + ".pbi");
        resource.FileIndices().Add(pbi);

        clericSet.ExternalResources().Add(resource);
        clericSet.Name(clericSet.TimeStampedName());

        const auto outputPrefix = outputFile.substr(0, outputFile.size() - 4);
        std::ofstream clericDSout(outputPrefix + ".consensusalignmentset.xml");
        clericSet.SaveToStream(clericDSout);

        out.reset(new BAM::BamWriter(outputFile, h));
    };

    // Convert and write to BAM
    for (auto read : *query) {
        if (!out) ProcessHeaderAndCreateBamWriter(read);
        std::string sourceStr = fromReferenceSequence_;
        std::string destStr = toReferenceSequence_;

        // Expand RLE cigar to flat vector
        std::string expandedCigarOps;
        for (const auto& c : read.CigarData(false))
            for (size_t i = 0; i < c.Length(); ++i)
                expandedCigarOps += c.Char();
        expandedCigarOps += "YZ";

        CigarOperation oldCigarState;  // UNKNOW_OP
        CigarOperation newCigarState;  // UNKNOW_OP

        bool foundStart = false;
        int posInRead = 0;
        int posInCigar = 0;
        int posInSourceRef = refToSourcePos_.at(read.ReferenceStart());

        Cigar newCigarTuple;

        int newSamStart = 0;
        int posInDestRef = 0;

        while (posInCigar < static_cast<int>(expandedCigarOps.size())) {
            char op = expandedCigarOps.at(posInCigar);

            CigarOperation newState;  // UNKNOWN_OP

            bool isFirstCigarAfterEnd = false;
            bool isSecondCigarAfterEnd = false;

            switch (op) {
                case 'M':
                case '=':
                case 'X':
                    if (!foundStart) {
                        if (sourceStr.at(posInSourceRef) == '-') {
                            // Dest:   A---AAA
                            // Source: AAA-AAA
                            // Read:     A-AAA
                            //            ^

                            ++posInSourceRef;
                            continue;
                        }

                        // don't have a start POS yet
                        if (sourceToRefPos_.find(posInSourceRef) != sourceToRefPos_.cend()) {
                            newSamStart = sourceToRefPos_.at(posInSourceRef);
                            // Dest:   ---AAA
                            // Source: AAAAAA
                            // Read:      AAA
                            //            ^

                            newState = newMatch_;
                            posInDestRef = posInSourceRef;
                            foundStart = true;
                            ++posInDestRef;
                        } else {
                            // Dest:   ----AA
                            // Source: AAAAAA
                            // Read:      AAA
                            //            ^

                            // left Clip
                            newState = newSoft_;
                        }

                        ++posInCigar;
                        ++posInRead;
                        ++posInSourceRef;
                    } else {
                        if (sourceStr.at(posInSourceRef) == '-') {
                            if (destStr.at(posInDestRef) == '-') {
                                // Dest:   AAA-AAA
                                // Source: AAA-AAA
                                // Read:   AAA-AAA
                                //            ^

                                ++posInSourceRef;
                                ++posInDestRef;
                                continue;
                            } else {
                                // Dest:   AAAAAAA
                                // Source: AAA-AAA
                                // Read:   AAA-AAA
                                //            ^

                                // Deletion
                                newState = newDel_;

                                ++posInSourceRef;
                                ++posInDestRef;
                            }
                        } else {
                            if (destStr.at(posInDestRef) == '-') {
                                // Dest:   AAA-AAA
                                // Source: AAAAAAA
                                // Read:   AAAAAAA
                                //            ^

                                // Insertion
                                newState = newIns_;

                                ++posInSourceRef;
                                ++posInDestRef;
                                ++posInCigar;
                                ++posInRead;
                            } else {
                                // Dest:   AAAAAAA
                                // Source: AAAAAAA
                                // Read:   AAAAAAA
                                //            ^

                                newState = newMatch_;

                                ++posInSourceRef;
                                ++posInDestRef;
                                ++posInCigar;
                                ++posInRead;
                            }
                        }
                    }
                    break;
                case 'I':
                    if (!foundStart) {
                        if (sourceStr.at(posInSourceRef) == '-') {
                            // Dest:   A---AAA
                            // Source: AAA-AAA
                            // Read:     AAAAA
                            //            ^

                            ++posInSourceRef;
                            continue;
                        }

                        // Dest:   -- AAA
                        // Source: AA AAA
                        // Read:    AGAAA
                        //           ^

                        // left Clip
                        newState = newSoft_;

                        ++posInCigar;
                        ++posInRead;
                    } else {
                        if (sourceStr.at(posInSourceRef) == '-') {
                            if (destStr.at(posInDestRef) == '-') {
                                // Dest:   AAA -AAA
                                // Source: AAA -AAA
                                // Read:   AAAA AAA
                                //            ^

                                ++posInSourceRef;
                                ++posInDestRef;
                                continue;
                            } else {
                                // Dest:   AAA AAAA
                                // Source: AAA -AAA
                                // Read:   AAAA AAA
                                //            ^

                                newState = newMatch_;

                                ++posInSourceRef;
                                ++posInDestRef;
                                ++posInCigar;
                                ++posInRead;
                            }
                        } else {
                            // Dest:   AAA -AAA
                            // Source: AAA AAAA
                            // Read:   AAAA AAA
                            //            ^
                            // OR
                            // Dest:   AAA AAA
                            // Source: AAA AAA
                            // Read:   AAAAAAA
                            //            ^

                            // Insertion
                            newState = newIns_;

                            ++posInCigar;
                            ++posInRead;
                        }
                    }
                    break;
                case 'N':
                case 'D':
                    if (!foundStart) {
                        if (sourceStr.at(posInSourceRef) == '-') {
                            // Dest:   A---AAA
                            // Source: AAA-AAA
                            // Read:     A--AA
                            //            ^

                            ++posInSourceRef;
                            continue;
                        }

                        // Dest:   ---AAA
                        // Source: AAAAAA
                        // Read:    A-AAA
                        //           ^

                        ++posInCigar;
                        ++posInSourceRef;
                        continue;
                    } else {
                        // have start POS
                        if (sourceStr.at(posInSourceRef) == '-') {
                            if (destStr.at(posInDestRef) == '-') {
                                // Dest:   AAA-AAA
                                // Source: AAA-AAA
                                // Read:   AAA-AAA
                                //            ^

                                ++posInSourceRef;
                                ++posInDestRef;
                                continue;
                            } else {
                                // Dest:   AAAAAAA
                                // Source: AAA-AAA
                                // Read:   AAA--AA
                                //            ^

                                // Deletion
                                newState = newDel_;

                                ++posInSourceRef;
                                ++posInDestRef;
                            }
                        } else {
                            if (destStr.at(posInDestRef) == '-') {
                                // Dest:   AAA-AAA
                                // Source: AAAAAAA
                                // Read:   AAA-AAA
                                //            ^

                                // Padded deletion
                                ++posInSourceRef;
                                ++posInDestRef;
                                ++posInCigar;

                                newState = newPad_;
                            } else {
                                // Dest:   AAAAAAA
                                // Source: AAAAAAA
                                // Read:   AAA-AAA
                                //            ^

                                // Deletion
                                newState = newDel_;

                                ++posInSourceRef;
                                ++posInDestRef;
                                ++posInCigar;
                            }
                        }
                    }
                    break;
                case 'S':
                    newState = newSoft_;

                    ++posInCigar;
                    ++posInRead;
                    break;
                case 'H':
                    newState = CigarOperation(CigarOperationType::HARD_CLIP, 1);

                    ++posInCigar;
                    break;
                case 'P':
                    if (foundStart) {
                        // Dest:   ---AAA
                        // Source: AAAAAA
                        // Read:    A-AAA
                        //           ^

                        ++posInCigar;
                        ++posInSourceRef;
                        continue;

                    } else {
                        // have start POS
                        if (sourceStr.at(posInSourceRef) == '-') {
                            if (destStr.at(posInDestRef) == '-') {
                                // Dest:   AAA-AAA
                                // Source: AAA-AAA
                                // Read:   AAA-AAA
                                //            ^

                                // Padded deletion
                                ++posInCigar;

                                newState = newPad_;
                            } else {
                                // Dest:   AAAAAAA
                                // Source: AAA-AAA
                                // Read:   AAA--AA
                                //            ^

                                // Deletion
                                newState = newDel_;

                                ++posInCigar;
                                ++posInSourceRef;
                                ++posInDestRef;
                            }
                        } else {
                            // Dest:   AAA--AAA
                            // Source: AAAAAAAA
                            // Read:   AAA-AAAA
                            //            ^
                            // OR
                            // Dest:   AAA AAAA
                            // Source: AAA AAAA
                            // Read:   AAA-AAAA
                            //            ^

                            // Padded deletion
                            ++posInCigar;

                            newState = newPad_;
                        }
                    }
                    break;
                case 'Y':
                    ++posInCigar;
                    isFirstCigarAfterEnd = true;
                    break;
                case 'Z':
                    ++posInCigar;
                    isSecondCigarAfterEnd = true;
                    break;
                default:
                    throw std::runtime_error("UNKNOWN CIGAR");
            }

            // If we reached Z, we have processed the CIGAR and can push the
            // lastest cigar operation.
            if (isSecondCigarAfterEnd) newCigarTuple.push_back(oldCigarState);

            if (newState.Type() != newCigarState.Type()) {
                // I ...... Y (end)
                if (newState.Type() == CigarOperationType::UNKNOWN_OP && isFirstCigarAfterEnd &&
                    newCigarState.Type() == CigarOperationType::INSERTION) {
                    newCigarState.Type(CigarOperationType::SOFT_CLIP);
                }

                // have to rewrite CIGAR tuples if, a D and I operations are adjacent
                // D + I
                if (oldCigarState.Type() == CigarOperationType::DELETION &&
                    newCigarState.Type() == CigarOperationType::INSERTION) {
                    const int numDel = oldCigarState.Length();
                    const int numInsert = newCigarState.Length();
                    const int numMatch = std::min(numDel, numInsert);

                    if (numDel == numInsert) {
                        // Dest:   GC AA-- TC      GC AA TC
                        // Read:   GC --AA TC  ->  GC AA TC
                        //            DDII            MM
                        oldCigarState = CigarOperation();
                        newCigarState =
                            CigarOperation(CigarOperationType::SEQUENCE_MATCH, numMatch);

                    } else if (numDel > numInsert) {
                        // Dest:   GC AAA-- TC      GC AAA TC
                        // Read:   GC ---AA TC  ->  GC -AA TC
                        //            DDDII            DMM
                        oldCigarState =
                            CigarOperation(CigarOperationType::DELETION, numDel - numMatch);
                        newCigarState =
                            CigarOperation(CigarOperationType::SEQUENCE_MATCH, numMatch);

                    } else {
                        // Dest:   GC AA--- TC      GC AA- TC
                        // Read:   GC --AAA TC  ->  GC AAA TC
                        //            DDIII            MMI
                        oldCigarState =
                            CigarOperation(CigarOperationType::SEQUENCE_MATCH, numMatch);
                        newCigarState =
                            CigarOperation(CigarOperationType::INSERTION, numInsert - numMatch);
                    }
                }

                // I + D
                if (oldCigarState.Type() == CigarOperationType::INSERTION &&
                    newCigarState.Type() == CigarOperationType::DELETION) {
                    const int numInsert = oldCigarState.Length();
                    const int numDel = newCigarState.Length();
                    const int numMatch = std::min(numDel, numInsert);

                    if (numDel == numInsert) {
                        // Dest:   GC --AA TC  ->  GC AA TC
                        // Read:   GC AA-- TC      GC AA TC
                        //            IIDD            MM
                        oldCigarState = CigarOperation();
                        newCigarState =
                            CigarOperation(CigarOperationType::SEQUENCE_MATCH, numMatch);

                    } else if (numDel > numInsert) {
                        // Dest:   GC --AAA TC  ->  GC AAA TC
                        // Read:   GC AA--- TC      GC AA- TC
                        //            IIDDD            MMD
                        oldCigarState =
                            CigarOperation(CigarOperationType::SEQUENCE_MATCH, numMatch);
                        newCigarState =
                            CigarOperation(CigarOperationType::DELETION, numDel - numMatch);

                    } else {
                        // Dest:   GC ---AA TC  ->  GC -AA TC
                        // Read:   GC AAA-- TC      GC AAA TC
                        //            IIIDD            IMM
                        oldCigarState =
                            CigarOperation(CigarOperationType::INSERTION, numInsert - numMatch);
                        newCigarState =
                            CigarOperation(CigarOperationType::SEQUENCE_MATCH, numMatch);
                    }
                }

                if ((oldCigarState.Type() != CigarOperationType::UNKNOWN_OP)) {
                    newCigarTuple.push_back(oldCigarState);
                }
                // swap old and new state
                oldCigarState = newCigarState;
                newCigarState = CigarOperation(newState.Type(), 1);
            } else {
                newCigarState.Length(newCigarState.Length() + 1);
            }
        }

        //////////////////////////////////////
        // POST-PROCESSING                  //
        //////////////////////////////////////
        // check left flanking region + merge M-M pairs
        int i = 0;
        while (i < static_cast<int>(newCigarTuple.size()) - 1) {
            CigarOperation leftOp = newCigarTuple.at(i);
            CigarOperation rightOp = newCigarTuple.at(i + 1);

            // clang-format off
            // M + M:
            if (leftOp.Type() == CigarOperationType::SEQUENCE_MATCH && rightOp.Type() == CigarOperationType::SEQUENCE_MATCH) {
                newCigarTuple[i] = CigarOperation(CigarOperationType::SEQUENCE_MATCH, leftOp.Length() + rightOp.Length());
                newCigarTuple.erase(newCigarTuple.begin() + i + 1);
            }
            // S + I:
            else if (leftOp.Type() == CigarOperationType::SOFT_CLIP && rightOp.Type() == CigarOperationType::INSERTION) {
                newCigarTuple[i] = CigarOperation(CigarOperationType::SOFT_CLIP, leftOp.Length() + rightOp.Length());
                newCigarTuple.erase(newCigarTuple.begin() + i + 1);
            }
            // S + D:
            else if (leftOp.Type() == CigarOperationType::SOFT_CLIP && rightOp.Type() == CigarOperationType::DELETION) {
                newCigarTuple[i] = CigarOperation(CigarOperationType::SOFT_CLIP, leftOp.Length());
                newCigarTuple.erase(newCigarTuple.begin() + i + 1);
            }
            // S + P:
            else if (leftOp.Type() == CigarOperationType::SOFT_CLIP && rightOp.Type() == CigarOperationType::PADDING) {
                newCigarTuple[i] = CigarOperation(CigarOperationType::SOFT_CLIP, leftOp.Length());
                newCigarTuple.erase(newCigarTuple.begin() + i + 1);
            }
            // H + I:
            else if (leftOp.Type() == CigarOperationType::HARD_CLIP && rightOp.Type() == CigarOperationType::INSERTION) {
                newCigarTuple[i + 1] = CigarOperation(CigarOperationType::SOFT_CLIP, rightOp.Length());
                ++i;
            }
            // H + D:
            else if (leftOp.Type() == CigarOperationType::HARD_CLIP && rightOp.Type() == CigarOperationType::DELETION) {
                newCigarTuple[i] = CigarOperation(CigarOperationType::HARD_CLIP, leftOp.Length());
                newCigarTuple.erase(newCigarTuple.begin() + i + 1);
            }
            // H + P:
            else if (leftOp.Type() == CigarOperationType::HARD_CLIP && rightOp.Type() == CigarOperationType::PADDING) {
                newCigarTuple[i] = CigarOperation(CigarOperationType::HARD_CLIP, leftOp.Length());
                newCigarTuple.erase(newCigarTuple.begin() + i + 1);
            // H + S:
            } else {
                ++i;
            }
            // clang-format on
        }

        // check right flanking region
        i = newCigarTuple.size() - 2;
        while (i >= 0) {

            CigarOperation leftOp = newCigarTuple.at(i);
            CigarOperation rightOp = newCigarTuple.at(i + 1);

            if (leftOp.Type() == CigarOperationType::SEQUENCE_MATCH) {
                // reached a match state, hence everything
                // before will be compliant
                break;
            }
            // I + S:
            if (leftOp.Type() == CigarOperationType::INSERTION &&
                rightOp.Type() == CigarOperationType::SOFT_CLIP) {
                newCigarTuple[i] = CigarOperation(CigarOperationType::SOFT_CLIP,
                                                  leftOp.Length() + rightOp.Length());
                newCigarTuple.erase(newCigarTuple.begin() + i + 1);
            }
            // D + S:
            else if (leftOp.Type() == CigarOperationType::DELETION &&
                     rightOp.Type() == CigarOperationType::SOFT_CLIP) {
                newCigarTuple[i] = CigarOperation(CigarOperationType::SOFT_CLIP, rightOp.Length());
                newCigarTuple.erase(newCigarTuple.begin() + i + 1);
            }
            // P + S:
            else if (leftOp.Type() == CigarOperationType::PADDING &&
                     rightOp.Type() == CigarOperationType::SOFT_CLIP) {
                newCigarTuple[i] = CigarOperation(CigarOperationType::SOFT_CLIP, rightOp.Length());
                newCigarTuple.erase(newCigarTuple.begin() + i + 1);
            }
            // I + H:
            else if (leftOp.Type() == CigarOperationType::INSERTION &&
                     rightOp.Type() == CigarOperationType::HARD_CLIP) {
                newCigarTuple[i] = CigarOperation(CigarOperationType::SOFT_CLIP, leftOp.Length());
            }
            // D + H:
            else if (leftOp.Type() == CigarOperationType::DELETION &&
                     rightOp.Type() == CigarOperationType::HARD_CLIP) {
                newCigarTuple[i] = CigarOperation(CigarOperationType::HARD_CLIP, rightOp.Length());
                newCigarTuple.erase(newCigarTuple.begin() + i + 1);
            }
            // P + H:
            else if (leftOp.Type() == CigarOperationType::PADDING &&
                     rightOp.Type() == CigarOperationType::HARD_CLIP) {
                newCigarTuple[i] = CigarOperation(CigarOperationType::HARD_CLIP, rightOp.Length());
                newCigarTuple.erase(newCigarTuple.begin() + i + 1);
            }
            --i;
        }

        std::string new_seq = read.Sequence(BAM::Orientation::GENOMIC);

        // calculate edit distance (and possibly replace match states)
        posInRead = 0;
        posInDestRef = newSamStart;
        int new_edit_distance = 0;
        Cigar replace_cigar_tuple;

        const auto match_state_det = [](char read_base, char genome_base) {
            if (read_base == genome_base)
                return CigarOperationType::SEQUENCE_MATCH;
            else
                return CigarOperationType::SEQUENCE_MISMATCH;
        };

        for (const auto& op : newCigarTuple) {
            const CigarOperationType cigarOp = op.Type();
            const int cigarOp_count = op.Length();

            if (cigarOp == CigarOperationType::SEQUENCE_MATCH) {
                auto old_state =
                    match_state_det(new_seq.at(posInRead), toReferenceGapless_.at(posInDestRef));
                int count = 1;
                for (int i = 1; i < cigarOp_count; ++i) {
                    const auto next_state = match_state_det(
                        new_seq.at(posInRead + i), toReferenceGapless_.at(posInDestRef + i));
                    if (old_state != next_state) {
                        if (old_state == CigarOperationType::SEQUENCE_MISMATCH)
                            new_edit_distance += count;
                        replace_cigar_tuple.emplace_back(old_state, count);
                        old_state = next_state;
                        count = 1;
                    } else {
                        ++count;
                    }
                }

                if (old_state == CigarOperationType::SEQUENCE_MISMATCH) new_edit_distance += count;
                replace_cigar_tuple.emplace_back(old_state, count);

                posInRead += cigarOp_count;
                posInDestRef += cigarOp_count;
            } else if (cigarOp == CigarOperationType::INSERTION) {
                new_edit_distance += cigarOp_count;
                replace_cigar_tuple.emplace_back(cigarOp, cigarOp_count);
                posInRead += cigarOp_count;
            } else if (cigarOp == CigarOperationType::DELETION) {
                new_edit_distance += cigarOp_count;
                replace_cigar_tuple.emplace_back(cigarOp, cigarOp_count);
                posInDestRef += cigarOp_count;
            } else if (cigarOp == CigarOperationType::SOFT_CLIP) {
                replace_cigar_tuple.emplace_back(cigarOp, cigarOp_count);
                posInRead += cigarOp_count;
            } else if (cigarOp == CigarOperationType::HARD_CLIP ||
                       cigarOp == CigarOperationType::PADDING) {
                replace_cigar_tuple.emplace_back(cigarOp, cigarOp_count);
            } else {
                throw std::runtime_error("STATE should not occur " +
                                         std::to_string(static_cast<int>(cigarOp)));
            }
        }
        read.Impl().CigarData(replace_cigar_tuple);
        read.Impl().Position(newSamStart);
        if (read.Impl().HasTag("NM"))
            read.Impl().EditTag("NM", new_edit_distance);
        else
            read.Impl().AddTag("NM", new_edit_distance);
        out->Write(read);
    }
    out.reset(nullptr);
    BAM::PbiFile::CreateFrom(outputFile);
}
}
}  // ::PacBio::Realign