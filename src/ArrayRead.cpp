// Copyright (c) 2011-2017, Pacific Biosciences of California, Inc.
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

#include <algorithm>
#include <array>
#include <iostream>
#include <string>
#include <vector>

#include <pbbam/BamRecord.h>

#include <pacbio/data/ArrayRead.h>

namespace PacBio {
namespace Data {

ArrayRead::ArrayRead(const int idx, const std::string& name) : idx_(idx), name_(name){};

BAMArrayRead::BAMArrayRead(const BAM::BamRecord& record, int idx)
    : ArrayRead(idx, record.FullName())
    , Record(record)  // Record(std::forward<BAM::BamRecord>(record))
{
    ArrayRead::referenceStart_ = record.ReferenceStart();
    ArrayRead::referenceEnd_ = record.ReferenceEnd();
    const auto seq = Record.Sequence(BAM::Orientation::GENOMIC, true, true);

    bool hasQualities = !Record.Qualities().empty();
    BAM::QualityValues qual;
    if (hasQualities) qual = Record.Qualities(BAM::Orientation::GENOMIC, true, true);

    std::string cigar;
    cigar.reserve(seq.size());
    for (const auto c : Record.CigarData(true))
        for (size_t i = 0; i < c.Length(); ++i)
            cigar += c.Char();

    BAM::QualityValues subQV;
    BAM::QualityValues delQV;
    BAM::QualityValues insQV;

    bool richQVs = Record.HasSubstitutionQV() && Record.HasDeletionQV() && Record.HasInsertionQV();
    if (richQVs) {
        subQV = Record.SubstitutionQV(BAM::Orientation::GENOMIC, true, true);
        delQV = Record.DeletionQV(BAM::Orientation::GENOMIC, true, true);
        insQV = Record.InsertionQV(BAM::Orientation::GENOMIC, true, true);
    }

    assert(cigar.size() == seq.size());

    if (hasQualities) {
        assert(seq.size() == qual.size());
    }

    bases_.reserve(cigar.length());
    if (richQVs)
        for (size_t i = 0; i < cigar.length(); ++i)
            bases_.emplace_back(cigar.at(i), seq.at(i), qual.at(i), subQV.at(i), delQV.at(i),
                                insQV.at(i));
    else if (hasQualities)
        for (size_t i = 0; i < cigar.length(); ++i)
            bases_.emplace_back(cigar.at(i), seq.at(i), qual.at(i));
    else
        for (size_t i = 0; i < cigar.length(); ++i)
            bases_.emplace_back(cigar.at(i), seq.at(i), 0);
}

ArrayBase::ArrayBase(char cigar, char nucleotide, uint8_t qualQV, uint8_t subQV, uint8_t delQV,
                     uint8_t insQV)
    : Cigar(cigar)
    , Nucleotide(nucleotide)
    , QualQV(qualQV)
    , DelQV(delQV)
    , SubQV(subQV)
    , InsQV(insQV)
    , ProbTrue(1 - pow(10, -1.0 * qualQV / 10.0))
    , ProbCorrectBase(1 - pow(10, -1.0 * subQV / 10.0))
    , ProbNoDeletion(1 - pow(10, -1.0 * delQV / 10.0))
    , ProbNoInsertion(1 - pow(10, -1.0 * insQV / 10.0))
{
}
ArrayBase::ArrayBase(char cigar, char nucleotide, uint8_t qualQV)
    : Cigar(cigar)
    , Nucleotide(nucleotide)
    , QualQV(qualQV)
    , ProbTrue(1 - pow(10, -1.0 * qualQV / 10.0))
{
}
ArrayBase::ArrayBase(char cigar, char nucleotide) : Cigar(cigar), Nucleotide(nucleotide) {}

}  // namespace Data
}  // namespace PacBio
