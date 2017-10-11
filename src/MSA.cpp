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

#include <limits>
#include <map>
#include <memory>
#include <numeric>
#include <string>
#include <vector>

#include <pacbio/data/ArrayRead.h>
#include <pacbio/data/FisherResult.h>
#include <pacbio/juliet/AminoAcidTable.h>

#include <pacbio/data/MSA.h>

namespace PacBio {
namespace Data {
MSAColumn::MSAColumn(int refPos) : refPos_(refPos) {}

MSAByColumn::MSAByColumn(const MSAByRow& msaRows)
{
    beginPos_ = msaRows.BeginPos() - 1;
    endPos_ = msaRows.EndPos() - 1;

    int pos = msaRows.BeginPos();
    int size = msaRows.EndPos() - msaRows.BeginPos();
    counts.reserve(size);
    for (int i = 0; i < size; ++i) {
        counts.emplace_back(pos);
        ++pos;
    }

    for (const auto& row : msaRows.Rows()) {
        int localPos = 0;
        for (const auto& c : row->Bases) {
            switch (c) {
                case 'A':
                case 'C':
                case 'G':
                case 'T':
                case '-':
                case 'N':
                    counts.at(localPos).IncCounts(c);
                    ++localPos;
                    break;
                case ' ':
                    ++localPos;
                    break;
                default:
                    throw std::runtime_error("Unexpected base " + std::string(1, c));
            }
        }
        for (const auto& ins : row->Insertions) {
            counts[ins.first].IncInsertion(ins.second);
        }
    }
}

MSAByRow::MSAByRow(const std::vector<std::shared_ptr<Data::ArrayRead>>& reads)
{
    for (const auto& r : reads)
        UpdateBoundaries(*r);

    for (const auto& r : reads) {
        auto row = AddRead(*r);
        row.Read = r;
        const auto x = std::make_shared<MSARow>(std::move(row));
        nameToRow_[r->Name()] = x;
        rows_.emplace_back(x);
    }

    ++beginPos_;
    ++endPos_;
}

MSAByRow::MSAByRow(const std::vector<Data::ArrayRead>& reads)
{
    for (const auto& r : reads)
        UpdateBoundaries(r);

    for (const auto& r : reads) {
        const auto x = std::make_shared<MSARow>(AddRead(r));
        nameToRow_[r.Name()] = x;
        rows_.emplace_back(x);
    }

    ++beginPos_;
    ++endPos_;
}

std::map<std::string, int> MSAByRow::CodonsAt(const int i) const
{
    std::map<std::string, int> codons;

    std::string codon;
    for (const auto& row : rows_) {
        if (row->CodingCodonAt(i, &codon)) ++codons[codon];
    }

    return codons;
}

void MSAByRow::UpdateBoundaries(const Data::ArrayRead& read)
{
    beginPos_ = std::min(beginPos_, read.ReferenceStart());
    endPos_ = std::max(endPos_, read.ReferenceEnd());
}

MSARow MSAByRow::AddRead(const Data::ArrayRead& read)
{
    MSARow row(endPos_ - beginPos_);

    int pos = read.ReferenceStart() - beginPos_;
    assert(pos >= 0);

    std::string insertion;
    auto CheckInsertion = [&insertion, &row, &pos]() {
        if (insertion.empty()) return;
        row.Insertions[pos] = insertion;
        insertion = "";
    };

    for (const auto& b : read.Bases()) {
        switch (b.Cigar) {
            case 'X':
            case '=':
                CheckInsertion();
                if (b.MeetQVThresholds(qvThresholds_))
                    row.Bases[pos++] = b.Nucleotide;
                else
                    row.Bases[pos++] = 'N';
                break;
            case 'D':
                CheckInsertion();
                row.Bases[pos++] = '-';
                break;
            case 'I':
                insertion += b.Nucleotide;
                break;
            case 'P':
                CheckInsertion();
                break;
            case 'S':
                CheckInsertion();
                break;
            default:
                throw std::runtime_error("Unexpected cigar " + std::to_string(b.Cigar));
        }
    }
    return row;
}

std::string MSARow::CodonAt(const int pos) const
{
    int length = Bases.size();
    std::string codon;
    for (int i = pos; i < pos + 3; ++i)
        if (i > 0 && i < length) codon += Bases.at(i);
    return codon;
}

bool MSARow::CodingCodonAt(const int winPos, std::string* codon) const
{
    using AAT = Juliet::AminoAcidTable;

    const auto CodonContains = [this, &winPos](const char x) {
        return (Bases.at(winPos + 0) == x || Bases.at(winPos + 1) == x ||
                Bases.at(winPos + 2) == x);
    };

    // Read does not cover codon
    if (winPos + 2 >= static_cast<int>(Bases.size()) || winPos < 0) return false;
    if (CodonContains(' ')) return false;

    // Read has a deletion
    if (CodonContains('-')) return false;

    auto proposedCodon = CodonAt(winPos);
    if (AAT::FromCodon.find(proposedCodon) == AAT::FromCodon.cend()) return false;

    *codon = proposedCodon;

    return true;
}

int MSAColumn::Coverage() const { return std::accumulate(counts_.cbegin(), counts_.cend(), 0); }

std::vector<std::string> MSAColumn::SignificantInsertions() const
{
    std::vector<std::string> results;
    for (const auto& bases_pvalue : insertionsPValues_)
        if (bases_pvalue.second < 0.01) results.push_back(bases_pvalue.first);
    return results;
}

void MSAColumn::AddFisherResult(const FisherResult& f)
{
    pValues_ = f.PValues;
    mask_ = f.Mask;
    hit_ = f.Hit;
    argMax_ = f.ArgMax;
}
}  // namespace Data
}  // namespace PacBio
