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

#include <numeric>

#include <pacbio/data/MSAColumn.h>

namespace PacBio {
namespace Data {
MSAColumn::MSAColumn(int refPos) : refPos_(refPos) {}

int MSAColumn::Coverage() const { return std::accumulate(counts_.cbegin(), counts_.cend(), 0); }

int MSAColumn::MaxElement() const
{
    return std::distance(counts_.begin(), std::max_element(counts_.begin(), counts_.end()));
}
char MSAColumn::MaxBase() const
{
    static const char bases[]{'A', 'C', 'G', 'T', '-'};
    int maxElement = MaxElement();
    if (maxElement == 5)
        return ' ';
    else
        return bases[maxElement];
}
int MSAColumn::Max() const { return counts_.at(MaxElement()); }

void MSAColumn::AddFisherResult(const FisherResult& f)
{
    pValues_ = f.PValues;
    mask_ = f.Mask;
    hit_ = f.Hit;
    argMax_ = f.ArgMax;
}

void MSAColumn::AddFisherResult(const std::map<std::string, double>& f) { insertionsPValues_ = f; }

std::ostream& MSAColumn::InDels(std::ostream& stream)
{
    stream << refPos_ << "\t";
    if (mask_.at(4) == 1) stream << "(-," << counts_.at(4) << "," << pValues_.at(4) << ")\t";
    for (const auto& bases_pvalue : insertionsPValues_)
        if (bases_pvalue.second < 0.01)
            stream << "(" << bases_pvalue.first << "," << insertions_.at(bases_pvalue.first) << ","
                   << bases_pvalue.second << ")\t";
    stream << std::endl;
    return stream;
}

std::vector<std::string> MSAColumn::SignificantInsertions() const
{
    std::vector<std::string> results;
    for (const auto& bases_pvalue : insertionsPValues_)
        if (bases_pvalue.second < 0.01) results.push_back(bases_pvalue.first);
    return results;
}

int MSAColumn::RefPos() const { return refPos_; }

void MSAColumn::IncInsertion(const std::string& seq) { insertions_[seq]++; }

const std::map<std::string, int>& MSAColumn::Insertions() const { return insertions_; }

double MSAColumn::PValue(const int i) const { return pValues_.at(i); }

std::ostream& operator<<(std::ostream& stream, const MSAColumn& r)
{
    for (int j = 0; j < 6; ++j)
        stream << r[j] << "\t" << r.PValue(j) << "\t";
    return stream;
}
}  // namespace Data
}  // namespace PacBio
