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

#include <pacbio/data/NucleotideConversion.h>

namespace PacBio {
namespace Data {
inline MSAColumn& MSAByColumn::operator[](int i) const
{
    return const_cast<MSAColumn&>(counts[i - beginPos_]);
}

inline bool MSAByColumn::has(int i) { return i >= beginPos_ && i < endPos_; }

inline MSAByColumn::MsaIt MSAByColumn::begin() { return counts.begin(); }
inline MSAByColumn::MsaIt MSAByColumn::end() { return counts.end(); }
inline MSAByColumn::MsaItConst MSAByColumn::begin() const { return counts.begin(); }
inline MSAByColumn::MsaItConst MSAByColumn::end() const { return counts.end(); }
inline MSAByColumn::MsaItConst MSAByColumn::cbegin() const { return counts.cbegin(); }
inline MSAByColumn::MsaItConst MSAByColumn::cend() const { return counts.cend(); }

inline double MSAColumn::Frequency(int i) const
{
    return (*this)[i] / static_cast<double>(Coverage());
}
inline double MSAColumn::Frequency(char c) const { return Frequency(NucleotideToTag(c)); }

inline int MSAColumn::operator[](char c) const { return counts_[NucleotideToTag(c)]; }

inline MSAColumn::operator int() { return Coverage(); }

inline void MSAColumn::IncCounts(const char c) { ++counts_[NucleotideToTag(c)]; }

inline int MSAColumn::MaxElement() const
{
    return std::distance(counts_.begin(), std::max_element(counts_.begin(), counts_.end()));
}
inline char MSAColumn::MaxBase() const
{
    static const char bases[]{'A', 'C', 'G', 'T', '-'};
    int maxElement = MaxElement();
    if (maxElement == 5)
        return ' ';
    else
        return bases[maxElement];
}

inline void MSAColumn::AddFisherResult(const std::map<std::string, double>& f)
{
    insertionsPValues_ = f;
}

inline std::ostream& MSAColumn::InDels(std::ostream& stream)
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

inline int MSAColumn::RefPos() const { return refPos_; }

inline void MSAColumn::IncInsertion(const std::string& seq) { insertions_[seq]++; }

inline const std::map<std::string, int>& MSAColumn::Insertions() const { return insertions_; }

inline double MSAColumn::PValue(const char c) const { return pValues_.at(NucleotideToTag(c)); }

inline std::ostream& operator<<(std::ostream& stream, const MSAColumn& r)
{
    for (const auto& b : {'A', 'C', 'G', 'T', '-'})
        stream << r[b] << "\t";
    return stream;
}
}  // namespace Data
}  // namespace PacBio
