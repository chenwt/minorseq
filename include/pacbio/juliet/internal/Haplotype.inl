// Copyright (c) 2017, Pacific Biosciences of California, Inc.
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

namespace PacBio {
namespace Juliet {

inline double Haplotype::Size() const { return readNames_.size() + softCollapses_; }

inline const std::vector<std::string>& Haplotype::ReadNames() const
{
    return const_cast<std::vector<std::string>&>(readNames_);
}

inline const std::string& Haplotype::Codon(const int i) { return codons_.at(i); }

inline size_t Haplotype::NumCodons() const { return numCodons_; }

inline int Haplotype::Flags() const { return flags_; }

inline std::string Haplotype::Name() { return name_; }

inline void Haplotype::AddFlag(const HaplotypeType& flag) { flags_ |= static_cast<int>(flag); }

inline void Haplotype::Frequency(const double& freq) { frequency_ = freq; }

inline void Haplotype::AddReadName(const std::string& name) { readNames_.push_back(name); }

inline void Haplotype::AddSoftReadCount(const double s) { softCollapses_ += s; }

inline void Haplotype::Name(const std::string& name) { name_ = name; }

inline std::ostream& operator<<(std::ostream& stream, const Haplotype& h)
{
    stream << h.Size() << "\t";
    for (const auto& c : h.codons_) {
        stream << " " << c;
    }
    return stream;
}
}
}
