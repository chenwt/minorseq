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

#pragma once

#include <pacbio/util/Termcolor.h>

#include <pbcopper/json/JSON.h>

namespace PacBio {
namespace Juliet {
enum class HaplotypeType : int
{
    REPORT = 0,
    WITH_GAP = 1,
    WITH_HETERODUPLEX = 2,
    PARTIAL = 4,
    LOW_COV = 8,
    OFFTARGET = 16
};

class Haplotype
{
public:
    Haplotype() = delete;
    Haplotype(const std::string readName, const std::vector<std::string>& codons,
              const HaplotypeType& flag)
        : readNames_({readName}), codons_(codons), numCodons_(codons_.size())
    {
        AddFlag(flag);
        SetFlagsByCodons();
    }
    Haplotype(const std::vector<std::string> readNames, std::vector<std::string>&& codons,
              const HaplotypeType flag)
        : readNames_(readNames)
        , codons_(std::forward<std::vector<std::string>>(codons))
        , numCodons_(codons_.size())
    {
        AddFlag(flag);
        SetFlagsByCodons();
    }

public:  // non-mod methods
    /// How many reads contributed to this haplotype
    double Size() const;
    /// Concat all codons to one string without seperator
    std::string ConcatCodons() const;
    /// Convert this to a JSON string
    JSON::Json ToJson() const;
    // All read names
    const std::vector<std::string>& ReadNames() const;
    // All codons
    const std::string& Codon(const int i);
    // Number of codons
    size_t NumCodons() const;
    // Combined flags
    int Flags() const;
    // Name of this haplotype
    std::string Name();

public:  // mod methods
    /// Set appropriate HaplotypeFlags from already stored codons
    void SetFlagsByCodons();
    /// Add HaplotypeFlag to this haplotype
    void AddFlag(const HaplotypeType& flag);
    /// Set the frequency of this
    void Frequency(const double& freq);
    /// Add additional read
    void AddReadName(const std::string& name);
    /// Add a fraction of reads as soft counts
    void AddSoftReadCount(const double s);
    /// Set name of this haplotype
    void Name(const std::string& name);

public:
    friend std::ostream& operator<<(std::ostream& stream, const Haplotype& h);

private:
    std::string name_;
    std::vector<std::string> readNames_;
    const std::vector<std::string> codons_;
    size_t numCodons_;
    double softCollapses_ = 0;
    double frequency_ = 0;
    int flags_ = 0;
};
}
}

#include "pacbio/juliet/internal/Haplotype.inl"