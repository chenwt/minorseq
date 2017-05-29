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

struct Haplotype
{
    std::string Name;
    std::vector<std::string> Names;
    std::vector<std::string> Codons;
    double SoftCollapses = 0;
    double GlobalFrequency = 0;
    int Flags = 0;

    double Size() const { return Names.size() + SoftCollapses; }

    std::string ConcatCodons() const
    {
        return std::accumulate(Codons.begin(), Codons.end(), std::string(""));
    }

    void SetCodons(std::vector<std::string>&& codons)
    {
        Codons = std::forward<std::vector<std::string>>(codons);
        for (const auto& c : Codons) {
            if (c.find('-') != std::string::npos)
                Flags |= static_cast<int>(HaplotypeType::WITH_GAP);
            if (c.find('N') != std::string::npos)
                Flags |= static_cast<int>(HaplotypeType::WITH_HETERODUPLEX);
            if (c.find(' ') != std::string::npos) Flags |= static_cast<int>(HaplotypeType::PARTIAL);
        }
    }

    friend std::ostream& operator<<(std::ostream& stream, const Haplotype& h)
    {
        stream << h.Size() << "\t";
        for (const auto& c : h.Codons) {
            if (c != "ATG" && c != "AAA" && c != "TAT" && c != "GGA" && c != "ACC")
                stream << termcolor::white;
            if (c == "TTG") stream << termcolor::red;
            if (c == "AGA") stream << termcolor::green;
            if (c == "TGT") stream << termcolor::blue;
            if (c == "GCA") stream << termcolor::blue;
            if (c == "TAC") stream << termcolor::yellow;
            stream << " " << c << termcolor::reset;
        }
        return stream;
    }

    JSON::Json ToJson() const
    {
        using namespace JSON;
        Json root;
        root["name"] = Name;
        root["reads_hard"] = Names.size();
        root["reads_soft"] = Size();
        root["frequency"] = GlobalFrequency;
        root["read_names"] = Names;
        root["codons"] = Codons;
        return root;
    }
};
}
}