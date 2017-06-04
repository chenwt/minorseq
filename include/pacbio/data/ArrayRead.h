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

#pragma once

#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

#include <boost/optional.hpp>

#include <pbbam/BamRecord.h>

#include <pacbio/data/NucleotideConversion.h>
#include <pacbio/data/QvThresholds.h>

namespace PacBio {
namespace Data {

struct ArrayBase;

/// A single array read that is "unrolled", as in an array of bases.
class ArrayRead
{
public:  // ctors
    ArrayRead(const int idx = -1, const std::string& name = "");

public:  // non-mod methods
    int ReferenceStart() const;
    int ReferenceEnd() const;
    const std::vector<ArrayBase>& Bases() const;
    const std::string& Name() const;
    virtual std::string SequencingChemistry() const;

public:
    friend std::ostream& operator<<(std::ostream& stream, const ArrayRead& r);

protected:
    std::vector<ArrayBase> bases_;
    const int idx_;
    const std::string name_;
    size_t referenceStart_;
    size_t referenceEnd_;
};

class BAMArrayRead : public ArrayRead
{
public:  // ctors
    /// Constructor that needs the BamRecord to be "unrolled" and a unique index
    BAMArrayRead(const BAM::BamRecord& record, int idx);

    virtual std::string SequencingChemistry() const override;

private:
    const BAM::BamRecord Record;
};

/// A single base in an ArrayRead with its associated qvs and cigar
struct ArrayBase
{
    ArrayBase(char cigar, char nucleotide, uint8_t qualQV, uint8_t subQV, uint8_t delQV,
              uint8_t insQV);
    ArrayBase(char cigar, char nucleotide, uint8_t qualQV);
    ArrayBase(char cigar, char nucleotide);

    bool MeetQVThresholds(const QvThresholds& qvs) const;
    bool MeetQualQVThreshold(boost::optional<uint8_t> threshold) const;
    bool MeetDelQVThreshold(boost::optional<uint8_t> threshold) const;
    bool MeetSubQVThreshold(boost::optional<uint8_t> threshold) const;
    bool MeetInsQVThreshold(boost::optional<uint8_t> threshold) const;

    char Cigar;
    char Nucleotide;
    boost::optional<uint8_t> QualQV;
    boost::optional<uint8_t> DelQV;
    boost::optional<uint8_t> SubQV;
    boost::optional<uint8_t> InsQV;
    double ProbTrue = 0;
    double ProbCorrectBase = 0;
    double ProbNoDeletion = 0;
    double ProbNoInsertion = 0;
};
}  // namespace Data
}  // namespace PacBio

#include "pacbio/data/internal/ArrayRead.inl"