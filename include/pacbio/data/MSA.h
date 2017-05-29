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

#include <pacbio/data/QvThresholds.h>

#include <array>
#include <map>
#include <string>
#include <vector>

namespace PacBio {
namespace Data {
class ArrayRead;
struct FisherResult;

class MSAByColumn;
class MSAColumn;
struct MSARow;

class MSAByRow
{
public:
    MSAByRow() = default;
    MSAByRow(const std::vector<std::shared_ptr<Data::ArrayRead>>& reads);
    MSAByRow(const std::vector<Data::ArrayRead>& reads);

public:
    void BeginEnd(const Data::ArrayRead& read);
    int BeginPos() const { return beginPos_; }
    int EndPos() const { return endPos_; }
    const std::vector<std::shared_ptr<MSARow>>& Rows() const { return rows_; }
    std::shared_ptr<MSARow> NameToRow(const std::string& name) const { return nameToRow_.at(name); }

private:
    std::vector<std::shared_ptr<MSARow>> rows_;
    std::map<std::string, std::shared_ptr<MSARow>> nameToRow_;
    const Data::QvThresholds qvThresholds_;
    int beginPos_ = std::numeric_limits<int>::max();
    int endPos_ = 0;

private:
    MSARow AddRead(const Data::ArrayRead& read);
};

struct MSARow
{
    MSARow(const int size) : Bases(size, ' ') {}
    std::vector<char> Bases;
    std::map<int, std::string> Insertions;
    std::shared_ptr<Data::ArrayRead> Read;
};

/// Multiple sequence alignment containing counts
class MSAByColumn
{
private:
    using MsaVec = std::vector<MSAColumn>;

public:
    using MsaIt = MsaVec::iterator;
    using MsaItConst = MsaVec::const_iterator;

public:
    MSAByColumn(const MSAByRow& nucMat);

public:
    /// Parameter is an index in ABSOLUTE reference space
    MSAColumn operator[](int i) const;
    /// Parameter is an index in ABSOLUTE reference space
    MSAColumn& operator[](int i);

    bool has(int i);

    MsaIt begin();
    MsaIt end();
    MsaItConst begin() const;
    MsaItConst end() const;
    MsaItConst cbegin() const;
    MsaItConst cend() const;

public:
    int BeginPos() const { return beginPos_; }
    int EndPos() const { return endPos_; }

private:
    MsaVec counts;
    int beginPos_ = std::numeric_limits<int>::max();
    int endPos_ = 0;

private:
    void BeginEnd(const Data::ArrayRead& read);
    void FillCounts(const ArrayRead& read, const QvThresholds& qvThresholds);
};

class MSAColumn
{
public:
    MSAColumn(int refPos);

public:
    // Relative per nucleotide abundance
    double Frequency(int i) const;
    double Frequency(char c) const;

    // Nucleotide counts_
    int operator[](int i) const;
    int& operator[](int i);
    int operator[](char c) const;
    int& operator[](char c);

    operator std::array<int, 6>();
    explicit operator int();

public:
    int Coverage() const;
    char MaxBase() const;
    int MaxElement() const;
    int Max() const;
    int RefPos() const;
    std::vector<std::string> SignificantInsertions() const;
    const std::map<std::string, int>& Insertions() const;
    double PValue(const int i) const;

public:
    void AddFisherResult(const FisherResult& f);
    void AddFisherResult(const std::map<std::string, double>& f);
    std::ostream& InDels(std::ostream& stream);
    void IncInsertion(const std::string& seq);

private:
    std::array<int, 6> counts_{{0, 0, 0, 0, 0, 0}};
    std::map<std::string, int> insertions_;
    std::map<std::string, double> insertionsPValues_;
    std::array<double, 6> pValues_{{1, 1, 1, 1, 1, 1}};
    std::array<double, 6> mask_{{0, 0, 0, 0, 0, 0}};
    bool hit_ = false;
    int argMax_ = 0;
    int refPos_ = -1;
};
}  // namespace Data
}  // namespace PacBio
