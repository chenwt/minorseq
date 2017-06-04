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
class MSAByRow;

/// Represents a multiple sequence alignment (MSA) via individual rows.
/// Insertions are omitted and saved a special variable of each row.
class MSAByRow
{
public:
    MSAByRow() = default;
    MSAByRow(const std::vector<std::shared_ptr<Data::ArrayRead>>& reads);
    MSAByRow(const std::vector<Data::ArrayRead>& reads);

public:
    /// The left-most position of all reads in the MSA.
    int BeginPos() const { return beginPos_; }
    /// The right-most position of all reads in the MSA.
    int EndPos() const { return endPos_; }
    /// The individual rows of the MSA.
    const std::vector<std::shared_ptr<MSARow>>& Rows() const { return rows_; }
    /// Easy way to access a row by its name.
    /// The name equivalent to BamRecord::FullName()
    std::shared_ptr<MSARow> NameToRow(const std::string& name) const { return nameToRow_.at(name); }

private:
    std::vector<std::shared_ptr<MSARow>> rows_;
    std::map<std::string, std::shared_ptr<MSARow>> nameToRow_;
    const Data::QvThresholds qvThresholds_;
    int beginPos_ = std::numeric_limits<int>::max();
    int endPos_ = 0;

private:
    /// Adds an ArrayRead to the MSA, by converting it to a MSARow object.
    MSARow AddRead(const Data::ArrayRead& read);

    /// Updates the begin and end positions of the MSA, if this read
    /// is extending the current boundaries.
    void UpdateBoundaries(const Data::ArrayRead& read);
};

/// A particular row of a MSA.
struct MSARow
{
public:
    MSARow(const int size) : Bases(size, ' ') {}
    /// Invididual bases with '-' as deletion.
    std::vector<char> Bases;
    /// Position to insertion string.
    std::map<int, std::string> Insertions;
    /// The underlying ArrayRead.
    std::shared_ptr<Data::ArrayRead> Read;

public:
    bool CodonAt(const int winPos, std::string* codon) const;
};

/// Represents a MSA by columns. Each column is a distribution of counts.
/// Offers iterators supporting a for each loop.
/// Index parameters are in ABSOLUTE reference space.
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
    MSAColumn& operator[](int i) const;
    /// Checks if the index is available.
    bool has(int i);

    MsaIt begin();
    MsaIt end();
    MsaItConst begin() const;
    MsaItConst end() const;
    MsaItConst cbegin() const;
    MsaItConst cend() const;

public:
    /// The left-most position of all reads in the MSA.
    int BeginPos() const { return beginPos_; }
    /// The right-most position of all reads in the MSA.
    int EndPos() const { return endPos_; }

private:
    MsaVec counts;
    int beginPos_ = std::numeric_limits<int>::max();
    int endPos_ = 0;
};

/// Represents a single MSA column with counts for each nucleotide, insertions,
/// and results of the Fisher's Exact test can be associated.
/// Nucleotide alphabet is {A, C, G, T, -, N}.
class MSAColumn
{
public:
    MSAColumn(int refPos);

public:
    /// Relative abundance for given nucleotide.
    double Frequency(char c) const;

    /// Access counts for a given nucleotide.
    int operator[](char c) const;

    // operator std::array<int, 6>();
    explicit operator int();

public:
    /// Coverage including deletions and Ns.
    int Coverage() const;
    /// Nucleotide with the highest counts.
    char MaxBase() const;
    /// Position in the absolute reference space.
    int RefPos() const;
    /// Insertions called significantly abundant
    std::vector<std::string> SignificantInsertions() const;
    /// Insertions and respective counts
    const std::map<std::string, int>& Insertions() const;
    /// P-value for given nucleotide.
    double PValue(const char c) const;

public:  // NOT YET USED
    void IncInsertion(const std::string& seq);
    void AddFisherResult(const FisherResult& f);
    void AddFisherResult(const std::map<std::string, double>& f);
    std::ostream& InDels(std::ostream& stream);

public:
    friend std::ostream& operator<<(std::ostream& stream, const MSAColumn& r);

private:
    std::array<int, 6> counts_{{0, 0, 0, 0, 0, 0}};
    std::map<std::string, int> insertions_;
    std::map<std::string, double> insertionsPValues_;
    std::array<double, 6> pValues_{{1, 1, 1, 1, 1, 1}};
    std::array<double, 6> mask_{{0, 0, 0, 0, 0, 0}};
    bool hit_ = false;
    int argMax_ = 0;
    int refPos_ = -1;

private:
    /// Relative per nucleotide abundance for given index. Index is wrt the
    /// underlying data structures.
    double Frequency(int i) const;
    // The maximal element of all nucleotide, as index.
    int MaxElement() const;
    // Increase counts for given nucleotide.
    void IncCounts(const char c);

public:
    friend MSAByColumn;
};
}  // namespace Data
}  // namespace PacBio

#include "pacbio/data/internal/MSA.inl"