// Copyright (c) 2016-2017, Pacific Biosciences of California, Inc.
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

#include <limits>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include <pacbio/data/ArrayRead.h>
#include <pacbio/data/MSARow.h>
#include <pacbio/data/QvThresholds.h>

namespace PacBio {
namespace Data {

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
}
}  // ::PacBio::Juliet