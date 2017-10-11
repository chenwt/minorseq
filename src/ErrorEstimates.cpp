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

#include <pacbio/juliet/ErrorEstimates.h>
#include <stdexcept>

namespace PacBio {
namespace Juliet {

ErrorEstimates::ErrorEstimates(const std::string& chemistry)
{
    Match = 0.9956844883;
    Substitution = 0.0005244257 / 3.0;
    Deletion = 0.003791086;
    Insertion = 0;
    if (chemistry == "P6-C4" || chemistry == "S/P1-C1/beta") {
        std::cerr << "+---------------------------------------------------+" << std::endl
                  << "|                     ATTENTION!                    |" << std::endl
                  << "| - - - - - - - - - - - - - - - - - - - - - - - - - |" << std::endl
                  << "|           This chemistry is unsupported.          |" << std::endl
                  << "|            Running in permissive mode.            |" << std::endl
                  << "|   Possibly increased type I and II error rates!   |" << std::endl
                  << "+---------------------------------------------------+" << std::endl;
    }
}

ErrorEstimates::ErrorEstimates(const double substitutionRate, const double deletionRate)
    : Match(1 - substitutionRate - deletionRate)
    , Substitution(substitutionRate / 3.0)
    , Deletion(deletionRate)
    , Insertion(0)
{
}
}
}