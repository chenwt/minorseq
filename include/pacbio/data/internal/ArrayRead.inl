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

namespace PacBio {
namespace Data {

inline int ArrayRead::ReferenceStart() const { return referenceStart_; }
inline int ArrayRead::ReferenceEnd() const { return referenceEnd_; }
inline const std::vector<ArrayBase>& ArrayRead::Bases() const
{
    return const_cast<std::vector<ArrayBase>&>(bases_);
}
inline const std::string& ArrayRead::Name() const { return name_; }

inline std::string ArrayRead::SequencingChemistry() const { return ""; }

inline std::string BAMArrayRead::SequencingChemistry() const
{
    return Record.ReadGroup().SequencingChemistry();
}

inline std::ostream& operator<<(std::ostream& stream, const ArrayRead& r)
{
    stream << r.ReferenceStart() << std::endl;
    for (const auto& b : r.bases_)
        stream << b.Cigar;
    stream << std::endl;
    for (const auto& b : r.bases_)
        stream << b.Nucleotide;
    return stream;
}

inline bool ArrayBase::MeetQVThresholds(const QvThresholds& qvs) const
{
    return MeetQualQVThreshold(qvs.QualQV) && MeetDelQVThreshold(qvs.DelQV) &&
           MeetSubQVThreshold(qvs.SubQV) && MeetInsQVThreshold(qvs.InsQV);
}
inline bool ArrayBase::MeetQualQVThreshold(boost::optional<uint8_t> threshold) const
{
    return !threshold || !QualQV || *QualQV >= *threshold;
}
inline bool ArrayBase::MeetDelQVThreshold(boost::optional<uint8_t> threshold) const
{
    return !threshold || !DelQV || *DelQV >= *threshold;
}
inline bool ArrayBase::MeetSubQVThreshold(boost::optional<uint8_t> threshold) const
{
    return !threshold || !SubQV || *SubQV >= *threshold;
}
inline bool ArrayBase::MeetInsQVThreshold(boost::optional<uint8_t> threshold) const
{
    return !threshold || !InsQV || *InsQV >= *threshold;
}
}  // namespace Data
}  // namespace PacBio
