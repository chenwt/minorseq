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

#include <memory>
#include <vector>

#include <pacbio/data/MSA.h>
#include <pacbio/juliet/ErrorEstimates.h>
#include <pacbio/juliet/Haplotype.h>
#include <pacbio/juliet/TargetConfig.h>
#include <pacbio/juliet/VariantGene.h>
#include <pbcopper/json/JSON.h>

namespace PacBio {
namespace Juliet {
struct JulietSettings;

/// Stores performance metrics of juliet
struct PerformanceMetrics
{
public:
    double TruePositives = 0;
    double FalsePositives = 0;
    double FalseNegative = 0;
    double TrueNegative = 0;
    const double NumberOfTests;
    const double NumExpectedMinors;

public:
    PerformanceMetrics(double numberOfTests, double numExpectedMinors)
        : NumberOfTests(numberOfTests), NumExpectedMinors(numExpectedMinors)
    {
    }

public:
    std::string ToJson() const;

public:
    double TruePositiveRate() const;
    double FalsePositiveRate() const;
    double Accuracy() const;
};

/// Contains a single codon, it's abundance, and the translated amino acid
struct MajorityCall
{
    std::string Codon;
    int Coverage = 0;
    char AA;
};

/// Given a MSA, target config, and noise model, compute variant amino acids
/// and generate machine-interpretable output.
class AminoAcidCaller
{
public:
    AminoAcidCaller(const std::vector<std::shared_ptr<Data::ArrayRead>>& reads,
                    const ErrorEstimates& error, const JulietSettings& settings);

public:
    /// Generate JSON output of variant amino acids
    JSON::Json JSON();

public:
    void PhaseVariants();

private:
    /// Finds the major codon given the codon map
    static MajorityCall FindMajorityCodon(const std::map<std::string, int>& codons);

private:
    static constexpr float alpha = 0.01;
    void CallVariants();

    /// Counts the number of tests that will be performed.
    /// This number can be used to bonferroni correct p-values.
    int CountNumberOfTests(const std::vector<TargetGene>& genes) const;

    /// Find those drugs associated with the current variant and generate
    /// a summary string.
    std::string FindDRMs(const std::string& geneName, const std::vector<TargetGene>& genes,
                         const DMutation curDRM) const;

    /// Compute the probability that the two strings generated each other
    /// via sequencing noise.
    double Probability(const std::string& a, const std::string& b);

    /// Compute if the current variant hits an expected minor and
    /// use it to measure the performance of juliet.
    bool MeasurePerformance(const TargetGene& tg, const std::string& codon,
                            const bool& variableSite, const int& aaPos, const double& p,
                            PerformanceMetrics* pm);

private:
    Data::MSAByRow msaByRow_;

public:
    Data::MSAByColumn msaByColumn_;

private:
    std::vector<VariantGene> variantGenes_;
    std::vector<Haplotype> reconstructedHaplotypes_;
    std::vector<Haplotype> filteredHaplotypes_;
    const ErrorEstimates error_;
    const TargetConfig targetConfig_;
    const bool verbose_;
    const bool debug_;
    const bool drmOnly_;
    const double minimalPerc_;
    const double maximalPerc_;

    int genCounts_ = 0;
    int margWithGap_ = 0;
    int margWithHetero_ = 0;
    int margPartial_ = 0;
    int lowCov_ = 0;
    int margOfftarget_ = 0;
};
}
}  // ::PacBio::Juliet

#include "pacbio/juliet/internal/AminoAcidCaller.inl"