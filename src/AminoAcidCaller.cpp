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

#include <array>
#include <cmath>
#include <cstdlib>
#include <exception>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <list>
#include <map>
#include <memory>
#include <numeric>
#include <set>
#include <sstream>
#include <unordered_map>
#include <vector>

#include <boost/optional.hpp>

#include <pacbio/data/ArrayRead.h>
#include <pacbio/juliet/AminoAcidCaller.h>
#include <pacbio/juliet/AminoAcidTable.h>
#include <pacbio/juliet/ErrorEstimates.h>
#include <pacbio/juliet/JulietSettings.h>
#include <pacbio/statistics/Fisher.h>
#include <pacbio/util/Termcolor.h>
#include <pbcopper/json/JSON.h>

namespace PacBio {
namespace Juliet {
using AAT = AminoAcidTable;

AminoAcidCaller::AminoAcidCaller(const std::vector<std::shared_ptr<Data::ArrayRead>>& reads,
                                 const ErrorEstimates& error, const JulietSettings& settings)
    : msaByRow_(reads)
    , msaByColumn_(msaByRow_)
    , error_(error)
    , targetConfig_(settings.TargetConfigUser)
    , verbose_(settings.Verbose)
    , debug_(settings.Debug)
    , drmOnly_(settings.DRMOnly)
    , minimalPerc_(settings.MinimalPerc)
    , maximalPerc_(settings.MaximalPerc)
{
    CallVariants();
}

int AminoAcidCaller::CountNumberOfTests(const std::vector<TargetGene>& genes) const
{
    int numberOfTests = 0;
    for (const auto& gene : genes) {
        for (int i = gene.begin; i < gene.end - 2; ++i) {
            // Relative to gene begin
            const int relPos = i - gene.begin;
            // Only work on beginnings of a codon
            if (relPos % 3 != 0) continue;
            // Relative to window begin
            const int winPos = i - msaByRow_.BeginPos();
            // Gather all observed codons and count number of different codons
            std::map<std::string, int> codons = msaByRow_.CodonsAt(winPos);
            numberOfTests += codons.size();
        }
    }
    return numberOfTests == 0 ? 1 : numberOfTests;
}

std::string AminoAcidCaller::FindDRMs(const std::string& geneName,
                                      const std::vector<TargetGene>& genes,
                                      const DMutation curDRM) const
{
    std::string drmSummary;
    for (const auto& gene : genes) {
        if (geneName == gene.name) {
            for (const auto& drms : gene.drms) {

                if (std::find(drms.positions.cbegin(), drms.positions.cend(), curDRM) !=
                    drms.positions.cend()) {
                    if (!drmSummary.empty()) drmSummary += " + ";
                    drmSummary += drms.name;
                }
            }
            break;
        }
    }
    return drmSummary;
};

void AminoAcidCaller::PhaseVariants()
{
    // Store variant positions by their absolute position
    std::vector<std::pair<int, std::shared_ptr<VariantGene::VariantPosition>>> variantPositions;
    for (const auto& vg : variantGenes_) {
        for (const auto& pos_vp : vg.relPositionToVariant)
            if (pos_vp.second->IsVariant())
                variantPositions.emplace_back(
                    std::make_pair(vg.geneOffset + pos_vp.first * 3, pos_vp.second));
    }

    // Print positions
    if (verbose_) {
        std::cerr << "Variant positions:";
        for (const auto& pos_var : variantPositions)
            std::cerr << " " << pos_var.first;
        std::cerr << std::endl;
    }

    // Storage for all observed haplotypes, for now
    std::vector<std::shared_ptr<Haplotype>> observations;

    // For each read
    for (const auto& row : msaByRow_.Rows()) {
        // Get all codons for this row
        std::vector<std::string> codons;
        HaplotypeType flag = HaplotypeType::REPORT;
        for (const auto& pos_var : variantPositions) {
            const std::string codon = row->CodonAt(pos_var.first - msaByRow_.BeginPos() - 3);

            // If this codon is not a variant, flag haplotype as off-target
            if (!pos_var.second->IsHit(codon)) flag = HaplotypeType::OFFTARGET;

            codons.emplace_back(std::move(codon));
        }

        // There are already haplotypes to compare against
        int miss = true;

        // Compare current row to existing haplotypes
        auto CompareHaplotypes = [&miss, &codons,
                                  &row](std::vector<std::shared_ptr<Haplotype>>& haplotypes) {
            for (auto& h : haplotypes) {
                // Don't trust if the number of codons differ.
                // That should only be the case if reads are not full-spanning.
                if (h->NumCodons() != codons.size()) {
                    continue;
                }
                bool same = true;
                for (size_t i = 0; i < codons.size(); ++i) {
                    if (h->Codon(i) != codons.at(i)) {
                        same = false;
                        break;
                    }
                }
                if (same) {
                    h->AddReadName(row->Read->Name());
                    miss = false;
                    break;
                }
            }
        };

        CompareHaplotypes(observations);

        // If row could not be collapsed into an existing haplotype
        if (miss) {
            observations.emplace_back(
                std::make_shared<Haplotype>(row->Read->Name(), std::move(codons), flag));
        }
    }

    // Generators are haplotypes that have been identified as on target
    std::vector<std::shared_ptr<Haplotype>> generators;
    // Filtered are all other haplotypes
    std::vector<std::shared_ptr<Haplotype>> filtered;
    for (auto& h : observations) {
        // Minimal evidence not reached
        if (h->Size() < 10) h->AddFlag(HaplotypeType::LOW_COV);
        if (h->Flags() == 0)
            generators.emplace_back(std::move(h));
        else
            filtered.emplace_back(std::move(h));
    }

    // Haplotype comparator, ascending
    auto HaplotypeComp = [](const std::shared_ptr<Haplotype>& a,
                            const std::shared_ptr<Haplotype>& b) { return a->Size() < b->Size(); };

    std::sort(generators.begin(), generators.end(), HaplotypeComp);
    std::sort(filtered.begin(), filtered.end(), HaplotypeComp);

    // TODO: matrix solution in J

    if (verbose_) std::cerr << "#Haplotypes: " << generators.size() << std::endl;
    double counts = 0;
    for (auto& hn : generators)
        counts += hn->Size();
    if (verbose_) std::cerr << "#Counts: " << counts << std::endl;

    // Sort generators descending
    std::stable_sort(generators.begin(), generators.end(),
                     [](const std::shared_ptr<Haplotype>& a, const std::shared_ptr<Haplotype>& b) {
                         return a->Size() >= b->Size();
                     });

    static constexpr int alphabetSize = 26;
    bool doubleName = generators.size() > alphabetSize;
    for (size_t genNumber = 0; genNumber < generators.size(); ++genNumber) {
        auto& hn = generators.at(genNumber);

        // Set frequency of this haplotype, among the generators
        hn->Frequency(hn->Size() / counts);

        // Assign each haplotype an unique name
        if (doubleName) {
            hn->Name(std::string(1, 'A' + genNumber / alphabetSize) +
                     std::string(1, 'a' + genNumber % alphabetSize));
        } else {
            hn->Name(std::string(1, 'A' + genNumber));
        }

        // Print
        if (verbose_) std::cerr << (hn->Size() / counts) << "\t" << hn->Size() << "\t";

        // For each variant position, loop through all amino acids and their
        // mutated codons and store which variant this haplotype hit
        size_t numCodons = hn->NumCodons();
        for (size_t i = 0; i < numCodons; ++i) {
            for (auto& kv : variantPositions.at(i).second->aminoAcidToCodons) {
                for (auto& vc : kv.second) {
                    bool hit = hn->Codon(i) == vc.codon;
                    vc.haplotypeHit.push_back(hit);
                    if (hit) {
                        std::cerr << termcolor::red;
                    }
                }
            }
            if (verbose_) std::cerr << hn->Codon(i) << termcolor::reset << " ";
        }
        if (verbose_) std::cerr << std::endl;

        // Store this haplotype
        reconstructedHaplotypes_.push_back(*hn);
    }
    std::cerr << termcolor::reset;

    // From here on only verbose output
    const auto PrintHaplotype = [&variantPositions, this](std::shared_ptr<Haplotype> h) {
        for (const auto& name : h->ReadNames()) {
            std::cerr << name << "\t";
            const auto& row = msaByRow_.NameToRow(name);
            for (const auto& pos_var : variantPositions)
                std::cerr << row->CodonAt(pos_var.first - msaByRow_.BeginPos() - 3) << "\t";
            std::cerr << std::endl;
        }
        std::cerr << std::endl;
    };

    if (verbose_) std::cerr << std::endl << "HAPLOTYPES" << std::endl;
    for (auto& hn : generators) {
        genCounts_ += hn->ReadNames().size();
        if (verbose_) std::cerr << "HAPLOTYPE: " << hn->Name() << std::endl;
        if (verbose_) PrintHaplotype(hn);
    }

    std::map<int, int> filteredCounts;

    if (verbose_) std::cerr << "FILTERED" << std::endl;
    for (auto& h : filtered) {
        filteredCounts[h->Flags()] += h->ReadNames().size();
        if (verbose_) PrintHaplotype(h);
        filteredHaplotypes_.emplace_back(*h);
    }

    int sumFiltered = 0;
    for (const auto& kv : filteredCounts) {
        sumFiltered += kv.second;
        if (kv.first & static_cast<int>(HaplotypeType::WITH_GAP)) margWithGap_ += kv.second;
        if (kv.first & static_cast<int>(HaplotypeType::WITH_HETERODUPLEX))
            margWithHetero_ += kv.second;
        if (kv.first & static_cast<int>(HaplotypeType::PARTIAL)) margPartial_ += kv.second;
        if (kv.first == static_cast<int>(HaplotypeType::LOW_COV)) lowCov_ += kv.second;
        if (kv.first & static_cast<int>(HaplotypeType::OFFTARGET)) margOfftarget_ += kv.second;
    }

    if (verbose_) {
        std::cerr << "HEALTHY, REPORTED\t\t: " << genCounts_ << std::endl;
        std::cerr << "HEALTHY, TOO LOW COVERAGE\t: " << lowCov_ << std::endl;
        std::cerr << "---" << std::endl;
        std::cerr << "ALL DAMAGED\t\t\t: " << margOfftarget_ << std::endl;
        std::cerr << "MARGINAL WITH GAPS\t\t: " << margWithGap_ << std::endl;
        std::cerr << "MARGINAL WITH HETERODUPLEXES\t: " << margWithHetero_ << std::endl;
        std::cerr << "MARGINAL PARTIAL READS\t\t: " << margPartial_ << std::endl;
        std::cerr << "---" << std::endl;
        std::cerr << "SUM\t\t\t: " << genCounts_ + sumFiltered << std::endl;
    }
}

double AminoAcidCaller::Probability(const std::string& a, const std::string& b)
{
    if (a.size() != b.size()) return 0.0;

    double p = 1;
    for (size_t i = 0; i < a.size(); ++i) {
        if (a[i] == '-' || b[i] == '-')
            p *= error_.Deletion;
        else if (a[i] != b[i])
            p *= error_.Substitution;
        else
            p *= error_.Match;
    }
    return p;
};

bool AminoAcidCaller::MeasurePerformance(const TargetGene& tg, const std::string& codon,
                                         const bool& variableSite, const int& aaPos,
                                         const double& p, PerformanceMetrics* pm)
{
    const char aminoacid = AAT::FromCodon.at(codon);
    auto Predictor = [&]() {
        if (!tg.minors.empty()) {
            for (const auto& minor : tg.minors) {
                if (aaPos == minor.position && aminoacid == minor.aminoacid[0] &&
                    codon == minor.codon) {
                    return true;
                }
            }
        }
        return false;
    };
    const bool predictor = Predictor();
    if (variableSite) {
        if (predictor) {
            if (p < alpha)
                ++pm->TruePositives;
            else
                ++pm->FalseNegative;
        } else {
            if (p < alpha)
                ++pm->FalsePositives;
            else
                ++pm->TrueNegative;
        }
    } else if (predictor) {
        if (p < alpha)
            ++pm->TruePositives;
        else
            ++pm->FalseNegative;
    }

    return predictor;
}

MajorityCall AminoAcidCaller::FindMajorityCodon(const std::map<std::string, int>& codons)
{
    MajorityCall mc;
    mc.Coverage = -1;
    for (const auto& codon_counts : codons) {
        if (codon_counts.second > mc.Coverage) {
            mc.Coverage = codon_counts.second;
            mc.Codon = codon_counts.first;
        }
    }
    if (AAT::FromCodon.find(mc.Codon) == AAT::FromCodon.cend()) {
        return MajorityCall();
    }
    mc.AA = AAT::FromCodon.at(mc.Codon);
    return mc;
}

void AminoAcidCaller::CallVariants()
{
    auto genes = targetConfig_.targetGenes;

    // If no user config has been provided, use complete input region
    if (genes.empty()) {
        TargetGene tg(msaByRow_.BeginPos(), msaByRow_.EndPos(), "Unnamed ORF", {});
        genes.emplace_back(tg);
    }

    const int numberOfTests = CountNumberOfTests(genes);
    PerformanceMetrics pm(numberOfTests, targetConfig_.NumExpectedMinors());

    const bool hasExpectedMinors = pm.NumExpectedMinors > 0;
    const bool hasReference = !targetConfig_.referenceSequence.empty();

    for (const auto& gene : genes) {
        VariantGene curVariantGene(gene.name, gene.begin);

        // For each codon in the gene
        for (int i = gene.begin; i < gene.end - 2; ++i) {
            // Absolute reference position
            const int absPos = i - 1;
            // Relative to gene begin
            const int relPos = i - curVariantGene.geneOffset;
            // Only work on beginnings of a codon
            if (relPos % 3 != 0) continue;
            // Relative to window begin
            const int winPos = i - msaByRow_.BeginPos();
            // Relative amino acid position
            const int aaPos = 1 + relPos / 3;

            // Each position is stored in the variant gene
            curVariantGene.relPositionToVariant.emplace(
                aaPos, std::make_shared<VariantGene::VariantPosition>());
            auto& curVariantPosition = curVariantGene.relPositionToVariant.at(aaPos);

            // Gather all observed codons and count actual coverage
            std::map<std::string, int> codons = msaByRow_.CodonsAt(winPos);
            int coverage = 0;
            for (const auto& codon_size : codons)
                coverage += codon_size.second;

            // Get the majority codon of the sample
            MajorityCall mc = FindMajorityCodon(codons);

            // In case a reference has been provided
            if (hasReference) {
                // Get the reference codon
                curVariantPosition->refCodon = targetConfig_.referenceSequence.substr(absPos, 3);
                if (AAT::FromCodon.find(curVariantPosition->refCodon) == AAT::FromCodon.cend()) {
                    continue;
                }
                // And the corresponding amino acid
                curVariantPosition->refAminoAcid = AAT::FromCodon.at(curVariantPosition->refCodon);

                // best alternative to the reference
                if (mc.Coverage == 0) continue;
                if (mc.Coverage * 100.0 / coverage > maximalPerc_) {
                    curVariantPosition->altRefCodon = mc.Codon;
                    curVariantPosition->altRefAminoAcid = mc.AA;
                }
            } else {  // In case no reference has been provided
                if (mc.Coverage == 0) continue;
                curVariantPosition->refCodon = mc.Codon;
                curVariantPosition->refAminoAcid = mc.AA;
            }

            for (const auto& codon_counts : codons) {
                // Skip if the codon of interest is the reference codon
                if (curVariantPosition->refCodon == codon_counts.first) continue;
                // Skip if an alternative reference codon is available and it
                // equals the codon of interest
                if (!curVariantPosition->altRefCodon.empty() &&
                    curVariantPosition->altRefCodon == codon_counts.first)
                    continue;

                // Compute expected counts for null hypothesis that the codon
                // of interest has been generated by the reference via
                // sequencing errors.
                auto expected =
                    coverage * Probability(curVariantPosition->refCodon, codon_counts.first);

                // Compute Fisher's Exact test
                double p =
                    (Statistics::Fisher::fisher_exact_tiss(
                         std::ceil(codon_counts.second), std::ceil(coverage - codon_counts.second),
                         std::ceil(expected), std::ceil(coverage - expected)) *
                     numberOfTests);

                // Handle possible overflows
                if (p > 1) p = 1;

                // Check if there is variability
                const double frequency = 1.0 * codon_counts.second / coverage;
                bool variableSite = frequency < 0.8;
                // Check if this site is a predictor for known minor variants,
                // annotated in the TargetConfig.
                bool predictorSite =
                    MeasurePerformance(gene, codon_counts.first, variableSite, aaPos, p, &pm);

                // Helper to store an actual variant at the current variant position
                auto StoreVariant = [&](const std::string& drmString = "") {
                    // Store if minimal percentage is reached or in debug mode
                    if (debug_ || frequency * 100 >= minimalPerc_) {
                        const char curAA = AAT::FromCodon.at(codon_counts.first);
                        VariantGene::VariantPosition::VariantCodon curVariantCodon;
                        curVariantCodon.codon = codon_counts.first;
                        curVariantCodon.frequency = frequency;
                        curVariantCodon.pValue = p;
                        if (!drmString.empty())
                            curVariantCodon.knownDRM = drmString;
                        else
                            curVariantCodon.knownDRM =
                                FindDRMs(gene.name, genes,
                                         DMutation(curVariantPosition->refAminoAcid, aaPos, curAA));

                        curVariantPosition->aminoAcidToCodons[curAA].push_back(curVariantCodon);
                    }
                };

                // In debug mode, store every codon candidate
                if (debug_) {
                    StoreVariant();
                } else if (p < alpha) {  // otherwise, if smaller than the p-value threshold
                    // In DRM-only mode, only store it, if we found DRMs at this position
                    if (drmOnly_) {
                        const std::string drmString = FindDRMs(
                            gene.name, genes, DMutation(curVariantPosition->refAminoAcid, aaPos,
                                                        AAT::FromCodon.at(codon_counts.first)));
                        if (!drmString.empty()) StoreVariant();
                    } else {  // If we are not in DRM-only mode
                        // In case this is a predictor site of a known variant
                        if (predictorSite) StoreVariant();
                        // If minors are expected and this is a variable site,
                        // not a major codon
                        else if (hasExpectedMinors && variableSite)
                            StoreVariant();
                        // If minors are not expected, normal mode
                        else if (!hasExpectedMinors)
                            StoreVariant();
                    }
                }
            }

            // Fill in the MSA counts of the surrounding positions.
            // Instead of using a special data structure, go directly to JSON
            if (!curVariantPosition->aminoAcidToCodons.empty()) {
                curVariantPosition->coverage = coverage;
                for (int j = -3; j < 6; ++j) {
                    if (i + j >= msaByRow_.BeginPos() && i + j < msaByRow_.EndPos()) {
                        int abs = absPos + j;
                        JSON::Json msaCounts;
                        msaCounts["rel_pos"] = j;
                        msaCounts["abs_pos"] = abs;
                        msaCounts["A"] = msaByColumn_[abs]['A'];
                        msaCounts["C"] = msaByColumn_[abs]['C'];
                        msaCounts["G"] = msaByColumn_[abs]['G'];
                        msaCounts["T"] = msaByColumn_[abs]['T'];
                        msaCounts["-"] = msaByColumn_[abs]['-'];
                        msaCounts["N"] = msaByColumn_[abs]['N'];
                        if (hasReference)
                            msaCounts["wt"] =
                                std::string(1, targetConfig_.referenceSequence.at(abs));
                        else
                            msaCounts["wt"] = std::string(1, msaByColumn_[abs].MaxBase());
                        curVariantPosition->msa.push_back(msaCounts);
                    }
                }
            }
        }
        // Store the gene
        variantGenes_.emplace_back(std::move(curVariantGene));
    }
    // If minors are expected generate performance metrics
    if (hasExpectedMinors) {
        std::ofstream outValJson("validation.json");
        outValJson << pm.ToJson();
        std::cerr << pm << std::endl;
    }
}

std::string PerformanceMetrics::ToJson() const
{
    std::stringstream ss;
    ss << "{\"true_positive_rate\":" << TruePositiveRate() << ",";
    ss << "\"false_positive_rate\":" << FalsePositiveRate() << ",";
    ss << "\"num_tests\":" << NumberOfTests << ",";
    ss << "\"num_false_positives\":" << FalsePositives << ",";
    ss << "\"accuracy\":" << Accuracy() << "}";
    return ss.str();
}

JSON::Json AminoAcidCaller::JSON()
{
    using JSON::Json;
    Json root;
    std::vector<Json> genes;
    for (const auto& v : variantGenes_) {
        Json j = v.ToJson();
        if (j.find("variant_positions") != j.cend()) genes.push_back(j);
    }
    root["genes"] = genes;
    auto HapsToJson = [](const std::vector<Haplotype>& haps) {
        std::vector<Json> haplotypes;
        for (const auto& h : haps) {
            haplotypes.push_back(h.ToJson());
        }
        return haplotypes;
    };
    root["haplotypes"] = HapsToJson(reconstructedHaplotypes_);
    // root["haplotypes_low_counts"] = HapsToJson(lowCountHaplotypes_);
    // root["haplotypes_filtered"] = HapsToJson(filteredHaplotypes_);
    Json counts;
    counts["healthy_reported"] = genCounts_;
    counts["healthy_low_coverage"] = lowCov_;
    counts["all_damaged"] = margOfftarget_;
    counts["marginal_with_gaps"] = margWithGap_;
    counts["marginal_with_heteroduplexes"] = margWithHetero_;
    counts["marginal_partial_reads"] = margPartial_;
    root["haplotype_read_counts"] = counts;
    return root;
}
}
}  // ::PacBio::Juliet
