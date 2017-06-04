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
#include <exception>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <numeric>
#include <vector>

#include <pbbam/BamReader.h>
#include <pbbam/BamRecord.h>
#include <pbbam/DataSet.h>

#include <pbcopper/json/JSON.h>
#include <pbcopper/utility/FileUtils.h>

#include <pacbio/data/ArrayRead.h>
#include <pacbio/data/MSA.h>
#include <pacbio/io/BamParser.h>
#include <pacbio/juliet/AminoAcidCaller.h>
#include <pacbio/juliet/JsonToHtml.h>
#include <pacbio/juliet/JulietSettings.h>
#include <pacbio/statistics/Fisher.h>
#include <pacbio/statistics/Tests.h>

#include <pacbio/juliet/JulietWorkflow.h>

namespace PacBio {
namespace Juliet {

std::ostream& JulietWorkflow::LogCI(const std::string& prefix)
{
    std::cout << std::setw(20) << std::left << prefix << ": ";
    return std::cout;
}

void JulietWorkflow::Run(const JulietSettings& settings)
{
    if (settings.Mode == AnalysisMode::AMINO || settings.Mode == AnalysisMode::PHASING) {
        AminoPhasing(settings);
    } else if (settings.Mode == AnalysisMode::ERROR) {
        Error(settings);
    }
}

void JulietWorkflow::AminoPhasing(const JulietSettings& settings)
{
    using BAM::DataSet;

    std::string outputHtml;
    std::string outputJson;
    std::string outputMsa;
    std::string bamInput;
    for (const auto& i : settings.InputFiles) {
        const auto fileExt = PacBio::Utility::FileExtension(i);
        if (fileExt == "json") {
            if (!outputJson.empty()) throw std::runtime_error("Only one json output file allowed");
            outputJson = i;
            continue;
        }
        if (fileExt == "html") {
            if (!outputHtml.empty()) throw std::runtime_error("Only one html output file allowed");
            outputHtml = i;
            continue;
        }
        if (fileExt == "msa") {
            if (!outputMsa.empty()) throw std::runtime_error("Only one msa output file allowed");
            outputMsa = i;
            continue;
        }
        DataSet ds(i);
        switch (ds.Type()) {
            case DataSet::TypeEnum::SUBREAD:
            case DataSet::TypeEnum::ALIGNMENT:
            case DataSet::TypeEnum::CONSENSUS_ALIGNMENT:
                bamInput = i;
                break;
            default:
                throw std::runtime_error("Unsupported input file: " + i + " of type " +
                                         DataSet::TypeToName(ds.Type()));
        }
    }

    if (bamInput.empty()) throw std::runtime_error("Missing input file!");
    if (outputHtml.empty() && outputJson.empty() && outputMsa.empty()) {
        const auto prefix = PacBio::Utility::FilePrefix(bamInput);
        outputHtml = prefix + ".html";
        outputJson = prefix + ".json";
    }

    auto sharedReads = IO::BamToArrayReads(bamInput, settings.RegionStart, settings.RegionEnd);

    if (sharedReads.empty()) {
        std::cerr << "Empty input." << std::endl;
        exit(1);
    }

    std::string chemistry = sharedReads.front()->SequencingChemistry();
    for (size_t i = 1; i < sharedReads.size(); ++i)
        if (chemistry != sharedReads.at(i)->SequencingChemistry())
            throw std::runtime_error("Mixed chemistries are not allowed");

    ErrorEstimates error;
    if (settings.SubstitutionRate != 0.0 && settings.DeletionRate != 0.0) {
        error = ErrorEstimates(settings.SubstitutionRate, settings.DeletionRate);
    } else {
        error = ErrorEstimates(chemistry);
    }

    // Call variants
    AminoAcidCaller aac(sharedReads, error, settings);
    if (settings.Mode == AnalysisMode::PHASING) aac.PhaseVariants();

    const auto json = aac.JSON();

    if (!outputJson.empty()) {
        std::ofstream jsonStream(outputJson);
        jsonStream << json.dump(2) << std::endl;
    }

    if (!outputHtml.empty()) {
        std::ofstream htmlStream(outputHtml);
        JsonToHtml::HTML(htmlStream, json, settings.TargetConfigUser, settings.DRMOnly, bamInput,
                         settings.CLI);
    }

    // Store msa + p-values
    if (!outputMsa.empty()) {
        std::ofstream msaStream(outputMsa);
        msaStream << "pos A C G T - N" << std::endl;
        int pos = aac.msaByColumn_.BeginPos();
        for (auto& column : aac.msaByColumn_) {
            ++pos;
            msaStream << pos;
            for (const auto& c : {'A', 'C', 'G', 'T', '-', 'N'})
                msaStream << " " << column[c];
            msaStream << std::endl;
        }
        msaStream.close();
    }
}
void JulietWorkflow::Error(const JulietSettings& settings)
{
    for (const auto& inputFile : settings.InputFiles) {
        auto reads = IO::BamToArrayReads(inputFile, settings.RegionStart, settings.RegionEnd);
        Data::MSAByColumn msa(reads);
        double sub = 0;
        double del = 0;
        int columnCount = 0;
        for (const auto& column : msa) {
            if (column.Coverage() > 100) {
                del += column.Frequency('-');
                sub += 1.0 - column.Frequency('-') - column.Frequency(column.MaxBase());
                ++columnCount;
            }
        }
        std::cout << inputFile << std::endl;
        std::cout << "sub: " << (sub / columnCount) << std::endl;
        std::cout << "del: " << (del / columnCount) << std::endl;
    }
}
}
}  // ::PacBio::Juliet