// Copyright (c) 2014-2015, Pacific Biosciences of California, Inc.
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

#include <thread>

#include <pacbio/Version.h>
#include <pacbio/data/PlainOption.h>
#include <boost/algorithm/string.hpp>

#include <pacbio/juliet/JulietSettings.h>
#include <pacbio/juliet/TargetConfig.h>

namespace PacBio {
namespace Juliet {
namespace OptionNames {
using PlainOption = Data::PlainOption;
// clang-format off
const PlainOption Region{
    "region",
    { "region", "r"},
    "Region of Interest",
    "Genomic region of interest, reads will be clipped to that region, empty means all reads.",
    CLI::Option::StringType("")
};
const PlainOption DRMOnly{
    "only_known_drms",
    { "drm-only", "k" },
    "Only Known DRMs",
    "Only report known DRM positions.",
    CLI::Option::BoolType()
};
const PlainOption Mode{
    "execution_mode",
    { "mode", "m" },
    "Execution mode",
    "Execution mode: amino, phasing, or error",
    CLI::Option::StringType("amino")
};
const PlainOption SubstitutionRate{
    "substitution_rate",
    { "sub", "s" },
    "Substitution Rate",
    "Substitution Rate, specify to override the learned rate",
    CLI::Option::FloatType(0)
};
const PlainOption DeletionRate{
    "deletion_rate",
    { "del", "d" },
    "DeletionRate",
    "Deletion Rate, specify to override the learned rate",
    CLI::Option::FloatType(0)
};
const PlainOption TargetConfig{
    "target_config",
    { "config", "c" },
    "Target config",
    "Predefined target config tag: [none, HIV, ABL1]",
    CLI::Option::StringType("none"),
    {"none", "HIV", "HIV-PB", "ABL1"}
};
const PlainOption TargetConfigFile{
    "target_config_file",
    { "config-file"},
    "Target config file",
    "Path to the target config JSON file",
    CLI::Option::StringType("")
};
const PlainOption TargetConfigJSON{
    "target_config_json",
    { "config-json"},
    "Target config",
    "Target config JSON string",
    CLI::Option::StringType("")
};
const PlainOption Verbose{
    "verbose",
    { "verbose" },
    "Verbose",
    "Verbose",
    CLI::Option::BoolType()
};
const PlainOption MergeOutliers{
    "merge_outliers",
    { "merge-outliers" },
    "Merge Outliers",
    "Merge outliers haplotypes.",
    CLI::Option::BoolType()
};
const PlainOption Debug{
    "debug",
    { "debug" },
    "Debug",
    "Debug returns all amino acids, irrelevant of their significance.",
    CLI::Option::BoolType()
};
// clang-format on
}  // namespace OptionNames

JulietSettings::JulietSettings(const PacBio::CLI::Results& options)
    : CLI(options.InputCommandLine())
    , InputFiles(options.PositionalArguments())
    , DRMOnly(options[OptionNames::DRMOnly])
    , MergeOutliers(options[OptionNames::MergeOutliers])
    , Verbose(options[OptionNames::Verbose])
    , Debug(options[OptionNames::Debug])
    , Mode(AnalysisModeFromString(options[OptionNames::Mode]))
    , SubstitutionRate(options[OptionNames::SubstitutionRate])
    , DeletionRate(options[OptionNames::DeletionRate])
{
    std::string targetConfigTag = options[OptionNames::TargetConfig];
    std::string targetConfigJSON = options[OptionNames::TargetConfigJSON];
    std::string targetConfigFile = options[OptionNames::TargetConfigFile];

    int targetCounter =
        !targetConfigTag.empty() + !targetConfigJSON.empty() + !targetConfigFile.empty();

    if (targetCounter > 1) {
        throw std::runtime_error("Target config options are mutually exclusive!");
    }
    if (targetCounter == 0 || targetConfigTag == "none") {
        TargetConfigUser = TargetConfig("{}");
    } else if (!targetConfigTag.empty()) {
        if (targetConfigTag.size() >= 3) {
            if (targetConfigTag.at(0) != '<') targetConfigTag = '<' + targetConfigTag;
            if (targetConfigTag.at(targetConfigTag.size() - 1) != '>') targetConfigTag += '>';
        }
        TargetConfigUser = targetConfigTag;
    } else if (!targetConfigJSON.empty()) {
        TargetConfigUser = targetConfigJSON;
    } else if (!targetConfigFile.empty()) {
        TargetConfigUser = targetConfigFile;
    }

    SplitRegion(options[OptionNames::Region], &RegionStart, &RegionEnd);
}

size_t JulietSettings::ThreadCount(int n)
{
    const int m = std::thread::hardware_concurrency();

    if (n < 1) return std::max(1, m + n);

    return std::min(m, n);
}

void JulietSettings::SplitRegion(const std::string& region, int* start, int* end)
{
    if (region.compare("") != 0) {
        std::vector<std::string> splitVec;
        boost::split(splitVec, region, boost::is_any_of("-"));
        *start = stoi(splitVec[0]);
        *end = stoi(splitVec[1]);
        if (*start <= 0 || *end <= 0) throw std::runtime_error("Indexing is 1-based");
    }
}

AnalysisMode JulietSettings::AnalysisModeFromString(const std::string& input)
{
    std::string s = input;
    std::transform(s.begin(), s.end(), s.begin(), ::tolower);
    if (s.find("amino") != std::string::npos || s.find("acid") != std::string::npos)
        return AnalysisMode::AMINO;
    else if (s.find("base") != std::string::npos || s.find("nuc") != std::string::npos)
        return AnalysisMode::BASE;
    else if (s.find("phas") != std::string::npos || s.find("hap") != std::string::npos)
        return AnalysisMode::PHASING;
    else if (s.find("error") != std::string::npos)
        return AnalysisMode::ERROR;
    else
        throw std::runtime_error("Unknown mode " + s);
}

PacBio::CLI::Interface JulietSettings::CreateCLI()
{
    using Option = PacBio::CLI::Option;
    using Task = PacBio::CLI::ToolContract::Task;

    PacBio::CLI::Interface i{
        "juliet",
        "Juliet, minimal minor variant calling "
        "software.\nAttention: Juliet is for research usage "
        "only. Predictions have not been validated.",
        PacBio::MinorseqVersion() + " (commit " + PacBio::MinorseqGitSha1() + ")"};

    i.AddHelpOption();     // use built-in help output
    i.AddVersionOption();  // use built-in version output

    // clang-format off
    i.AddPositionalArguments({
        {"source", "Source BAM or DataSet XML file.", "FILE"}
    });

    i.AddOptions(
    {
        OptionNames::Mode,
        OptionNames::Region,
        OptionNames::DRMOnly,
        OptionNames::MergeOutliers,
        OptionNames::SubstitutionRate,
        OptionNames::DeletionRate,
        OptionNames::Verbose,
        OptionNames::Debug
    });

    i.AddGroup("Target Configuration (mutually exclusive)",
    {
        OptionNames::TargetConfig,
        OptionNames::TargetConfigJSON,
        OptionNames::TargetConfigFile
    });

    const std::string id = "minorseq.tasks.juliet";
    Task tcTask(id);
    tcTask.AddOption(OptionNames::Mode);
    tcTask.AddOption(OptionNames::Region);
    tcTask.AddOption(OptionNames::DRMOnly);
    tcTask.AddOption(OptionNames::TargetConfig);
    tcTask.AddOption(OptionNames::TargetConfigJSON);
    tcTask.AddOption(OptionNames::TargetConfigFile);
    tcTask.AddOption(OptionNames::MergeOutliers);
    tcTask.AddOption(OptionNames::SubstitutionRate);
    tcTask.AddOption(OptionNames::DeletionRate);
    tcTask.AddOption(OptionNames::Debug);

        tcTask.InputFileTypes({
        {
            "alignment_set",
            "AlignmentSet",
            "Consensus (CCS) Alignment DataSet or aligned .bam file",
            "PacBio.DataSet.ConsensusAlignmentSet"
        }
    });

    tcTask.OutputFileTypes({
        {
            "html_report",
            "HTML Report",
            "Human-readable HTML report generated by juliet",
            "PacBio.FileTypes.html",
            "juliet_report"
        },
        {
            "json_report",
            "JSON Report",
            "JSON report generated by juliet",
            "PacBio.FileTypes.json",
            "juliet_report"
        }
    });

    CLI::ToolContract::Config tcConfig(tcTask);
    i.EnableToolContract(tcConfig);
    // clang-format on

    return i;
}
}
}  // ::PacBio::CCS
