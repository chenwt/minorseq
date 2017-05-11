// Copyright (c) 2014-2017, Pacific Biosciences of California, Inc.
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
const PlainOption Phasing{
    "mode_phasing",
    { "mode-phasing", "p" },
    "Phase Variants",
    "Phase variants and cluster haplotypes.",
    CLI::Option::BoolType()
};
const PlainOption Error{
    "mode_error",
    { "mode-error" },
    "Alignment Error Rates",
    "Compute alignment error rates.",
    CLI::Option::BoolType(),
    JSON::Json(nullptr),
    CLI::OptionFlags::HIDE_FROM_HELP
};
const PlainOption SubstitutionRate{
    "substitution_rate",
    { "sub", "s" },
    "Substitution Rate",
    "Substitution Rate, specify to override the learned rate.",
    CLI::Option::FloatType(0)
};
const PlainOption DeletionRate{
    "deletion_rate",
    { "del", "d" },
    "Deletion Rate",
    "Deletion Rate, specify to override the learned rate.",
    CLI::Option::FloatType(0)
};
const PlainOption MinimalPerc{
    "minimal_percentage",
    { "min-perc", "m" },
    "Minimal Variant Percentage.",
    "Minimal variant percentage to report.",
    CLI::Option::FloatType(0)
};
const PlainOption TargetConfigTC{
    "target_config",
    { "target-config-tc" },
    "Target",
    "Predefined target config tag, one of \"none\" or \"HIV_HXB2\".",
    CLI::Option::StringType("none"),
    {"none", "HIV_HXB2"},
    CLI::OptionFlags::HIDE_FROM_HELP
};
const PlainOption TargetConfigCLI{
    "target_config_universal",
    { "config", "c" },
    "Target Config",
    "Path to the target config JSON file, predefined target config tag, or the JSON string.",
    CLI::Option::StringType(""),
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
    "Merge outlier haplotypes.",
    CLI::Option::BoolType(),
    JSON::Json(nullptr),
    CLI::OptionFlags::HIDE_FROM_HELP
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
    , Mode(AnalysisModeFromOptions(options))
    , SubstitutionRate(options[OptionNames::SubstitutionRate])
    , DeletionRate(options[OptionNames::DeletionRate])
    , MinimalPerc(options[OptionNames::MinimalPerc])
{
    const std::string targetConfigTC = options[OptionNames::TargetConfigTC];
    const std::string targetConfigCLI = options[OptionNames::TargetConfigCLI];

    if (targetConfigTC != "none")
        TargetConfigUser = targetConfigTC;
    else
        TargetConfigUser = targetConfigCLI;

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

AnalysisMode JulietSettings::AnalysisModeFromOptions(const PacBio::CLI::Results& options)
{
    bool phasing = options[OptionNames::Phasing];
    bool error = options[OptionNames::Error];
    int counter = phasing + error;
    if (counter > 1) throw std::runtime_error("Overriding mode is mutually exclusive!");

    if (!phasing && !error)
        return AnalysisMode::AMINO;
    else if (phasing)
        return AnalysisMode::PHASING;
    else if (error)
        return AnalysisMode::ERROR;
    else
        throw std::runtime_error("Cannot execute mode, undefined behaviour!");
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
        OptionNames::Verbose,
        OptionNames::Debug,
        OptionNames::MergeOutliers,
        OptionNames::TargetConfigTC,
        OptionNames::Error
    });

    i.AddGroup("Configuration",
    {
        OptionNames::TargetConfigCLI,
        OptionNames::Phasing
    });

    i.AddGroup("Restrictions",
    {
        OptionNames::Region,
        OptionNames::DRMOnly,
        OptionNames::MinimalPerc
    });

    i.AddGroup("Chemistry override (specify both)",
    {
        OptionNames::SubstitutionRate,
        OptionNames::DeletionRate
    });

    const std::string id = "minorseq.tasks.juliet";
    Task tcTask(id);
    tcTask.AddOption(OptionNames::Phasing);
    tcTask.AddOption(OptionNames::Region);
    tcTask.AddOption(OptionNames::DRMOnly);
    tcTask.AddOption(OptionNames::TargetConfigTC);
    tcTask.AddOption(OptionNames::TargetConfigCLI);
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
