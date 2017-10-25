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

#include <map>
#include <regex>
#include <string>
#include <unordered_map>

#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/replace.hpp>

#include <pbbam/DataSet.h>

#include <pacbio/Version.h>
#include <pacbio/juliet/JsonToHtml.h>

namespace PacBio {
namespace Juliet {

void JsonToHtml::DRMView(std::ostream& out, const JSON::Json& j, const TargetConfig& config,
                         bool onlyKnownDRMs)
{
    std::map<std::string, std::string> geneStart;
    for (const auto& targetGene : config.targetGenes) {
        geneStart[targetGene.name] = std::to_string(targetGene.begin);
    }

    struct VariantDRM
    {
        std::string refAA;
        std::string refPos;
        std::string curAA;
        double frequency;
    };
    std::map<std::string, std::map<std::string, std::vector<VariantDRM>>> drmsWithVariants;
    for (const auto& gene : j["genes"]) {
        for (auto& variantPosition : gene["variant_positions"]) {
            for (auto& variant_amino_acid : variantPosition["variant_amino_acids"]) {
                for (auto& variant_codon : variant_amino_acid["variant_codons"]) {
                    std::string knownDRM = variant_codon["known_drm"];
                    if (knownDRM.empty()) continue;
                    std::vector<std::string> singleDrms;
                    boost::algorithm::split(singleDrms, knownDRM, boost::is_any_of("+"));
                    for (auto& d : singleDrms) {
                        d = std::regex_replace(d, std::regex("^ +| +$|( ) +"), "$1");
                        VariantDRM v;
                        v.refAA = variantPosition["ref_amino_acid"];
                        v.refPos =
                            std::to_string(static_cast<int>(variantPosition["ref_position"]));
                        v.curAA += variant_amino_acid["amino_acid"].get_ref<const std::string&>();
                        v.frequency = variant_codon["frequency"];
                        std::string key = geneStart[gene["name"]] + "|";
                        key += gene["name"].get_ref<const std::string&>();
                        drmsWithVariants[d][key].emplace_back(v);
                    }
                }
            }
        }
    }
    if (drmsWithVariants.empty()) {
        out << "No known drug-resistance mutations present." << std::endl;
        return;
    }
    size_t geneWidth = 0;
    for (const auto& targetGene : config.targetGenes) {
        geneWidth = std::max(targetGene.name.size(), geneWidth);
    }
    geneWidth *= 8;
    size_t drugWidth = 0;
    for (const auto& a : drmsWithVariants) {
        drugWidth = std::max(a.first.size(), drugWidth);
    }
    drugWidth *= 10;
    out << "<table class=\"drmview\">" << std::endl;
    // clang-format off
    out << "<tr><th colspan=2></th><th colspan=2 style=\"border-right: 1px dashed black\">Reference</th><th colspan=2>Sample</th></tr>";
    out << "<tr><th>Drug</th><th>Gene</th><th>AA</th><th style=\"border-right: 1px dashed black\">Pos</th><th>AA</th><th>%</th></tr>";
    // clang-format on
    for (const auto& a : drmsWithVariants) {
        int numRow = 0;
        for (const auto& b : a.second) {
            for (const auto& c : b.second) {
                ++numRow;
            }
        }
        out << "<tr><td rowspan=\"" << numRow << "\" class=\"gene\" style=\"width:" << drugWidth
            << "px\">" << a.first << "</td>" << std::endl;
        bool firstInDrug = true;
        for (const auto& b : a.second) {

            std::string drugSuffix;
            if (!firstInDrug)
                out << "<tr>";
            else
                drugSuffix = "First";

            std::string geneName;
            const auto idx = b.first.find_first_of('|');
            if (std::string::npos != idx) geneName = b.first.substr(idx + 1);

            out << "<td rowspan=\"" << b.second.size() << "\" class=\"drug" << drugSuffix
                << "\" style=\"width:" << geneWidth << "px\">" << geneName << "</td>" << std::endl;
            bool firstInGene = true;
            for (const auto& c : b.second) {
                std::string classSuffix;
                if (!firstInGene) out << "<tr>";
                if (firstInDrug)
                    classSuffix = "FirstDrug";
                else if (firstInGene)
                    classSuffix = "FirstGene";
                out << "<td class=\"refaa" << classSuffix << "\">" << c.refAA << "</td>";
                out << "<td class=\"refpos" << classSuffix << "\">" << c.refPos << "</td>";
                out << "<td class=\"curaa" << classSuffix << "\">" << c.curAA << "</td>";
                double fOrig = c.frequency;
                double fTmp;
                int exp = 0;
                do {
                    fTmp = fOrig * std::pow(10, ++exp);
                } while (static_cast<int>(fTmp) < 10);
                fOrig = static_cast<int>(fOrig * std::pow(10, exp));
                fOrig /= std::pow(10, exp - 2);
                out << "<td class=\"freq" << classSuffix << "\">" << fOrig << "</td>";
                out << "</tr>";
                firstInDrug = false;
                firstInGene = false;
            }
            firstInDrug = false;
        }
    }
    out << "</table>" << std::endl;
}

void JsonToHtml::Escape(std::string& data)
{
    std::string buffer;
    buffer.reserve(data.size());
    for (size_t pos = 0; pos != data.size(); ++pos) {
        switch (data[pos]) {
            case '&':
                buffer.append("&amp;");
                break;
            case '\"':
                buffer.append("&quot;");
                break;
            case '\'':
                buffer.append("&apos;");
                break;
            case '<':
                buffer.append("&lt;");
                break;
            case '>':
                buffer.append("&gt;");
                break;
            default:
                buffer.append(&data[pos], 1);
                break;
        }
    }
    data.swap(buffer);
}

std::string JsonToHtml::Strip(const std::string& input)
{
    std::string s = input;
    s.erase(std::remove(s.begin(), s.end(), '\"'), s.end());
    return s;
};

void JsonToHtml::HTML(std::ostream& out, const JSON::Json& j, const TargetConfig& config,
                      bool onlyKnownDRMs, std::string filename, std::string parameters)
{
    Escape(filename);
    Escape(parameters);
    // Count number of haplotypes
    auto CountNumHaplotypes = [&j]() {
        int i = -1;
        for (const auto& gene : j["genes"]) {
            for (auto& variantPosition : gene["variant_positions"]) {
                for (auto& variant_amino_acid : variantPosition["variant_amino_acids"]) {
                    for (auto& variant_codons : variant_amino_acid["variant_codons"]) {
                        const int tmp = variant_codons["haplotype_hit"].size();
                        if (i == -1)
                            i = tmp;
                        else if (i != tmp)
                            throw std::runtime_error("Different number of haplotypes.");
                    }
                }
            }
        }
        return i;
    };

    int numHaplotypes = CountNumHaplotypes();

    out << "<!-- Juliet Minor Variant Summary by Dr. Armin Toepfer (Pacific Biosciences) -->"
        << std::endl
        << "<html>" << std::endl
        << "<head>" << std::endl
        << R"(
            <script src="http://ajax.googleapis.com/ajax/libs/jquery/1.11.1/jquery.min.js"></script>
            <script type="text/javascript">
            $(document).ready(function() {
                $(".var").bind( "click", function( event ) {
                    $(this).next().slideToggle(0);
            });
            });
            </script>)"
        << std::endl
        << "<style>" << std::endl
        << R"(
        *,
        *:before,
        *:after {
            -moz-box-sizing: border-box;
            -webkit-box-sizing: border-box;
            box-sizing: border-box;
        }

        html {
            font-family: Helvetica, Arial, sans-serif;
            font-size: 100%;
            background: #fff;
            -webkit-font-smoothing: antialiased;
        }

        details {
            border-radius: 5px;
            border-left: 2px solid black;
        }

        summary {
            border-radius: 3px;
            padding: 5px 10px;
            outline: none;
            font-weight: bold;
        }

        /* Tooltip container */
        .tooltip {
            position: relative;
            display: inline-block;
        }

        /* Tooltip text */
        .tooltip .tooltiptext, .tooltip .tooltiptextlarge {
            visibility: hidden;
            width: 50px;
            border: 1px dotted #2d2d2d;
            color: black;
            text-align: center;
            padding: 5px 0;
            border-radius: 6px;
            background-color: white;

            bottom: 100%;
            left: 50%;
            margin-left: -25px;

            /* Position the tooltip text - see examples below! */
            position: absolute;
            z-index: 1;
        }

        /* Show the tooltip text when you mouse over the tooltip container */
        .tooltip:hover .tooltiptext, .tooltip:hover .tooltiptextlarge {
            visibility: visible;
        }

        .tooltip .tooltiptextlarge {
            width: 380px;
            margin-left: -190px;
            align:center;
        }

        table.hapcounts {
            margin-left: auto;
            margin-right: auto;
        }

        table.discovery {
            border-collapse: collapse;
            margin-bottom: 5px;
        }

        table.discovery>tbody>tr:nth-child(1),
        table.msacounts>tbody>tr:nth-child(1) {
            background-color: #3d3d3d;
            color: white;
        }

        table.discovery tr:nth-child(3):not(.msa) th {
            padding: 5px 5px 5px 5px;
            text-align: center;
            border-bottom: 1px solid #2d2d2d;
        }

        table.discovery tr:nth-child(2):not(.msa) th:nth-child(2) {
            border-left: 1px dashed black;
        }

        table.discovery tr:nth-child(3):not(.msa) th:nth-child(3) {
            border-right: 1px dashed black;
        })";
    if (numHaplotypes)
        out << R"(
        table.discovery tr:nth-child(2):not(.msa) th:nth-child(2) {
            border-right: 1px dashed black;
        })";

    out << R"(
        table.discovery tr:nth-child(3):not(.msa) th:nth-child(9) {
            border-left: 1px dashed black;
        }

        table.discovery tr.var td {
            padding: 15px 5px 15px 5px;
            text-align: center;
            border-bottom: 1px solid white;
        }

        table.discovery tr.var td:nth-child(1) {
            background-color: #ddd;
            border-right: 1px solid #eee;
        }

        table.discovery tr.var td:nth-child(2) {
            background-color: #eee;
            border-right: 1px solid #ddd;
        }

        table.discovery tr.var td:nth-child(3) {
            background-color: #fff;
            border-right: 1px solid #ddd;
            font-weight: bold;
        }

        table.discovery tr.var td:nth-child(4) {
            background-color: #eee;
            border-right: 1px dashed #ccc;
        }

        table.discovery tr.var td:nth-child(5) {
            background-color: #ddd;
            border-right: 1px dashed #bbb;
        }

        table.discovery tr.var td:nth-child(6) {
            background-color: #ccc;
            border-right: 1px dashed #aaa;
        }

        table.discovery tr.var td:nth-child(7) {
            background-color: #bbb;
        }

        table.discovery tr.var td:nth-child(8) {
            background-color: #aaa;
            color: white
        }

        table.discovery tr.var td:nth-child(8) {
            border-right: 1px solid white;
        }

        table.discovery tr.var td:nth-child(n+9) {
            border-left: 1px dotted white;
        }

        table.discovery tr.var td:nth-child(n+9) {
            background-color: #4a4a4a;
        }

        table.discovery tr.var:hover td {
            background-color: white;
        }

        table.discovery tr.var:hover td:nth-child(8) {
            color: purple;
        }

        table.discovery tr.msa table tr:hover td {
            background-color: gray;
            color: white;
        }

        table.discovery tr.msa table {
            border-collapse: collapse;
            background-color: white;
            border: 0;
        }

        table.discovery tr.msa table td {
            background-color: white;
            text-align: center;
            padding: 15px 5px 15px 5px;
            border: 0;
            border-bottom: 1px solid gray;
            font-weight: normal
        }

        table.discovery tr.msa table tr {
            border: 0;
        }

        table.discovery tr.msa table th {
            border: 0;
        }

        .msa {
            display: none;
        })";
    out << R"(
        table.drmview {
            border-collapse: collapse;
            margin-bottom: 20px;
            min-width: 200px;
            text-align: center;
            margin-left: 20px;
        }

        table.drmview td {
            padding: 15px;
        }

        table.drmview td.gene {
            border-top: 3px solid white;
            background-color: #b50937;
            color:white;
            vertical-align:top;
            font-weight: bold;
        }

        table.drmview td.drug, table.drmview td.drugFirst {
            border-top: 1px dashed white;
            background-color: #2d2d2d;
            color:white;
            vertical-align:top;
        }

        table.drmview td.drugFirst {
            border-top: 3px solid white;
        }

        table.drmview td.refaa, table.drmview td.refaaFirstDrug, table.drmview td.refaaFirstGene {
            background-color: #bbb;
            border-right: 1px solid #ddd;
        }
        table.drmview td.refaaFirstDrug {
            border-top: 3px solid white;
        }
        table.drmview td.refaaFirstGene {
            border-top: 1px dashed white;
        }

        table.drmview td.refpos, table.drmview td.refposFirstDrug, table.drmview td.refposFirstGene {
            background-color: #ccc;
            border-right: 1px solid #ddd;
        }
        table.drmview td.refposFirstDrug {
            border-top: 3px solid white;
        }
        table.drmview td.refposFirstGene {
            border-top: 1px dashed white;
        }

        table.drmview td.curaa, table.drmview td.curaaFirstDrug, table.drmview td.curaaFirstGene {
            background-color: #ddd;
            border-right: 1px dashed #ccc;
        }
        table.drmview td.curaaFirstDrug {
            border-top: 3px solid white;
        }
        table.drmview td.curaaFirstGene {
            border-top: 1px dashed white;
        }

        table.drmview td.freq, table.drmview td.freqFirstDrug, table.drmview td.freqFirstGene {
            background-color: #eee;
        }
        table.drmview td.freqFirstDrug {
            border-top: 3px solid white;
        }
        table.drmview td.freqFirstGene {
            border-top: 1px dashed white;
        })"
        << std::endl
        << "</style>" << std::endl
        << "</head>" << std::endl
        << R"(<body>
            <h1 style="margin-top:5px">Minor Variants Summary (Juliet)</h1>
            <details style="margin-bottom: 20px">
            <summary>Input data</summary>
            <div style="margin-left:20px; padding-top: 10px">)";
    out << "<table>";
    out << "<tr><td>Timestamp:</td><td><code>" << BAM::ToIso8601(std::chrono::system_clock::now())
        << "</code></td></tr>";
    out << "<tr><td>Input File:</td><td><code>" << filename << "</td></tr>";
    out << R"(<tr><td>Command Line Call:</td><td><code>)";
    if (parameters.empty())
        out << "Invoked from SMRTLink, please check SMRTLink logs for parameters";
    else
        out << parameters;
    out << "</td></tr>"
        << R"(<tr><td>Juliet Version:</td><td><code>)" << PacBio::MinorseqVersion() << " (commit "
        << PacBio::MinorseqGitSha1() << ")"
        << "</td></tr>";
    out << "</table>";
    out << "</div></details>" << std::endl;

    out << R"(
            <details style="margin-bottom: 20px;margin-top:10px">
            <summary>Target config</summary>
            <div style="padding-left:20px;padding-top:10px">)";
    out << "<table>";
    const std::string version = config.version.empty() ? "NA" : config.version;
    const std::string referenceName = config.referenceName.empty() ? "NA" : config.referenceName;
    const std::string referenceSequenceLength =
        config.referenceSequence.empty() ? "NA" : std::to_string(config.referenceSequence.size());
    out << "<tr><td>Config Version:</td><td><code>" << version << "</code></td></tr>";
    out << "<tr><td>Reference Name:</td><td><code>" << referenceName << "</code></td></tr>";
    out << "<tr><td>Reference Length:</td><td><code>" << referenceSequenceLength
        << "</code></td></tr>";
    if (config.targetGenes.empty()) out << "<tr><td>Genes:</td><td><code>NA</code></td></tr>";
    out << "</table>";
    if (!config.targetGenes.empty()) {
        out << "<span style=\"padding-left:3px\">Genes:</span><ul style=\"margin-top:0px\">"
            << std::endl;
        for (const auto& gene : config.targetGenes) {
            out << "<li style=\"margin-top:5px;\">"
                << "<b>" << gene.name << "</b>"
                << " (" << gene.begin << "-" << gene.end << ")";
            if (!gene.drms.empty()) {
                out << "<ul>";
                for (const auto& drm : gene.drms) {
                    out << "<li>"
                        << "<code>" << drm.name << ":";

                    for (const auto& pos : drm.positions) {
                        out << " " << std::string(pos);
                    }
                    out << "</code>"
                        << "</li>";
                }
                out << "</ul>" << std::endl;
            }
            out << "</li>" << std::endl;
        }
        out << "</ul>" << std::endl;
    }
    out << "</div></details>" << std::endl;

    out << R"(<details open style="margin-bottom: 20px">
            <summary>Variant Discovery</summary>
            <div style="margin-left:20px; padding-top:10px">)";
    Discovery(out, j, config, onlyKnownDRMs, numHaplotypes);
    out << "</div></details>" << std::endl;
    out << R"(<details style="margin-bottom: 20px">
            <summary>Drug Summaries</summary>)";
    DRMView(out, j, config, onlyKnownDRMs);
    out << "</details>" << std::endl;
    out << std::endl << "</body>" << std::endl << "</html>" << std::endl;
}

void JsonToHtml::Discovery(std::ostream& out, const JSON::Json& j, const TargetConfig& config,
                           bool onlyKnownDRMs, int numHaplotypes)
{
    bool hasConf = !config.referenceName.empty() && !config.referenceSequence.empty();

    static std::vector<std::string> colors = {"#ea3c1c", "#f48e00", "#ebff0a", "#56e400",
                                              "#51c6ff", "#4a80ff", "#ae37ff", "#db005f"};
    const std::string referenceName = config.referenceName;

    if (j.find("genes") == j.cend() || j["genes"].is_null()) return;
    for (const auto& gene : j["genes"]) {
        out << "<table class=\"discovery\">" << std::endl
            << R"(
                <col width="40px"/>
                <col width="40px"/>
                <col width="40px"/>
                <col width="40px"/>
                <col width="40px"/>
                <col width="60px"/>
                <col width="60px"/>
                <col width="180px"/>)";
        for (int hap = 0; hap < numHaplotypes; ++hap) {
            out << R"(<col width="40"/>)";
        }
        out << R"(<tr>
                <th colspan=")"
            << 8 << R"(">)" << Strip(gene["name"]) << "</th>";
        for (int hap = 0; hap < numHaplotypes; ++hap) {
            out << "<th style=\"color:" << colors.at(hap % colors.size()) << "\">"
                << Strip(j["haplotypes"][hap]["name"]);
            out << "</th>";
        }

        out << R"(</tr><tr>
                <th colspan="3">)";
        if (referenceName.empty()) out << "Majority Call";
        if (referenceName.size() > 11)
            out << referenceName.substr(0, 11) << "...";
        else
            out << referenceName;
        out << R"(</th>
                <th colspan="5">Sample Variants</th>)";
        if (numHaplotypes > 0) {
            out << R"(<th colspan=")" << (numHaplotypes) << R"("><div class="tooltip">)";
            out << "<span class=\"tooltiptextlarge\">";
            out << R"(<table class="hapcounts"><col width="280px" /><col width="60px" />)";
            out << "<tr><td>"
                << "<b>Haplotype Category</b>"
                << "</td><td>"
                << "<b>#Reads</b>"
                << "</td</tr>" << std::endl;
            out << "<tr><td>"
                << "Reported"
                << "</td><td>" << j["haplotype_read_counts"]["healthy_reported"] << "</td</tr>"
                << std::endl;
            out << "<tr><td>"
                << "Insufficient Coverage (unreported)"
                << "</td><td>" << j["haplotype_read_counts"]["healthy_low_coverage"] << "</td</tr>"
                << std::endl;
            out << "<tr><td>"
                << "Overall Damaged (unreported)"
                << "</td><td>" << j["haplotype_read_counts"]["all_damaged"] << "</td</tr>"
                << std::endl;
            out << "<tr><td>"
                << R"(<span style="padding-left:10px">- Marginal Gaps</span>)"
                << "</td><td>" << j["haplotype_read_counts"]["marginal_with_gaps"] << "</td</tr>"
                << std::endl;
            out << "<tr><td>"
                << R"(<span style="padding-left:10px">- Marginal Heteroduplexes</span>)"
                << "</td><td>" << j["haplotype_read_counts"]["marginal_with_heteroduplexes"]
                << "</td</tr>" << std::endl;
            out << "<tr><td>"
                << R"(<span style="padding-left:10px">- Marginal Partial</span>)"
                << "</td><td>" << j["haplotype_read_counts"]["marginal_partial_reads"]
                << "</td</tr>" << std::endl;
            out << "</table>";
            out << "</span>"
                << "Haplotypes %</div></th>";
        }
        out << R"(
                </tr>
                <tr>
                <th>Codon</th>
                <th>AA</th>
                <th>Pos</th>
                <th>AA</th>
                <th>Codon</th>
                <th>%</th>
                <th>Coverage</th>
                <th>Affected Drugs)";
        if (!config.dbVersion.empty()) out << "<sup>*</sup>";
        out << "</th>";
        for (int hap = 0; hap < numHaplotypes; ++hap) {
            out << R"(<th><div class="tooltip">)"
                << std::round(1000 * static_cast<double>(j["haplotypes"][hap]["frequency"])) / 10.0;
            out << "<span class=\"tooltiptext\">" << j["haplotypes"][hap]["reads_hard"]
                << "</span>";
            out << "</div></th>";
        }
        out << R"(</tr>)" << std::endl;

        for (auto& variantPosition : gene["variant_positions"]) {
            std::stringstream line;
            const std::string refCodon = Strip(variantPosition["ref_codon"]);
            line << "<tr class=\"var\">\n"
                 << "<td>" << refCodon[0] << " " << refCodon[1] << " " << refCodon[2] << "</td>\n"
                 << "<td>" << Strip(variantPosition["ref_amino_acid"]) << "</td>\n"
                 << "<td>" << variantPosition["ref_position"] << "</td>";
            std::string prefix = line.str();
            line.str("");
            bool first = true;
            for (auto& variant_amino_acid : variantPosition["variant_amino_acids"]) {
                for (auto& variant_codons : variant_amino_acid["variant_codons"]) {
                    bool mutated[]{refCodon[0] != Strip(variant_codons["codon"])[0],
                                   refCodon[1] != Strip(variant_codons["codon"])[1],
                                   refCodon[2] != Strip(variant_codons["codon"])[2]};
                    line << "<td>" << Strip(variant_amino_acid["amino_acid"]) << "</td>";
                    line << "<td>";
                    for (int j = 0; j < 3; ++j) {
                        if (mutated[j]) line << "<b style=\"color:#E90032; font-weight:normal\">";
                        line << Strip(variant_codons["codon"])[j] << " ";
                        if (mutated[j]) line << "</b>";
                    }

                    double fOrig = variant_codons["frequency"];
                    double fTmp;
                    int exp = 0;
                    do {
                        fTmp = fOrig * std::pow(10, ++exp);
                    } while (static_cast<int>(fTmp) < 10);
                    fOrig = static_cast<int>(fOrig * std::pow(10, exp));
                    fOrig /= std::pow(10, exp - 2);
                    line << "<td>" << fOrig << "</td>";
                    if (first) {
                        out << prefix << line.str();
                        out << "<td>" << variantPosition["coverage"] << "</td>";
                        first = false;
                    } else {
                        out << "<tr class=\"var\"><td></td><td></td><td></td>" << line.str()
                            << "<td></td>";
                    }
                    out << "<td>" << Strip(variant_codons["known_drm"]) << "</td>";
                    int col = 0;
                    for (auto& hit : variant_codons["haplotype_hit"]) {
                        if (hit)
                            out << "<td style=\"background-color:" << colors.at(col % colors.size())
                                << "\"></td>";
                        else
                            out << "<td></td>";
                        ++col;
                    }
                    out << "</tr>" << std::endl;
                    line.str("");

                    out << R"(
                        <tr class="msa">
                        <td colspan=3 style="background-color: white"></td>
                        <td colspan=14 style="padding:0; margin:0">
                        <table style="padding:0; margin:0" class="msacounts">
                        <col width="50px" />
                        <col width="67px" />
                        <col width="67px" />
                        <col width="67px" />
                        <col width="67px" />
                        <col width="67px" />
                        <col width="67px" />
                        <tr style="padding:0">
                        <th style="padding:2px 0 0px 0">Pos</th>
                        <th style="padding:2px 0 0px 0">A</th>
                        <th style="padding:2px 0 0px 0">C</th>
                        <th style="padding:2px 0 0px 0">G</th>
                        <th style="padding:2px 0 0px 0">T</th>
                        <th style="padding:2px 0 0px 0">-</th>
                        <th style="padding:2px 0 0px 0">N</th>
                        </tr>
                        )";

                    for (auto& column : variantPosition["msa"]) {
                        int relPos = column["rel_pos"];
                        out << "<tr><td>" << relPos << "</td>" << std::endl;
                        for (int j = 0; j < 6; ++j) {
                            out << "<td style=\"";
                            if (relPos >= 0 && relPos < 3 &&
                                j ==
                                    Data::NucleotideToTag(Strip(variant_codons["codon"])[relPos])) {
                                out << "color:#B50A36;";
                            }
                            if (j == Data::NucleotideToTag(Strip(column["wt"])[0]))
                                out << "font-weight:bold;";
                            out << "\">" << column[std::string(1, Data::TagToNucleotide(j))]
                                << "</td>" << std::endl;
                        }
                        out << "</tr>" << std::endl;
                    }
                    out << "</table></tr>" << std::endl;
                }
            }
        }
    }
    out << "</table>" << std::endl;

    if (!config.dbVersion.empty()) out << "<b><sup>*</sup>" << config.dbVersion << "</b>";
    out << R"(
            <details style="margin-bottom: 20px;margin-top:15px">
            <summary>Legend</summary>
            <div style="padding-left:20px">)";

    out << "<p>General:<br/><ul>" << std::endl;
    if (hasConf) {
        out << "<li>Every table represents a gene.</li>" << std::endl;
        out << "<li>Positions are relative to the current gene.</li>" << std::endl;
    } else {
        out << "<li>There is at maximum one table with an \"Unnamed ORF\"</li>" << std::endl;
        out << "<li>Reading frame starts at the first position of the reference used for "
               "alignment.</li>"
            << std::endl;
        out << "<li>The left side of the table shows major codons / AAs observed in this "
               "sample.</li>"
            << std::endl;
    }
    out << "<li>Each row stands for a mutated amino acid.</li>" << std::endl;
    out << "<li>Positions without significant mutations are omitted.</li>" << std::endl;
    out << "<li>All coordinates are in reference space.</li>" << std::endl;
    out << "<li>The mutated nucleotide is highlighted in the codon.</li>" << std::endl;
    out << "<li>Percentage is per codon.</li>" << std::endl;
    out << "<li>Coverage includes deletions.</li>" << std::endl;
    out << "<li>Drugs affected by known drug resistance mutations are listed in the corresponding "
           "column.</li>"
        << std::endl;
    out << "</ul>" << std::endl;
    out << "<p>Alignment Details:</p>" << std::endl;
    out << "<ul>" << std::endl;
    out << "<li>Clicking on a row unfolds the counts of the multiple sequence alignment of the "
           "codon position and up to ±3 surrounding positions.</li>"
        << std::endl;
    out << "<li>Nucleotides of this codon are in red and wild type in bold.</li>" << std::endl;
    out << "</ul>" << std::endl;
    out << "<p>Limitations:</p>" << std::endl;
    out << "<ul>" << std::endl;
    out << "<li>Deletions and insertions are being ignored in this version.</li>" << std::endl;
    out << "</ul>" << std::endl;
    if (numHaplotypes > 0) {
        out << "<p>Haplotypes:</p>" << std::endl;
        out << "<ul>" << std::endl;
        out << "<li>The row-wise variant calls are \"transposed\" onto the per column "
               "haplotypes.</li>"
            << std::endl;
        out << "<li>For each variant, the haplotype shows a colored box, wild type is represented "
               "by plain dark gray.</li>"
            << std::endl;
        out << "<li>A color gradiant helps to distinguish between columns. Colors are purely for "
               "the visualization.</li>"
            << std::endl;
        out << "<li>Haplotypes are sorted in descending order by their relative abundance in "
               "percent.</li>"
            << std::endl;
        out << "<li>Haplotypes are assigned a single or combination of letters for documentation "
               "purposes.</li>"
            << std::endl;
        if (hasConf) out << "<li>Haplotypes are phased across genes.</li>" << std::endl;
        out << "</ul>" << std::endl;
    }
    out << "<p>This software is for research only and has not been clinically "
           "validated!</p></div></details>"
        << std::endl;
}
}
}  //::PacBio::Juliet