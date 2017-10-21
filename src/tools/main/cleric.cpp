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

#include <exception>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include <pbcopper/cli/CLI.h>
#include <pbcopper/utility/FileUtils.h>

#include <pbbam/BamReader.h>
#include <pbbam/DataSet.h>
#include <pbbam/FastaReader.h>

#include <pacbio/cleric/Cleric.h>
#include <pacbio/cleric/ClericSettings.h>

namespace PacBio {
namespace Cleric {
static void ParsePositionalArgs(const std::vector<std::string>& args, std::string* bamPath,
                                std::string* fromReference, std::string* fromReferenceName,
                                std::string* toReference, std::string* toReferenceName,
                                std::string* outputFile)
{
    using BAM::DataSet;

    auto SetBamInput = [&bamPath, &fromReferenceName](const DataSet& ds) {
        if (!bamPath->empty()) throw std::runtime_error("Only one BAM input is allowed!");

        const auto bamfiles = ds.BamFiles();
        if (bamfiles.size() != 1) throw std::runtime_error("Only one bam file is allowed!");

        const auto header = bamfiles.front().Header();
        *bamPath = bamfiles.front().Filename();
        if (header.Sequences().empty())
            throw std::runtime_error("Could not find reference sequence name");
        *fromReferenceName = header.Sequences().begin()->Name();
    };

    std::vector<std::string> fastaPaths;
    for (const auto& i : args) {
        const bool fileExist = PacBio::Utility::FileExists(i);
        if (!fileExist) {
            if (!outputFile->empty())
                throw std::runtime_error(
                    "Only one output file allowed. Following files do not exist: " + *outputFile +
                    " and " + i);
            *outputFile = i;
            continue;
        }
        DataSet ds(i);

        switch (ds.Type()) {
            case BAM::DataSet::TypeEnum::SUBREAD:
            case BAM::DataSet::TypeEnum::ALIGNMENT:
            case BAM::DataSet::TypeEnum::CONSENSUS_ALIGNMENT:
                SetBamInput(ds);
                break;
            case BAM::DataSet::TypeEnum::REFERENCE:
                fastaPaths.push_back(i);
                break;
            default:
                throw std::runtime_error("Unsupported input file: " + i + " of type " +
                                         DataSet::TypeToName(ds.Type()));
        }
    }

    for (const auto& fasta : fastaPaths) {
        DataSet ds(fasta);
        const auto fastaFiles = ds.FastaFiles();
        if (fastaFiles.size() != 1)
            throw std::runtime_error("Only one fasta file allowed per dataset: " + fasta);
        BAM::FastaReader msaReader(fastaFiles.front());

        BAM::FastaSequence f;
        while (msaReader.GetNext(f)) {
            if (f.Name() == *fromReferenceName) {
                if (fromReference->empty()) {
                    *fromReference = f.Bases();
                    std::transform(fromReference->begin(), fromReference->end(),
                                   fromReference->begin(), ::toupper);
                } else
                    throw std::runtime_error("Multiple original references provided!");
            } else if (toReference->empty()) {
                *toReference = f.Bases();
                std::transform(toReference->begin(), toReference->end(), toReference->begin(),
                               ::toupper);
                *toReferenceName = f.Name();
            } else {
                throw std::runtime_error("Multiple target references provided!");
            }
        }
    }
}

static int Runner(const PacBio::CLI::Results& options)
{
    // Check args size, as pbcopper does not enforce the correct number
    if (options.PositionalArguments().empty()) {
        std::cerr << "ERROR: Please provide BAM input, see --help" << std::endl;
        return EXIT_FAILURE;
    }
    if (options.PositionalArguments().size() < 2 || options.PositionalArguments().size() >= 5) {
        std::cerr << "ERROR: Please provide _one_ BAM input, maximal _two_ "
                     "FASTA files, and _one_ output file. See --help"
                  << std::endl;
        return EXIT_FAILURE;
    }

    // Parse options
    ClericSettings settings(options);

    std::string bamPath;
    std::string fromReference;
    std::string fromReferenceName;
    std::string toReference;
    std::string toReferenceName;
    std::string outputFile;
    bool alreadyAligned = false;

    ParsePositionalArgs(settings.InputFiles, &bamPath, &fromReference, &fromReferenceName,
                        &toReference, &toReferenceName, &outputFile);

    // parse pre-aligned FASTA file
    if (options.PositionalArguments().size() == 2) {
        if (settings.PrealignedFile.empty()) {
            std::cerr << "ERROR: You need to provide a pre-aligned FASTA file with --aln"
                      << std::endl;
            return EXIT_FAILURE;
        }

        if (!PacBio::Utility::FileExists(settings.PrealignedFile)) {
            std::cerr << "ERROR: The pre-aligned FASTA file '" << settings.PrealignedFile
                      << "' does not exist" << std::endl;
            return EXIT_FAILURE;
        }

        // Load all pre-aligned sequence at once
        std::vector<BAM::FastaSequence> allPrealignedSeqs{
            BAM::FastaReader::ReadAll(settings.PrealignedFile)};

        if (allPrealignedSeqs.size() != 2) {
            std::cerr << "ERROR: The pre-aligned FASTA file '" << settings.PrealignedFile
                      << "' has to contain _exactly_ 2 sequences (contains "
                      << allPrealignedSeqs.size() << ')' << std::endl;
            return EXIT_FAILURE;
        }

        if (allPrealignedSeqs.front().Name() == fromReferenceName) {
            fromReference = allPrealignedSeqs.front().Bases();

            toReferenceName = allPrealignedSeqs.back().Name();
            toReference = allPrealignedSeqs.back().Bases();
        } else {
            if (allPrealignedSeqs.back().Name() == fromReferenceName) {
                fromReference = allPrealignedSeqs.back().Bases();

                toReferenceName = allPrealignedSeqs.front().Name();
                toReference = allPrealignedSeqs.front().Bases();
            } else {
                std::cerr << "ERROR: The pre-aligned FASTA file '" << settings.PrealignedFile
                          << "' does not contain a sequence with name '" << fromReferenceName
                          << '\'' << std::endl;
                return EXIT_FAILURE;
            }
        }

        if (fromReference.size() != toReference.size()) {
            std::cerr << "ERROR: The reference sequence '" << fromReferenceName
                      << "' and the query sequence '" << toReferenceName
                      << "' have different lengths (" << fromReference.size() << " vs "
                      << toReference.size() << ')' << std::endl;
            return EXIT_FAILURE;
        }

        alreadyAligned = true;
    }

    if (outputFile.empty()) outputFile = PacBio::Utility::FilePrefix(bamPath) + "_cleric";

    Cleric cleric(bamPath, outputFile, fromReference, fromReferenceName, toReference,
                  toReferenceName, alreadyAligned);

    return EXIT_SUCCESS;
}
}
};

// Entry point
int main(int argc, char* argv[])
{
    return PacBio::CLI::Run(argc, argv, PacBio::Cleric::ClericSettings::CreateCLI(),
                            &PacBio::Cleric::Runner);
}