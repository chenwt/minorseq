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
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <pbbam/DataSet.h>
#include <pbcopper/cli/CLI.h>
#include <pbcopper/utility/FileUtils.h>

#include <pacbio/fuse/Fuse.h>
#include <pacbio/fuse/FuseSettings.h>

namespace PacBio {
namespace Fuse {

static int Runner(const PacBio::CLI::Results& options)
{
    // Check args size, as pbcopper does not enforce the correct number
    if (options.PositionalArguments().empty()) {
        std::cerr << "ERROR: Please provide BAM input, see --help" << std::endl;
        return EXIT_FAILURE;
    }

    // Parse options
    FuseSettings settings(options);

    Fuse fuse(settings.InputFile, settings.MinCoverage);

    auto outputFile = settings.OutputFile;
    const bool isXml = Utility::FileExtension(outputFile) == "xml";
    if (isXml) boost::ireplace_all(outputFile, ".referenceset.xml", ".fasta");

    std::ofstream outputFastaStream(outputFile);
    const auto consensus = fuse.ConsensusSequence();
    outputFastaStream << ">CONSENSUS" << std::endl;
    outputFastaStream << consensus << std::endl;

#if 0
    // Write Dataset
    using BAM::DataSet;
    const std::string metatype = "PacBio.ReferenceFile.ReferenceFastaFile";
    DataSet fuseSet(DataSet::TypeEnum::REFERENCE);
    BAM::ExternalResource resource(metatype, outputFile);
    fuseSet.ExternalResources().Add(resource);
    fuseSet.Name(fuseSet.TimeStampedName());

    const auto outputPrefix = outputFile.substr(0, outputFile.size() - 6);
    std::ofstream fuseDSout(outputPrefix + ".referenceset.xml");
    fuseSet.SaveToStream(fuseDSout);
#endif

    return EXIT_SUCCESS;
}
}
};

// Entry point
int main(int argc, char* argv[])
{
    return PacBio::CLI::Run(argc, argv, PacBio::Fuse::FuseSettings::CreateCLI(),
                            &PacBio::Fuse::Runner);
}