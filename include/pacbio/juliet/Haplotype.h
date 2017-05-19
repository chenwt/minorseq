#pragma once

#include <pacbio/juliet/HaplotypeType.h>
#include <pacbio/util/Termcolor.h>

#include <pbcopper/json/JSON.h>

namespace PacBio {
namespace Juliet {
struct Haplotype
{
    std::string Name;
    std::vector<std::string> Names;
    std::vector<std::string> Codons;
    double SoftCollapses = 0;
    double GlobalFrequency = 0;
    int Flags = 0;

    double Size() const { return Names.size() + SoftCollapses; }

    std::string ConcatCodons() const
    {
        return std::accumulate(Codons.begin(), Codons.end(), std::string(""));
    }

    void SetCodons(std::vector<std::string>&& codons)
    {
        Codons = std::forward<std::vector<std::string>>(codons);
        for (const auto& c : Codons) {
            if (c.find('-') != std::string::npos)
                Flags |= static_cast<int>(HaplotypeType::WITH_GAP);
            if (c.find('N') != std::string::npos)
                Flags |= static_cast<int>(HaplotypeType::WITH_HETERODUPLEX);
            if (c.find(' ') != std::string::npos) Flags |= static_cast<int>(HaplotypeType::PARTIAL);
        }
    }

    friend std::ostream& operator<<(std::ostream& stream, const Haplotype& h)
    {
        stream << h.Size() << "\t";
        for (const auto& c : h.Codons) {
            if (c != "ATG" && c != "AAA" && c != "TAT" && c != "GGA" && c != "ACC")
                stream << termcolor::white;
            if (c == "TTG") stream << termcolor::red;
            if (c == "AGA") stream << termcolor::green;
            if (c == "TGT") stream << termcolor::blue;
            if (c == "GCA") stream << termcolor::blue;
            if (c == "TAC") stream << termcolor::yellow;
            stream << " " << c << termcolor::reset;
        }
        return stream;
    }

    JSON::Json ToJson() const
    {
        using namespace JSON;
        Json root;
        root["name"] = Name;
        root["reads_hard"] = Names.size();
        root["reads_soft"] = Size();
        root["frequency"] = GlobalFrequency;
        root["read_names"] = Names;
        root["codons"] = Codons;
        return root;
    }
};
}
}