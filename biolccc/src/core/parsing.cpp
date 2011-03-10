#include <cmath>
#include "parsing.h"

namespace BioLCCC
{
ParsingException::ParsingException(std::string message):
        BioLCCCException(message) 
{
};

std::vector<double> calculateMonomerEnergyProfile(
    const std::vector<ChemicalGroup> &parsedSequence,
    const ChemicalBasis & chemBasis,
    const double secondSolventConcentration,
    const double columnRelativeStrength, 
    const double temperature) throw (BioLCCCException)
{
    if (parsedSequence.size() < 3)
    {
        throw BioLCCCException(
            "The parsed sequence contains too little chemical groups.");
    }

    if (columnRelativeStrength == 0.0)
    {
        return std::vector<double> (parsedSequence.size()-2, 0.0);
    }

    // Due to the preliminary scaling the binding energy of water equals zero.
    double Q = exp(columnRelativeStrength * 
                   (chemBasis.secondSolventBindEnergy() - 0.0) *
                   293.0 / temperature);
    // Nb = (DensityB * %B / MB) / 
    //      (DensityB * %B / MB + DensityA * %A / MA)
    // where DensityA and DensityB are the corresponding densities and
    // MA and MB are the corresponding molecular weights.
    double Nb =
        secondSolventConcentration * chemBasis.secondSolventDensity()
        / chemBasis.secondSolventAverageMass() 
        / ( secondSolventConcentration * chemBasis.secondSolventDensity() 
                / chemBasis.secondSolventAverageMass()
            + (100.0 - secondSolventConcentration) 
              * chemBasis.firstSolventDensity()
              / chemBasis.firstSolventAverageMass());

    double Eab = 0.0;
    if (chemBasis.snyderApproximation()) 
    {
        Eab = Nb * chemBasis.secondSolventBindEnergy();
    }
    else
    {
        Eab = 0.0 + 1.0 / columnRelativeStrength * log( 1.0 - Nb + Nb * Q );
    }

    std::vector<double> monomerEnergyProfile;
    for (std::vector<BioLCCC::ChemicalGroup>::const_iterator residue =
                ++(parsedSequence.begin());
            residue != --(parsedSequence.end());
            ++residue)
    {
        double residueEnergy = residue->bindEnergy();
        double residueArea = residue->bindArea();

        // Adding the energy of the N-terminal group to the first residue.
        if (residue == ++(parsedSequence.begin()))
        {
            residueEnergy += parsedSequence.begin()->bindEnergy();
            residueArea += parsedSequence.begin()->bindArea();
        }

        // Adding the energy of the C-terminal group to the last residue.
        else if (residue == --(--(parsedSequence.end())))
        {
            residueEnergy += (--(parsedSequence.end()))->bindEnergy();
            residueArea += (--(parsedSequence.end()))->bindArea();
        }

        monomerEnergyProfile.push_back(
            columnRelativeStrength * 
                //(residueEnergy - Eab) * 293.0 / temperature);
                (residueEnergy - residueArea * Eab) * 293.0 / temperature);
    }
    return monomerEnergyProfile;
}

std::vector<double> calculateSegmentEnergyProfile(
    const std::vector<double> &monomerEnergyProfile,
    const double monomerLength,
    const double kuhnLength)
{
    std::vector<double> segmentEnergyProfile; 
    std::vector<double>::const_iterator monomerEnergyIt =
        monomerEnergyProfile.begin();
    double kuhnLeftBorder  = 0;
    double monomerLeftBorder  = 0;
    double sumEnergy = 0.0;
    bool kuhnSegmentOpen = true;

    while (monomerEnergyIt != monomerEnergyProfile.end()) 
    {
        if ((kuhnLeftBorder + kuhnLength) >= 
                (monomerLeftBorder + monomerLength))
        {
            sumEnergy += *(monomerEnergyIt) * 
                (monomerLeftBorder + monomerLength - 
                    std::max(monomerLeftBorder, kuhnLeftBorder)) / 
                monomerLength;
            kuhnSegmentOpen = true;
            monomerLeftBorder += monomerLength;
            ++monomerEnergyIt;
        }
        else 
        {
            sumEnergy += *(monomerEnergyIt) * 
                (kuhnLeftBorder + kuhnLength - 
                    std::max(monomerLeftBorder, kuhnLeftBorder)) / 
                monomerLength;
            segmentEnergyProfile.push_back(sumEnergy);
            sumEnergy = 0.0;
            kuhnSegmentOpen = false;
            kuhnLeftBorder += kuhnLength;
        }
    }

    if (kuhnSegmentOpen)
    {
        segmentEnergyProfile.push_back(sumEnergy);
    }

    return segmentEnergyProfile;
}

std::vector<ChemicalGroup> parseSequence(
    const std::string &source,
    const ChemicalBasis &chemBasis) throw(BioLCCCException)
{
    std::vector<ChemicalGroup> parsedSequence;
    ChemicalGroup NTerminus;
    ChemicalGroup CTerminus;

    // At first we need to strip the sequence from adjacent amino acids.
    // If a source sequence contains them, it should contain two dots, so we
    // need the part of sequence between them.
    std::size_t firstDotPosition = 0;
    std::size_t secondDotPosition = 0;

    // We'll use this variable to contain peptide sequence without adjacent
    // amino acids.
    std::string strippedSource = source;

    firstDotPosition = source.find(".");

    if (firstDotPosition != std::string::npos)
    {
        secondDotPosition = source.find(".", firstDotPosition+1);
        if (secondDotPosition != std::string::npos)
        {

            // If a source sequence contains more that two dots, it's broken.
            if (source.find(".", secondDotPosition+1) != std::string::npos)
            {
                throw ParsingException(
                    "The sequence " + source +" contains more than two dots.");
            }
            else
            {
                strippedSource = source.substr(firstDotPosition+1,
                    secondDotPosition - firstDotPosition - 1);
            }
        }
        // If a source sequence contains only one dot, it's broken.
        else
        {
            throw ParsingException(
                "The sequence " + source + " contains only one dot.");
        }
    }

    // Than goes parsing.
    std::size_t NTerminusPosition = 0;

    // First we need to check the stripped source sequence for
    // the N-Terminal group.
    NTerminus = chemBasis.defaultNTerminus();
    for (std::map<std::string,ChemicalGroup>::const_iterator
            NTerminusIterator = chemBasis.chemicalGroups().begin();
            NTerminusIterator != chemBasis.chemicalGroups().end();
            NTerminusIterator++)
    {
        if (NTerminusIterator->second.isNTerminal())
        {
            if (strippedSource.find(NTerminusIterator->second.label()) ==
                    (size_t)0)
            {
                NTerminus = NTerminusIterator->second;
                NTerminusPosition = NTerminusIterator->second.label().size();
            }
        }
    }

    // Then we need to found the location of the C-Terminus.
    CTerminus = chemBasis.defaultCTerminus();
    std::size_t CTerminusPosition;
    CTerminusPosition = strippedSource.find("-", NTerminusPosition);
    if (CTerminusPosition != std::string::npos)
    {
        // The sequence should not contain hyphens after C-terminal group.
        if (strippedSource.find("-", CTerminusPosition+1) != std::string::npos)
        {
            throw ParsingException(
                "The sequence " + source +
                " contains hyphen after C-terminal group.");
        }

        // Searching for known C-terminal groups.
        bool CTerminusFound = false;
        for (std::map<std::string,ChemicalGroup>::const_iterator
                CTerminusIterator = chemBasis.chemicalGroups().begin();
                CTerminusIterator != chemBasis.chemicalGroups().end();
                CTerminusIterator++)
        {
            if (CTerminusIterator->second.isCTerminal())
            {
                if (strippedSource.find(CTerminusIterator->second.label(),
                    CTerminusPosition) != std::string::npos)
                {
                    CTerminus = CTerminusIterator->second;
                    CTerminusFound = true;
                }
            }
        }
        if (!CTerminusFound)
        {
            throw ParsingException(
                "The sequence " + source +
                " contains unknown C-terminal group\"" + 
                source.substr(CTerminusPosition) + "\".");
        }
    }
    else
    {
        CTerminusPosition = strippedSource.size();
    }

    // At this step we obtain the sequence of a peptide without adjacent
    // amino acids and terminal groups.
    strippedSource = strippedSource.substr(
        NTerminusPosition, CTerminusPosition-NTerminusPosition);

    // We need to check whether it contains any non-letter characters.
    for (std::size_t i=0; i<strippedSource.size(); i++)
    {
        if (!(((int(strippedSource[i]) >= int('a')) &&
                (int(strippedSource[i]) <= int('z'))) ||
                ((int(strippedSource[i]) >= int('A')) &&
                 (int(strippedSource[i]) <= int('Z')))))
        {
            throw ParsingException(
                "The sequence " + source +
                " contains a non-letter character.");
        }
    }

    // Then we divide the whole sequence into aminoacids.
    bool aminoAcidFound;
    size_t curPos = 0;
    while (curPos < strippedSource.size())
    {
        aminoAcidFound = false;
        for (std::map<std::string,ChemicalGroup>::const_iterator
                aminoAcidIterator = chemBasis.chemicalGroups().begin();
                aminoAcidIterator != chemBasis.chemicalGroups().end();
                aminoAcidIterator++)
        {
            if (strippedSource.compare(curPos,
                aminoAcidIterator->second.label().size(),
                aminoAcidIterator->second.label()) == 0)
            {
                curPos += aminoAcidIterator->second.label().size();
                parsedSequence.push_back(aminoAcidIterator->second);
                aminoAcidFound = true;
                break;
            }
        }

        if (!aminoAcidFound)
        {
            throw ParsingException(
                "The sequence " + source + " contains unknown amino acid \"" + 
                source.substr(curPos, 1) + "\".");
        }
    }
    parsedSequence.insert(parsedSequence.begin(), NTerminus);
    parsedSequence.push_back(CTerminus);
    return parsedSequence;
}
}

