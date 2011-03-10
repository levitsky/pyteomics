#ifndef PARSING_H
#define PARSING_H

#include "biolcccexception.h"
#include "chemicalbasis.h"

namespace BioLCCC
{

//! This exception is raised when a parsing process cannot be completed.
class ParsingException : public BioLCCCException
{
public:
    //! Constructs an instance ParsingException with the given \a message.
    ParsingException(std::string message);
};

//! Parses a given peptide sequence.
/*!
    Parses a given peptide sequence \a source using \a chemBasis.
    Returns a vector with chemical groups.

    Throws ParsingException if the peptide is not parseable.
*/
std::vector<ChemicalGroup> parseSequence(
    const std::string &source,
    const ChemicalBasis &chemBasis) throw(BioLCCCException);

//! Calculates the effective energy profile of monomers of the polymer chain.
/*!
    Each element in this profile contains the adsorption energy of a single
    monomer.

    E_{effective} = alpha * ( E_{residue} - E_{ab} ) * 293.0 / temperature,
    and Eab is the effective energy of binary solvent binding to the
    stationary phase.
    Eab = Ea + 1 / alpha * ln( 1 + Nb + Nb * exp( alpha * ( Ea - Eb ) ) )
*/
std::vector<double> calculateMonomerEnergyProfile(
    const std::vector<ChemicalGroup> &parsedSequence,
    const ChemicalBasis & chemBasis,
    const double secondSolventConcentration,
    const double columnRelativeStrength, 
    const double temperature) throw (BioLCCCException);

//! Calculates the effective energy profile of segments of the polymer chain.
/*!
    The energy of a single segment equals to sum of the energies of the whole
    monomers containing in the segment PLUS the proportional parts of the
    energies of monomers which cross the borders of the segment. 
*/
std::vector<double> calculateSegmentEnergyProfile(
    const std::vector<double> &monomerEnergyProfile,
    const double monomerLength,
    const double kuhnLength);
}

#endif
