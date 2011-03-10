#ifndef CHAIN_MODEL_H
#define CHAIN_MODEL_H

#include "biolcccexception.h"
#include "chemicalbasis.h"

namespace BioLCCC
{

//! Converts energy profile of a peptide/protein to a profile of probabilities.
/*!
    This function converts the energy profile of a peptide/protein to 
    a profile of distribution probabilities. Probability = exp(E_effective).
*/
std::vector<double> calculateBoltzmannFactorProfile(
    const std::vector<double> &effectiveEnergyProfile);

//! Calculates coefficient of distribution of a polymer using the chain model.
double calculateKdChain(
    const std::vector<ChemicalGroup> &parsedSequence,
    const double secondSolventConcentration,
    const ChemicalBasis &chemBasis,
    const double columnPoreSize,
    const double columnRelativeStrength,
    const double temperature) throw (BioLCCCException);
}

#endif
