#ifndef BIOLCCC_H
#define BIOLCCC_H

#include "auxiliary.h"
#include "biolcccexception.h"
#include "chemicalbasis.h"
#include "chromoconditions.h"
#include "parsing.h"
#include "chain_model.h"
#include "rod_model.h"

//! Apart from classes, BioLCCC contains calculation methods and constants.
namespace BioLCCC
{

//! A ChromoConditions instance with the standard chromatographic conditions.
const ChromoConditions standardChromoConditions = ChromoConditions();

//! A ChemicalBasis instance of predefined RP_ACN_TFA_CHAIN.
const ChemicalBasis rpAcnTfaChain =
    ChemicalBasis(RP_ACN_TFA_CHAIN);

//! A ChemicalBasis instance of predefined RP_ACN_FA_ROD.
const ChemicalBasis rpAcnFaRod =
    ChemicalBasis(RP_ACN_FA_ROD);

//! Calculates the retention time of a peptide.
/*!
    Calculates the retention time of a peptide with given \a sequence
    using the given description of chromatographic conditions \a conditions and
    set of physicochemical constants \a chemBasis.
    
    If \a continueGradient is true, than the last section of a gradient is
    prolonged.

    If \a backwardCompatibility is true, than the calculated RT will be 
    a multiple of dV/flow rate. This type of calculations was used in the 
    version 1.1.0 and earlier.
*/
double calculateRT(const std::string &sequence,
                   const ChemicalBasis & chemBasis,
                   const ChromoConditions & conditions =
                       standardChromoConditions,
                   const int numInterpolationPoints = 0,
                   const bool continueGradient = true,
                   const bool backwardCompatibility = false) 
                   throw(BioLCCCException);

//! Calculates the average (molar) mass of a peptide.
/*!
    Calculates the average (molar) mass of a peptide with given
    \a sequence using the given set of physicochemical constants \a chemBasis.
*/
double calculateAverageMass(const std::string &sequence,
                            const ChemicalBasis &chemBasis)
                            throw(BioLCCCException);
//! Calculates the monoisotopic mass of a peptide.
/*!
    Calculates the monoisotopic mass of a peptide with given
    \a sequence using the given set of physicochemical constants \a chemBasis.
*/
double calculateMonoisotopicMass(const std::string &sequence,
                                 const ChemicalBasis &chemBasis)
                                 throw(BioLCCCException);

//! Calculates the coefficient of distribution Kd for the given peptide.
/*!
    Calculates the coefficient of distribution Kd (i.e. the ratio of
    concentrations of a peptide in the pores and in the interstitial volume).
    
    \param sequence The sequence of a peptide.
    \param secondSolventConcentration The concentration of the second solvent in
    the liquid phase
    \param chemBasis The set of the physicochemical constants.
    \param columnPoreSize The size of adsorbent pores.
    \param columnRelativeStrength The relative strength of adsorption.
    \param temperature Temperature of the column.
*/
double calculateKd(const std::string &sequence,
                   const double secondSolventConcentration,
                   const ChemicalBasis &chemBasis,
                   const double columnPoreSize = 100.0,
                   const double columnRelativeStrength = 1.0,
                   const double temperature = 293.0)
                   throw(BioLCCCException);
}
#endif
