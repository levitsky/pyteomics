#ifndef ROD_MODEL_H
#define ROD_MODEL_H

#include "biolcccexception.h"
#include "chemicalbasis.h"

namespace BioLCCC
{

//! Calculates a term in the special expression for the partition function.
double partitionFunctionRodPartiallySubmergedTermSpecial(
    double segmentLength, double slitWidth, double layerWidth,
    int N, int n1);

//! Calculates a term in the general expression for the partition function.
double partitionFunctionRodPartiallySubmergedTermGeneral(
    double segmentLength, double slitWidth, double layerWidth,
    int N, int n1, int n2);

//! Calculates the adsorption energy of the first n segments of a rod.
double rodAdsorptionEnergy(const std::vector<double> & rodEnergyProfile,
                           int n1, int n2) throw(BioLCCCException);

//! Calculates the partition function of a rod in a slit with impenetrable walls.
double partitionFunctionRodFreeSlit(double rodLength,
                                    double slitWidth);

//! Calculates Z of a rod partially submerged into the adsorbing layer.
/*!
    This function employs the general algorithm of partion function calculation
    which is valid for all arguments' values.
 */
double partitionFunctionRodPartiallySubmergedGeneral(
    double segmentLength,
    double slitWidth,
    double layerWidth,
    const std::vector<double> & rodEnergyProfile,
    bool reversed = false) throw(BioLCCCException);

//! Calculates Z of the rod partially submerged into an adsorbing layer.
/*!
    This function employs the special algorithm of partion function calculation
    which is valid only when rodWidth < slitWidth - 2 * layerWidth.
 */
double partitionFunctionRodPartiallySubmergedSpecial(
    double segmentLength,
    double slitWidth,
    double layerWidth,
    const std::vector<double> & rodEnergyProfile,
    bool reversed = false) throw(BioLCCCException);

//! Calculates the partition function of a rod in a slit without walls.
double partitionFunctionRodFreeVolume(double rodLength,
                                      double slitWidth);

//! Calculates the coefficient of distribution of a polymer using the rod model.
double calculateKdRod(
    const std::vector<ChemicalGroup> &parsedSequence,
    const double secondSolventConcentration,
    const ChemicalBasis &chemBasis,
    const double columnPoreSize,
    const double columnRelativeStrength,
    const double temperature
    ) throw(BioLCCCException);
}

#endif
