#include <cmath>
#include <stdlib.h>
#include "chain_model.h"
#include "parsing.h"

namespace BioLCCC
{

std::vector<double> calculateBoltzmannFactorProfile(
    const std::vector<double> &effectiveEnergyProfile)
{
    std::vector<double> boltzmannFactorProfile; 

    for (std::vector<double>::const_iterator energy =
                effectiveEnergyProfile.begin();
            energy != effectiveEnergyProfile.end();
            ++energy)
    {
        boltzmannFactorProfile.push_back(exp(*energy));
    }

    return boltzmannFactorProfile;
}

double calculateKdChain(
						const std::vector<ChemicalGroup> &parsedSequence,
						const double secondSolventConcentration,
						const ChemicalBasis &chemBasis,
						const double columnPoreSize,
						const double columnRelativeStrength,
						const double temperature
						) throw (BioLCCCException)
{
    // At first, we need to convert the energy profile to a profile of 
    // distribution probabilities. Probability = exp(E_effective),
    // where E_effective = E_of_residue - Eab,
    // and Eab is an energy of binding for a binary solvent
    // and Eab = Ea + ln ( 1 + Nb + Nb*exp (Ea - Eb) )
    // also corrections of energies due to temperature (energies in exponents
    // are scaled to the new temperature) and the relative strength 
    // are introduced.
    std::vector<std::vector<double> > boltzmannFactorProfiles;
    for (std::vector<double>::const_iterator adsorptionLayerFactor = 
            chemBasis.adsorptionLayerFactors().begin();
         adsorptionLayerFactor != chemBasis.adsorptionLayerFactors().end();
         ++adsorptionLayerFactor) 
    {
        boltzmannFactorProfiles.push_back(
            calculateBoltzmannFactorProfile(
                calculateSegmentEnergyProfile(
                    calculateMonomerEnergyProfile(
                        parsedSequence,
                        chemBasis,
                        secondSolventConcentration,
                        columnRelativeStrength * (*adsorptionLayerFactor),
                        temperature),
                    chemBasis.monomerLength(),
                    chemBasis.kuhnLength())));
    }

    // The size of the lattice must be greater than 
    // (number of adsorbing layers) * 2.
    // double round (double x) {return floor(x+0.5);}
    unsigned int latticeSize = 
        floor(columnPoreSize / chemBasis.kuhnLength() + 0.5);

    // If we want to neglect the partially desorbed states, we need to insert
    // two impenetrable layers right after the near-wall ones.
    if (chemBasis.neglectPartiallyDesorbedStates())
    {
        boltzmannFactorProfiles.push_back(
            std::vector<double>(
                boltzmannFactorProfiles.back().size(), 0.0));
        latticeSize += 2;
    }

    if (latticeSize < boltzmannFactorProfiles.size() * 2)
    {
        throw BioLCCCException(
            "The pore size is too small for the given number of adsorbing "
            "layers.");
    }

    // The density vector correspond to a probability of n-th residue to be in
    // a certain layer between pore walls.
    double *density;

    // The transition matrix used to calculate a density vector of n-th
    // particle from a density vector of (n-1)-th particle.
    double *transitionMatrix;

    // The density buffer vector is used during matrix calculations.
    double *densityBuffer;

    // Memory managment.
    try
    {
        density = new double[latticeSize];
        densityBuffer = new double[latticeSize];
        transitionMatrix = new double[latticeSize*latticeSize];
    }
    catch (...)
    {
        throw BioLCCCException("Cannot allocate memory for calculations");
    }

    // Constructing a density vector for the first amino acid residue.
    // A density is distributed uniformly over all non-adsorbing layers of 
    // the lattice.
    // The density in adsorbing layers is multiplied by a Boltzmann factor of
    // the first segment.

    for (unsigned int i = 0; i < latticeSize; i++)
    {
        density[i] = 1.0;
    }

    for (unsigned int i = 0; i < boltzmannFactorProfiles.size(); ++i) 
    {
        density[i] = boltzmannFactorProfiles[i][0];
        density[latticeSize - i - 1] = boltzmannFactorProfiles[i][0];
    }

    // Debugging facilities.
    //for (unsigned int i = 0; i < latticeSize; i++)
    //{
    //    std::cout << density[i] << " ";
    //}
    //std::cout << std::endl;
    //std::cout << std::endl;

    // Than we construct a basis for the transition matrix. The basis is
    // a diagonal matrix with 4.0/6.0 on the main diagonal and 1.0/6.0 on
    // the side diagonals.

    // Filling the matrix.
    for (unsigned int i = 0; i < latticeSize; i++)
    {
        for (unsigned int j = 0; j < latticeSize; j++)
        {
            switch (abs( j - i ))
            {
            case 0:
            {
                transitionMatrix[j + latticeSize * i] = 4.0/6.0;
                break;
            }
            case 1:
            {
                transitionMatrix[j + latticeSize * i] = 1.0/6.0;
                break;
            }
            default:
                transitionMatrix[j + latticeSize * i] = 0.0;
            }
        }
    }

    // On each step we obtain the density vector for the n-th segment
    // multiplying the transition matrix and the density vector of the 
    // (n-1)th residue.
    // The multiplication starts from the second segment.
    for (unsigned int segmentIndex = 1; 
         segmentIndex < boltzmannFactorProfiles[0].size();
         ++segmentIndex) 
    {
        // Filling the matrix elements that describe the adsorption.
        for (unsigned int layerIndex = 0;
             layerIndex < boltzmannFactorProfiles.size();
             ++layerIndex) 
        {
            int indexShift = layerIndex * ( latticeSize + 1 );
            double boltzmannFactor = 
                boltzmannFactorProfiles[layerIndex][segmentIndex];

            transitionMatrix[indexShift + 0] = 4.0/6.0 * boltzmannFactor;
            transitionMatrix[indexShift + 1] = 1.0/6.0 * boltzmannFactor;
            transitionMatrix[latticeSize*latticeSize - 1 - indexShift - 0] =
                    4.0 / 6.0 * boltzmannFactor;
            transitionMatrix[latticeSize*latticeSize - 1 - indexShift - 1] =
                    1.0 / 6.0 * boltzmannFactor;

            // A segment may enter the second and further adsorption layers from
            // the inner layer (i.e. the layer lying closer to the walls).
            if (layerIndex > 0) 
            {
                transitionMatrix[indexShift - 1] = 1.0/6.0 * boltzmannFactor;
                transitionMatrix[latticeSize*latticeSize - 1 - indexShift + 1] =
                        1.0 / 6.0 * boltzmannFactor;
            }
        }

        // Zeroing the calculation buffer.
        for (unsigned int i = 0; i < latticeSize; i++)
        {
            densityBuffer[i] = 0.0;
        }

        // Multiplying the transition matrix by the density vector. The result
        // is stored in the buffer vector.
        for (unsigned int i = 0; i < latticeSize; i++)
        {
            for (unsigned int j = 0; j < latticeSize; j++)
            {
                densityBuffer[i] = densityBuffer[i] + density[j] *
                                   transitionMatrix[j + i * latticeSize];
            }
        }

        // Transferring the results from the density vector.
        for (unsigned int i = 0; i < latticeSize; i++)
        {
            density[i] = densityBuffer[i];
        }
    }

    // Kd is a sum of elements of the density vector, normalized to the size 
    // of the lattice.
    double Kd=0;
    for (unsigned int i=0; i < latticeSize; i++)
    {
        Kd += density[i];
    }

    if (chemBasis.neglectPartiallyDesorbedStates())
    {
        Kd = Kd / (double)(latticeSize - 2);
    }
    else
    {
        Kd = Kd / (double)(latticeSize);
    }

    // Cleaning memory.
    try
    {
        delete[] density;
        delete[] densityBuffer;
        delete[] transitionMatrix;
    }
    catch (...)
    {
        throw BioLCCCException("Cannot allocate memory for calculations");
    }

    return Kd;
}

}
