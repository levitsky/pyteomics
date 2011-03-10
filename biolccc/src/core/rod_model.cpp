#include <cmath>
#include <cfloat>
#include <algorithm>
#include <functional>
#include <numeric>
#include <vector>
#include <utility>
#include "rod_model.h"
#include "parsing.h"

#include <iostream>
#include <iomanip>

#define PI 3.14159265

namespace BioLCCC
{

double rodAdsorptionEnergy(const std::vector<double> & rodEnergyProfile,
                           int n1, int n2) throw(BioLCCCException)
{
    if ((n1 < 0) || (n1 > rodEnergyProfile.size())
        || (n2 < 0) || (n2 > rodEnergyProfile.size()))
    {
        throw BioLCCCException("Index is out of range.");
    }

    double init = 0;
    return std::accumulate(rodEnergyProfile.begin(),
                           rodEnergyProfile.begin()+n1,
                           init)
           + std::accumulate(rodEnergyProfile.end()-n2,
                           rodEnergyProfile.end(),
                           init);
}

double partitionFunctionRodFreeSlit(double rodLength,
                                    double slitWidth)
{
    // the equation works only if the slit is wider than the rod
    if (rodLength <= slitWidth)
    {
        // full volume without exclusion caused by walls
        return (4 * PI * slitWidth * rodLength * rodLength -
                // minus volume excluded by walls
                2 * PI * rodLength * rodLength * rodLength);
    }
    else
    {
        // This case is considered in the paper as Zc.
        return (2 * PI * slitWidth * slitWidth * rodLength);
    }

}

namespace {
    std::pair<double, double> intersection(
        std::pair<double, double> first_line,
        std::pair<double, double> second_line)
    {
        std::pair<double, double> output;
        if (first_line.first != second_line.first)
        {
            double x = (second_line.second - first_line.second) 
                       / (first_line.first - second_line.first);
            x = ceil(x * 1.0e10) / 1.0e10;
            output = std::make_pair(x, first_line.first*x + first_line.second);
        }
        else
        {
            if (first_line.first >= 0.0)
            {
                output = std::make_pair(-DBL_MAX, -DBL_MAX);
            }
            else
            {
                output = std::make_pair(-DBL_MAX, DBL_MAX);
            }
        }
        return output;
    }
}

double partitionFunctionRodPartiallySubmergedTermGeneral(
    double segmentLength, double slitWidth, double layerWidth,
    int N, int n1, int n2)
{
    const double rodLength = (N - 1) * segmentLength;
    std::vector<std::pair<double, double> > lowerLimits;
    std::vector<std::pair<double, double> > upperLimits;

    // Calculating the limits of integration over theta.
    // A limit is the minimal or maximal value of cos(theta) at which 
    // a rod conformation is still acceptable. 

    // The maximal angle at which the bead next to the n1-th one is still out
    // of the adsorbing layer.
    upperLimits.push_back(std::make_pair(
        -1.0 / n1 / segmentLength,
        layerWidth / n1 / segmentLength));

    // The minimal angle theta at which the bead next to the n2-th one is still
    // out of the adsorbing layer.
    lowerLimits.push_back(std::make_pair(
        -1.0 / (N - n2 - 1) / segmentLength,
        (slitWidth - layerWidth) / (N - n2 - 1) / segmentLength));

    if (n1 > 1)
    {
        // The minimal angle at which the n1-th bead is still in the adsorbing
        // layer.
        lowerLimits.push_back(std::make_pair(
            -1.0 / (n1 - 1) / segmentLength,
            layerWidth / (n1 - 1) / segmentLength));
    }
    if (n2 != 0)
    {
        // The minimal angle at which the last bead does not touch the wall.  
        lowerLimits.push_back(std::make_pair(
            -1.0 / (N - 1) / segmentLength,
            slitWidth / (N - 1) / segmentLength));
        // The maximal angle at which the n2-th bead is still in the adsorbing
        // layer.  
        upperLimits.push_back(std::make_pair(
            -1.0 / (N - n2) / segmentLength,
            (slitWidth - layerWidth) / (N - n2) / segmentLength));
    }

    std::vector<std::pair<double, double> > allLimits;
    allLimits.insert(allLimits.begin(), lowerLimits.begin(), lowerLimits.end());
    allLimits.insert(allLimits.begin(), upperLimits.begin(), upperLimits.end());

    // Finding the point at which the dependency of the limits of integration
    // over theta on z may change.
    std::vector<double> critPoints;
    for (std::vector<std::pair<double, double> >::const_iterator lineIter1 =
            allLimits.begin();
         lineIter1 < allLimits.end();
         ++lineIter1)
    {
        for (std::vector<std::pair<double, double> >::const_iterator lineIter2 =
                lineIter1 + 1;
             lineIter2 < allLimits.end();
             ++lineIter2)
        {
            critPoints.push_back(intersection(*lineIter1, *lineIter2).first);
        }
        critPoints.push_back(
            intersection(*lineIter1, std::make_pair(0.0, 1.0)).first);
    }
    critPoints.push_back(layerWidth);
    critPoints.push_back(0.0);

    // Removing points lying out the range [0; layerWidth].
    critPoints.erase(std::remove_if(critPoints.begin(), critPoints.end(),
        bind2nd(std::less <double>(), 0)), critPoints.end());

    critPoints.erase(std::remove_if(critPoints.begin(), critPoints.end(),
        bind2nd(std::greater <double>(), layerWidth)), critPoints.end());

    // Excluding repeating points.
    std::sort(critPoints.begin(), critPoints.end());
    critPoints.erase(
        std::unique(critPoints.begin(), critPoints.end()), critPoints.end());

    // Finding the exact dependencies for each subrange of z.
    std::vector<std::pair<double, double> > terms;
    for (int i = 0; i < critPoints.size() - 1; i++) 
    {
        std::pair<double, double> lowerLimit = *(lowerLimits.begin());
        for (std::vector<std::pair<double, double> >::const_iterator lineIter =
                lowerLimits.begin() + 1;
             lineIter < lowerLimits.end();
             ++lineIter)
        {
            if ((lineIter->first * (critPoints[i] + critPoints[i+1]) / 2.0
                   + lineIter->second) <
                (lowerLimit.first * (critPoints[i] + critPoints[i+1]) / 2.0
                   + lowerLimit.second))
            {
                lowerLimit = *(lineIter);
            }

        }
        if ((lowerLimit.first * (critPoints[i] + critPoints[i+1]) / 2.0
               + lowerLimit.second) > 1.0)
        {
            lowerLimit = std::make_pair(0.0, 1.0);
        }

        std::pair<double, double> upperLimit = *(upperLimits.begin());
        for (std::vector<std::pair<double, double> >::const_iterator lineIter =
                upperLimits.begin() + 1;
             lineIter < upperLimits.end();
             ++lineIter)
        {
            if ((lineIter->first * (critPoints[i] + critPoints[i+1]) / 2.0
                   + lineIter->second) >
                (upperLimit.first * (critPoints[i] + critPoints[i+1]) / 2.0
                   + upperLimit.second))
            {
                upperLimit = *(lineIter);
            }
            
        }
        if ((upperLimit.first * (critPoints[i] + critPoints[i+1]) / 2.0
               + upperLimit.second) > 1.0)
        {
            upperLimit = std::make_pair(0.0, 1.0);
        }

        terms.push_back(
            std::make_pair(
                lowerLimit.first - upperLimit.first,
                lowerLimit.second - upperLimit.second));
    }

    // Calculating the integral.
    double integral = 0.0;
    for (int i = 0; i < critPoints.size() - 1; i++) 
    {
        double term = 
            (terms[i].first * (critPoints[i+1] + critPoints[i]) / 2.0
                + terms[i].second)
            * (critPoints[i+1] - critPoints[i]);
        if (term > 0.0)
        {
            integral += term;
        }
    }

    return 2.0 * PI * rodLength * rodLength * integral;
}

double partitionFunctionRodPartiallySubmergedGeneral(
    double segmentLength,
    double slitWidth,
    double layerWidth,
    const std::vector<double> & rodEnergyProfile,
    bool reversed) throw(BioLCCCException)
{
    double partitionFunction = 0.0;
    const double N = rodEnergyProfile.size();
    for (int n1 = 1; n1 < N; n1++)
    {
        for (int n2 = 0; n2 <= N-n1; n2++)
        {
            if (reversed)
            {
                partitionFunction +=
                    ( (n2 == 0 ) ? 1.0 : 0.5 ) *
                    partitionFunctionRodPartiallySubmergedTermGeneral(
                        segmentLength, slitWidth, layerWidth, N, n1, n2)
                    * exp(rodAdsorptionEnergy(rodEnergyProfile, n2, n1));
            }
            else
            {
                partitionFunction +=
                    ( (n2 == 0 ) ? 1.0 : 0.5 ) *
                    partitionFunctionRodPartiallySubmergedTermGeneral(
                        segmentLength, slitWidth, layerWidth, N, n1, n2)
                    * exp(rodAdsorptionEnergy(rodEnergyProfile, n1, n2));
            }
        }
    }
    return partitionFunction;
}

double partitionFunctionRodPartiallySubmergedTermSpecial(
    double segmentLength, double slitWidth, double layerWidth,
    int N, int n1)
{
    double output;
    double rodLength = (N-1) * segmentLength;
    if (layerWidth >= segmentLength * double(n1))
    {
        output = 2.0 * PI * rodLength * rodLength *
                             segmentLength / 2.0;
    }
    else if ((segmentLength * (double)(n1-1) < layerWidth) &&
             (layerWidth < segmentLength * double(n1)))
    {
        output = 2.0 * PI * rodLength * rodLength *
                             (layerWidth
                              - segmentLength * (double)(n1-1) / 2.0
                              - layerWidth * layerWidth / 2.0 / 
                                  (double)n1 / segmentLength);
    }
    else
    {
        output = 2.0 * PI * rodLength * rodLength 
                 * layerWidth * layerWidth / 2.0 / double(n1)
                 / double(n1-1) / segmentLength;
    }
    return output;
}

double partitionFunctionRodPartiallySubmergedSpecial(
    double segmentLength,
    double slitWidth,
    double layerWidth,
    const std::vector<double> & rodEnergyProfile,
    bool reversed) throw(BioLCCCException)
{
    double partitionFunction = 0.0;
    for (unsigned int n1 = 1; n1 < rodEnergyProfile.size(); ++n1)
    {
        if (reversed)
        {
            partitionFunction +=
                partitionFunctionRodPartiallySubmergedTermSpecial(
                    segmentLength, slitWidth, layerWidth, 
                    rodEnergyProfile.size(), n1)
                * exp(rodAdsorptionEnergy(rodEnergyProfile, 0, n1));
        }
        else
        {
            partitionFunction += 
                partitionFunctionRodPartiallySubmergedTermSpecial(
                    segmentLength, slitWidth, layerWidth, 
                    rodEnergyProfile.size(), n1)
                * exp(rodAdsorptionEnergy(rodEnergyProfile, n1, 0));
        }
    }
    return partitionFunction;
}

double partitionFunctionRodFreeVolume(double rodLength,
                                      double slitWidth)
{
    return (4 * PI * slitWidth * rodLength * rodLength);
}

double calculateKdRod(
					  const std::vector<ChemicalGroup> &parsedSequence,
					  const double secondSolventConcentration,
					  const ChemicalBasis &chemBasis,
					  const double columnPoreSize,
					  const double columnRelativeStrength,
					  const double temperature
					  ) throw(BioLCCCException)
{
    if (parsedSequence.size() == 0)
    {
        return 0.0;
    }

    std::vector<double> segmentEnergyProfile = 
        calculateSegmentEnergyProfile(
            calculateMonomerEnergyProfile(
                parsedSequence,
                chemBasis,
                secondSolventConcentration,
                columnRelativeStrength,
                temperature),
            chemBasis.monomerLength(),
            chemBasis.kuhnLength());

    double rodLength = chemBasis.kuhnLength() *
        (segmentEnergyProfile.size() - 1);

    double Kd =
            partitionFunctionRodFreeSlit(
                rodLength,
                columnPoreSize - 2.0 * chemBasis.adsorptionLayerWidth())

            + 2.0 * partitionFunctionRodFreeSlit(
                rodLength,
                chemBasis.adsorptionLayerWidth())
            * exp(rodAdsorptionEnergy(
                segmentEnergyProfile,
                segmentEnergyProfile.size(), 
                0));

    if (!chemBasis.neglectPartiallyDesorbedStates())
    {
        if (chemBasis.specialRodModel())
        {
            Kd += 2.0 * partitionFunctionRodPartiallySubmergedSpecial(
                    chemBasis.kuhnLength(),
                    columnPoreSize,
                    chemBasis.adsorptionLayerWidth(),
                    segmentEnergyProfile,
                    false)

                  + 2.0 * partitionFunctionRodPartiallySubmergedSpecial(
                      chemBasis.kuhnLength(),
                      columnPoreSize,
                      chemBasis.adsorptionLayerWidth(),
                      segmentEnergyProfile,
                      true);
        }
        else
        {
            Kd += 2.0 * partitionFunctionRodPartiallySubmergedGeneral(
                    chemBasis.kuhnLength(),
                    columnPoreSize,
                    chemBasis.adsorptionLayerWidth(),
                    segmentEnergyProfile,
                    false)

                  + 2.0 * partitionFunctionRodPartiallySubmergedGeneral(
                      chemBasis.kuhnLength(),
                      columnPoreSize,
                      chemBasis.adsorptionLayerWidth(),
                      segmentEnergyProfile,
                      true);
        }
    }

    Kd /= partitionFunctionRodFreeVolume(
            rodLength,
            columnPoreSize);

    return Kd;
}

}
