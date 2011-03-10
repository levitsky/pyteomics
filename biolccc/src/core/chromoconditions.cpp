#include "chromoconditions.h"

namespace BioLCCC
{
ChromoConditionsException::ChromoConditionsException(std::string message):
        BioLCCCException(message) {};

ChromoConditions::ChromoConditions(double iColumnLength,
                                   double iColumnDiameter,
                                   double iColumnPoreSize,
                                   Gradient iGradient,
                                   double iSecondSolventConcentrationA,
                                   double iSecondSolventConcentrationB,
                                   double iDelayTime,
                                   double iFlowRate,
                                   double iDV,
                                   double iColumnRelativeStrength,
                                   double iColumnVpToVtot,
                                   double iColumnPorosity,
                                   double iTemperature)
                                   throw(ChromoConditionsException)
{
    setColumnLength(iColumnLength);
    setColumnDiameter(iColumnDiameter);
    setColumnPoreSize(iColumnPoreSize);
    setGradient(iGradient);
    setColumnVpToVtot(iColumnVpToVtot);
    setColumnPorosity(iColumnPorosity);
    setTemperature(iTemperature);
    setColumnRelativeStrength(iColumnRelativeStrength);
    setFlowRate(iFlowRate);
    setDV(iDV);
    setDelayTime(iDelayTime);
    setSecondSolventConcentrationA(iSecondSolventConcentrationA);
    setSecondSolventConcentrationB(iSecondSolventConcentrationB);
}

double ChromoConditions::columnLength() const
{
    return mColumnLength;
}

void ChromoConditions::setColumnLength(double newColumnLength)
    throw(ChromoConditionsException)
{
    if (newColumnLength < 0.0)
    {
        throw(ChromoConditionsException("The new column length is negative."));
    }
    mColumnLength = newColumnLength;
}

double ChromoConditions::columnDiameter() const
{
    return mColumnDiameter;
}

void ChromoConditions::setColumnDiameter(double newColumnDiameter)
    throw(ChromoConditionsException)
{
    if (newColumnDiameter < 0.0)
    {
        throw(ChromoConditionsException(
            "The new column diameter is negative."));
    }
    mColumnDiameter = newColumnDiameter;
}

double ChromoConditions::columnPoreSize() const
{
    return mColumnPoreSize;
}

void ChromoConditions::setColumnPoreSize(double newColumnPoreSize)
    throw(ChromoConditionsException)
{
    if (newColumnPoreSize < 0.0)
    {
        throw(ChromoConditionsException(
            "The new column pore size is negative."));
    }
    mColumnPoreSize = newColumnPoreSize;
}

double ChromoConditions::columnVpToVtot() const
{
    return mColumnVpToVtot;
}

void ChromoConditions::setColumnVpToVtot(double newColumnVpToVtot)
    throw(ChromoConditionsException)
{
    if (newColumnVpToVtot < 0.0)
    {
        throw(ChromoConditionsException(
            "The new column VpToVtot is negative."));
    }
    if (newColumnVpToVtot > 1.0)
    {
        throw(ChromoConditionsException(
            "The new column VpToVtot is greater than 1.0."));
    }
    mColumnVpToVtot = newColumnVpToVtot;
}

double ChromoConditions::columnPorosity() const
{
    return mColumnPorosity;
}

void ChromoConditions::setColumnPorosity(double newColumnPorosity)
    throw(ChromoConditionsException)
{
    if (newColumnPorosity < 0.0)
    {
        throw(ChromoConditionsException(
            "The new column porosity is negative."));
    }
    if (newColumnPorosity > 1.0)
    {
        throw(ChromoConditionsException(
            "The new column porosity is greater than 1.0"));
    }
    mColumnPorosity = newColumnPorosity;
}

double ChromoConditions::temperature() const
{
    return mTemperature;
}

void ChromoConditions::setTemperature(double newTemperature)
    throw(ChromoConditionsException)
{
    if (newTemperature < 0.0)
    {
        throw(ChromoConditionsException("The new temperature is negative."));
    }
    mTemperature = newTemperature;
}

double ChromoConditions::columnRelativeStrength() const
{
    return mColumnRelativeStrength;
}

void ChromoConditions::setColumnRelativeStrength(
    double newColumnRelativeStrength)
{
    mColumnRelativeStrength = newColumnRelativeStrength;
}

double ChromoConditions::flowRate() const
{
    return mFlowRate;
}

void ChromoConditions::setFlowRate(double newFlowRate)
    throw(ChromoConditionsException)
{
    if (newFlowRate < 0.0)
    {
        throw(ChromoConditionsException("The new flow rate is negative."));
    }
    mFlowRate = newFlowRate;
}

double ChromoConditions::dV() const
{
    return mDV;
}

void ChromoConditions::setDV(double newDV)
    throw(ChromoConditionsException)
{
    if (newDV < 0.0)
    {
        throw(ChromoConditionsException("The new dV is negative."));
    }
    mDV = newDV;
}

double ChromoConditions::delayTime() const
{
    return mDelayTime;
}

void ChromoConditions::setDelayTime(double newDelayTime)
{
    mDelayTime = newDelayTime;
}

double ChromoConditions::secondSolventConcentrationA() const
{
    return mSecondSolventConcentrationA;
}

void ChromoConditions::setSecondSolventConcentrationA(
    double newSecondSolventConcentrationA)
    throw(ChromoConditionsException)
{
    if (newSecondSolventConcentrationA < 0.0)
    {
        throw(ChromoConditionsException(
            "The new concentration of second solvent in the component A "
            "is negative."));
    }
    if (newSecondSolventConcentrationA > 100.0)
    {
        throw(ChromoConditionsException(
            "The new concentration of second solvent in the component A "
            "is greater than 100%."));
    }
    mSecondSolventConcentrationA = newSecondSolventConcentrationA;
}

double ChromoConditions::secondSolventConcentrationB() const
{
    return mSecondSolventConcentrationB;
}

void ChromoConditions::setSecondSolventConcentrationB(
    double newSecondSolventConcentrationB)
    throw(ChromoConditionsException)
{
    if (newSecondSolventConcentrationB < 0.0)
    {
        throw(ChromoConditionsException(
            "The new concentration of second solvent in the component A "
            "is negative."));
    }
    if (newSecondSolventConcentrationB > 100.0)
    {
        throw(ChromoConditionsException(
            "The new concentration of second solvent in the component A "
            "is greater than 100%."));
    }
    mSecondSolventConcentrationB = newSecondSolventConcentrationB;
}

Gradient ChromoConditions::gradient() const
{
    return mGradient;
}

void ChromoConditions::setGradient(Gradient newGradient)
    throw(ChromoConditionsException)
{
    if (newGradient.size() < 2)
    {
        throw ChromoConditionsException(
            "The gradient must contain at least two points.");
    }
    mGradient = newGradient;
}
}

