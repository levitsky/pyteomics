#ifndef CHROMOCONDITION_H
#define CHROMOCONDITION_H

#include <string>
#include <utility>
#include <vector>

#include "biolcccexception.h"
#include "gradient.h"

namespace BioLCCC
{

//! This exception is raised when something goes wrong with a ChromoConditions.
class ChromoConditionsException : public BioLCCCException
{
public:
    //! Constructs a ChromoConditionsException instance with a given message.
    ChromoConditionsException(std::string message);
};

//! A ChromoConditions instance describes conditions of chromatography.
/*!
    An instance of ChromoConditions manages all the parameters of
    chromatographic equipment. It contains:
        - The geometry of the column.
        - The properties of the adsorbent: average size of the pores, porosity
          (i.e. percentage of volume not filled with the solid phase),
          (volume of pores)/(total volume of column) ratio, relative strength of
          adsorption.
        - Elution parameters: the shape of the gradient, the composition of
          components, flow rate, delay time.
        - Temperature of a column (EXPERIMENTAL).
        - The step of integration.

 */
class ChromoConditions
{
public:

    //! Constructs a custom ChromoConditions object.
    /*!
        Default values are the same as for standardChromoConditions. 
     */
    ChromoConditions(double iColumnLength = 150.0,
                     double iColumnDiameter = 0.075,
                     double iColumnPoreSize = 100.0,
                     Gradient iGradient = Gradient(0.0, 50.0, 60.0),
                     double iSecondSolventConcentrationA = 2.0,
                     double iSecondSolventConcentrationB = 80.0,
                     double iDelayTime = 0.0,
                     double iFlowRate = 0.0003,
                     double iDV = 0.0,
                     double iColumnRelativeStrength = 1.0,
                     double iColumnVpToVtot = 0.5,
                     double iColumnPorosity = 0.9,
                     double iTemperature = 293.0)
                     throw(ChromoConditionsException);

    //! Returns the length of the column in mm.
    /*!
        Note that it is the length of the area filled with an adsorbent.
     */
    double columnLength() const;

    //!  Sets the length of the column in mm.
    /*!
        Note that it is the length of the area filled with an adsorbent.
     */
    void setColumnLength(double newColumnLength)
        throw(ChromoConditionsException);

    //! Returns the internal diameter of the column in mm.
    double columnDiameter() const;

    //! Sets the internal diameter of the column in mm.
    void setColumnDiameter(double newColumnDiameter)
        throw(ChromoConditionsException);

    //! Returns the size of the pores in angstroms.
    double columnPoreSize() const;

    //! Sets the size of the pores in angstroms.
    void setColumnPoreSize(double newColumnPoreSize)
        throw(ChromoConditionsException);

    //! Returns the ratio of the volume of pores to the total column volume.
    double columnVpToVtot() const;

    //! Sets the ratio of the volume of pores to the total column volume.
    void setColumnVpToVtot(double newColumnVpToVtot)
        throw(ChromoConditionsException);

    //! Returns the porosity of a column.
    /*!
        Porosity of a column describes which part of the column in not filled
        with a solid phase. This part is made up by pores and interstitial
        volume.
     */
    double columnPorosity() const;

    //! Sets the porosity of a column.
    /*!
        Porosity of a column describes which part of the column in not filled
        with a solid phase. This part is made up by pores and interstitial
        volume.
     */
    void setColumnPorosity(double newColumnPorosity)
        throw(ChromoConditionsException);

    //! Returns the temperature of the column in kelvin degrees.
    double temperature() const;

    //! Sets the temperature of the column in kelvin degrees.
    void setTemperature(double newTemperature)
        throw(ChromoConditionsException);

    //! Returns the relative strength of the adsorbent.
    /*!
        All the adsorption energies are multiplied by the relative strength of
        the column. Please, check the BioLCCC theory for the further details of
        implementation.
     */
    double columnRelativeStrength() const;

    //! Sets the relative strength of the adsorbent.
    /*!
        All the adsorption energies are multiplied by the relative strength of
        the column. Please, check the BioLCCC theory for the further details of
        implementation.
     */
    void setColumnRelativeStrength(double newColumnRelativeStrength);

    //! Returns the flow rate in ml/min.
    double flowRate() const;

    //! Sets the flow rate in ml/min.
    void setFlowRate(double newFlowRate)
        throw(ChromoConditionsException);

    //! Returns the step of integration over volume in ml.
    /*!
        The main equation of chromatography includes the integration over
        retention volume. This parameter describes the step of this
        intergration.
        The physical interpretation could be a volume of the pump mixer.
        
        Note that if the dV is set to zero, than it is assumed to be equal to
        flowRate*1 min/20.
     */
    double dV() const;

    //! Sets the step of integration over volume in ml.
    /*!
        The main equation of chromatography includes the integration over
        retention volume. This parameter describes the step of this
        intergration.
        The physical interpretation could be a volume of the pump mixer.
        
        Note that if the dV is set to zero, than it is assumed to be equal to
        flowRate*1 min/20.
     */
    void setDV(double newDV)
        throw(ChromoConditionsException);

    //! Returns the delay time.
    /*!
        Delay time is simply added to a calculated retention time.
     */
    double delayTime() const;

    //! Sets the delay time.
    /*!
        Delay time is simply added to a calculated retention time.
     */
    void setDelayTime(double newDelayTime);

    //! Returns the concentration of the second solvent in component A.
    double secondSolventConcentrationA() const;

    //!  Sets the concentration of the second solvent in component A.
    void setSecondSolventConcentrationA(
        double newSecondSolventConcentrationA)
        throw(ChromoConditionsException);

    //! Returns the concentration of the second solvent in component B.
    double secondSolventConcentrationB() const;

    //! Sets the concentration of the second solvent in component B.
    void setSecondSolventConcentrationB(
        double newSecondSolventConcentrationB)
        throw(ChromoConditionsException);

    //! Returns the elution gradient.
    Gradient gradient() const;

    //! Sets the elution gradient.
    void setGradient(Gradient newGradient)
        throw(ChromoConditionsException);

private:
    double mColumnLength;
    double mColumnDiameter;
    double mColumnPoreSize;
    double mColumnVpToVtot;
    double mColumnPorosity;
    double mTemperature;
    double mColumnRelativeStrength;
    double mFlowRate;
    double mDV;
    double mDelayTime;
    Gradient mGradient;
    double mSecondSolventConcentrationA;
    double mSecondSolventConcentrationB;
};

}

#endif
