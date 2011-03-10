#ifndef GRADIENTPOINT_H
#define GRADIENTPOINT_H

#include "biolcccexception.h"

namespace BioLCCC
{

//! This exception is raised when something goes wrong with a GradientPoint.
class GradientPointException : public BioLCCCException
{
public:
    //! Constructs an instance of GradientPointException with the given message.
    GradientPointException(std::string message);
};

//! An instance of GradientPoint keeps the properties of a point of a gradient.
/*!
    Briefly, an instance of GradientPoint is a pair of values. The first one
    is time measured in minutes from the start of the gradient. The second is 
    the concentration of component B in binary solvent.  
 */

class GradientPoint
{
public:
    //! Default constructor.
    /*!
        Constructs a point of a gradient with at time = \a iTime and 
        with the concentration of component B of binary solution = \a
        iConcentrationB.

        Valid values are \a iTime >= 0.0 and 0.0 <= \a iConcentrationB <= 100.0
    */
    GradientPoint(double iTime = 0.0,
                  double iConcentrationB = 0.0)
                  throw (GradientPointException);

    //! Returns the time of the point in minutes.
    double time() const;

    //! Returns the concentration of component B at the point in percents.
    double concentrationB() const;

    //! Sets the time of the point in minutes.
    void setTime(double newTime)
        throw (GradientPointException);

    //! Sets the concentration of component B at the point in percents.
    void setConcentrationB(double newConcentrationB)
        throw (GradientPointException);

private:
    double mTime;
    double mConcentrationB;
};

}

#endif

