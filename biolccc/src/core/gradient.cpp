#include "gradient.h"

namespace BioLCCC
{

GradientException::GradientException(std::string message):
        BioLCCCException(message) {};

Gradient::Gradient()
{
}

Gradient::Gradient(double initialConcentrationB,
                   double finalConcentrationB,
                   double time)
                   throw(GradientException, GradientPointException)
{
    this->addPoint(GradientPoint(0.0, initialConcentrationB));
    this->addPoint(GradientPoint(time, finalConcentrationB));
}

Gradient Gradient::addPoint(GradientPoint iPoint)
    throw(GradientException, GradientPointException)
{
    // The gradient should start at 0.0 min.
    if ((this->size() == 0) && (iPoint.time() != 0.0))
    {
        throw GradientException("The gradient doesn't start at 0.0 min.");
    }

    // Each new point should be later than the previous one.
    if ((this->size() > 0) && (iPoint.time() < this->back().time()))
    {
        throw GradientException(
            "The time of the last point is less than the time "
            "of the previous one");
    }

    this->push_back(iPoint);
    return (*this);
}

Gradient Gradient::addPoint(double iTime, double iConcentrationB)
    throw(GradientException, GradientPointException)
{
    this->addPoint(GradientPoint(iTime, iConcentrationB));
    return (*this);
}
}

