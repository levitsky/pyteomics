#include "biolcccexception.h"
#include <string.h>

namespace BioLCCC
{

BioLCCCException::BioLCCCException(std::string message)
{
    mMessage = message;
};

BioLCCCException::~BioLCCCException() throw() {};

const char* BioLCCCException::what() const throw()
{
    return mMessage.c_str();
};
}

