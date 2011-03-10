#ifndef BIOLCCCEXCEPTION_H
#define BIOLCCCEXCEPTION_H

#include <string>

namespace BioLCCC
{

//! Base class for all BioLCCC exceptions. Can be used by itself.
class BioLCCCException : public std::exception
{
public:
    //! A default constructor.
    BioLCCCException(std::string message);

    //! A default destructor.
    ~BioLCCCException() throw();

    //! Returns a message.
    virtual const char* what() const throw();
private:
    std::string mMessage;
};
}

#endif
