#include <iostream>
#include <string>
#include <vector>
#include <biolccc.h>

int main()
{
    std::string peptide("aaaa"); 

    try 
    {
        std::vector<BioLCCC::ChemicalGroup> parsedSequence =
            BioLCCC::parseSequence(peptide, BioLCCC::rpAcnFaRod);
        std::cout << "No exception has been caught." << std::endl;
    }
    catch (BioLCCC::ParsingException & e)
    {
        std::cout << "Parsing exception has been caught." << std::endl;
    }
    catch (BioLCCC::BioLCCCException & e)
    {
        std::cout << "BioLCCC exception has been caught." << std::endl;
    }
    catch (std::exception & e)
    {
        std::cout << "std::exception has been caught." << std::endl;
    }
    catch (...)
    {
        std::cout << "An exception of other type has been caught." << std::endl;
    }

    return 0;
}
