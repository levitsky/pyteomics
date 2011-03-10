#include <iostream>
#include <string>
#include <biolccc.h>

int main() 
{
    std::string peptide("Ac-PEPTIDE-NH2");
    double RT = BioLCCC::calculateRT(peptide,
        BioLCCC::rpAcnFaRod,
        BioLCCC::standardChromoConditions);
    std::cout << "The retention time of " 
            << peptide << " is " << RT << std::endl;
    return 0;
}
