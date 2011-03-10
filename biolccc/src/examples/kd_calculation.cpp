#include <iostream>
#include <string>
#include <biolccc.h>

int main() 
{
    std::string peptide("Ac-PEPTIDE-NH2");

    double kd = BioLCCC::calculateKd(
        peptide, // the peptide sequence 
        15.0,    // the concentration of the second solvent, %
        BioLCCC::rpAcnFaRod, // the chemical basis
        100.0);  // the size of the pores, angstroms

    std::cout << "The coefficient of distribution of " << peptide 
              << " is " << kd << std::endl;
    return 0;
}

