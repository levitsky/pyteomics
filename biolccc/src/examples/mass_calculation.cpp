#include <iostream>
#include <string>
#include <biolccc.h>

int main()
{
    std::string peptide("Ac-PEPTIDE-NH2");

    double averageMass = BioLCCC::calculateAverageMass(
        peptide, BioLCCC::rpAcnFaRod);
    double monoisotopicMass = BioLCCC::calculateMonoisotopicMass(
        peptide, BioLCCC::rpAcnFaRod);

    std::cout << "Average mass of " << peptide << " is " 
        << averageMass << " Da" << std::endl;
    std::cout << "Monoisotopic mass of " << peptide << " is " 
        << monoisotopicMass << " Da" << std::endl;

    return 0;
}
