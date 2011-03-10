#include <iostream>
#include <string>
#include <biolccc.h>

int main()
{
    std::string peptide("Ac-PEPTIDE-NH2");
    BioLCCC::ChromoConditions myChromoConditions;
    BioLCCC::Gradient myGradient;
    myGradient.addPoint(0.0, 5.0);
    myGradient.addPoint(20.0, 5.0);
    myGradient.addPoint(60.0, 45.0);
    myGradient.addPoint(65.0, 100.0);
    myGradient.addPoint(85.0, 100.0);
    myChromoConditions.setGradient(myGradient);

    double RT = BioLCCC::calculateRT(peptide,
        BioLCCC::rpAcnFaRod,
        myChromoConditions);
    std::cout << "The retention time of " 
              << peptide << " in the custom gradient is " 
              << RT << std::endl;
    return 0;
}
