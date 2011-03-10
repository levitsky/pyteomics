#include <iostream>
#include <string>
#include <biolccc.h>

int main() 
{
    BioLCCC::ChromoConditions myChromoConditions;

    // The column length in mm.
    myChromoConditions.setColumnLength(100.0);

    // The internal column diameter in mm.
    myChromoConditions.setColumnDiameter(0.1);

    // The average pore size in A.
    myChromoConditions.setColumnPoreSize(300.0);

    // The concentration of the eluting solvent (ACN for the reversed
    // phase) in component A in %.
    myChromoConditions.setSecondSolventConcentrationA(5.0);

    // The concentration of the eluting solvent (ACN for the reversed
    // phase) in component B in %.
    myChromoConditions.setSecondSolventConcentrationB(80.0);

    // The shape of the gradient. The example is a linear gradient
    // from 0% to 90% of component B over 60 minutes.
    myChromoConditions.setGradient(
        BioLCCC::Gradient(0.0, 90.0, 60.0));

    // The flow rate in ml/min. 
    myChromoConditions.setFlowRate(0.0005);

    std::string peptide("Ac-PEPTIDE-NH2");
    double RT = BioLCCC::calculateRT(peptide,
        BioLCCC::rpAcnFaRod,
        myChromoConditions);
    std::cout << "The retention time of " 
        << peptide << " is " << RT << std::endl;
    return 0;
}
