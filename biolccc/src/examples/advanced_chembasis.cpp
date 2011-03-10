#include <iostream>
#include <string>
#include <biolccc.h>

int main() 
{
    // Deriving a new ChemicalBasis instance from a predefined one.
    BioLCCC::ChemicalBasis
        myChemicalBasis(BioLCCC::RP_ACN_FA_ROD);

    // Changing the bind energy of a chemical group.
    myChemicalBasis.chemicalGroups()["E"].setBindEnergy(0.0);
    myChemicalBasis.chemicalGroups()["-NH2"].setBindEnergy(0.0);

    std::cout << "The bind energy of E is "
        << myChemicalBasis.chemicalGroups()["E"].bindEnergy()
        << std::endl;
    std::cout << "The bind energy of -NH2 is "
        << myChemicalBasis.chemicalGroups()["-NH2"].bindEnergy()
        << std::endl;

    // Adding a new chemical group. The energy is not valid.
    myChemicalBasis.addChemicalGroup(
        BioLCCC::ChemicalGroup(
            "Hydroxyproline",      // full name
            "hoP",                 // label
            0.40,                  // hypothetical bind energy
            97.1167+15.9994,       // average mass
            97.05276+15.9994915)); // monoisotopic mass

    // Setting a new type of model. Without a massive recalibration
    // it will ruin the accuracy of prediction.
    myChemicalBasis.setPolymerModel(BioLCCC::CHAIN);

    std::string peptide("Ac-PEhoPTIDE-NH2");
    double RT = BioLCCC::calculateRT(peptide,
        myChemicalBasis,
        BioLCCC::standardChromoConditions);

    double monoisotopicMass = BioLCCC::calculateMonoisotopicMass(
        peptide, myChemicalBasis);

    std::cout << "The retention time of " 
              << peptide << " is " << RT << std::endl;
    std::cout << "The monoisotopic mass of " << peptide << " is " 
              << monoisotopicMass << " Da" << std::endl;

    return 0;
}
