#include <iostream>
#include <biolccc.h>

int main()
{
    for (std::map<std::string, BioLCCC::ChemicalGroup>::const_iterator
           chemicalGroupIt = BioLCCC::rpAcnFaRod.chemicalGroups().begin();
       chemicalGroupIt != BioLCCC::rpAcnFaRod.chemicalGroups().end();
       ++chemicalGroupIt)
    {
        std::cout << "Name: " << chemicalGroupIt->second.name() 
            << std::endl;
        std::cout << "Label: " << chemicalGroupIt->second.label()
            << std::endl;
        std::cout << "Bind energy: " << 
            chemicalGroupIt->second.bindEnergy() << std::endl;
        std::cout << "Average mass: "
            << chemicalGroupIt->second.averageMass() << std::endl;
        std::cout << "Monoisotopic mass: "
            << chemicalGroupIt->second.monoisotopicMass()
            << std::endl;
        std::cout << std::endl;
    }

    return 0;
}
