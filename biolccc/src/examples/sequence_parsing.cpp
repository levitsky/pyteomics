#include <iostream>
#include <string>
#include <vector>
#include <biolccc.h>

int main()
{
    std::string peptide("PEPTIDE"); 

    std::vector<BioLCCC::ChemicalGroup> parsedSequence =
        BioLCCC::parseSequence(peptide, BioLCCC::rpAcnFaRod);

    for (unsigned int i = 0;
         i < parsedSequence.size();
         ++i)
    {
        std::cout << parsedSequence[i].name() << std::endl;
    }

    return 0;
}
