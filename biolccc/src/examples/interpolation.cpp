#include <iostream>
#include <string>
#include <biolccc.h>

int main()
{
    std::string peptide("QWERTYIPASDFGHKLCVNM");

    double RT_old = BioLCCC::calculateRT(peptide,
        BioLCCC::rpAcnTfaChain,
        BioLCCC::standardChromoConditions);
    double RT_fast = BioLCCC::calculateRT(peptide,
        BioLCCC::rpAcnTfaChain,
        BioLCCC::standardChromoConditions, 21);
    std::cout << "For peptide " << peptide
              << " the old algorithm predicts RT=" << RT_old
              << " while the new fast algorithm gives RT=" << RT_fast
              << std::endl;
    return 0;
}
