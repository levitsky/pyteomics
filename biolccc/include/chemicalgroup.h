#ifndef AMINOACID_H
#define AMINOACID_H

#include <string>

namespace BioLCCC
{

//! A ChemicalGroup instance contains the properties of a group of atoms.
/*!
    An instance of ChemicalGroup contains the physical properties of a group of
    atoms inside a protein molecule.

    This group can be an amino acid residue or a terminal group, depending on 
    its label.

    Please, see the "Peptide sequence notation" for the following help on
    labels.
*/
class ChemicalGroup
{

public:
    //! Constructs a chemical group with the given parameters.
    /*!
        \param name The full name of the chemical group.
        \param label The label of the chemical group used in sequence notation.
        \param bindEnergy The energy of binding to the surface of a solid phase
        measured in kT.
        \param averageMass The average mass of the group in Da.
        \param monoisotopicMass The monoisotopic mass of the group in Da.
        \param bindArea The area of the contact between the group and the surface.
     */
    ChemicalGroup(std::string name = "",
                  std::string label = "",
                  double bindEnergy = 0.0,
                  double averageMass = 0.0,
                  double monoisotopicMass = 0.0,
                  double bindArea = 1.0
                 );

    //! Returns the full name of the chemical group.
    /*!
        Examples:
        - Methionine
        - Phosphorylated threonine
        - C-Terminal amidation
    */
    std::string name() const;

    //! Returns the brief code of the group used in peptide sequence notation.
    /*!
        Examples:
        - M
        - pT
        - Ac-
    */
    std::string label() const;

    //! Returns the average mass of the chemical group.
    /*!
        The average mass of an amino acid is measured for R-CH(NH-)-CO-
        structure WITHOUT terminal H- and -OH (equals to the average mass
        of a whole amino acid molecule minus 18.01528).
    */
    double averageMass() const;

    //! Returns the monoisotopic mass of the chemical group.
    /*!
        The monoisotopic mass of an amino acid is measured for R-CH(NH-)-CO-
        structure WITHOUT terminal H- and -OH (equals to the monoisotopic
        mass of a whole amino acid molecule minus 18.010565).
    */
    double monoisotopicMass() const;

    //! Returns the energy of binding to the surface of a solid phase.
    /*!
        The bind energy of water is zero and unit is kT.
        The bind energy of a terminal group is added to the binding group of
        the corresponding terminal amino acid.
    */
    double bindEnergy() const;

    //! Returns the area of the contact with the surface of a solid phase.
    /*!
        The unit of the bind area is the area of ACN contact.
        The bind area of a terminal group is added to the area of
        the corresponding terminal amino acid.
    */
    double bindArea() const;


    //! Shows whether the group is N-Terminal.
    bool isNTerminal() const;

    //! Shows whether the group is C-Terminal.
    bool isCTerminal() const;

    //! Shows whether the group is an amino acid.
    bool isAminoAcid() const;

    //! Sets the full name of the chemical group.
    void setName(std::string newName);

    ////! Sets the brief code of the group used in peptide sequence notation.
    //void setLabel(std::string newLabel);

    //! Sets the binding energy value.
    void setBindEnergy(double newBindEnergy);

    //! Sets the binding energy value.
    void setBindArea(double newBindArea);

    //! Sets the average mass of the chemical group.
    void setAverageMass(double newAverageMass);

    //! Sets the monoisotopic mass of the chemical group.
    void setMonoisotopicMass(double newMonoisotopicMass);

private:
    std::string mName;
    std::string mLabel;
    double mBindEnergy;
    double mBindArea;
    double mAverageMass;
    double mMonoisotopicMass;
};

}

#endif
