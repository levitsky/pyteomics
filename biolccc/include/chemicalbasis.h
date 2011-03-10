#ifndef CHEMICALBASIS_H
#define CHEMICALBASIS_H

#include <map>
#include <vector>
#include "biolcccexception.h"
#include "chemicalgroup.h"

namespace BioLCCC
{

//! This exception is raised when something goes wrong with a ChemicalBasis.
class ChemicalBasisException : public BioLCCCException
{
public:
    //! Constructs an instance of ChemicalBasisException with the given message.
    ChemicalBasisException(std::string message);
};

//! The model of polymer being used in calculations.
/*!
    There are different representations of a polymer model, each suitable for
    different substances and described by different set of equations. At this
    step, the BioLCCC theory describes a polymer molecule as a free-joint chain
    or a rigid rod.
 */

enum PolymerModel
{
    /*! The model in which a polymer molecule is assumed to be absolutely rigid.
        This model works better for short molecules, e.g. peptides. The
        advantage of this model is that it uses explicit expressions for Kd.
        The equations are valid only for molecules, which are shorter than the
        size of a pore. More precisely:
        Length < PoreSize - 2 * AdsorptionLayerWidth 
        */
    ROD,

    /*! The free-joint chain model of a polymer. This model relies heavily on
        the matrix equations.
        */
    CHAIN
};

//! This enum describes the predefined sets of physicochemical constants.
/*!
    The BioLCCC library contains several predifined sets of physicochemical
    constants. Please note that usually changing only one parameter in a whole
    set of constants deteriorate the quality of prediction.
 */
enum PredefinedChemicalBasis
{
    //! Reversed phase, ACN, trifluoracetic acid, CHAIN model.
    /*! A ChemicalBasis calibrated for reversed phase, ACN as a second solvent,
        0.1% TFA and CHAIN type of BioLCCC model. The data was 
        obtained in Guo et al, Journal of Chromatography, 359 (1986) 449-517.
        */
    RP_ACN_TFA_CHAIN, 
    //! Reversed phase, ACN, formic acid, ROD model.
    /*! A ChemicalBasis calibrated for reversed phase, ACN as a second solvent,
        0.1% FA and ROD type of BioLCCC model. The data was obtained
        in the joint research of Harvard University and Institute for Energy 
        Problems for Chemical Physics, Russian Academy of Science.
        */
    RP_ACN_FA_ROD 
};

//! An instance of ChemicalBasis contains a set of BioLCCC constants.
/*!
    An instance of ChemicalBasis manages all the physicochemical constants,
    which are used in the calculations. Currently it contains:
        - the list of amino acids and peptide terminal groups;
        - which terminal groups are set by default (cannon be changed);
        - the chemical properties of solvents: densities, molar mass and
          adsorption energies (adsorption energy of the first solvent always
          equals zero);
        - the type of BioLCCC model being used in calculations and
          approximations used in the equations;
        - peptide geometry: the length of amino acid and the Kuhn length;
        - the width of the adsorbing layer on the walls.
       
    Note that the constants are highly interconnected. That leads to the fact,
    that a change in a single constant, like the width of adsorbing layer or 
    type of BioLCCC model, deteriorates the accuracy of RT prediction. 
 */
class ChemicalBasis
{
public:
    //! Constructs an empty ChemicalBasis instance.
    ChemicalBasis();

    //! Constructs a ChemicalBasis instance with a predefined set of constants.
    ChemicalBasis(PredefinedChemicalBasis predefinedChemicalBasisId);

    //! Returns the map of all chemical groups. Non-constant version.
    /*!
        A chemical group can be retrieved from the map by its label.
     */
    std::map<std::string, ChemicalGroup> & chemicalGroups();

    //! Returns the map of all chemical groups. Constant version.
    /*!
        A chemical group can be retrieved from the map by its label.
     */
    const std::map<std::string, ChemicalGroup> & chemicalGroups() const;

    //! Returns the default N-terminal group.
    const ChemicalGroup & defaultNTerminus() const
        throw(ChemicalBasisException);

    //! Returns the default C-terminal group.
    const ChemicalGroup & defaultCTerminus() const
        throw(ChemicalBasisException);

    //! Adds \a newChemicalGroup to the ChemicalBasis.
    /*!
        If an instance of ChemicalBasis already contains a chemical group with
        the same label than it is overwritten.
     */
    void addChemicalGroup(ChemicalGroup newChemicalGroup);

    //! Removes a chemical group with the given \a label;
    /*!
        Throws ChemicalBasisException if a chemical group with the given label
        is not found.
    */
    void removeChemicalGroup(std::string label)
        throw(ChemicalBasisException);

    //! Removes all chemical groups in a ChemicalBasis.
    void clearChemicalGroups();

    ////!Sets \a newBindEnergy as the binding energy of chemical group \a label.
    ///*!
    //    Throws ChemicalBasisException if the chemical group is not found.
    //    \param label The label of the chemical group to be modified.
    //    \param newBindEnergy The new value of the bind energy.
    //*/
    //void setChemicalGroupBindEnergy(std::string label, double newBindEnergy);

    //! Returns the bind energy of the second solvent. 
    /*! 
        Note that the bind energy of water is zero and the unit is kT.
    */
    double secondSolventBindEnergy() const;

    //! Sets \a newEnergy as the bind energy of the second solvent. 
    /*!
        Note that the bind energy of water is zero and the unit is kT.
    */
    void setSecondSolventBindEnergy(double newEnergy);

    //! Sets the model of a polymer (CHAIN or ROD).
    void setPolymerModel(PolymerModel newModel);

    //! Returns the model of a polymer (CHAIN or ROD).
    const PolymerModel polymerModel() const;

    //! Returns the length of a single monomer in a polymer chain in angstroms.
    /*!
        Due to the complex geometry of peptide molecule, this length is defined
        only approximately. The definition is the average length of an amino
        acid residue along backbone. In other terms, it is the length of 
        a backbone divided by the number of amino acid residues.
     */
    double monomerLength() const;

    //! Sets the length of a single monomer in a polymer chain in angstroms.
    /*!
        Due to the complex geometry of peptide molecule, this length is defined
        only approximately. The definition is the average length of an amino
        acid residue along backbone. In other terms, it is the length of 
        a backbone divided by the number of amino acid residues.
     */
    void setMonomerLength(double newMonomerLength)
        throw(ChemicalBasisException);

    //! Returns the Kuhn length of a polymer molecule in angstroms.
    /*!
        A polymer molecule can be modelled as a chain of equal-sized rigid rods,
        freely joined with each other. In this case, the rods would be called 
        Kuhn segments, and the length of a segment would be the Kuhn length.
        The effective adsorption energy of a Kuhn segment equals to the total
        adsorption energy of all monomers that contains in this segment. If only
        a part of monomer contains in a segment than its energy is taken
        proportionally.

        There are no joints in the ROD model, the whole molecule is assumed to
        be shorter than a single Kuhn segment. However, in ROD model kuhnLength
        is still used to calculate the energy profile of a rod. The whole rod is
        divided into segments of kuhnLength and each segment transforms into an
        adsorbing bead. The effective energy of adsorption equals to the total
        effective energy of a segment, with the same expression as in the CHAIN
        model.
    */
    double kuhnLength() const;

    //! Sets the Kuhn length of a molecule in angstroms.
    /*!
        A polymer molecule can be modelled as a chain of equal-sized rigid rods,
        freely joined with each other. In this case, the rods would be called 
        Kuhn segments, and the length of a segment would be the Kuhn length.
        The effective adsorption energy of a Kuhn segment equals to the total
        adsorption energy of all monomers that contains in this segment. If only
        a part of monomer contains in a segment than its energy is taken
        proportionally.

        There are no joints in the ROD model, the whole molecule is assumed to
        be shorter than a single Kuhn segment. However, in ROD model kuhnLength
        is still used to calculate the energy profile of a rod. The whole rod is
        divided into segments of kuhnLength and each segment transforms into an
        adsorbing bead. The effective energy of adsorption equals to the total
        effective energy of a segment, with the same expression as in the CHAIN
        model.
    */
    void setKuhnLength(double newKuhnLength)
        throw(ChemicalBasisException);

    //! Returns the width of a solid phase adsorption layer in ROD model.
    /*!
        The width of a solid phase adsorption layer can be defined as a
        characteristic distance of interaction between an amino acid residue and
        the surface of a solid phase. 

        This value is used only in the ROD model.
    */
    double adsorptionLayerWidth() const;

    //! Sets the width of a solid phase adsorption layer in ROD model.
    /*!
        The width of a solid phase adsorption layer can be defined as a
        characteristic distance of interaction between an amino acid residue and
        the surface of a solid phase. 

        This value is used only in the ROD model.
    */
    void setAdsorptionLayerWidth(double newAdsorptionLayerWidth)
        throw(ChemicalBasisException);

    //! Returns the absorption factors of the near-wall layers in CHAIN model.
    /*! 
        The standard CHAIN BioLCCC model assumes that adsorption occurs only in
        one layer, which is closest to the wall. This assumption can be
        generalized to the case when several near-wall layers adsorb segments 
        of a polymer chain. adsorptionLayerFactors() vector contains the 
        relative adsorption
        strengths of the near-wall layers. This adsorption strength have the 
        same
        meaning as the relative adsorption strength of a column and is
        multiplied by it. The first element of the vector corresponds to the
        layer closest to the wall, second to the next and so on.

        This value is used only in the CHAIN model.
     */
    const std::vector<double> & adsorptionLayerFactors() const;

    //! Sets the absorption factors of the near-wall layers in CHAIN model.
    /*! 
        The standard CHAIN BioLCCC model assumes that adsorption occurs only in
        one layer, which is closest to the wall. This assumption can be
        generalized to the case when several near-wall layers adsorb segments 
        of a polymer chain. adsorptionLayerFactors() vector contains the 
        relative adsorption
        strengths of the near-wall layers. This adsorption strength have the 
        same
        meaning as the relative adsorption strength of a column and is
        multiplied by it. The first element of the vector corresponds to the
        layer closest to the wall, second to the next and so on.

        This value is used only in the CHAIN model.
     */
    void setAdsorptionLayerFactors(
        std::vector<double> newAdsorptionLayerFactors);

    //! Returns true if the energy of binary solvent is linearly fitted.
    /*!
        If the value of snyderApproximation is true then the energy of binary
        solvent is expressed by:
        E_{ab} = secondSolventBindEnergy * Nb
     */
    bool snyderApproximation() const;

    //! Enables the linear approximation of the energy of binary solvent.
    /*!
        If the value of snyderApproximation is true then the energy of binary
        solvent is expressed by:
        E_{ab} = secondSolventBindEnergy * Nb
     */
    void setSnyderApproximation(bool flag);

    //! Returns true if equations for the special case of rod model are used.
    bool specialRodModel() const;
    
    //! Enables the usage of the special equation for the rod model.
    /*!
      The special equation of the rod model do not account for the bridge
      conformations of a polymer, in which both its ends are adsorbed to
      the opposite walls of a pore. Therefore, it is valid only when
      polymerLength < slitWidth - 2 * layerWidth.

      If specialRodModel is set to false, than the general equation is
      used. It is valid for polymers of any length, but requires more
      computational resources.
     */
	void setSpecialRodModel(bool flag);
	
	//! Returns true if partially desorbed states are neglected.
	bool neglectPartiallyDesorbedStates() const;

	//! Excludes partially desorbed states from calculation. FOR MODEL STUDY ONLY.
    /*!  
      In BioLCCC model, the distribution coefficient of a polymer
      is calculated by integration over all its possible
      conformations. These could be differentiated into three distinct
      groups: totally adsorbed, totally desorbed and partially
      desorbed states. The latter are of particular interest, because
      these conformations give rise to sequence specificity of BioLCCC
      model. Using setNeglectPartiallyDesorbedStates function, these
      conformations can be excluded from calculation. This option is
      intended for study of BioLCCC properties, and should not be used
      in routine applications.
     */
	void setNeglectPartiallyDesorbedStates(bool flag);

    //! Returns the density of the first solvent in kg/m^3.
    double firstSolventDensity() const;

    //! Sets the density of the first solvent in kg/m^3.
    void setFirstSolventDensity(double newFirstSolventDensity)
        throw(ChemicalBasisException);

    //! Returns the density of the second solvent in kg/m^3.
    double secondSolventDensity() const;

    //! Sets the density of the second solvent in kg/m^3.
    void setSecondSolventDensity(double newSecondSolventDensity)
        throw(ChemicalBasisException);

    //! Returns the molecular mass of the first solvent in g/mol.
    double firstSolventAverageMass() const;

    //! Sets the molecular mass of the first solvent in g/mol.
    void setFirstSolventAverageMass(double newFirstSolventAverageMass)
        throw(ChemicalBasisException);

    //! Returns the molecular mass of the second solvent in g/mol.
    double secondSolventAverageMass() const;

    //! Sets the molecular mass of the second solvent in g/mol.
    void setSecondSolventAverageMass(double newSecondSolventAverageMass)
        throw(ChemicalBasisException);

    //! Sets one of predefined chemical basis.
    ChemicalBasis setPredefinedChemicalBasis(
        PredefinedChemicalBasis predefinedChemicalBasisId);

private:
    std::map<std::string,ChemicalGroup> mChemicalGroups;
    double mSecondSolventBindEnergy;
    double mMonomerLength;
    double mKuhnLength;
    double mAdsorptionLayerWidth;
    std::vector<double> mAdsorptionLayerFactors;
    PolymerModel mPolymerModel;
    double mFirstSolventDensity;
    double mSecondSolventDensity;
    double mFirstSolventAverageMass;
    double mSecondSolventAverageMass;
    bool mSnyderApproximation;
    bool mSpecialRodModel;
	bool mNeglectPartiallyDesorbedStates;
};

}

#endif
