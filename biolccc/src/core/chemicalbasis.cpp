#include "chemicalbasis.h"

namespace BioLCCC
{
ChemicalBasisException::ChemicalBasisException(std::string message):
        BioLCCCException(message) {};

ChemicalBasis::ChemicalBasis()
{
    setPolymerModel(CHAIN);
    setFirstSolventDensity(1e-10);
    setFirstSolventAverageMass(1e-10);
    setSecondSolventDensity(1e-19);
    setSecondSolventAverageMass(1e-10);
    setSecondSolventBindEnergy(1.0e-10);
    setMonomerLength(1.0e-10);
    setKuhnLength(1.0e-10);
    setAdsorptionLayerWidth(1.0e-10);
    setAdsorptionLayerFactors(std::vector<double>(1, 1.0));
    setSnyderApproximation(false);
    setSpecialRodModel(true);
	setNeglectPartiallyDesorbedStates(false);
}

ChemicalBasis::ChemicalBasis(PredefinedChemicalBasis predefinedChemicalBasisId)
{
    setPredefinedChemicalBasis(predefinedChemicalBasisId);
}

std::map<std::string,ChemicalGroup> & 
    ChemicalBasis::chemicalGroups()
{
    return mChemicalGroups;
}

const std::map<std::string,ChemicalGroup> & 
    ChemicalBasis::chemicalGroups() const
{
    return mChemicalGroups;
}

const ChemicalGroup & ChemicalBasis::defaultNTerminus() const
    throw(ChemicalBasisException)
{
    std::map<std::string, ChemicalGroup>::const_iterator NTerminusIterator =
        mChemicalGroups.find("H-");
    if (NTerminusIterator == mChemicalGroups.end())
    {
        throw ChemicalBasisException("The default H- N-terminus not found.");
    }
    return NTerminusIterator->second;
}

const ChemicalGroup & ChemicalBasis::defaultCTerminus() const
    throw(ChemicalBasisException)
{
    std::map<std::string, ChemicalGroup>::const_iterator CTerminusIterator =
        mChemicalGroups.find("-OH");
    if (CTerminusIterator == mChemicalGroups.end())
    {
        throw ChemicalBasisException("The default -OH C-terminus not found.");
    }
    return CTerminusIterator->second;
}

double ChemicalBasis::secondSolventBindEnergy() const
{
    return mSecondSolventBindEnergy;
}

void ChemicalBasis::setSecondSolventBindEnergy(double newEnergy)
{
    mSecondSolventBindEnergy = newEnergy;
}

double ChemicalBasis::monomerLength() const
{
    return mMonomerLength;
}

void ChemicalBasis::setMonomerLength(double newMonomerLength)
    throw(ChemicalBasisException)
{
    if (newMonomerLength <= 0.0)
    {
        throw ChemicalBasisException(
            "The new length of a monomer is not positive.");
    }

    mMonomerLength = newMonomerLength;
}

double ChemicalBasis::kuhnLength() const
{
    return mKuhnLength;
}

void ChemicalBasis::setKuhnLength(double newKuhnLength)
    throw(ChemicalBasisException)
{
    if (newKuhnLength <= 0)
    {
        throw ChemicalBasisException(
            "The new persistent length is not positive.");
    }
    mKuhnLength = newKuhnLength;
}

double ChemicalBasis::adsorptionLayerWidth() const
{
    return mAdsorptionLayerWidth;
}

void ChemicalBasis::setAdsorptionLayerWidth(double newAdsorptionLayerWidth)
    throw(ChemicalBasisException)
{
    if (newAdsorptionLayerWidth < 0.0)
    {
        throw ChemicalBasisException(
            "The new adsorption layer width is negative.");
    }
    mAdsorptionLayerWidth = newAdsorptionLayerWidth;
}

const std::vector<double> & ChemicalBasis::adsorptionLayerFactors() const
{
    return mAdsorptionLayerFactors;
}

void ChemicalBasis::setAdsorptionLayerFactors(
    std::vector<double> newAdsorptionLayerFactors)
{
    mAdsorptionLayerFactors = newAdsorptionLayerFactors;
}

void ChemicalBasis::addChemicalGroup(ChemicalGroup newChemicalGroup)
{
    mChemicalGroups[newChemicalGroup.label()] = newChemicalGroup;
}

void ChemicalBasis::removeChemicalGroup(std::string label)
    throw(ChemicalBasisException)
{
    if (mChemicalGroups.erase(label) == (std::size_t)0)
    {
        throw ChemicalBasisException(
            "The chemical group " + label + " is not found.");
    }
}

void ChemicalBasis::clearChemicalGroups()
{
    mChemicalGroups.clear();
}

//void ChemicalBasis::setChemicalGroupBindEnergy(std::string label,
//        double newBindEnergy)
//{
//    std::map<std::string,ChemicalGroup>::iterator it =
//        mChemicalGroups.find(label);
//    if (it == mChemicalGroups.end())
//    {
//        throw ChemicalBasisException(
//            "The chemical group " + label + " is not found.");
//    }
//    it->second.setBindEnergy(newBindEnergy);
//}

const PolymerModel ChemicalBasis::polymerModel() const
{
    return mPolymerModel;
}

void ChemicalBasis::setPolymerModel(PolymerModel newModel)
{
    mPolymerModel = newModel;
}

bool ChemicalBasis::snyderApproximation() const 
{
    return mSnyderApproximation;
}

void ChemicalBasis::setSnyderApproximation(bool flag) 
{
    mSnyderApproximation = flag;
}

bool ChemicalBasis::specialRodModel() const 
{
    return mSpecialRodModel;
}

void ChemicalBasis::setSpecialRodModel(bool flag) 
{
    mSpecialRodModel = flag;
}

bool ChemicalBasis::neglectPartiallyDesorbedStates() const
{
	return mNeglectPartiallyDesorbedStates;
}

void ChemicalBasis::setNeglectPartiallyDesorbedStates(bool flag)
{
	mNeglectPartiallyDesorbedStates = flag;
}

double ChemicalBasis::firstSolventDensity() const
{
    return mFirstSolventDensity;
}

void ChemicalBasis::setFirstSolventDensity(double newFirstSolventDensity)
    throw(ChemicalBasisException)
{
    if (newFirstSolventDensity < 0.0)
    {
        throw ChemicalBasisException(
            "The density must have a non-negative value.");
    }
    mFirstSolventDensity = newFirstSolventDensity;
}

double ChemicalBasis::secondSolventDensity() const
{
    return mSecondSolventDensity;
}

void ChemicalBasis::setSecondSolventDensity(double newSecondSolventDensity)
    throw(ChemicalBasisException)
{
    if (newSecondSolventDensity < 0.0)
    {
        throw ChemicalBasisException(
            "The density must have a not-negative value.");
    }
    mSecondSolventDensity = newSecondSolventDensity;
}

double ChemicalBasis::firstSolventAverageMass() const
{
    return mFirstSolventAverageMass;
}

void ChemicalBasis::setFirstSolventAverageMass(
    double newFirstSolventAverageMass)
    throw(ChemicalBasisException)
{
    if (newFirstSolventAverageMass < 0.0)
    {
        throw ChemicalBasisException(
            "The average mass must have a non-negative value.");
    }
    mFirstSolventAverageMass = newFirstSolventAverageMass;
}

double ChemicalBasis::secondSolventAverageMass() const
{
    return mSecondSolventAverageMass;
}

void ChemicalBasis::setSecondSolventAverageMass(
    double newSecondSolventAverageMass)
    throw(ChemicalBasisException)
{
    if (newSecondSolventAverageMass < 0.0)
    {
        throw ChemicalBasisException(
            "The average mass must have a not-negative value.");
    }
    mSecondSolventAverageMass = newSecondSolventAverageMass;
}

ChemicalBasis ChemicalBasis::setPredefinedChemicalBasis(
    PredefinedChemicalBasis predefinedChemicalBasisId)
{
    switch ( predefinedChemicalBasisId ) 
    {
        case RP_ACN_TFA_CHAIN:
        {
            setPolymerModel(CHAIN);
            // Water as the first solvent.
            setFirstSolventDensity(1000.0);
            setFirstSolventAverageMass(18.02);
            //setFirstSolventDensity(5.56);
            //setFirstSolventAverageMass(1.0);
            // Acetonitrile as the second solvent.
            setSecondSolventDensity(786.0);
            setSecondSolventAverageMass(41.05);
            //setSecondSolventDensity(1.91);
            //setSecondSolventAverageMass(1.0);
            setSecondSolventBindEnergy(2.4);
            setAdsorptionLayerWidth(15.0);
            setAdsorptionLayerFactors(std::vector<double>(1, 1.0));
            setMonomerLength(10.0);
            setKuhnLength(10.0);
            setSnyderApproximation(false);
            setSpecialRodModel(true);
			setNeglectPartiallyDesorbedStates(false);

            clearChemicalGroups();
            addChemicalGroup(ChemicalGroup ("Alanine",
                                            "A",
                                            1.1425,
                                            71.0788,
                                            71.03711));
            addChemicalGroup(ChemicalGroup ("Cysteine",
                                            "C",
                                            1.2955,
                                            103.1388,
                                            103.00919));
            addChemicalGroup(ChemicalGroup ("Carboxyamidomethylated cysteine",
                                            "camC",
                                            0.77,
                                            160.1901,
                                            160.03065));
            addChemicalGroup(ChemicalGroup ("Aspartic acid",
                                            "D",
                                            0.7805,
                                            115.0886,
                                            115.02694));
            addChemicalGroup(ChemicalGroup ("Glutamic acid",
                                            "E",
                                            0.9835,
                                            129.1155,
                                            129.04259));
            addChemicalGroup(ChemicalGroup ("Phenylalanine",
                                            "F",
                                            2.3185,
                                            147.1766,
                                            147.06841));
            addChemicalGroup(ChemicalGroup ("Glycine",
                                            "G",
                                            0.6555,
                                            57.0519,
                                            57.02146));
            addChemicalGroup(ChemicalGroup ("Histidine",
                                            "H",
                                            0.3855,
                                            137.1411,
                                            137.05891));
            addChemicalGroup(ChemicalGroup ("Isoleucine",
                                            "I",
                                            2.1555,
                                            113.1594,
                                            113.08406));
            addChemicalGroup(ChemicalGroup ("Lysine",
                                            "K",
                                            0.2655,
                                            128.1741,
                                            128.09496));
            addChemicalGroup(ChemicalGroup ("Leucine",
                                            "L",
                                            2.2975,
                                            113.1594,
                                            113.08406));
            addChemicalGroup(ChemicalGroup ("Methionine",
                                            "M",
                                            1.8215,
                                            131.1926,
                                            131.04049));
            addChemicalGroup(ChemicalGroup ("Oxidated methionine",
                                            "oxM",
                                            1.8215,
                                            131.1926 + 15.9994,
                                            131.04049 + 15.994915));
            addChemicalGroup(ChemicalGroup ("Asparagine",
                                            "N",
                                            0.6135,
                                            114.1038,
                                            114.04293));
            addChemicalGroup(ChemicalGroup ("Proline",
                                            "P",
                                            1.1425,
                                            97.1167,
                                            97.05276));
            addChemicalGroup(ChemicalGroup ("Glutamine",
                                            "Q",
                                            0.7455,
                                            128.1307,
                                            128.05858));
            addChemicalGroup(ChemicalGroup ("Arginine",
                                            "R",
                                            0.5155,
                                            156.1875,
                                            156.10111));
            addChemicalGroup(ChemicalGroup ("Serine",
                                            "S",
                                            0.6975,
                                            87.0782,
                                            87.03203));
            addChemicalGroup(ChemicalGroup ("Phosphorylated serine",
                                            "pS",
                                            0.45,
                                            167.0581,
                                            166.99836));
            addChemicalGroup(ChemicalGroup ("Threonine",
                                            "T",
                                            0.8755,
                                            101.1051,
                                            101.04768));
            addChemicalGroup(ChemicalGroup ("Phosphorylated threonine",
                                            "pT",
                                            0.74,
                                            181.085,
                                            181.01401));
            addChemicalGroup(ChemicalGroup ("Valine",
                                            "V",
                                            1.7505,
                                            99.1326,
                                            99.06841));
            addChemicalGroup(ChemicalGroup ("Tryptophan",
                                            "W",
                                            2.4355,
                                            186.2132,
                                            186.07931));
            addChemicalGroup(ChemicalGroup ("Tyrosine",
                                            "Y",
                                            1.6855,
                                            163.176,
                                            163.06333));
            addChemicalGroup(ChemicalGroup ("Phosphorylated tyrosine",
                                            "pY",
                                            1.32,
                                            243.1559,
                                            243.02966));
            addChemicalGroup(ChemicalGroup("N-terminal hydrogen",
                                           "H-",
                                           -1.69,
                                           1.0079,
                                           1.00782,
                                           0.0));
            addChemicalGroup(ChemicalGroup("N-terminal acetyl",
                                           "Ac-",
                                           0.0,
                                           43.0452,
                                           43.01839,
                                           0.0));
            addChemicalGroup(ChemicalGroup("C-terminal carboxyl group",
                                           "-OH",
                                           -0.03,
                                           17.0073,
                                           17.00274,
                                           0.0));
            addChemicalGroup(ChemicalGroup("C-terminal amide",
                                           "-NH2",
                                           0.0,
                                           16.0226,
                                           16.01872,
                                           0.0));
            break;
        }

        case RP_ACN_FA_ROD:
        {
            setPolymerModel(ROD);
            setAdsorptionLayerWidth(16.0);
            setAdsorptionLayerFactors(std::vector<double>(1, 1.0));
            setMonomerLength(4.0);
            setKuhnLength(4.0);
            // Water as the first solvent.
            setFirstSolventDensity(1000.0);
            setFirstSolventAverageMass(18.02);
            // Acetonitrile as the second solvent.
            setSecondSolventDensity(786.0);
            setSecondSolventAverageMass(41.05);
            setSecondSolventBindEnergy(2.4);
            setSnyderApproximation(false);
            setSpecialRodModel(true);
			setNeglectPartiallyDesorbedStates(false);

            clearChemicalGroups();
            addChemicalGroup(ChemicalGroup ("Alanine",
                                            "A",
                                            0.81,
                                            71.0788,
                                            71.03711));
            addChemicalGroup(ChemicalGroup ("Cysteine",
                                            "C",
                                            0.90,
                                            103.1388,
                                            103.00919));
            addChemicalGroup(ChemicalGroup ("Carboxyamidomethylated cysteine",
                                            "camC",
                                            0.24,
                                            160.1901,
                                            160.03065));
            addChemicalGroup(ChemicalGroup ("Aspartic acid",
                                            "D",
                                            0.60,
                                            115.0886,
                                            115.02694));
            addChemicalGroup(ChemicalGroup ("Glutamic acid",
                                            "E",
                                            0.62,
                                            129.1155,
                                            129.04259));
            addChemicalGroup(ChemicalGroup ("Phenylalanine",
                                            "F",
                                            2.65,
                                            147.1766,
                                            147.06841));
            addChemicalGroup(ChemicalGroup ("Glycine",
                                            "G",
                                            0.47,
                                            57.0519,
                                            57.02146));
            addChemicalGroup(ChemicalGroup ("Histidine",
                                            "H",
                                            -0.77,
                                            137.1411,
                                            137.05891));
            addChemicalGroup(ChemicalGroup ("Isoleucine",
                                            "I",
                                            2.20,
                                            113.1594,
                                            113.08406));
            addChemicalGroup(ChemicalGroup ("Lysine",
                                            "K",
                                            -0.64,
                                            128.1741,
                                            128.09496));
            addChemicalGroup(ChemicalGroup ("Leucine",
                                            "L",
                                            2.38,
                                            113.1594,
                                            113.08406));
            addChemicalGroup(ChemicalGroup ("Methionine",
                                            "M",
                                            1.73,
                                            131.1926,
                                            131.04049));
            addChemicalGroup(ChemicalGroup ("Oxidated methionine",
                                            "oxM",
                                            1.73,
                                            131.1926 + 15.9994,
                                            131.04049 + 15.994915));
            addChemicalGroup(ChemicalGroup ("Asparagine",
                                            "N",
                                            0.28,
                                            114.1038,
                                            114.04293));
            addChemicalGroup(ChemicalGroup ("Proline",
                                            "P",
                                            0.62,
                                            97.1167,
                                            97.05276));
            addChemicalGroup(ChemicalGroup ("Glutamine",
                                            "Q",
                                            0.47,
                                            128.1307,
                                            128.05858));
            addChemicalGroup(ChemicalGroup ("Arginine",
                                            "R",
                                            -0.59,
                                            156.1875,
                                            156.10111));
            addChemicalGroup(ChemicalGroup ("Serine",
                                            "S",
                                            0.50,
                                            87.0782,
                                            87.03203));
            addChemicalGroup(ChemicalGroup ("Phosphorylated serine",
                                            "pS",
                                            0.55,
                                            167.0581,
                                            166.99836));
            addChemicalGroup(ChemicalGroup ("Threonine",
                                            "T",
                                            0.69,
                                            101.1051,
                                            101.04768));
            addChemicalGroup(ChemicalGroup ("Phosphorylated threonine",
                                            "pT",
                                            0.59,
                                            181.085,
                                            181.01401));
            addChemicalGroup(ChemicalGroup ("Valine",
                                            "V",
                                            1.49,
                                            99.1326,
                                            99.06841));
            addChemicalGroup(ChemicalGroup ("Tryptophan",
                                            "W",
                                            2.87,
                                            186.2132,
                                            186.07931));
            addChemicalGroup(ChemicalGroup ("Tyrosine",
                                            "Y",
                                            1.40,
                                            163.176,
                                            163.06333));
            addChemicalGroup(ChemicalGroup ("Phosphorylated tyrosine",
                                            "pY",
                                            0.92,
                                            243.1559,
                                            243.02966));
            addChemicalGroup(ChemicalGroup("N-terminal hydrogen",
                                           "H-",
                                           -2.39,
                                           1.0079,
                                           1.00782,
                                           0.0));
            addChemicalGroup(ChemicalGroup("N-terminal acetyl",
                                           "Ac-",
                                           -0.05,
                                           43.0452,
                                           43.01839,
                                           0.0));
            addChemicalGroup(ChemicalGroup("C-terminal carboxyl group",
                                           "-OH",
                                           1.02,
                                           17.0073,
                                           17.00274,
                                           0.0));
            addChemicalGroup(ChemicalGroup("C-terminal amide",
                                           "-NH2",
                                           0.55,
                                           16.0226,
                                           16.01872,
                                           0.0));
            break;
        }
    }

    return (*this);
}
}
