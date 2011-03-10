#include <cmath>
#include <vector>
#include <time.h>
#include <gtest/gtest.h>
#include "biolccc.h"

// The fixture for testing BioLCCC functions.
class BioLCCCTest: public ::testing::Test
{
protected:
    BioLCCCTest()
    {
        // The set-up work for each test.
        testChemBasis = BioLCCC::ChemicalBasis();
        testChemBasis.setPolymerModel(BioLCCC::CHAIN);
        testChemBasis.setFirstSolventDensity(1.0);
        testChemBasis.setFirstSolventAverageMass(1.0);
        testChemBasis.setSecondSolventDensity(1.0);
        testChemBasis.setSecondSolventAverageMass(1.0);
        testChemBasis.setSecondSolventBindEnergy(0.0);
        testChemBasis.setMonomerLength(10.0);
        testChemBasis.setKuhnLength(10.0);
        testChemBasis.setAdsorptionLayerWidth(1.0e-10);
        testChemBasis.setAdsorptionLayerFactors(std::vector<double>(1, 1.0));
        testChemBasis.setSnyderApproximation(false);
        testChemBasis.addChemicalGroup(BioLCCC::ChemicalGroup(
            "Blank amino acid",
            "B",
            0.0,
            0.0,
            0.0));
        testChemBasis.addChemicalGroup(BioLCCC::ChemicalGroup(
            "Hydrophobic amino acid",
            "O",
            1.0,
            0.0,
            0.0));
        testChemBasis.addChemicalGroup(BioLCCC::ChemicalGroup(
            "Hydrophilic amino acid",
            "Z",
            -1.0,
            0.0,
            0.0));
        testChemBasis.addChemicalGroup(BioLCCC::ChemicalGroup(
            "empty N-terminal group",
            "H-",
            0.0,
            0.0,
            0.0, 
            0.0));
        testChemBasis.addChemicalGroup(BioLCCC::ChemicalGroup(
            "empty C-Terminal group",
            "-OH",
            0.0,
            0.0,
            0.0,
            0.0));

    }

    virtual ~BioLCCCTest()
    {
        // You can do clean-up work that doesn't throw exceptions here.
    }

    // If the constructor and destructor are not enough for setting up
    // and cleaning up each test, you can define the following methods:

    virtual void SetUp()
    {
        // Code here will be called immediately after the constructor (right
        // before each test).
    }

    virtual void TearDown()
    {
        // Code here will be called immediately after each test (right
        // before the destructor).
    }

    void parsingRulesTester(BioLCCC::ChemicalBasis basisToTest)
    {
        std::string peptideSource = "QWERTYIPASDFGHKLCVNM";

        std::vector<BioLCCC::ChemicalGroup> parsedSequence = 
            BioLCCC::parseSequence("foo.QWERTYIPASDFGHKLCVNM.bar", basisToTest);

        for (unsigned int i=0; i<peptideSource.size(); i++)
        {
            ASSERT_EQ(parsedSequence[i+1].label(),
                      peptideSource.substr(i, 1));
        }

        peptideSource = "Ac-QWERTYIPASDFGHKLCVNM-NH2";

        ASSERT_THROW(
            BioLCCC::parseSequence("ac-QWERTY-H", basisToTest),
            BioLCCC::ParsingException);

        ASSERT_THROW(
            BioLCCC::parseSequence("QWERTY-oH", basisToTest),
            BioLCCC::ParsingException);

        ASSERT_THROW(
            BioLCCC::parseSequence("QWERaTY", basisToTest),
            BioLCCC::ParsingException);

        ASSERT_THROW(
            BioLCCC::parseSequence("Ac-QWERaTY-NH2", basisToTest),
            BioLCCC::ParsingException);

        ASSERT_THROW(
            BioLCCC::parseSequence("Ac-QWERaTY-NH-NH2", basisToTest),
            BioLCCC::ParsingException);
        
        ASSERT_THROW(
            BioLCCC::parseSequence("A.Ac-QWERaTY-NH-NH2.K.K", basisToTest),
            BioLCCC::ParsingException);
    }

    void parsesStandardAminoacidsTester(BioLCCC::ChemicalBasis basisToTest)
    {
        std::string peptideSource = "QWERTYIPASDFGHKLCVNM";

        std::vector<BioLCCC::ChemicalGroup> parsedSequence = 
            BioLCCC::parseSequence(peptideSource, basisToTest);

        for (unsigned int i=0; i<peptideSource.size(); i++)
        {
            ASSERT_EQ(parsedSequence[i+1].label(),
                      peptideSource.substr(i, 1));
        }
    }

    void parsesPhosphoAminoacidsTester(BioLCCC::ChemicalBasis basisToTest)
    {
        std::string peptideSource = "pSpTpY";

        std::vector<BioLCCC::ChemicalGroup> parsedSequence =
            BioLCCC::parseSequence(peptideSource, basisToTest);

        for (unsigned int i=0; i<peptideSource.size()/2; i++)
        {
            ASSERT_EQ(parsedSequence[i+1].label(),
                      peptideSource.substr(i*2, 2));
        }
    }

    void parsesStandardTerminalGroupsTester(BioLCCC::ChemicalBasis basisToTest)
    {
        std::vector<BioLCCC::ChemicalGroup> parsedSequence = 
            BioLCCC::parseSequence("GGGG", basisToTest);

        ASSERT_EQ(parsedSequence.begin()->label(), "H-");
        ASSERT_EQ((--(parsedSequence.end()))->label(), "-OH");

        parsedSequence = 
            BioLCCC::parseSequence("H-GGGG-OH", basisToTest);

        ASSERT_EQ(parsedSequence.begin()->label(), "H-");
        ASSERT_EQ((--(parsedSequence.end()))->label(), "-OH");

        parsedSequence = 
            BioLCCC::parseSequence("Ac-GGGG-NH2", basisToTest);

        ASSERT_EQ(parsedSequence.begin()->label(), "Ac-");
        ASSERT_EQ((--(parsedSequence.end()))->label(), "-NH2");
    }

    void calculatesMonoisotopicMassTester(BioLCCC::ChemicalBasis basisToTest)
    {
        ASSERT_LE(
            fabs(BioLCCC::calculateMonoisotopicMass("QWERTYIPASDFGHKLCVNM",
                basisToTest) - 2394.1248),
            0.0001);
        ASSERT_LE(
            fabs(BioLCCC::calculateMonoisotopicMass("Ac-QWERTYIPASDFGHKLCVNM",
                basisToTest) - 2436.1354), 0.0001);
        ASSERT_LE(
            fabs(BioLCCC::calculateMonoisotopicMass("QWERTYIPASDFGHKLCVNM-NH2",
                basisToTest) - 2393.1408), 0.0001);
    }

    void calculatesKdTester(BioLCCC::ChemicalBasis basisToTest)
    {
        ASSERT_GT(BioLCCC::calculateKd("QWERTYIPASDFGHKLCVNM", 0.0,
            basisToTest), 0.0);
        ASSERT_GT(BioLCCC::calculateKd("QWERTYIPASDFGHKLCVNM", 50.0,
            basisToTest), 0.0);
        ASSERT_GT(BioLCCC::calculateKd("QWERTYIPASDFGHKLCVNM", 100.0,
            basisToTest), 0.0);
    }

    void calculatesRTTester(BioLCCC::ChemicalBasis basisToTest)
    {
        ASSERT_GT(BioLCCC::calculateRT("QWERTYIPASDFGHKLCVNM",
            basisToTest), 0.0);
    }

    BioLCCC::ChemicalBasis testChemBasis;
};

TEST_F(BioLCCCTest, parsingRules)
{
    parsingRulesTester(BioLCCC::rpAcnTfaChain);
    parsingRulesTester(BioLCCC::rpAcnFaRod);
}

TEST_F(BioLCCCTest, parsesStandardAminoacids)
{
    parsesStandardAminoacidsTester(BioLCCC::rpAcnTfaChain);
    parsesStandardAminoacidsTester(BioLCCC::rpAcnFaRod);
}

TEST_F(BioLCCCTest, parsesPhosphoAminoacids)
{
    parsesPhosphoAminoacidsTester(BioLCCC::rpAcnTfaChain);
    parsesPhosphoAminoacidsTester(BioLCCC::rpAcnFaRod);
}

TEST_F(BioLCCCTest, parsesStandardTerminalGroups)
{
    parsesStandardTerminalGroupsTester(BioLCCC::rpAcnTfaChain);
    parsesStandardTerminalGroupsTester(BioLCCC::rpAcnFaRod);
}

TEST_F(BioLCCCTest, calculatesMonoisotopicMass)
{
    calculatesMonoisotopicMassTester(BioLCCC::rpAcnTfaChain);
    calculatesMonoisotopicMassTester(BioLCCC::rpAcnFaRod);
}

TEST_F(BioLCCCTest, calculatesKd)
{
    calculatesKdTester(BioLCCC::rpAcnTfaChain);
    calculatesKdTester(BioLCCC::rpAcnFaRod);
}

TEST_F(BioLCCCTest, calculatesRT)
{
    calculatesRTTester(BioLCCC::rpAcnTfaChain);
    calculatesRTTester(BioLCCC::rpAcnFaRod);
}

TEST_F(BioLCCCTest, matrixTest)
{
    for (float i = 3.0; i < 11.0; ++i)
    {
        ASSERT_EQ(
            BioLCCC::calculateKd("H-B-OH", 0.0,
                testChemBasis, i * 10.0),
            1.0);
        ASSERT_LT(
            abs(BioLCCC::calculateKd("H-BB-OH", 0.0, testChemBasis, i * 10.0) -
               (1.0 - 2.0/6.0/i)),
            1e-5);
    }

    for (float i = 3.0; i < 11.0; ++i)
    {
        ASSERT_LT(
            abs(BioLCCC::calculateKd("H-O-OH", 0.0, testChemBasis, i * 10.0) 
                - (1.0 + (2.0 * (exp(1.0) - 1.0) / i))),
            1e-5);
        ASSERT_LT(
            abs(BioLCCC::calculateKd("H-OO-OH", 0.0, testChemBasis, i * 10.0) 
                - (1.0 + (4.0 * exp(1.0) * exp(1.0) + 2.0 * exp(1.0)  - 7.0) 
                         / 3.0 / i)),
            1e-5);
    }

    std::vector<double> adsorptionLayerFactors;
    adsorptionLayerFactors.push_back(0.0);
    adsorptionLayerFactors.push_back(1.0);
    testChemBasis.setAdsorptionLayerFactors(adsorptionLayerFactors);

    for (float i = 6.0; i < 11.0; ++i)
    {
        ASSERT_LT(
            abs(BioLCCC::calculateKd("H-O-OH", 0.0, testChemBasis, i * 10.0) 
                - (1.0 + (2.0 * (exp(1.0) - 1.0) / i))),
            1e-5);

        ASSERT_LT(
            abs(BioLCCC::calculateKd("H-OO-OH", 0.0, testChemBasis, i * 10.0) 
                - (1.0 + (4.0 * exp(1.0) * exp(1.0) + 4.0 * exp(1.0) - 9.0) 
                         / 3.0 / i)),
            1e-5);
    }
}

TEST_F(BioLCCCTest, adsorptionStrengthTest)
{
    testChemBasis.setSecondSolventBindEnergy(1.0);
    double kd1 = BioLCCC::calculateKd("H-OBO-OH", 100.0, testChemBasis, 100.0);

    testChemBasis.setSecondSolventBindEnergy(2.0);
    testChemBasis.chemicalGroups()["O"].setBindEnergy(2.0);
    double kd2 = 
        BioLCCC::calculateKd("H-OBO-OH", 100.0, testChemBasis, 100.0, 0.5);
    ASSERT_EQ(kd1, kd2);

    double kd3 = 
        BioLCCC::calculateKd("H-OBO-OH", 100.0, testChemBasis, 100.0, 0.6);
    ASSERT_NE(kd1, kd3);

    testChemBasis.setSecondSolventBindEnergy(1.0);
    testChemBasis.chemicalGroups()["O"].setBindEnergy(1.0);
    testChemBasis.setAdsorptionLayerFactors(std::vector<double>(1, 0.5));
    double kd4 = 
        BioLCCC::calculateKd("H-OBO-OH", 100.0, testChemBasis, 100.0, 2.0);
    ASSERT_EQ(kd1, kd4);

    // Check if zero column adsorption strength works as a non-adsorbed peptide 
    // in a non-adsorbed second solvent.
    double kd5 = 
        BioLCCC::calculateKd("H-OOO-OH", 100.0, testChemBasis, 100.0, 0.0);
    double kd6 = 
        BioLCCC::calculateKd("H-BBB-OH", 0.0, testChemBasis, 100.0, 1.0);
    ASSERT_EQ(kd5, kd6);
}

TEST_F(BioLCCCTest, monomerEnergyProfileTest)
{
    testChemBasis.setSecondSolventBindEnergy(1.0);
    double kd1 = BioLCCC::calculateKd("H-OOO-OH", 100.0, testChemBasis, 100.0);

    testChemBasis.setSecondSolventBindEnergy(0.0);
    double kd2 = BioLCCC::calculateKd("H-BBB-OH", 100.0, testChemBasis, 100.0);
    ASSERT_EQ(kd1, kd2);

    srand( time(NULL));
    testChemBasis.chemicalGroups()["H-"].setBindEnergy(
        (float)(rand() % 1000) / 1000.0);
    testChemBasis.chemicalGroups()["-OH"].setBindEnergy(
        (float)(rand() % 1000) / 1000.0);
    testChemBasis.addChemicalGroup(BioLCCC::ChemicalGroup(
        "N-terminal hydrophylic amino acid",
        "nO",
        1.0 - testChemBasis.chemicalGroups()["H-"].bindEnergy(),
        0.0,
        0.0));
    testChemBasis.addChemicalGroup(BioLCCC::ChemicalGroup(
        "C-terminal hydrophylic amino acid",
        "cO",
        1.0 - testChemBasis.chemicalGroups()["-OH"].bindEnergy(),
        0.0,
        0.0));

    testChemBasis.setSecondSolventBindEnergy(1.0);
    double kd3 = BioLCCC::calculateKd("H-nOOcO-OH", 100.0, testChemBasis, 
        100.0);
    ASSERT_EQ(kd1, kd3);
}

TEST_F(BioLCCCTest, segmentEnergyProfileTest)
{
    testChemBasis.setSecondSolventBindEnergy(1.0);
    testChemBasis.addChemicalGroup(BioLCCC::ChemicalGroup(
        "Hydrophobic amino acid with half energy",
        "hO",
        0.5,
        0.0,
        0.0));
    testChemBasis.addChemicalGroup(BioLCCC::ChemicalGroup(
        "Hydrophobic amino acid with double energy",
        "dO",
        2.0,
        0.0,
        0.0));

    double kd1 = BioLCCC::calculateKd("H-OOOO-OH", 100.0, testChemBasis, 100.0);

    testChemBasis.setKuhnLength(20.0);
    testChemBasis.setSecondSolventBindEnergy(0.5);
    double kd2 = 
        BioLCCC::calculateKd("H-hOhOhOhOhOhOhOhO-OH", 100.0, 
            testChemBasis, 200.0);
    ASSERT_EQ(kd1, kd2);

    testChemBasis.setKuhnLength(5.0);
    testChemBasis.setSecondSolventBindEnergy(2.0);
    double kd3 = BioLCCC::calculateKd("H-dOdO-OH", 100.0, testChemBasis, 50.0);
    ASSERT_EQ(kd1, kd3);
}

TEST_F(BioLCCCTest, backwardCalculationCompatibility)
{
    BioLCCC::ChemicalBasis chemBasisChain(BioLCCC::RP_ACN_TFA_CHAIN);
    chemBasisChain.setFirstSolventDensity(5.56);
    chemBasisChain.setFirstSolventAverageMass(1.0);
    chemBasisChain.setSecondSolventDensity(1.91);
    chemBasisChain.setSecondSolventAverageMass(1.0);

    ASSERT_LT(
        abs(BioLCCC::calculateKd("QWERTYIPASDFGHKLCVNM", 20.0, chemBasisChain)
            - 104.76633),
        1e-5);

    ASSERT_LT(
        abs(BioLCCC::calculateRT("QWERTYIPASDFGHKLCVNM", chemBasisChain,
            BioLCCC::standardChromoConditions, 0, true, true) 
            - 43.48354),
        1e-5);

    BioLCCC::ChromoConditions accurate_chromoconditions;
    accurate_chromoconditions.setDV(
        accurate_chromoconditions.flowRate() / 1000.0);
    ASSERT_LT(
        abs(BioLCCC::calculateRT("QWERTYIPASDFGHKLCVNM", chemBasisChain,
            BioLCCC::standardChromoConditions, 0, true, true) 
            - BioLCCC::calculateRT("QWERTYIPASDFGHKLCVNM", chemBasisChain,
            BioLCCC::standardChromoConditions, 0, true, false)),
        1e-5);

    BioLCCC::ChemicalBasis chemBasisRod(BioLCCC::RP_ACN_FA_ROD);
    chemBasisRod.setFirstSolventDensity(5.56);
    chemBasisRod.setFirstSolventAverageMass(1.0);
    chemBasisRod.setSecondSolventDensity(1.91);
    chemBasisRod.setSecondSolventAverageMass(1.0);

    ASSERT_LT(
        abs(BioLCCC::calculateKd("QWERTYIPASDFGHKLCVNM", 20.0, chemBasisRod)
            - 22.64802),
        1e-5);

    ASSERT_LT(
        abs(BioLCCC::calculateRT("QWERTYIPASDFGHKLCVNM", chemBasisRod,
            BioLCCC::standardChromoConditions, 0, true, true) 
            - 36.53354),
        1e-5);

    ASSERT_LT(
        abs(BioLCCC::calculateRT("QWERTYIPASDFGHKLCVNM", chemBasisRod,
            BioLCCC::standardChromoConditions, 0, true, true) 
            - BioLCCC::calculateRT("QWERTYIPASDFGHKLCVNM", chemBasisRod,
            BioLCCC::standardChromoConditions, 0, true, false)),
        1e-5);
}

TEST_F(BioLCCCTest, assignsChemicalBasis) 
{
    BioLCCC::ChemicalBasis myChemicalBasis;
    srand( time(NULL) );
    double testValue;

    myChemicalBasis.addChemicalGroup(BioLCCC::ChemicalGroup(
        "Test1",
        "tT",
        1234.0,
        1235.0,
        1236.0));
    ASSERT_EQ(myChemicalBasis.chemicalGroups()["tT"].name(), "Test1");
    ASSERT_EQ(myChemicalBasis.chemicalGroups()["tT"].label(), "tT");
    ASSERT_EQ(myChemicalBasis.chemicalGroups()["tT"].bindEnergy(), 1234.0);
    ASSERT_EQ(myChemicalBasis.chemicalGroups()["tT"].averageMass(), 1235.0);
    ASSERT_EQ(myChemicalBasis.chemicalGroups()["tT"].monoisotopicMass(), 
        1236.0);

    myChemicalBasis.chemicalGroups()["tT"].setName("Test");
    ASSERT_EQ(myChemicalBasis.chemicalGroups()["tT"].name(), "Test");

    testValue = (float)(rand() % 1000);
    myChemicalBasis.chemicalGroups()["tT"].setBindEnergy(testValue);
    ASSERT_EQ(myChemicalBasis.chemicalGroups()["tT"].bindEnergy(), testValue);

    testValue = (float)(rand() % 1000);
    myChemicalBasis.chemicalGroups()["tT"].setAverageMass(testValue);
    ASSERT_EQ(myChemicalBasis.chemicalGroups()["tT"].averageMass(), testValue);

    testValue = (float)(rand() % 1000);
    myChemicalBasis.chemicalGroups()["tT"].setMonoisotopicMass(testValue);
    ASSERT_EQ(myChemicalBasis.chemicalGroups()["tT"].monoisotopicMass(),
              testValue);

    testValue = (float)(rand() % 1000);
    myChemicalBasis.setFirstSolventDensity(testValue);
    ASSERT_EQ(myChemicalBasis.firstSolventDensity(), testValue);

    testValue = (float)(rand() % 1000);
    myChemicalBasis.setSecondSolventDensity(testValue);
    ASSERT_EQ(myChemicalBasis.secondSolventDensity(), testValue);

    testValue = (float)(rand() % 1000);
    myChemicalBasis.setFirstSolventAverageMass(testValue);
    ASSERT_EQ(myChemicalBasis.firstSolventAverageMass(), testValue);

    testValue = (float)(rand() % 1000);
    myChemicalBasis.setSecondSolventAverageMass(testValue);
    ASSERT_EQ(myChemicalBasis.secondSolventAverageMass(), testValue);

    testValue = (float)(rand() % 1000);
    myChemicalBasis.setSecondSolventBindEnergy(testValue);
    ASSERT_EQ(myChemicalBasis.secondSolventBindEnergy(), testValue);

    testValue = (float)(rand() % 1000);
    myChemicalBasis.setAdsorptionLayerWidth(testValue);
    ASSERT_EQ(myChemicalBasis.adsorptionLayerWidth(), testValue);

    testValue = (float)(rand() % 1000);
    myChemicalBasis.setKuhnLength(testValue);
    ASSERT_EQ(myChemicalBasis.kuhnLength(), testValue);

    testValue = (float)(rand() % 1000);
    myChemicalBasis.setMonomerLength(testValue);
    ASSERT_EQ(myChemicalBasis.monomerLength(), testValue);
}

TEST_F(BioLCCCTest, fitsSpline) 
{
    double * x = new double[11];
    double * y = new double[11];
    double * y2 = new double[11];

    for (int i=0; i<11; i++)
    {
        x[i] = i;
        y[i] = 1.0 / (x[i] + 1.0);
    }
    BioLCCC::fitSpline(x, y, 11, y2);

    for (int i=0; i<11; i++)
    {
        ASSERT_LE(
            fabs(BioLCCC::calculateSpline(x, y, y2, 11, i/2.0) * (i / 2.0 + 1.0)
                 - 1.0),
            0.1);
    }
    delete[] x, y, y2;
}

TEST_F(BioLCCCTest, partialPolynomialInterpolation) 
{
    double * x = new double[11];
    double * y = new double[11];

    for (int i=0; i<11; i++)
    {
        x[i] = i;
        y[i] = 1.0 / (x[i] + 1.0);
    }

    for (int i=0; i<11; i++)
    {
        ASSERT_LE(
            fabs(BioLCCC::partPolInterpolate(x, y, 11, 3, i/2.0) 
                 * (i / 2.0 + 1.0)
                 - 1.0),
            0.05);
    }
    delete[] x, y;
}

TEST_F(BioLCCCTest, calculatesRTwithInterpolation) 
{
    std::vector<std::string> peptideStandard;
    peptideStandard.push_back("KYIPGTK");
    peptideStandard.push_back("YIPGTK");
    peptideStandard.push_back("IFVQK");
    peptideStandard.push_back("KTGQAPGFSYTDANK");
    peptideStandard.push_back("TGQAPGFSYTDANK");
    peptideStandard.push_back("GEREDLIAYLKK");
    peptideStandard.push_back("TGPNLHGLFGR");
    peptideStandard.push_back("MIFAGIK");
    peptideStandard.push_back("EDLIAYLK");
    peptideStandard.push_back("IFVQKCAQCHTVEK");
    peptideStandard.push_back("GITWGEETLMEYLENPKK");
    peptideStandard.push_back("GITWGEETLMEYLENPK");
    peptideStandard.push_back(std::string(20, 'T'));
    peptideStandard.push_back(std::string(30, 'T'));
    peptideStandard.push_back(std::string(40, 'T'));
    peptideStandard.push_back(std::string(50, 'T'));

    for (std::vector<std::string>::const_iterator peptide =
            peptideStandard.begin();
        peptide != peptideStandard.end();
        peptide++)
    {
        ASSERT_LE(
            fabs(BioLCCC::calculateRT(*peptide, BioLCCC::rpAcnTfaChain,
                                      BioLCCC::standardChromoConditions, 21) 
                 / BioLCCC::calculateRT(*peptide, BioLCCC::rpAcnTfaChain,
                                      BioLCCC::standardChromoConditions, 0) 
                 - 1.0),
            0.01);
    }
}

TEST_F(BioLCCCTest, rodModelSpecial) 
{
    BioLCCC::ChemicalBasis chemBasis1(BioLCCC::RP_ACN_FA_ROD);
    BioLCCC::ChemicalBasis chemBasis2(BioLCCC::RP_ACN_FA_ROD);
    chemBasis1.setMonomerLength(10.0);
    chemBasis1.setKuhnLength(10.0);
    chemBasis1.setAdsorptionLayerWidth(20.0);
    chemBasis2.setSpecialRodModel(true);
    chemBasis2.setMonomerLength(10.0);
    chemBasis2.setKuhnLength(10.0);
    chemBasis2.setAdsorptionLayerWidth(1.0e-10);
    chemBasis2.setSpecialRodModel(true);
    std::string sequence1("QWERTYIPAS");
    std::string sequence2("DFGHKLCVNM");
    
    double kd1 = BioLCCC::calculateKd(sequence1, 0.0, chemBasis1, 300.0, 0.0);
    double kd2 = BioLCCC::calculateKd(sequence2, 0.0, chemBasis2, 300.0, 0.0);
    ASSERT_DOUBLE_EQ(kd1, kd2);
}

TEST_F(BioLCCCTest, rodModelGeneral) 
{
    // Test if the width of adsorbing layer does not affect Kd when no
    // adsorption occurs.
    BioLCCC::ChemicalBasis chemBasis1(BioLCCC::RP_ACN_FA_ROD);
    BioLCCC::ChemicalBasis chemBasis2(BioLCCC::RP_ACN_FA_ROD);
    chemBasis1.setSpecialRodModel(false);
    chemBasis1.setMonomerLength(10.0);
    chemBasis1.setKuhnLength(10.0);
    chemBasis1.setAdsorptionLayerWidth(30.0);
    chemBasis2.setSpecialRodModel(false);
    chemBasis2.setMonomerLength(10.0);
    chemBasis2.setKuhnLength(10.0);
    chemBasis2.setAdsorptionLayerWidth(1.0e-10);
    std::string sequence(30, 'W');
    
    double kd1 = BioLCCC::calculateKd(sequence, 0.0, chemBasis1, 300.0, 0.0);
    double kd2 = BioLCCC::calculateKd(sequence, 0.0, chemBasis2, 300.0, 0.0);
    ASSERT_DOUBLE_EQ(kd1, kd2);

    // Check each term in general model.
    double segmentLength = 5.0;
    double slitWidth = 100.0;
    double layerWidth = 20.0;

    int Ns[] = {10, 20, 30, 40};
    
    for (int i=0; i<4; i++)
    {
        int N = Ns[i];
        for (int n1 = 1; n1 < N; ++n1)
        {
            for (int n2 = 1; n2 < N - n1; ++n2)
            {
                // Checking whether the general term remains the same when n1
                // and n2 are swapped.
                double term1 =
                    BioLCCC::partitionFunctionRodPartiallySubmergedTermGeneral(
                        segmentLength, slitWidth, layerWidth, N, n1, n2);
                double term2 = 
                    BioLCCC::partitionFunctionRodPartiallySubmergedTermGeneral(
                        segmentLength, slitWidth, layerWidth, N, n2, n1);
                term1 = (term1 > 1e-20) ? term1 : 0.0;
                term2 = (term2 > 1e-20) ? term2 : 0.0;
                ASSERT_FLOAT_EQ(term1, term2);

                // If the length of the part of the rod between the adsorbed 
                // parts of the rod is less then the width of the free slit.
                if ((N - n1 - n2 + 1) * segmentLength 
                       < slitWidth - 2.0*layerWidth)
                {
                    ASSERT_EQ(term1, 0.0);
                }

            }

            // Checking whether the general algorithm gives the same number
            // for the special case, when the rod do not reach the opposite
            // layer.
            bool normalRodTouchesOppositeLayer = (
                (N - n1) * segmentLength > slitWidth - 2.0 * layerWidth);
            if (!normalRodTouchesOppositeLayer) 
            {
                double pft1 =
                    BioLCCC::partitionFunctionRodPartiallySubmergedTermSpecial(
                        segmentLength, slitWidth, layerWidth, N, n1);
                double pft2 =
                    BioLCCC::partitionFunctionRodPartiallySubmergedTermGeneral(
                        segmentLength, slitWidth, layerWidth, N, n1, 0);
                ASSERT_FLOAT_EQ(pft1, pft2);
            }
            
            // Checking whether the general algorithm gives a different number
            // when the rod do reach the opposite layer.
            bool rodCanBeNormalToWall = 
                ((n1 - 1) * segmentLength <= layerWidth);
            bool inclinedRodTouchesOppositeLayer = 
                ((N - 1) / (n1 -1) > (slitWidth - layerWidth) / layerWidth);
            if ((rodCanBeNormalToWall && normalRodTouchesOppositeLayer) ||
                ((!rodCanBeNormalToWall) && inclinedRodTouchesOppositeLayer))
            {
                double pft1 =
                    BioLCCC::partitionFunctionRodPartiallySubmergedTermSpecial(
                        segmentLength, slitWidth, layerWidth, N, n1);
                double pft2 =
                    BioLCCC::partitionFunctionRodPartiallySubmergedTermGeneral(
                        segmentLength, slitWidth, layerWidth, N, n1, 0);
                SCOPED_TRACE(N);
                SCOPED_TRACE(n1);
                EXPECT_GE(fabs(pft1 - pft2) / pft1, 1e-2);
                //EXPECT_PRED_FORMAT2(::testing::FloatLE, 
                //    (pft1+pft2) * 1.0e-5, fabs(pft1 - pft2));
            }
        }
    }

    std::vector<double> rodEnergyProfile;
    for (int i = 0; i < 10; ++i)
    {
        rodEnergyProfile.push_back(0.0);
    }
    double pf1 = BioLCCC::partitionFunctionRodPartiallySubmergedSpecial(
        10.0,   //segmentLength,
        140.0, //slitWidth,
        20.0,  //layerWidth,
        rodEnergyProfile);
    double pf2 = BioLCCC::partitionFunctionRodPartiallySubmergedGeneral(
        10.0,   //segmentLength,
        140.0, //slitWidth,
        20.0,  //layerWidth,
        rodEnergyProfile);
    ASSERT_FLOAT_EQ(pf1, pf2);
}

TEST_F(BioLCCCTest, rodModelGeneralVsSpecial) 
{
    // Both models must give the same results when there is a strong
    // exclusion from the adsorbing layer. Testing the case where the
    // bridge conformations occurs.
    BioLCCC::ChemicalBasis chemBasis1(BioLCCC::RP_ACN_FA_ROD);
    BioLCCC::ChemicalBasis chemBasis2(BioLCCC::RP_ACN_FA_ROD);
    chemBasis1.setSpecialRodModel(false);
    chemBasis1.setMonomerLength(10.0);
    chemBasis1.setKuhnLength(10.0);
    chemBasis1.setAdsorptionLayerWidth(30.0);
    chemBasis2.setSpecialRodModel(true);
    chemBasis2.setMonomerLength(10.0);
    chemBasis2.setKuhnLength(10.0);
    chemBasis2.setAdsorptionLayerWidth(30.0);
    std::string sequence(30, 'W');
    
    double kd1 = BioLCCC::calculateKd(sequence, 0.0, chemBasis1, 300.0, -100.0);
    double kd2 = BioLCCC::calculateKd(sequence, 0.0, chemBasis2, 300.0, -100.0);
    ASSERT_DOUBLE_EQ(kd1, kd2);
}

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
