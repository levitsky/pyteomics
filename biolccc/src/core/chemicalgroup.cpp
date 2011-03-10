#include "chemicalgroup.h"

namespace BioLCCC
{

ChemicalGroup::ChemicalGroup(std::string name,
                             std::string label,
                             double bindEnergy,
                             double averageMass,
                             double monoisotopicMass,
                             double bindArea
                            )
{
    mName = name;
    mLabel = label;
    mAverageMass = averageMass;
    mMonoisotopicMass = monoisotopicMass;
    mBindEnergy = bindEnergy;
    mBindArea = bindArea;
}

std::string ChemicalGroup::name() const
{
    return mName;
}

std::string ChemicalGroup::label() const
{
    return mLabel;
}

double ChemicalGroup::bindEnergy() const
{
    return mBindEnergy;
}

double ChemicalGroup::bindArea() const
{
    return mBindArea;
}

double ChemicalGroup::averageMass() const
{
    return mAverageMass;
}

double ChemicalGroup::monoisotopicMass() const
{
    return mMonoisotopicMass;
}

void ChemicalGroup::setBindEnergy(double newBindEnergy)
{
    mBindEnergy = newBindEnergy;
}

void ChemicalGroup::setBindArea(double newBindArea)
{
    mBindArea = newBindArea;
}

void ChemicalGroup::setName(std::string newName)
{
    mName = newName;
}

void ChemicalGroup::setAverageMass(double newAverageMass)
{
    mAverageMass = newAverageMass;
}

void ChemicalGroup::setMonoisotopicMass(double newMonoisotopicMass)
{
    mMonoisotopicMass = newMonoisotopicMass;
}

bool ChemicalGroup::isNTerminal() const
{
    return (mLabel.find("-") == (size_t)(mLabel.size()-1));
}

bool ChemicalGroup::isCTerminal() const
{
    return (mLabel.find("-") == (size_t)(0));
}

bool ChemicalGroup::isAminoAcid() const
{
    return (!(isCTerminal() || isNTerminal()));
}
}
