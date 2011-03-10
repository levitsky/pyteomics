// biolccc.i - SWIG interface
%module biolccc 

%feature("autodoc", "0");

%{
#include "biolcccexception.h"
#include "chemicalgroup.h"
#include "chemicalbasis.h"
#include "gradientpoint.h"
#include "gradient.h"
#include "chromoconditions.h"
#include "auxiliary.h"
#include "parsing.h"
#include "rod_model.h"
#include "chain_model.h"
#include "biolccc.h"
%}

%include "std_string.i"
%include "std_map.i"
%include "std_vector.i"
%include "std_list.i"
%include "carrays.i"
%include "pyabc.i"
%include "exception.i"

%array_class(double, doubleArray);
%pythonabc(ChemicalGroup, collections.MutableMapping);
%pythonabc(ChemicalBasis, collections.MutableMapping);
%pythonabc(ChromoConditions, collections.MutableMapping);
%pythonabc(GradientPoint, collections.MutableMapping);
%template(GradientPointVector) std::vector<BioLCCC::GradientPoint>;
%template(ChemicalGroupVector) std::vector<BioLCCC::ChemicalGroup>;
%template(DoubleVector) std::vector<double>;
%template(StringChemicalGroupMap) std::map<std::string,BioLCCC::ChemicalGroup>;
%template(StringChemicalGroupPtrMap) std::map<std::string,BioLCCC::ChemicalGroup *> ;
%rename(__chemicalGroups__) chemicalGroups();
%ignore StringChemicalGroupMap::operator[] const;
%ignore BioLCCC::ChemicalBasis::chemicalGroups() const;

// Parse the original header file
%include "biolcccexception.h"
%include "chemicalgroup.h"
%include "chemicalbasis.h"
%include "gradientpoint.h"
%include "gradient.h"
%include "chromoconditions.h"
%include "auxiliary.h"
%include "parsing.h"
%include "chain_model.h"
%include "rod_model.h"
%include "biolccc.h"

%extend std::map<std::string,BioLCCC::ChemicalGroup>{
    %insert("python") %{
        def __str__(self):
            return str(dict(self))

        def __repr__(self):
            return str(dict(self))

        def __eq__(self, other):
            return dict(self) == dict(other)
    %}
}

%extend std::map<std::string,BioLCCC::ChemicalGroup *>{
    %insert("python") %{
        def __str__(self):
            return str(dict(self))

        def __repr__(self):
            return str(dict(self))

        def __eq__(self, other):
            return dict(self) == dict(other)
    %}
}

%extend BioLCCC::Gradient {
    %insert("python") %{
        def __str__(self):
            return str(list(self))

        def __repr__(self):
            return str(list(self))

        def __eq__(self, other):
            return list(self) == list(other)
    %}
}

%extend BioLCCC::ChemicalGroup{
    %insert("python") %{
        def __eq__(self, other):
            return dict(self) == dict(other)

        def __str__(self):
            return str(dict(self))

        def __repr__(self):
            return str(dict(self))

        def __len__(self):
            return len(self._keys);

        def __iter__(self):
            return iter(self._keys)

        def __contains__(self, key):
            return key in self._keys 

        def __getitem__(self, key):
            return {
                'name' : self.name,
                'label': self.label,
                'bindEnergy': self.bindEnergy,
                'bindArea': self.bindArea,
                'averageMass': self.averageMass,
                'monoisotopicMass': self.monoisotopicMass,
            }[key]()

        def __raiseLabelException__(self, value):
            raise RuntimeError('Label cannot be set')

        def __setitem__(self, key, value):
            {
                'name' : self.setName,
                'label': self.__raiseLabelException__,
                'bindEnergy': self.setBindEnergy,
                'bindArea': self.setBindArea,
                'averageMass': self.setAverageMass,
                'monoisotopicMass': self.setMonoisotopicMass,
            }[key](value)

        def __delitem__(self, key):
            pass

        _keys = ['name', 'label', 'bindEnergy', 'bindArea', 'averageMass',
                'monoisotopicMass']

        def keys(self):
            return self._keys
    %}
}

%extend BioLCCC::ChemicalBasis {
    std::map<std::string, BioLCCC::ChemicalGroup *> 
        __ptrChemicalGroups__()
    {
        std::map<std::string, BioLCCC::ChemicalGroup * > 
            ptrChemicalGroupsMap;
        for (
            std::map<std::string, BioLCCC::ChemicalGroup>::iterator it = 
                $self->BioLCCC::ChemicalBasis::chemicalGroups().begin();
            it != $self->BioLCCC::ChemicalBasis::chemicalGroups().end();
            it++)
        {
            ptrChemicalGroupsMap[it->first] = &(it->second);
        }
        return ptrChemicalGroupsMap;
    }
}

%extend BioLCCC::ChemicalBasis {
    %insert("python") %{
        def __eq__(self, other):
            return dict(self) == dict(other)

        def __repr__(self):
            return str(dict(self))

        def __str__(self):
            return str(dict(self))

        def __len__(self):
            return len(self._keys);

        def __iter__(self):
            return iter(self._keys)

        def __contains__(self, key):
            return key in self._keys 

        def chemicalGroups(self):
            return self.__ptrChemicalGroups__()

        def __getitem__(self, key):
            return {
                'chemicalGroups': self.__ptrChemicalGroups__,
                'firstSolventDensity': self.firstSolventDensity,
                'firstSolventAverageMass': self.firstSolventAverageMass,
                'secondSolventDensity': self.secondSolventDensity,
                'secondSolventAverageMass': self.secondSolventAverageMass,
                'secondSolventBindEnergy': self.secondSolventBindEnergy,
                'adsorptionLayerWidth': self.adsorptionLayerWidth,
                'adsorptionLayerFactors': self.adsorptionLayerFactors,
                'kuhnLength': self.kuhnLength,
                'monomerLength': self.monomerLength,
                'polymerModel': self.polymerModel,
                'snyderApproximation': self.snyderApproximation,
                'specialRodModel': self.specialRodModel,
            }[key]()

        def __setitem__(self, key, value):
            {
                'chemicalGroups' : self.setChemicalGroups,
                'firstSolventDensity': self.setFirstSolventDensity,
                'firstSolventAverageMass': self.setFirstSolventAverageMass,
                'secondSolventDensity': self.setSecondSolventDensity,
                'secondSolventAverageMass': self.setSecondSolventAverageMass,
                'secondSolventBindEnergy': self.setSecondSolventBindEnergy,
                'adsorptionLayerWidth': self.setAdsorptionLayerWidth,
                'adsorptionLayerFactors': self.setAdsorptionLayerFactors,
                'kuhnLength': self.setKuhnLength,
                'monomerLength': self.setMonomerLength,
                'polymerModel': self.setPolymerModel,
                'snyderApproximation': self.setSnyderApproximation,
                'specialRodModel': self.setSpecialRodModel,
            }[key](value)

        def __delitem__(self, key):
            pass

        _keys = ['chemicalGroups', 'firstSolventDensity', 
                 'firstSolventAverageMass',
                 'secondSolventDensity', 'secondSolventAverageMass',
                 'secondSolventBindEnergy', 'adsorptionLayerWidth', 
                 'adsorptionLayerFactors', 'kuhnLength', 'monomerLength', 
                 'polymerModel', 'snyderApproximation', 'specialRodModel']

        def keys(self):
            return self._keys

        def __getstate__(self):
            state_dict = {}
            for key in self:
                if key != 'chemicalGroups':
                    state_dict[key] = self[key]
            state_dict['chemicalGroups'] = {}
            for key in self['chemicalGroups']:
                state_dict['chemicalGroups'][key] = (
                    dict(self['chemicalGroups'][key]))
            return state_dict

        def __setstate__(self, state_dict):
            for key in state_dict:
                self[key] = state_dict[key]

        def __reduce__(self):
            return (ChemicalBasis, (), self.__getstate__(),)

        def setChemicalGroups(self, chemicalGroupsDict):
            self.clearChemicalGroups()
            for key, value in chemicalGroupsDict.items():
                if type(value).__name__ == 'dict':
                    self.addChemicalGroup(
                        ChemicalGroup(
                            value['name'],
                            value['label'],
                            value['bindEnergy'],
                            value['averageMass'],
                            value['monoisotopicMass'],
                            value['bindArea']))
                elif type(value).__name__ == 'biolccc.ChemicalBasis':
                    self.addChemicalGroup(value)
                else:
                    raise Exception('biolccc', 'wrong type for ChemicalGroup')
    %}
};

%extend BioLCCC::GradientPoint{
    %insert("python") %{
        def __eq__(self, other):
            return dict(self) == dict(other)

        def __str__(self):
            return str(dict(self))

        def __repr__(self):
            return str(dict(self))

        def __len__(self):
            return len(self._keys);

        def __iter__(self):
            return iter(self._keys)

        def __contains__(self, key):
            return key in self._keys 

        def __getitem__(self, key):
            return {
                'time' : self.time,
                'concentrationB': self.concentrationB,
            }[key]()

        def __setitem__(self, key, value):
            {
                'time' : self.setName,
                'concentrationB': self.setConcentrationB,
            }[key](value)

        def __delitem__(self, key):
            pass

        _keys = ['time', 'concentrationB']

        def keys(self):
            return self._keys
    %}
};

%extend BioLCCC::ChromoConditions{
    %insert("python") %{
        def __eq__(self, other):
            return dict(self) == dict(other)

        def __str__(self):
            return str(dict(self))

        def __repr__(self):
            return str(dict(self))

        def __len__(self):
            return len(self._keys);

        def __iter__(self):
            return iter(self._keys)

        def __contains__(self, key):
            return key in self._keys 

        def __getitem__(self, key):
            return {
                'columnLength': self.columnLength,
                'columnDiameter': self.columnDiameter,
                'columnPoreSize': self.columnPoreSize,
                'columnRelativeStrength': self.columnRelativeStrength,
                'columnVpToVtot': self.columnVpToVtot,
                'columnPorosity': self.columnPorosity,
                'delayTime': self.delayTime,
                'dV': self.dV,
                'flowRate': self.flowRate,
                'gradient': self.gradient,
                'secondSolventConcentrationA':
                    self.secondSolventConcentrationA,
                'secondSolventConcentrationB':
                    self.secondSolventConcentrationB,
                'temperature': self.temperature,
            }[key]()

        def __setitem__(self, key, value):
            return {
                'columnLength': self.setColumnLength,
                'columnDiameter': self.setColumnDiameter,
                'columnPoreSize': self.setColumnPoreSize,
                'columnRelativeStrength': self.setColumnRelativeStrength,
                'columnVpToVtot': self.setColumnVpToVtot,
                'columnPorosity': self.setColumnPorosity,
                'delayTime': self.setDelayTime,
                'dV': self.setDV,
                'flowRate': self.setFlowRate,
                'gradient': self.setGradient,
                'secondSolventConcentrationA':
                    self.setSecondSolventConcentrationA,
                'secondSolventConcentrationB':
                    self.setSecondSolventConcentrationB,
                'temperature': self.setTemperature,
            }[key](value)

        def __delitem__(self, key):
            pass

        _keys = ['columnLength', 'columnDiameter', 'columnPoreSize', 
                 'gradient', 'secondSolventConcentrationA',
                 'secondSolventConcentrationB', 'delayTime', 'flowRate',
                 'dV', 'columnRelativeStrength', 'columnVpToVtot',
                 'columnPorosity', 'temperature']

        def keys(self):
            return self._keys

        def __getstate__(self):
            state_dict = {}
            for key in self:
                if key != 'gradient':
                    state_dict[key] = self[key]
            state_dict['gradient'] = []
            for point in self['gradient']:
                state_dict['gradient'].append(dict(point))
            return state_dict

        def __setstate__(self, state_dict):
            for key in state_dict:
                if key != 'gradient':
                    self[key] = state_dict[key]
            gradient = Gradient()
            for point in state_dict['gradient']:
                if type(point).__name__ == 'dict':
                    gradient.addPoint(
                        point['time'], point['concentrationB'])
                elif type(point).__name__ == 'biolccc.GradientPoint':
                    gradient.addPoint(point)
                else:
                    raise Exception('biolccc', 'wrong type for GradientPoint')
            self['gradient'] = gradient

        def __reduce__(self):
            return (ChromoConditions, (), self.__getstate__(),)
    %}
};

// Instantiate some templates

