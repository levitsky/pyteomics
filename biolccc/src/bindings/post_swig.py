import sys
import os
current_dir = os.path.dirname(__file__)
modified_file = []

with open(os.path.join(current_dir, 'biolccc.py'), 'r') as f:
    for line in f:
        modified_line = line
        modified_line = modified_line.replace(
            'class ChemicalGroup(_object)',
            'class ChemicalGroup(collections.MutableMapping, _object)')
        modified_line = modified_line.replace(
            'class ChemicalBasis(_object)',
            'class ChemicalBasis(collections.MutableMapping, _object)')
        modified_line = modified_line.replace(
            'class ChromoConditions(_object)',
            'class ChromoConditions(collections.MutableMapping, _object)')
        modified_line = modified_line.replace(
            'class GradientPoint(_object)',
            'class GradientPoint(collections.MutableMapping, _object)')
        modified_file.append(modified_line)

    newlines = f.newlines

if modified_file[0].count('import collections') == 0:
    modified_file.insert(0, 'import collections' + (newlines or '\n'))

with open(os.path.join(current_dir, 'biolccc.py'), 'w') as f:
    f.writelines(modified_file)
    
