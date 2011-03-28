import timeit

import os, sys
sys.path.insert(0, os.path.abspath('../'))

NUM_RUNS = 10000

stmt_calc_mass = """
sequence=''.join(
    [random.choice(mass.STD_AMINO_ACID_MASS.keys()) for i in range(20)])
x = mass.calculate_mass(sequence=sequence)
"""
print 'mass.calculate_mass: %.6f usec/pass' % (
    timeit.timeit(stmt=stmt_calc_mass,
                  setup='import mass, random; random.seed()',
                  number=NUM_RUNS) / NUM_RUNS, )

stmt_parsing_only = """
sequence=''.join(
    [random.choice(mass.STD_AMINO_ACID_MASS.keys()) for i in range(20)])
x = sum([mass.STD_AMINO_ACID_MASS[i] for i in modx.parse_sequence(sequence)])
"""
print 'modx.parse_sequence + dict.get: %.6f usec/pass' % (
    timeit.timeit(stmt=stmt_parsing_only,
                  setup='import mass, modx, random; random.seed()',

                  number=NUM_RUNS) / NUM_RUNS, )

stmt_fast_mass = """
sequence=''.join(
    [random.choice(mass.STD_AMINO_ACID_MASS.keys()) for i in range(20)])
x = mass.fast_mass(sequence=sequence, ion_type='M')
"""
print 'fast_mass: %.6f usec/pass' % (
    timeit.timeit(stmt=stmt_fast_mass,
                  setup='import mass, random; random.seed()',
                  number=NUM_RUNS) / NUM_RUNS, )


stmt_no_parsing_ion_comp = """
sequence=''.join(
    [random.choice(mass.STD_AMINO_ACID_MASS.keys()) for i in range(20)])
x = sum([mass.STD_AMINO_ACID_MASS[i] for i in sequence])
x += mass.calculate_mass(composition=mass.STD_ION_COMP['M'])
"""
print 'str + dict.get + mass.calculate_mass for ion type: %.6f usec/pass' % (
    timeit.timeit(stmt=stmt_no_parsing_ion_comp,
                  setup='import mass, random; random.seed()',
                  number=NUM_RUNS) / NUM_RUNS, )

stmt_no_parsing_ion_mass_dict = """
sequence=''.join(
    [random.choice(mass.STD_AMINO_ACID_MASS.keys()) for i in range(20)])
x = sum([mass.STD_AMINO_ACID_MASS[i] for i in sequence])
x += mass.STD_AMINO_ACID_MASS['A'] # emulate dict.get
"""
print 'str + dict.get + dict.get for ion type: %.6f usec/pass' % (
    timeit.timeit(stmt=stmt_no_parsing_ion_mass_dict,
                  setup='import mass, random; random.seed()',
                  number=NUM_RUNS) / NUM_RUNS, )

