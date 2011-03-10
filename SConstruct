import os
import inspect
import distutils.sysconfig

# Configuring the build.
#========================

Help("""
pyteomics build script. 
Usage: scons -Y <path_to_repository> [OPTIONS] [TARGETS]
Note that building in a source directory is strictly prohibited.

Options:
    buildtype=<buildtype>    Build type. Possible values are release and debug.

Targets:
    all                      All targets.

    pyteomics                Build prerequisites for pyteomics package. 
                             The package itself is built using setup.py.
                             Compiled by default.
    tests                    A test suite for pyteomics.
    doc                      Documentation.

    libbiolccc_shared        Shared BioLCCC library. Compiled by default.
    libbiolccc_static        Static BioLCCC library. 
    libbiolccc_examples      Examples for libbiolccc.
    libbiolccc_tests         A test suite for libbiolccc.
    libgtest_static          Static Google Test library.
""")

# Get the mode flag from the command line.
# Default to 'release' if the user didn't specify.
buildtype = ARGUMENTS.get('buildtype', 'release')
if buildtype not in ['debug', 'release']:
   print "Error: expected 'debug' or 'release', found: %s" % (buildtype,)
   Exit(1)

# It is strictly prohibited to build the application in the source directory.
build_in_repository = not bool(GetOption('repository'))
if not build_in_repository:
    for dir in GetOption('repository'):
        if Dir(dir).abspath == Dir('.').abspath:
            build_in_repository = True
if build_in_repository:
    print 'Error:'
    print 'Avoid building the application in the source directory.'
    print 'Print \'scons -Y source_directory\' in the build directory.'
    Exit(1)

VariantDir('.', GetOption('repository'), duplicate=True)

# Setting the platform specific options.
platform = ARGUMENTS.get('OS', Platform())
version = open(str(File('./VERSION').srcnode())).readline().strip()

ccflags = ' -DVERSION=\"%s\"' % version
if platform.name in ['posix', 'linux', 'unix']:
    if buildtype=='release':
        ccflags += ' -O2'
    else:
        ccflags += ' -Wall -g'
    tools = 'default'

if platform.name in ['win32', 'windows']:
    if buildtype=='release':
        ccflags += ' -O2'
    else:
        ccflags += ' -Wall -g'
    tools = 'mingw'

env = Environment(
    PLATFORM=platform.name,
    CPPPATH=[Dir(os.path.join('biolccc', 'include')).abspath,
             distutils.sysconfig.get_python_inc()],
    CCFLAGS=ccflags,
    tools=[tools],
    BUILDTYPE=buildtype,
    LIBPATH=[os.path.join('#lib', 'static', platform.name, buildtype),
             os.path.join('#lib', 'shared', platform.name, buildtype)],
    ROOTBUILDDIR=Dir('.').abspath,
    )

# Building targets.
#===================

# pyteomics package.
#-------------------
pyteomics = SConscript(
    os.path.join('biolccc', 'src', 'bindings', 'SConscript'),
    exports = {'env':env},
    variant_dir=os.path.join(
        'build', platform.name, env['BUILDTYPE'], 'bindings'), 
    duplicate=True,
    )

# Copying source files required for the python source package.
env.AddPostAction(
   pyteomics,
   Copy(os.path.join('build', platform.name, env['BUILDTYPE'], 'bindings',
                     'post_swig.py'),
        os.path.join(Dir('#.').abspath, 'biolccc', 'src', 
                     'bindings', 'post_swig.py')))
env.AddPostAction(
   pyteomics, 
   'python ' + os.path.join('build', platform.name, env['BUILDTYPE'],
                            'bindings', 'post_swig.py') + ' ' + version)
env.AddPostAction(
   pyteomics, 
   Copy(os.path.join(Dir('#.').abspath, 'biolccc', 'src', 
                     'bindings', 'biolccc.py'),
        os.path.join('build', platform.name, env['BUILDTYPE'], 'bindings',
                     'biolccc.py')))
env.AddPostAction(
   pyteomics, 
   Copy(os.path.join(Dir('#.').abspath, 'biolccc', 'src',
                     'bindings', 'biolccc_wrap.cc'),
        os.path.join('build', platform.name, env['BUILDTYPE'], 
                     'bindings', 'biolccc_wrap.cc')))

env.AddPostAction(
   pyteomics,
   Copy(os.path.join('build', platform.name, env['BUILDTYPE'], 'bindings',
                     '__init__.py'),
        os.path.join(Dir('#.').abspath, 'biolccc', 'src', 
                     'bindings', '__init__.py')))

env.AddPostAction(
   pyteomics, 
   Touch(os.path.join(Dir('#.').abspath, '__init__.py')))

# Copying the remaining files for pyteomics.
Depends(pyteomics, Glob(os.path.join('biolccc', 'src', 'core', '*.cpp')))
Depends(pyteomics, os.path.join('biolccc', 'src', 'bindings', 'post_swig.py'))
Depends(pyteomics, os.path.join('biolccc', 'src', 'bindings', '__init__.py'))
Depends(pyteomics, 'setup.py')
Depends(pyteomics, 'MANIFEST.in')
# Copying the documentation to the build dir.
Depends(pyteomics, 'VERSION')
Depends(pyteomics, 'README')
Alias('pyteomics', pyteomics)

# A test suite for pyteomics.
#----------------------------
tests = env.Install(
   os.path.join('build', platform.name, env['BUILDTYPE'], 'bindings'),
   os.path.join(Dir('#.').abspath, 'biolccc', 'src', 
                'bindings', 'test_biolccc.py'))
Depends(tests, 
        os.path.join(Dir('#.').abspath, 'biolccc', 'src', 
                     'bindings', 'test_biolccc.py'))
Alias('tests', tests)

# Shared BioLCCC library.
#-------------------------
libbiolccc_shared = SConscript(
    os.path.join('biolccc', 'src', 'core', 'SConscript'),
    exports = {'env':env, 'libtype':'shared'},
    variant_dir=os.path.join(
        'build', platform.name, env['BUILDTYPE'], 'core', 'shared'), 
    duplicate=True
    )
Alias('libbiolccc_shared', libbiolccc_shared)

# Static BioLCCC library.
#-------------------------
libbiolccc_static = SConscript(
    os.path.join('biolccc', 'src', 'core', 'SConscript'),
    exports = {'env':env, 'libtype':'static'},
    variant_dir=os.path.join(
        'build', platform.name, env['BUILDTYPE'], 'core', 'static'), 
    duplicate=True
    )
Alias('libbiolccc_static', libbiolccc_static)

# Google test library.
#----------------------
libgtest_static = SConscript(
    os.path.join('biolccc', 'src', 'gtest', 'SConscript'),
    exports = {'env':env},
    variant_dir=os.path.join(
        'build', platform.name, env['BUILDTYPE'], 'gtest'), 
    duplicate=True,
    )
Alias('libgtest_static', libgtest_static)

# A test suite for libbiolccc.
#-----------------------------
libbiolccc_tests = SConscript(
    os.path.join('biolccc', 'src', 'apps', 'SConscript'),
    exports={'env':env},
    variant_dir=os.path.join(
        'build', platform.name, env['BUILDTYPE'], 'apps'), 
    duplicate=True,
    )

Requires(tests, libbiolccc_static)
Requires(tests, libgtest_static)
Alias('libbiolccc_tests', libbiolccc_tests)

# Examples of libbiolccc usage.
#------------------------------
libbiolccc_examples = SConscript(
    os.path.join('biolccc', 'src', 'examples', 'SConscript'),
    exports={'env':env},
    variant_dir=os.path.join(
        'build', platform.name, env['BUILDTYPE'], 'examples'), 
    duplicate=True,
    )

Requires(libbiolccc_examples, libbiolccc_static)
Alias('libbiolccc_examples', libbiolccc_examples)

# Doxygen documentation for libbiolccc.
#-------------------------------------
libbiolccc_doc = env.Command('libbiolccc_doc', 'doc/Doxyfile', 
    [Mkdir('doc'), 'doxygen $SOURCE'])
Depends(libbiolccc_doc, 'doc/Doxyfile')
# Source code needs to be copied.
Depends(libbiolccc_doc, libbiolccc_shared)

# Sphinx documentation for pyteomics.
#------------------------------------
pyteomics_doc = env.Command('pyteomics_doc', '', 
    'mkdir ./doc/pyteomics/source/examples; '
    'cp src/examples/*.py ./doc/pyteomics/source/examples; '
    'cp src/examples/*.cpp ./doc/pyteomics/source/examples; '
    'cd ./doc/pyteomics; make html; cd ../')
Depends(pyteomics_doc, Glob(os.path.join('doc', 'pyteomics', '*')))
Depends(pyteomics_doc, Glob(os.path.join('doc', 'pyteomics', 'source',' *')))
Depends(pyteomics_doc, Glob(os.path.join('doc', 'pyteomics', 
                                         'source', '_static', '*')))
Depends(pyteomics_doc, Glob(os.path.join('doc', 'pyteomics', 'source', 
                                         '_sphinxext', '*.py')))
Depends(pyteomics_doc, Glob(os.path.join('doc', 'pyteomics', 'source',
                                         '_templates', '*')))
Depends(pyteomics_doc, Glob(os.path.join('src', 'examples', '*.py')))
Depends(pyteomics_doc, Glob(os.path.join('src', 'examples', '*.cpp')))
Depends(pyteomics_doc, 'VERSION')
Depends(pyteomics_doc, 'README')
Depends(pyteomics_doc, 'INSTALL')
Depends(pyteomics_doc, 'CHANGELOG')

# Complete documentation.
#------------------------
docs = env.Command('docs', '',
    #[''])
    ['mv doc/sphinx/build/html/* doc/',
     'mkdir doc/build',
     'mkdir doc/API',
     'mv doc/doxygen/html/* doc/API',
     'rm -r doc/doxygen',
     'rm -r doc/pyteomics'])
Requires('docs', libbiolccc_doc)
Requires('docs', pyteomics_doc)
# For some funny reason name 'doc' doesn't work, so we need to use an alias.
Alias('doc', docs)

# Final configuration of the build.
#===================================
env.Default([pyteomics])
Alias('all', 
      [libbiolccc_static, libbiolccc_shared, libgtest_static, pyteomics,
       tests, libbiolccc_tests, libbiolccc_examples, docs])

