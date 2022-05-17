#!/bin/bash
if [ $# -eq 0 ]; then
    echo "Installing Python versions..."
    xargs -a python-versions.txt -I{} pyenv install -s {}
    echo "Creating environments..."
    xargs -a python-versions.txt -I{} pyenv virtualenv {} tmp-pyteomics-{}
    echo "Installing dependencies..."
    while read pv; do
        eval "$(pyenv root)/versions/tmp-pyteomics-${pv}/bin/pip" install -U pip wheel
        eval "$(pyenv root)/versions/tmp-pyteomics-${pv}/bin/pip" install -r ../test-requirements.txt
        eval "$(pyenv root)/versions/tmp-pyteomics-${pv}/bin/pip" install -U ..
    done < python-versions.txt
    echo "Running tests..."
    find . -name 'test_*.py' -exec bash -c 'echo '"$(pyenv root)"'/versions/tmp-pyteomics-${0}/bin/python $0' {} \;
    echo "Deleting environments..."
    xargs -a python-versions.txt -I{} pyenv virtualenv-delete tmp-pyteomics-{}
else
    for f; do
        while read pv; do
            echo "Creating environment..."
            pyenv virtualenv "$pv" "tmp-pyteomics-single-${pv}"
            echo "Installing dependencies..."
            eval "$(pyenv root)/versions/tmp-pyteomics-single-${pv}/bin/pip" install -U pip wheel
            eval "$(pyenv root)/versions/tmp-pyteomics-single-${pv}/bin/pip" install -r ../test-requirements.txt
            eval "$(pyenv root)/versions/tmp-pyteomics-single-${pv}/bin/pip" install -U ..
            echo "Running test..."
            if [ -f "$f" ]; then
                fname="$f"
            elif [ -f "test_${f}.py" ]; then
                fname="test_${f}.py"
            fi
            echo "Executing $(pyenv root)/versions/tmp-pyteomics-single-${pv}/bin/python ${fname}"
            eval "$(pyenv root)/versions/tmp-pyteomics-single-${pv}/bin/python" "${fname}"
            echo "Deleting environment..."
            echo pyenv virtualenv-delete "tmp-pyteomics-single-${pv}"
            eval pyenv virtualenv-delete "tmp-pyteomics-single-${pv}"
        done < python-versions.txt
    done
fi
