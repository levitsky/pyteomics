#!/bin/bash
export PYTHONPATH=".."
if [ $# -eq 0 ]; then
    find . -name 'test_*.py' -exec bash -c 'declare -a versions=(2.7 3.3 3.4 3.5 3.6 3.7 3.8 3.9 3.10); for v in "${versions[@]}"; do command -v "python${v}" > /dev/null 2>&1 && { echo "Executing python${v}" "$0"; eval "python${v}" "$0"; }; done' {} \;
else
    for f; do
        for v in 2.7 3.3 3.4 3.5 3.6 3.7 3.8 3.9 3.10; do
            command -v "python${v}" >/dev/null 2>&1 && {
                if [ -f "$f" ]; then
                    fname="$f"
                elif [ -f "test_${f}.py" ]; then
                    fname="test_${f}.py"
                fi
                echo "Executing python${v}" "$fname"; eval "python${v}" "$fname"
            }
        done
    done
fi
