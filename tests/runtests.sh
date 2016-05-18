#!/bin/bash
export PYTHONPATH=".."
if [ $# -eq 0 ]; then
    find . -name 'test*.py' -exec echo "Executing python" {} \; -exec python {} \; -exec echo "Executing python2" {} \; -exec python2 {} \;
else
    for f; do
        echo "Executing python2" "$f"
        python2 "$f"
        echo "Executing python" "$f"
        python "$f"
    done
fi
