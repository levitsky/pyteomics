#!/bin/bash
export PYTHONPATH=".."
if [ $# -eq 0 ]; then
    find . -name 'test*.py' -exec echo "Executing" {} \; -exec python {} \; -exec python2 {} \;
else
    for f; do
        echo "Executing" "$f"
        python2 "$f"
        python "$f"
    done
fi
