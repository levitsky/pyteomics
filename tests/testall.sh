#!/bin/bash
export PYTHONPATH=".."
find . -name 'test*.py' -exec echo "Executing" {} \; -exec python {} \; -exec python2 {} \;
