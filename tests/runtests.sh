#!/bin/bash

export ROOT="$(pyenv root)"
export MODE="single"
export NODEPS=0
export PREFIX="test_pyteomics_"
export CLEAN=0
JOBLOG="parallel.log"
consequtive=("unimod")
declare -A exitcodes

rm -f "$JOBLOG"

for i in "$@" ; do
    if [[ "$i" == "--all" ]] ; then
        MODE="all"
    elif [[ "$i" == "--no-run" ]] ; then
        MODE="no-op"
    fi
    if [[ "$i" == "--clean" ]] ; then
        CLEAN=1
    fi
    if [[ "$i" == "--no-deps" ]] ; then
        NODEPS=1
    fi
done

create() {
    if [[ $CLEAN == 1 ]] || [ ! -d "${ROOT}/versions/${PREFIX}${1}" ]; then
        echo "Creating environment ${PREFIX}${1} ..."
        echo pyenv virtualenv -f "$1" "${PREFIX}${1}"
        pyenv virtualenv -f "$1" "${PREFIX}${1}"
    fi
}

deps() {
    eval "${ROOT}/versions/${PREFIX}${1}/bin/pip" install -U pip wheel
    eval "${ROOT}/versions/${PREFIX}${1}/bin/pip" install -U -r ../test-requirements.txt
    eval "${ROOT}/versions/${PREFIX}${1}/bin/pip" install -U ..
}

run() {
    if [ -f "$2" ]; then
        fname="$2"
    elif [ -f "test_${2}.py" ]; then
        fname="test_${2}.py"
    else
        return 0
    fi
    echo "${ROOT}/versions/${PREFIX}${1}/bin/python" "$fname"
    eval "${ROOT}/versions/${PREFIX}${1}/bin/python" "$fname"
}

export -f create deps run

echo "Checking installed Python versions..."
parallel pyenv install -s < python-versions.txt
if [[ $CLEAN == 1 ]] ; then
    echo "Deleting old environments..."
    xargs -a python-versions.txt -I{} echo pyenv virtualenv-delete -f "${PREFIX}"{}
    xargs -a python-versions.txt -I{} pyenv virtualenv-delete -f "${PREFIX}"{}
fi

parallel create :::: python-versions.txt

if [[ $NODEPS == 0 ]] ; then
    echo "Installing dependencies..."
    parallel deps :::: python-versions.txt
fi

if [[ $MODE == "all" ]] ; then
    files=$(find . -name 'test_*.py' -print)
elif [[ "$MODE" == "single" ]] ; then
    files=$@
fi

if [[ "$MODE" != "no-op" ]] ; then
    parallel_files=()
    consecutive_files=()
    for f in ${files[@]}; do
        if [[ "$f" =~ ^-- ]] ; then
            continue
        fi
        cons_flag=0
        for check in ${consequtive[@]}; do
            if [[ "$f" =~ "$check" ]] ; then
                cons_flag=1
                break
            fi
        done
        if [[ $cons_flag == 1 ]]; then
            consecutive_files+=("$f")
        else
            parallel_files+=("$f")
        fi
    done

    if [[ ${#consecutive_files[@]} != 0 ]] ; then
        echo "Running tests consecutively:" ${consecutive_files[@]}
        for f in ${consecutive_files[@]}; do
            while read pv; do
                run "$pv" "$f"
                exitcodes["run ${pv} ${f}"]=$?
            done < python-versions.txt
        done
    fi

    if [[ ${#parallel_files[@]} != 0 ]] ; then
        echo "Running tests in parallel:" "${parallel_files[@]}"
        parallel --joblog "$JOBLOG" run :::: python-versions.txt ::: ${parallel_files[@]}
    fi

    total=${#exitcodes[@]}
    nerrors=0
    if [ -f "$JOBLOG" ] ; then
        total=$(( $total + $( wc -l < "$JOBLOG" ) - 1 ))

        errors=$(awk -F '\t' '
            NR==1   { next }
            $7 != 0 { sum += 1; printf "exit code: %d; job %s\t(parallel)\n", $7, $9 }' "$JOBLOG")
        echo "$errors"
        if [[ "$errors" != "" ]]; then
            nerrors=$(echo "$errors" | wc -l)
        fi
    fi
    conseq_errors=0
    for key in "${!exitcodes[@]}"; do
        if [[ ${exitcodes["$key"]} != 0 ]]; then
            (( conseq_errors++ ))
            printf "exit code: %d; job %s\t(consequtive)\n" ${exitcodes["$key"]} "$key"
        fi
    done
    nerrors=$(( $nerrors + $conseq_errors ))
    if [[ $nerrors == 0 ]] ; then
        echo "All ${total} test runs complete, there were no errors."
    else
        echo "All ${total} test runs complete, there were ${nerrors} errors (see above)."
    fi
fi
