#!/bin/bash

export ROOT="$(pyenv root)"
export MODE="single"
export NODEPS=0
export PREFIX="test_pyteomics_"
export CLEAN=0
export PYTEOMICS_TEST_POSTGRES_URL="${PYTEOMICS_TEST_POSTGRES_URL:-}"
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
}

detect_postgres() {
    local probe_py
    probe_py="${ROOT}/versions/${PREFIX}$(head -n1 python-versions.txt)/bin/python"

    if [[ -z "$PYTEOMICS_TEST_POSTGRES_URL" ]] && command -v psql >/dev/null 2>&1; then
        if psql -d postgres -c 'select 1' >/dev/null 2>&1; then
            export PYTEOMICS_TEST_POSTGRES_URL="postgresql+psycopg:///postgres"
            echo "Detected local PostgreSQL via psql; enabling live Unimod PostgreSQL test."
        fi
    fi

    if [[ -z "$PYTEOMICS_TEST_POSTGRES_URL" ]]; then
        echo "PYTEOMICS_TEST_POSTGRES_URL is not set; live PostgreSQL Unimod test will be skipped."
        return 0
    fi

    if [[ ! -x "$probe_py" ]]; then
        echo "Could not find probe interpreter: $probe_py"
        return 0
    fi

    if ! "$probe_py" -c "import sys; from sqlalchemy import create_engine, text; engine = create_engine(sys.argv[1]);
with engine.connect() as conn:
    conn.execute(text('SELECT 1'))
engine.dispose()" "$PYTEOMICS_TEST_POSTGRES_URL"
    then
        echo "Could not connect to PYTEOMICS_TEST_POSTGRES_URL; disabling live PostgreSQL test."
        unset PYTEOMICS_TEST_POSTGRES_URL
    else
        echo "Live PostgreSQL Unimod test enabled via PYTEOMICS_TEST_POSTGRES_URL=$PYTEOMICS_TEST_POSTGRES_URL"
    fi
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
    while read pv; do
        echo "${ROOT}/versions/${PREFIX}${pv}/bin/pip" install -U ..
        eval "${ROOT}/versions/${PREFIX}${pv}/bin/pip" install -U ..
    done < python-versions.txt
fi

detect_postgres

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
