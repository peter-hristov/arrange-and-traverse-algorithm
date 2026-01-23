#!/bin/bash

cd ./build
make -j 8
cd ..

GREEN="\033[0;32m"
RED="\033[0;31m"
RESET="\033[0m"


run_check() {
    if "$@" > /dev/null 2>&1; then
        echo -e "${GREEN}OK${RESET} : $*"
    else
        echo -e "${RED}NOT MATCHING${RESET} : $*"
    fi
}

# Now call run_check with your test commands
run_check ./build/fv99 -f ./data/three-sheet-toy.txt -u
run_check ./build/fv99 -f ~/Projects/data/reeb-space-test-data/data.vtu -u

# ---- Loop over all files in smallNeighbourhoods ----
for file in ~/Projects/data/reeb-space-test-data/smallNeighbourhoods/*; do
    #echo $file
    run_check ./build/fv99 -f "$file" -u
done

run_check ./build/fv99 -e 0.0 -f ~/Projects/data/reeb-space-test-data/ttk/downsample-20-300.vtu -u
run_check ./build/fv99 -e 0.0 -f ~/Projects/data/reeb-space-test-data/ttk/downsample-18-400.vtu -u
run_check ./build/fv99 -e 0.0 -f ~/Projects/data/reeb-space-test-data/ttk/downsample-15-875.vtu -u
run_check ./build/fv99 -e 0.0 -f ~/Projects/data/reeb-space-test-data/ttk/downsample-12-1440.vtu -u
run_check ./build/fv99 -e 0.0 -f ~/Projects/data/reeb-space-test-data/ttk/downsample-10-2800.vtu -u
# run_check ./build/fv99 -e 0.0 -f ~/Projects/data/reeb-space-test-data/ttk/downsample-8-5800.vtu -u

run_check ./build/fv99 -e 0.1 -f ~/Projects/data/reeb-space-test-data/torus/torus-factor-40-tets-625.vtu  -u
run_check ./build/fv99 -e 0.1 -f ~/Projects/data/reeb-space-test-data/torus/torus-factor-30-tets-1080.vtu -u
run_check ./build/fv99 -e 0.1 -f ~/Projects/data/reeb-space-test-data/torus/torus-factor-24-tets-2560.vtu -u
run_check ./build/fv99 -e 0.1 -f ~/Projects/data/reeb-space-test-data/torus/torus-factor-20-tets-2560.vtu -u
run_check ./build/fv99 -e 0.1 -f ~/Projects/data/reeb-space-test-data/torus/torus-factor-15-tets-10985.vtu -u
run_check ./build/fv99 -e 0.1 -f ~/Projects/data/reeb-space-test-data/torus/torus-factor-12-tets-20480.vtu -u
#run_check ./build/fv99 -e 0.1 -f ~/Projects/data/reeb-space-test-data/torus/torus-factor-10-tets-40000.vtu -u

run_check ./build/fv99 -f ~/Projects/data/reeb-space-test-data/isabel/isabel1.40.90.vtu -u
run_check ./build/fv99 -f ~/Projects/data/reeb-space-test-data/isabel/isabel1.30.240.vtu -u
run_check ./build/fv99 -f ~/Projects/data/reeb-space-test-data/isabel/isabel1.20.720.vtu -u

