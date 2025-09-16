# Compile

Here's how to compile the code with cmake, using CGAL version 6.0.1 and VTK 9.4.1

Build version
```
cmake -DCGAL_DIR=/home/peter/Projects/libraries/cgal/build -DCMAKE_PREFIX_PATH="/home/peter/Projects/libraries/VTK-9.4.1/install" -DCMAKE_EXPORT_COMPILE_COMMANDS=On ..

```

Release version
```
cmake -DCMAKE_PREFIX_PATH="/home/peter/Projects/libraries/cgal/install;/home/peter/Projects/libraries/VTK-9.4.1/install" -DCMAKE_EXPORT_COMPILE_COMMANDS=On -DCMAKE_BUILD_TYPE=Release ..
```


make


# Test

make -j 8 && ./fv99 -f ../data/three-sheet-toy.txt
make -j 8 && ./fv99 -f ~/Projects/data/reeb-space-test-data/data.vtu
make -j 8 && ./fv99 -e 0.1 -f ~/Projects/data/reeb-space-test-data/torus/torus-factor-50-tets-320.vtu
make -j 8 && ./fv99 -e 0.0 -f ~/Projects/data/reeb-space-test-data/ttk/downsample-20-300.vtu


# Big Unit Test
bash ./unitTest.sh ./build/fv99 -f ./data/three-sheet-toy.txt -u
bash ./unitTest.sh ./build/fv99 -f ~/Projects/data/reeb-space-test-data/data.vtu -u

bash ./unitTest.sh ./build/fv99 -e 0.0 -f ~/Projects/data/reeb-space-test-data/ttk/downsample-20-300.vtu -u
bash ./unitTest.sh ./build/fv99 -e 0.0 -f ~/Projects/data/reeb-space-test-data/ttk/downsample-18-400.vtu -u
bash ./unitTest.sh ./build/fv99 -e 0.0 -f ~/Projects/data/reeb-space-test-data/ttk/downsample-15-875.vtu -u
bash ./unitTest.sh ./build/fv99 -e 0.0 -f ~/Projects/data/reeb-space-test-data/ttk/downsample-12-1440.vtu -u
bash ./unitTest.sh ./build/fv99 -e 0.0 -f ~/Projects/data/reeb-space-test-data/ttk/downsample-10-2800.vtu -u

bash ./unitTest.sh ./build/fv99 -e 0.1 -f ~/Projects/data/reeb-space-test-data/torus/torus-factor-40-tets-625.vtu  -u
bash ./unitTest.sh ./build/fv99 -e 0.1 -f ~/Projects/data/reeb-space-test-data/torus/torus-factor-30-tets-1080.vtu -u
bash ./unitTest.sh ./build/fv99 -e 0.1 -f ~/Projects/data/reeb-space-test-data/torus/torus-factor-20-tets-2560.vtu -u
bash ./unitTest.sh ./build/fv99 -e 0.1 -f ~/Projects/data/reeb-space-test-data/torus/torus-factor-15-tets-10985.vtu -u

bash ./unitTest.sh ./build/fv99 -f ~/Projects/data/reeb-space-test-data/isabel/isabel1.40.90.vtu -u
bash ./unitTest.sh ./build/fv99 -f ~/Projects/data/reeb-space-test-data/isabel/isabel1.30.240.vtu -u
bash ./unitTest.sh ./build/fv99 -f ~/Projects/data/reeb-space-test-data/isabel/isabel1.20.720.vtu -u





