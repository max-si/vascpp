#!/bin/zsh

# BASE_DIR=${pwd}
# BUILD_DIR=${BASE_DIR}/build
DIR="/Users/notmax/Desktop/VascularFiles/vascpp/build"

if [ -d "$DIR" ]; then 
    cd build
else
    echo "Directory 'build' not found"
    echo "Creating directory 'build' "
    mkdir build && cd build
fi

cmake ..
echo "Building . . ."
cmake --build .
echo "Executing build . . ."
mpiexec -n 2 vascpp 20 2
# h5dump build/test.h5
# python3 plot_file.py