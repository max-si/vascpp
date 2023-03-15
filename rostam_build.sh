#!/bin/zsh

#SBATCH --job-name=vascpp
#SBATCH --output=vascpp_build.out
#SBATCH --partition=medusa
#SBATCH --nodes=4
#SBATCH --time=30:00

# set xe

# load modules
module load gcc openmpi hdf5 python

# load spack env and packages
cd /work/maxwell/spack
. share/spack/setup-env.sh
spack load trilinos
spack load metis
spack load parmetis
cd /work/maxwell/vascpp

# check if build dir exists, if not then create it
DIR="/work/maxwell/vascpp/build"
if [ -d "$DIR" ]; then 
    cd build
else
    echo "Directory 'build' not found"
    echo "Creating directory 'build' "
    mkdir build && cd build
fi

METIS_DIR=$(spack find -p metis | egrep "/work.*" -o)
CPATH=$CPATH:$METIS_DIR/include

# build with cmake
cmake ..
echo "Building . . ."
cmake --build .

# run program
srun app/vessel_gen 4 2