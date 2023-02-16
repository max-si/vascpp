#!/bin/zsh

#SBATCH --job-name=vascpp
#SBATCH --output=vascpp_build.out
#SBATCH --partition=medusa
#SBATCH --nodes=2
#SBATCH --time=30:00

# set xe

# load modules
module load gcc
module load openmpi/4.1.4
module load hdf5

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
srun app/vessel_gen 3 2 


#valgrind srun vascpp 3 2
#mpiexec -n 2 vascpp 25 2
#python3 plot_file.py