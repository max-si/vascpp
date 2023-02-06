#!/bin/zsh

#SBATCH --job-name=vascpp
#SBATCH --output=vascpp_build.out
#SBATCH --partition=medusa
#SBATCH --nodes=1
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

# build with cmake
cmake ..
#cmake .. -DHIGHFIVE_PARALLEL_HDF5=ON -DHIGHFIVE_EXAMPLES=OFF -DHIGHFIVE_USE_BOOST=OFF
echo "Building . . ."
cmake --build .

# run program
srun vascpp 10 2


#valgrind srun vascpp 3 2
#mpiexec -n 2 vascpp 25 2
#python3 plot_file.py