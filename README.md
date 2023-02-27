# vascpp

This project is a collection of algorithms implemented in C++ to build large-scale vascular networks 
and analyze blood flow and radiation damage.

## Contact
Author: Maxwell Cole | Louisiana State University, Department of Physics & Astronomy, Medical and Health Physics Program

Date: 27 February 2023

# Usage
See rostam_build.sh for an example HPC build and launch script. 

### Dependencies
Dependencies can be installed with spack package manager (https://github.com/spack/spack) and linked in the master CMakeLists.txt file. The required packages are:

* OpenMPI
* HDF5
* ParMETIS
* Trilinos

### Build
To build the project, create a bash script containing the following information:
1. Find and load dependencies in spack or similar environment module
2. Create local build directory
3. Build project with cmake
Make any necessary adjustments in the CMakeLists.txt files.

### Launch
Command Line Inputs:
1. Step/Algorithm
2. Number of Generations
3. Number of Dimensions
4. Output File Name
Example: "srun app/vessel_gen 4 2 test.h5" will create a 2-dimensional network of 4 generations (30 vessels) and save the output data to "geometry_4_2.h5"

The executable can be run sequentially or distributed, however MPI is required either way (this may be changed in future implementations). 
The project was developed on HPC clusters utilizing the SLURM Workload Manager. For clusters using SLURM, submit batch scripts with "sbatch" or use "srun -N x ...", where "x" is the number of compute nodes. For personal computers, use the MPI process launcher from the command line (mpiexec or mpirun) with specified launch instructions.

The execution command can be added to the end of the build script to build and run from a single command.

### Output
The data is stored in HDF5 files (by default stored in the build directory). The save path can be changed in the command line inputs.

### Plot
2-dimensional plots can be generated from the HDF5 output file with the python script "plot_file.py" in the python_scripts directory. The file path can be changed with the "filename" variable.
![alt text](https://github.com/max-si/vascpp/blob/main/Plots/5-2-geometry.png?raw=true)

# Funding & Acknowledgments
The development of this software was supported by ...

## License
Boost Software License 1.0