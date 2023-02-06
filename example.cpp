
#include <mpi.h>
#include <string>
#include <iostream>
#include <vector>
#include <cmath>

#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5File.hpp>
#include <highfive/H5PropertyList.hpp>
#include <highfive/H5Attribute.hpp>

void ExportInitial(std::string filename, long long totalNumToExport) {
	int mpiSize = 1;
    int mpiRank = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);

    long long initialNumToExport = 2;

	// populate vector (contents don't matter)
	std::vector<std::vector<double>> initialVector(initialNumToExport, std::vector<double>(7, 0));
	for(int i = 0; i < initialNumToExport; i++) {
		for(int j = 0; j < 7; j++) {
			initialVector[i][j] = 3.14159265;
		}
	}

	using namespace HighFive;
	try {
		// only export from rank 0
		if(!mpiRank) {
			// create or overwrite file with filename
			File file(filename, File::Truncate);

			// create dataspace with initial and max shape
			DataSpace d_space = DataSpace({std::size_t(initialNumToExport), 7}, {std::size_t(totalNumToExport), 7});

			// add chunking
			DataSetCreateProps props;
			props.add(Chunking(std::vector<hsize_t>{std::size_t(1), 7}));
		
			// create datasets
			DataSet dataset = file.createDataSet("example_dset", d_space, create_datatype<double>(), props);

			// write initial part of datasets -- root vessels
			dataset.select({0,0}, {std::size_t(initialNumToExport), 7}).write(initialVector);

			file.flush();
		}
	} catch (Exception& err) {
        // catch and print any HDF5 error
        std::cerr << err.what() << std::endl;
    }
	initialVector.clear();
}


void ExportParallel(std::string filename, long long totalNumToExport) {
	int mpiSize = 1;
    int mpiRank = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);

    long long initialNumToExport = 2;
	long long numRowsOnRank = (totalNumToExport - initialNumToExport) / 2;

	if(!mpiRank) {
		std::cout << "totalNumToExport = " << totalNumToExport << std::endl;
		std::cout << "numRowsOnRank = " << numRowsOnRank << std::endl;
	}

	using namespace HighFive;
	{
		// initiate collective IO
		FileAccessProps fapl;
        fapl.add(MPIOFileAccess{MPI_COMM_WORLD, MPI_INFO_NULL});
		fapl.add(MPIOCollectiveMetadata{});
		auto xfer_props = DataTransferProps{};
        xfer_props.add(UseCollectiveIO{});

		// reopen file
		File file(filename, File::ReadWrite, fapl);

		// retrieve dataset and resize
		DataSet dset = file.getDataSet("example_dset");
		dset.resize({std::size_t(totalNumToExport), 7});

		// populate new vector on each MPI Rank
		std::vector<std::vector<double>> newData(numRowsOnRank, std::vector<double>(7, 0));
		for(int i = 0; i < numRowsOnRank; i++) {
			for(int j = 0; j < 7; j++) {
				newData[i][j] = 3.14159265;
			}
		}
		
		// define offset and count for each MPI Rank
		std::vector<size_t> offset{std::size_t(initialNumToExport + mpiRank * numRowsOnRank), 0ul};
        std::vector<size_t> count{std::size_t(numRowsOnRank), 7ul};

		dset.select(offset, count).write(newData, xfer_props);

		// free memory from vector
		newData.clear();
		file.flush();
	}
}

int main(int argc, char* argv[]) {

	MPI_Init(NULL, NULL);
	int mpiSize = 1;
    int mpiRank = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);

	std::string filename = "test.h5";
	long long totalNumToExport = pow(2, 20) - 2;

	ExportInitial(filename, totalNumToExport);

	if(!mpiRank) {
		std::cout << "Initial Export Complete" << std::endl;
	}
	
	ExportParallel(filename, totalNumToExport);

	if(!mpiRank) {
		std::cout << "Parallel Export Complete. Finalizing MPI." << std::endl;
	}

	MPI_Finalize();

	return 0;
}