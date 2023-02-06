
#include <cmath>
#include <vector>
#include <iostream>
#include <mpi.h>
#include <string>

#include "vessel.h"
#include "coordinate.h"
#include "IOPrep.h"

#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5File.hpp>
#include <highfive/H5PropertyList.hpp>
#include <highfive/H5Attribute.hpp>

typedef std::vector<Vessel> VesselVector;

void ConvertVesselToGeometryRow_HF(Vessel& vessel, std::vector<std::vector<double>>& geomArray, long long& geomTableIndex) 
{
	const Coordinate& start = vessel.get_startingPoint();
	const Coordinate& end = vessel.get_endingPoint();

	geomArray[geomTableIndex][0] = start.x;
    geomArray[geomTableIndex][1] = start.y;
    geomArray[geomTableIndex][2] = start.z;
    geomArray[geomTableIndex][3] = end.x;
    geomArray[geomTableIndex][4] = end.y;
    geomArray[geomTableIndex][5] = end.z;
    geomArray[geomTableIndex][6] = vessel.get_radius();
}

void ConvertVesselToNodeRow_HF(Vessel& vessel, std::vector<std::vector<int>>& nodeArray,long long& nodeTableIndex)
{
    nodeArray[nodeTableIndex][0] = vessel.get_nodes()[0];
    nodeArray[nodeTableIndex][1] = vessel.get_nodes()[1];
}

// Function that serializes the start pt, end pt, and radius of a vessel object
void ConvertVesselVectorsToGeometryArray_HF(long long numVesselsToWrite, VesselVector& vessels, std::vector<std::vector<double>>& geomArray)
{
    long long rowIndex = 0, geomTableIndex = 0;
    VesselVector::size_type i;

    for (i = 0; i < numVesselsToWrite; i++)
    {
        // geomTableIndex = i * 7;
		// ConvertVectorToGeometryRow(vessels[i], geomArray, geomTableIndex);
        ConvertVesselToGeometryRow_HF(vessels[i], geomArray, rowIndex);
		rowIndex += 1;
    } 
}

// Function to serialize node identifiers in vessel objects
void ConvertVesselVectorsToNodeArray_HF(long long numVesselsToWrite,
    VesselVector& arteries, VesselVector& veins, std::vector<std::vector<int>>& nodeArray)
{
    long long rowIndex = 0, nodeTableIndex = 0;
    VesselVector::size_type i;

    for (i = 0; i < numVesselsToWrite; i++)
    {
        //nodeTableIndex = i * 2;
        ConvertVesselToNodeRow_HF(arteries[i], nodeArray, rowIndex);
		rowIndex += 1;
    }

    for (i = 0; i < numVesselsToWrite; i++)
    {
        //nodeTableIndex = (numVesselsToWrite + i) * 2;
        ConvertVesselToNodeRow_HF(veins[i], nodeArray, rowIndex);
		rowIndex += 1;
    }
}

// Function to serialize vessel data for root vessels in parallel algorithm
void ConvertRootVesselVectorsToDataTables_HF(long long numVesselsToWrite, VesselVector& arteries, 
	VesselVector& veins, std::vector<std::vector<double>>& geomArray, std::vector<std::vector<int>>& nodeArray) 
{
	long long rowIndex = 0;
    VesselVector::size_type i;
	
    for (i = 0; i < numVesselsToWrite; i++)
    {
		ConvertVesselToGeometryRow_HF(arteries[i], geomArray, rowIndex);
        ConvertVesselToNodeRow_HF(arteries[i], nodeArray, rowIndex);
		rowIndex += 1;
    }
	for (i = 0; i < numVesselsToWrite; i++)
    {
		ConvertVesselToGeometryRow_HF(veins[i], geomArray, rowIndex);
        ConvertVesselToNodeRow_HF(veins[i], nodeArray, rowIndex);
		rowIndex += 1;
    }
}

// Function to serialize vessel data for root vessels in sequential algorithm
void ConvertSequentialVesselVectorsToDataTables_HF(long long numVesselsToWrite, VesselVector& vessels, 
			std::vector<std::vector<double>>& geomArray, std::vector<std::vector<int>>& nodeArray) 
{
	long long rowIndex = 0;
    VesselVector::size_type i;
	
    for (i = 0; i < numVesselsToWrite; i++)
    {
		ConvertVesselToGeometryRow_HF(vessels[i], geomArray, rowIndex);
        ConvertVesselToNodeRow_HF(vessels[i], nodeArray, rowIndex);
		rowIndex += 1;
    }
}

//*****
long long PopulateNumConnectedVesselsArray_HF(VesselVector& arteries, VesselVector& veins, 
	std::vector<short>& numConnectedVesselsArray, long long numVessels) 
{
	long long rowIndex = 0, tableIndex = 0;
    long long totalNumberOfConnectedVessels = 0;
    VesselVector::size_type i;
    for (i = 0; i < numVessels; i++, tableIndex++)
    {
        short numConnectedToThisVessel = arteries[i].get_numConnectedVessels();
        numConnectedVesselsArray[tableIndex] = numConnectedToThisVessel;
        totalNumberOfConnectedVessels += numConnectedToThisVessel;
    }
	//!for (i = 0; i < numVessels; i++, tableIndex++)
	for (i = 0; i < veins.size(); i++, tableIndex++)
    {
        short numConnectedToThisVessel = veins[i].get_numConnectedVessels();
        numConnectedVesselsArray[tableIndex] = numConnectedToThisVessel;
        totalNumberOfConnectedVessels += numConnectedToThisVessel;
    }
	return totalNumberOfConnectedVessels;
}

//****
void PopulateConnectedVesselArray_HF(long long numVesselsToWrite, VesselVector& arteries,
    VesselVector& veins, std::vector<long long>& connectedVesselArray)
{
    unsigned long long tableIndex = 0;
    VesselVector::size_type i;
    for (i = 0; i < numVesselsToWrite; i++) {
        const std::vector<long long>& connectedVessels = arteries[i].get_connectedVessels();
        for (std::vector<long long>::const_iterator itr = connectedVessels.begin(); itr != connectedVessels.end(); itr++) {
            connectedVesselArray[tableIndex++] = *itr;
        }
    }
    //!for (i = 0; i < veins.size(); i++) 
	for (i = 0; i < numVesselsToWrite; i++)
	{
        const std::vector<long long>& connectedVessels = veins[i].get_connectedVessels();
        for (std::vector<long long>::const_iterator itr = connectedVessels.begin(); itr != connectedVessels.end(); itr++) {
            connectedVesselArray[tableIndex++] = *itr;
        }
    }
}

//**** 
long long PopulateRootNumConnectedVesselsArray_HF(VesselVector& arteries, VesselVector& veins, 
	short* numConnectedVesselsArray, long long numVessels)
{
    long long rowIndex = 0, tableIndex = 0;
    long long totalNumberOfConnectedVessels = 0;
    VesselVector::size_type i;
    for (i = 0; i < numVessels; i++, tableIndex++)
    {
        short numConnectedToThisVessel = arteries[i].get_numConnectedVessels();
        numConnectedVesselsArray[tableIndex] = numConnectedToThisVessel;
        totalNumberOfConnectedVessels += numConnectedToThisVessel;
    }
    //!for (i = numVessels + 1; i < veins.size(); i++, tableIndex++)
	for (i = 0; i < numVessels; i++, tableIndex++)
    {
        short numConnectedToThisVessel = veins[i].get_numConnectedVessels();
        numConnectedVesselsArray[tableIndex] = numConnectedToThisVessel;
        totalNumberOfConnectedVessels += numConnectedToThisVessel;
    }
	return totalNumberOfConnectedVessels;
}

//**** PopulateRootConnectedVesselsArray
void PopulateRootConnectedVesselArray_HF(long long numVesselsToWrite, VesselVector& arteries,
    VesselVector& veins, long long* connectedVesselArray)
{
    unsigned long long tableIndex = 0;
    VesselVector::size_type i;
    for (i = 0; i < numVesselsToWrite; i++) {
        const std::vector<long long>& connectedVessels = arteries[i].get_connectedVessels();
        for (std::vector<long long>::const_iterator itr = connectedVessels.begin(); itr != connectedVessels.end(); itr++) {
            connectedVesselArray[tableIndex++] = *itr;
        }
    }
	//! Donahue has wrong indices here!
    for (i = 0; i < numVesselsToWrite; i++) {
        const std::vector<long long>& connectedVessels = veins[i].get_connectedVessels();
        for (std::vector<long long>::const_iterator itr = connectedVessels.begin(); itr != connectedVessels.end(); itr++) {
            connectedVesselArray[tableIndex++] = *itr;
        }
	}
}

// *****
void WriteSubConnectedVesselDataToFile_HF(std::string filename, int totalLevels, long long numVesselsExport,
    long long rootNumVessels, long long totalNumConnectedVessels,
    std::vector<short>& numConnectedVesselsArray, std::vector<long long>& connectedVesselArray)
{
	int mpiSize = 1;
    int mpiRank = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);

	long long connectedVesselsSize = pow(2, totalLevels) * 7 - 12;
	long long totalNumVessels = pow(2, totalLevels + 1) - 2;

	using namespace HighFive; {
		// initiate collective IO
		FileAccessProps fapl;
        fapl.add(MPIOFileAccess{MPI_COMM_WORLD, MPI_INFO_NULL});
		fapl.add(MPIOCollectiveMetadata{});
		auto xfer_props = DataTransferProps{};
        xfer_props.add(UseCollectiveIO{});

		// reopen file
		File file(filename, File::ReadWrite, fapl);

		// retreive connected vessels datasets and resize
		DataSet numConnections_dset = file.getDataSet("GEOM_GROUP/NUM_CONNECTED_VESSELS_ARRAY");
		numConnections_dset.resize({std::size_t(totalNumVessels), 1});
		DataSet connections_dset = file.getDataSet("GEOM_GROUP/CONNECTED_VESSELS_ARRAY");
		connections_dset.resize({std::size_t(connectedVesselsSize), 1});

		// write sub-vessels numConnectedVessels data to h5 file
		std::vector<size_t> offset{std::size_t(rootNumVessels + mpiRank * numVesselsExport), 0ul};
        std::vector<size_t> count{std::size_t(numVesselsExport), 1ul};
		numConnections_dset.select(offset, count).write(numConnectedVesselsArray, xfer_props);

		// determine number of blocks for successful write -- 4GB per write
		const long maxNumVesselsPerWrite = 1024 * 1024;	// from Donahue - 1,048,576
		const long long remainingNumberOfVesselsToExport = totalNumConnectedVessels % maxNumVesselsPerWrite;
		const int numLoopsNeeded = totalNumConnectedVessels / maxNumVesselsPerWrite;

		// write sub-vessels connectedVessels data to h5 file
		offset[0] = std::size_t(rootNumVessels * 4 - 4 + totalNumConnectedVessels * mpiRank);
		count[0] = std::size_t(totalNumConnectedVessels);
		connections_dset.select(offset, count).write(connectedVesselArray, xfer_props);
		
		// write excess first
		if (remainingNumberOfVesselsToExport > 0) {
			count[0] = std::size_t(remainingNumberOfVesselsToExport);
			connections_dset.select(offset, count).write(connectedVesselArray, xfer_props);
		}

		// loop for maximum write size
		count[0] = std::size_t(maxNumVesselsPerWrite);
		for (int i = 0; i < numLoopsNeeded; i++) {
			offset[0] = std::size_t(rootNumVessels * 4 - 4 + totalNumConnectedVessels * mpiRank +
            		i * count[0] + remainingNumberOfVesselsToExport);
			connections_dset.select(offset, count).write(connectedVesselArray, xfer_props);
		}

		file.flush();
	}
}

// ******
void WriteSubVesselsToFile_HF(std::string filename, long long totalNumVessels, long long subArterialTreeNumberOfVessels,
    long long rootArteriealTreeNumberOfVessels,
    VesselVector& arteries, VesselVector& veins)
{
	long long totalNumberOfVesselsToExport = 2 * subArterialTreeNumberOfVessels;
    long long totalNumberOfRootVessels = 2 * rootArteriealTreeNumberOfVessels;

    int mpiSize = 1;
    int mpiRank = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);

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

		// retrieve geom dataset, resize, and populate geom array
		DataSet geom_dset = file.getDataSet("GEOM_GROUP/GEOM_ARRAY");
		geom_dset.resize({std::size_t(totalNumVessels), 7});
		std::vector<std::vector<double>> geomVector(subArterialTreeNumberOfVessels, std::vector<double>(7, 0));
    	ConvertVesselVectorsToGeometryArray_HF(subArterialTreeNumberOfVessels, arteries, geomVector);

		// determine number of blocks for successful write -- 4GB per write
		const long maxNumVesselsPerWrite = 1024 * 1024;	// from Donahue - 1,048,576
		const long long remainingNumberOfVesselsToExport = subArterialTreeNumberOfVessels % maxNumVesselsPerWrite;
		const int numLoopsNeeded = subArterialTreeNumberOfVessels / maxNumVesselsPerWrite;

		// std::cout << "\nRemaining Number Of Vessels To Export = " << remainingNumberOfVesselsToExport << std::endl;
		// std::cout << "totalNumberOfVesselsToExport = " << totalNumberOfVesselsToExport << std::endl;
		// std::cout << "totalNumberOfRootVessels = " << totalNumberOfRootVessels << std::endl;
		// std::cout << "numLoopsNeeded = " << numLoopsNeeded << std::endl;

		// write sub-arteries geom data to h5 file
		// offset - starting location for the hyperslab
		// count - # of elements in the hyperslab selection
		std::vector<size_t> offset{std::size_t(totalNumberOfRootVessels + mpiRank * totalNumberOfVesselsToExport), 0ul};
        std::vector<size_t> count{std::size_t(remainingNumberOfVesselsToExport), 7ul};

		//! geomVector is too large, needs to be broken up
		//! may need to use regular hdf5 or have seperate dataset for Sub Vessels

		// arteries - write extra vessels first
		if(remainingNumberOfVesselsToExport > 0) {
			count[0] = remainingNumberOfVesselsToExport;
			offset[0] = std::size_t(totalNumberOfRootVessels + mpiRank * totalNumberOfVesselsToExport);
			geom_dset.select(offset, count).write(geomVector, xfer_props);
		}

		//std::cout << "---------CHECKPOINT--------" << std::endl;

		// arteries - loop for maximum write size
		count[0] = std::size_t(maxNumVesselsPerWrite);
		for (int i = 0; i < numLoopsNeeded; i++) {
			offset[0] = std::size_t(totalNumberOfRootVessels + mpiRank * totalNumberOfVesselsToExport + 
					i * count[0] + remainingNumberOfVesselsToExport);

			geom_dset.select(offset, count).write(geomVector, xfer_props);
		}

		//! repeat for sub-veins
		ConvertVesselVectorsToGeometryArray_HF(subArterialTreeNumberOfVessels, veins, geomVector);
		offset[0] = std::size_t(totalNumberOfRootVessels + mpiRank * totalNumberOfVesselsToExport + subArterialTreeNumberOfVessels);
		if (remainingNumberOfVesselsToExport > 0) {
			count[0] = std::size_t(remainingNumberOfVesselsToExport);
			offset[0] = std::size_t(totalNumberOfRootVessels + mpiRank * totalNumberOfVesselsToExport + subArterialTreeNumberOfVessels);
			geom_dset.select(offset, count).write(geomVector, xfer_props);
		}

		//! veins - loop for max write size 
		count[0] = std::size_t(maxNumVesselsPerWrite);
		for (int i = 0; i < numLoopsNeeded; i++) {
			offset[0] = std::size_t(totalNumberOfRootVessels + mpiRank * totalNumberOfVesselsToExport + 
					i * count[0] + remainingNumberOfVesselsToExport + subArterialTreeNumberOfVessels);
			geom_dset.select(offset, count).write(geomVector, xfer_props);
		}

		// free memory from geometry vector
		geomVector.clear();

		// retrieve node dataset and populate node vector
		DataSet node_dset = file.getDataSet("GEOM_GROUP/NODE_ARRAY");
		node_dset.resize({std::size_t(totalNumVessels), 2});
		std::vector<std::vector<int>> nodeVector(totalNumberOfVesselsToExport, std::vector<int>(2, 0));
		ConvertVesselVectorsToNodeArray_HF(totalNumberOfVesselsToExport / 2, arteries, veins, nodeVector);
		
		// write node data to file
		offset[0] = std::size_t(totalNumberOfRootVessels + mpiRank * totalNumberOfVesselsToExport);
		count[0] = totalNumberOfVesselsToExport;
		count[1] = 2ul;
		node_dset.select(offset, count).write(nodeVector, xfer_props);
		nodeVector.clear();
		
		file.flush();
	} // namespace HighFive;
}

// *****
//? combine with WriteSubVesselsToFile?
void ExportSubVesselsToHdf5_HF(std::string filename, int totalLevels, int subLevels,
    VesselVector& arteries, VesselVector& veins) 
{
	int mpiSize = 1;
	int mpiRank = 0;
	MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
	MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);

	long long totalNumVessels = pow(2, totalLevels + 1) - 2;
    long long subNumVessels = pow(2, subLevels + 1) - 1;
    long long rootNumVessels = pow(2, totalLevels - subLevels - 1) - 1;

	// write sub vessels to file
	WriteSubVesselsToFile_HF(filename, totalNumVessels, subNumVessels, rootNumVessels, arteries, veins);

	// populate num connected vessels array
	std::vector<short> numConnectedVesselsArray(2 * subNumVessels);
	long long totalNumConnectedVessels = PopulateNumConnectedVesselsArray_HF(
        arteries, veins, numConnectedVesselsArray, subNumVessels);

	// populate connected vessels array
	std::vector<long long> connectedVesselArray(totalNumConnectedVessels);
	PopulateConnectedVesselArray_HF(subNumVessels, arteries, veins, connectedVesselArray);
	
	// Write Sub Connected Vessel Data To File
	//WriteSubConnectedVesselDataToFile(filename, totalLevels, 2 * subNumVessels, 2 * rootNumVessels, totalNumConnectedVessels, numConnectedVesselsArray, connectedVesselArray);
}

// *****
void ExportRootVesselsToHdf5_HF(std::string filename, int totalLevels, int rootLevels,
    VesselVector& arteries,
    VesselVector& veins)
{
	int mpiSize = 1;
    int mpiRank = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    
	long long totalNumVessels = pow(2, totalLevels + 1) - 2;
    long long rootNumVessels = pow(2, rootLevels - 1) - 1;

    // populate arrays
	std::vector<std::vector<double>> geomArray(2*rootNumVessels, std::vector<double>(7, 0));
	std::vector<std::vector<int>> nodeArray(2*rootNumVessels, std::vector<int>(2, 0));
    ConvertRootVesselVectorsToDataTables_HF(rootNumVessels, arteries, veins, geomArray, nodeArray);

	// write data to HDF5
	using namespace HighFive;
	try {
		if(!mpiRank) {
			File file(filename, File::ReadWrite);

			// create dataspaces with initial and max shape
			DataSpace geom_space = DataSpace({std::size_t(2*rootNumVessels), 7}, {std::size_t(totalNumVessels), 7});
			DataSpace node_space = DataSpace({std::size_t(2*rootNumVessels), 2}, {std::size_t(totalNumVessels), 2});

			// add chunking
			DataSetCreateProps gprops, nprops;
			gprops.add(Chunking(std::vector<hsize_t>{std::size_t(totalNumVessels), 7}));
			nprops.add(Chunking(std::vector<hsize_t>{std::size_t(totalNumVessels), 2}));
			// gprops.add(Chunking(std::vector<hsize_t>{std::size_t(1), 7}));
			// nprops.add(Chunking(std::vector<hsize_t>{std::size_t(1), 2}));
		
			// create datasets
			DataSet geom_set = file.createDataSet("GEOM_GROUP/GEOM_ARRAY", geom_space, create_datatype<double>(), gprops);
			DataSet node_set = file.createDataSet("GEOM_GROUP/NODE_ARRAY", node_space, create_datatype<int>(), nprops);

			// write initial part of datasets -- root vessels
			geom_set.select({0,0}, {std::size_t(2*rootNumVessels), 7}).write(geomArray);
			node_set.select({0,0}, {std::size_t(2*rootNumVessels), 2}).write(nodeArray);

			file.flush();
		}
	} catch (Exception& err) {
        // catch and print any HDF5 error
        std::cerr << err.what() << std::endl;
    }

	nodeArray.clear();
	geomArray.clear();

	// populate array for root number of connected vessels
	short* numConnectedVesselsArray = new short[2 * rootNumVessels];
	long long totalNumConnectedVessels = PopulateRootNumConnectedVesselsArray_HF(arteries, veins, numConnectedVesselsArray, rootNumVessels);

	// populate array for root connected vessels
	long long* connectedVesselArray = new long long[totalNumConnectedVessels];
	PopulateRootConnectedVesselArray_HF(rootNumVessels, arteries, veins, connectedVesselArray);

	// Write Root Connected Vessel Data To File
	using namespace HighFive;
	try {
		if (!mpiRank) {
			File file(filename, File::ReadWrite);

			long long connectedVesselsSize = pow(2, totalLevels) * 7 - 12;

			// create dataspaces with initial and max shape
			DataSpace numCon_space = DataSpace({std::size_t(2*rootNumVessels), 1}, {std::size_t(totalNumVessels), 1});
			DataSpace con_space = DataSpace({std::size_t(totalNumConnectedVessels), 1}, {std::size_t(connectedVesselsSize), 1});

			// add chunking
			DataSetCreateProps ncprops, cprops;
			ncprops.add(Chunking(std::vector<hsize_t>{std::size_t(totalNumVessels), 1}));
			cprops.add(Chunking(std::vector<hsize_t>{std::size_t(connectedVesselsSize), 1}));
			// ncprops.add(Chunking(std::vector<hsize_t>{std::size_t(1), 1}));
			// cprops.add(Chunking(std::vector<hsize_t>{std::size_t(1), 1}));

			// create datasets
			DataSet numCon_set = file.createDataSet("GEOM_GROUP/NUM_CONNECTED_VESSELS_ARRAY", numCon_space, create_datatype<short>(), ncprops);
			DataSet con_set = file.createDataSet("GEOM_GROUP/CONNECTED_VESSELS_ARRAY", con_space, create_datatype<long long>(), cprops);

			// write initial part of datasets
			numCon_set.select({0,0}, {std::size_t(2*rootNumVessels), 1}).write(numConnectedVesselsArray);
			con_set.select({0,0}, {std::size_t(totalNumConnectedVessels), 1}).write(connectedVesselArray);

			file.flush();
		}
	} catch (Exception& err) {
        // catch and print any HDF5 error
        std::cerr << err.what() << std::endl;
    }

	delete[] connectedVesselArray;
    delete[] numConnectedVesselsArray;
}

//******
void ExportSequential_HF(std::string filename, int totalLevels, VesselVector& vessels) {
	long long totalNumVessels = pow(2, totalLevels + 1) - 2;

	std::vector<std::vector<double>> geomArray(totalNumVessels, std::vector<double>(7, 0));
	std::vector<std::vector<int>> nodeArray(totalNumVessels, std::vector<int>(2, 0));
	ConvertSequentialVesselVectorsToDataTables_HF(totalNumVessels, vessels, geomArray, nodeArray);

	// write data to HDF5
	using namespace HighFive;
	try {
		File file(filename, File::Truncate);
		std::vector<size_t> geom_dims{std::size_t(totalNumVessels), 7};
		std::vector<size_t> node_dims{std::size_t(totalNumVessels), 2};
		DataSet geom_array = file.createDataSet<double> ("GEOM_GROUP/GEOM_ARRAY", DataSpace(geom_dims));
		DataSet node_array = file.createDataSet<int> ("GEOM_GROUP/NODE_ARRAY", DataSpace(node_dims));
		geom_array.write(geomArray);
		node_array.write(nodeArray);
	} catch (Exception& err) {
        // catch and print any HDF5 error
        std::cerr << err.what() << std::endl;
    }

	nodeArray.clear();
	geomArray.clear();
}

//******
// prep file
void prepFile_HF(std::string filename, int numLevels, std::vector<double>& bbox, double output_bc) {
	// for (int i=0; i < bbox.size(); i++) {
	// 	std::cout << bbox[i] << std::endl;
	// }
	using namespace HighFive;
	{
		File file(filename, File::Truncate);

		// export bounding box
		std::vector<size_t> bbox_dims {6, 1};
		DataSet bbox_dset = file.createDataSet<double>("/boundingBox", DataSpace(bbox_dims));
		bbox_dset.write(bbox);

		// export number of levels
		file.createDataSet("/numLevels", numLevels);

		// create partition atributes
		file.createAttribute("Num_Partitions", 1);
		file.createAttribute("Is_Partitioned", 0);

		// export boundary conditions
		std::vector<long long> boundaryConditionNodes {0, (long long)(3 * pow(2, numLevels - 1) - 1)};
		std::vector<double> boundaryConditions {1, output_bc};
		std::vector<size_t> bc_dims{2ul, 1ul};
		short bcType[2] = {-1, 1};
		DataSet bc_nodes_dset = file.createDataSet<long long>("FLOW_GROUP/BOUNDARY_CONDITION_NODES", DataSpace(bc_dims));
		bc_nodes_dset.write(boundaryConditionNodes);
		DataSet bc_value_dset = file.createDataSet<double>("FLOW_GROUP/BOUNDARY_CONDITION_VALUE", DataSpace(bc_dims));
		bc_value_dset.write(boundaryConditions);
		DataSet bc_type_dset = file.createDataSet<short>("FLOW_GROUP/BOUNDARY_CONDITION_TYPE", DataSpace(bc_dims));
		bc_type_dset.write(bcType);

		file.flush();
	}
}