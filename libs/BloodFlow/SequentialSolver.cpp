
#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <stdexcept>
#include <sstream>
#include <hdf5.h>
#include <mpi.h>

// Trilinos Epetra headers
#include <Epetra_config.h>
#include <Epetra_Version.h>
#include <Epetra_MpiComm.h>
#include <Epetra_Map.h>
#include <Epetra_Vector.h>
#include <Epetra_CrsMatrix.h>
#include <EpetraExt_RowMatrixOut.h>
#include "Epetra_LinearProblem.h"
#include "Epetra_Time.h"
#include "AztecOO.h"
#include "EpetraExt_VectorOut.h"

// The ML include file required when working with Epetra objects.
#include "ml_include.h"
#include "ml_epetra_preconditioner.h"
#include "Teuchos_ParameterList.hpp"
#include "Trilinos_Util_CrsMatrixGallery.h"

using namespace Teuchos;
using namespace Trilinos_Util;

#include "SequentialSolver.h"
#include "VesselGenerator/vessel.h"
#include "VesselGenerator/ImportExportCommon.h"
#include "VesselGenerator/Hdf5FileStructureNames.h"
#include "VesselGenerator/ImportExportHdfBlock.h"


void AssembleMatrixSequential(Network& network, std::string filename, int numLevels) {
	//std::cout << std::endl << "ASSEMBLE MATRIX: " << std::endl;
	MPI_Comm yourComm = MPI_COMM_WORLD;
	Epetra_MpiComm comm (yourComm);
	const int myRank = comm.MyPID ();
  	const int numProcs = comm.NumProc ();

	int numVessels = pow(2, numLevels+1) - 2;

	// determine number of rows per node
	int count = numVessels / numProcs;
	int remainder = numVessels % numProcs;
	int start = myRank * count + std::min(myRank, remainder);
	int end = (myRank + 1) * count + std::min(myRank + 1, remainder);
	int numRows = end-start;

	// Import Arrays
	std::vector<long long> nodeArray;
	std::vector<double> conductanceArray;
	ImportNodesAndConductances(filename, start, numRows, nodeArray, conductanceArray);

	typedef int global_ordinal_type;

	// num rows and columns in matrix
	const global_ordinal_type numGlobalElements = network.getNumNodes();

	// construct a map that puts approx. same num equations per proc
	const global_ordinal_type indexBase = 0;
	Epetra_Map map (numGlobalElements, indexBase, comm);

	// get list of global indices this process owns
	const int numMyElements = map.NumMyElements();

	global_ordinal_type* myGlobalElements = NULL;

	myGlobalElements = map.MyGlobalElements();

	// create Epetra sparse matrix whose rows have distribution given by the map
	// with max num entries per row = 6
	Epetra_CrsMatrix A (Copy, map, 6);

	// local error code for use below
	int lclerr = 0;

	// Fill the sparse matrix one row at a time. InsertGlobalValues add entries
	// to the sparse matrix, using global column indices.
	// It changes both the graph structure and the values
	double tempVals[2];
	global_ordinal_type tempGblInds[2];
	int index = 0;
	for (int i = 0; i < numVessels; i++) {
		double conductance = conductanceArray[i];
		const int node1 = nodeArray[index];
		const int node2 = nodeArray[index+1];

		if (node1 == 0) {
			tempVals[0] = 1.0;
			tempGblInds[0] = node1;
			if (lclerr == 0) {
				lclerr = A.InsertGlobalValues(node1, 1, tempVals, tempGblInds);
			}
			if (lclerr!=0) {
				break;
			}
		} else {
			tempVals[0] = -conductance;
			tempVals[1] = conductance;
			tempGblInds[0] = node1;
			tempGblInds[1] = node2;
			if (lclerr == 0) {
				lclerr = A.InsertGlobalValues(node1, 2, tempVals, tempGblInds);
			}
			if (lclerr!=0) {
				break;
			}
		}
		if (node2 == 0) {
			tempVals[0] = 1.0;
			tempGblInds[0] = node2;
			if (lclerr == 0) {
				lclerr = A.InsertGlobalValues(node2, 1, tempVals, tempGblInds);
			}
			if (lclerr!=0) {
				break;
			}
		} else {
			tempVals[0] = -conductance;
			tempVals[1] = conductance;
			tempGblInds[0] = node2;
			tempGblInds[1] = node1;
			if (lclerr == 0) {
				lclerr = A.InsertGlobalValues(node2, 2, tempVals, tempGblInds);
			}
			if (lclerr!=0) {
				break;
			}
		}
		index += 2;
	}

	// if any process failed to insert at least 1 entry, throw
	int gblerr = 0;
	(void) comm.MaxAll (&lclerr, &gblerr, 1);
	if (gblerr !=0 ) {
		throw std::runtime_error ("Some process failed to insert an entry.");
	}

	// tell the sparse matrix we are done adding entries to it
	gblerr = A.FillComplete();
	if(gblerr != 0) {
		std::ostringstream os;
		os << "A.FillComplete() failed with error code " << gblerr << ".";
		throw std::runtime_error (os.str ());
	}

	// clear arrays
	nodeArray.clear();
	conductanceArray.clear();

	// output sparse matrix to txt file
	// char file[16] = "test_matrix.txt";
	// int mat = EpetraExt::RowMatrixToMatrixMarketFile(file, A);

	// create vectors x and b
	Epetra_Vector x (map);
	Epetra_Vector b (map);

	// create pressure solution vector
	std::cout << std::endl << "Creating Pressure Solution Vector." << std::endl;
	for (int i = 0; i < network.getNumNodes(); i++) {
		if (i == 0) {
			//int 	ReplaceGlobalValues (int NumEntries, const double *Values, const int *Indices)
			tempVals[0] = 1.0;
			tempGblInds[0] = i;
			if (lclerr == 0) {
				lclerr = b.ReplaceGlobalValues(1, tempVals, tempGblInds);
			}
			if (lclerr!=0) {
				break;
			}
		} else if (i == network.getNumNodes()-1) {
			tempVals[0] = network.getOutputFlowRate();
			tempGblInds[0] = i;
			if (lclerr == 0) {
				lclerr = b.ReplaceGlobalValues(1, tempVals, tempGblInds);
			}
			if (lclerr!=0) {
				break;
			}
		} 
	}

	// output flow solutions to txt file
	// char psol_file[16] = "psol_vector.txt";
	// int psol = EpetraExt::VectorToMatrixMarketFile(psol_file, b);

	// create linear problem
	std::cout << std::endl << "Initializing Linear Problem and Setting Preconditioner" << std::endl;
	Epetra_LinearProblem problem(&A, &x, &b);

	// create AztecOO instance
	AztecOO solver(problem);

	// Create the preconditioner object and compute the multilevel hierarchy.
	Teuchos::ParameterList MLList;

	// set prec values
	ML_Epetra::SetDefaults("DD-ML", MLList);
	MLList.set("max levels", 15);
	MLList.set("increasing or decreasing", "increasing");
	//MLList.set("aggregation: type", "METIS");
	MLList.set("smoother: type", "Aztec");
	MLList.set("smoother: Aztec as solver", false);
	MLList.set("coarse: type", "Amesos-UMFPACK");
  	ML_Epetra::MultiLevelPreconditioner * MLPrec = new ML_Epetra::MultiLevelPreconditioner(A, MLList, true);

	// Tell AztecOO to use this preconditioner.
  	solver.SetPrecOperator(MLPrec);

	// Tell AztecOO to use GMRES to solve the problem.
	solver.SetAztecOption(AZ_solver, AZ_gmres);
	solver.SetAztecOption(AZ_kspace, 35);

	// Tell AztecOO to output status information every iteration 
  	// (hence the 1, which is the output frequency in terms of 
  	// number of iterations).
  	solver.SetAztecOption(AZ_output, 1);

	// Maximum number of iterations to try.
  	int Niters = 150; 
  	// Convergence tolerance.
  	double tol = 1e-8;

	// Solve the linear problem.
	std::cout << std::endl << "Solving Linear Problem" << std::endl;
  	solver.Iterate (Niters, tol);

	// Print out some information about the preconditioner
  	if (myRank == 0) { std::cout << MLPrec->GetOutputList(); }

	// deallocate preconditioner
  	delete MLPrec;

	std::cout << "Solver performed " << solver.NumIters() << " iterations." << std::endl
		<< "Norm of true residual = " << solver.TrueResidual() << std::endl;

	// output pressure vector to txt file
	char pressure_file[20] = "seq_pressures.txt";
	int vect = EpetraExt::VectorToMatrixMarketFile(pressure_file, x);


	// Solve For Flows
	std::cout << std::endl << "Found Pressures. Setting up Flow Solutions" << std::endl;

	// create poiseuille matrix
	std::cout << std::endl << "Creating Poiseuille Matrix" << std::endl;
	const global_ordinal_type elements = numVessels;
	Epetra_Map map2 (elements, indexBase, comm);

	Epetra_CrsMatrix poiseuille (Copy, map2, 2);
	index = 0;
	for (int i = 0; i < numVessels; i++) {
		double conductance = conductanceArray[i];
		const int node1 = nodeArray[index];
		const int node2 = nodeArray[index+1];
		tempVals[0] = conductance;
		tempVals[1] = - conductance;
		tempGblInds[0] = node1;
		tempGblInds[1] = node2;
		if (lclerr == 0) {
			lclerr = poiseuille.InsertGlobalValues(i, 2, tempVals, tempGblInds);
		}
		if (lclerr!=0) {
			break;
		}
		index += 2;
	}

	// if any process failed to insert at least 1 entry, throw
	gblerr = 0;
	(void) comm.MaxAll (&lclerr, &gblerr, 1);
	if (gblerr !=0 ) {
		throw std::runtime_error ("Some process failed to insert an entry.");
	}

	// tell the sparse matrix we are done adding entries to it
	gblerr = poiseuille.FillComplete();
	if(gblerr != 0) {
		std::ostringstream os;
		os << "poiseuille.FillComplete() failed with error code " << gblerr << ".";
		throw std::runtime_error (os.str ());
	}

	// output poiseuille matrix to txt file
	char pois_file[20] = "poiseuille_mat.txt";
	int pois = EpetraExt::RowMatrixToMatrixMarketFile(pois_file, poiseuille);

	std::cout << std::endl << "Performing Matrix Vector Multiplication" << std::endl;

	Epetra_Vector y (map2);
	bool trans = false;
	(void) poiseuille.Multiply(trans, x, y);

	// check flows are correct
	bool correct = CheckFlowsAreCorrect(network, y);
	if (correct) { std::cout << std::endl << "All Healthy Flows Calculated Correctly." << std::endl; }
	else { std::cout << std::endl << "ERROR: Flows Not Inverse Power of 2." << std::endl; }

	// output flow solutions to txt file
	char flow_file[16] = "flow_mat.txt";
	int flow = EpetraExt::VectorToMatrixMarketFile(flow_file, y);

	std::cout << std::endl << "Flow Solution Vector Exported to " << flow_file << std::endl;
}


void MLAztecOO () {
	Epetra_MpiComm Comm (MPI_COMM_WORLD);
	Epetra_Time Time(Comm);

	// Initialize a Gallery object, for generating a 3-D Laplacian
  	// matrix distributed over the given communicator Comm.
  	CrsMatrixGallery Gallery("laplace_3d", Comm);

	Gallery.Set("problem_size", 1000);

	// Get pointers to the generated matrix and a test linear problem.
  	Epetra_RowMatrix* A = Gallery.GetMatrix();

	Epetra_LinearProblem* Problem = Gallery.GetLinearProblem();

	// Construct an AztecOO solver object for this problem.
  	AztecOO solver (*Problem);

  	// Create the preconditioner object and compute the multilevel hierarchy.
  	ML_Epetra::MultiLevelPreconditioner * MLPrec = new ML_Epetra::MultiLevelPreconditioner(*A, true);

	// Tell AztecOO to use this preconditioner.
  	solver.SetPrecOperator(MLPrec);

  	// Tell AztecOO to use CG to solve the problem.
  	solver.SetAztecOption(AZ_solver, AZ_cg);

	// Tell AztecOO to output status information every iteration 
  	// (hence the 1, which is the output frequency in terms of 
  	// number of iterations).
  	solver.SetAztecOption(AZ_output, 1);

  	// Maximum number of iterations to try.
  	int Niters = 150; 
  	// Convergence tolerance.
  	double tol = 1e-10;

	// Solve the linear problem.
  	solver.Iterate (Niters, tol);

  	// Print out some information about the preconditioner
  	if (Comm.MyPID() == 0) {std::cout << MLPrec->GetOutputList();}

	// We're done with the preconditioner now, so we can deallocate it.
  	delete MLPrec;

	// Verify the solution by computing the residual explicitly.
  	double residual = 0.0;
  	double diff = 0.0;
  	Gallery.ComputeResidual (&residual);
  	Gallery.ComputeDiffBetweenStartingAndExactSolutions (&diff);

	// The Epetra_Time object has been keeping track of elapsed time
  	// locally (on this MPI process).  Take the min and max globally
  	// to find the min and max elapsed time over all MPI processes.
  	double myElapsedTime = Time.ElapsedTime ();
  	double minElapsedTime = 0.0;
  	double maxElapsedTime = 0.0;
  	(void) Comm.MinAll (&myElapsedTime, &minElapsedTime, 1);
  	(void) Comm.MaxAll (&myElapsedTime, &maxElapsedTime, 1);

	if (Comm.MyPID()==0) {
    	const int numProcs = Comm.NumProc ();
		std::cout << "||b-Ax||_2 = " << residual << std::endl
			<< "||x_exact - x||_2 = " << diff << std::endl
			<< "Min total time (s) over " << numProcs << " processes: " 
			<< minElapsedTime << std::endl
			<< "Max total time (s) over " << numProcs << " processes: "
			<< maxElapsedTime << std::endl;
  	}
}


bool CheckFlowsAreCorrect(Network& network, Epetra_Vector& flows) {
	bool isCorrect = true;
	int index = 0;
	for (int i = 0; i < network.getNumLevels(); i++) {
		int numVesselsInLevel = pow(2, i);
		for (int j = 0; j < numVesselsInLevel; j++) {
			double symmetricFlow = network.getOutputFlowRate() / pow(2, i);
			double healthyFlow = flows[index];
			if ((fabs(symmetricFlow - healthyFlow) / symmetricFlow) > 0.001) {
				return false;
			}
			index++;
		}
	}
	return isCorrect;
}

void ImportNodesAndConductances(std::string filename, int start, int numRows, std::vector<long long>& nodeArray, std::vector<double>& conductanceArray) {
 	hid_t fileId = OpenHdfFile(filename);

 	hid_t nodeDatasetId = OpenHdfDataset(fileId, GEOM_GROUP_NAME, NODE_DATASET);
 	CollectiveHdfBlockImport<long long, 2, 2>(nodeDatasetId, start, numRows, nodeArray);
 	H5Dclose(nodeDatasetId);

 	hid_t conductanceDatasetId = OpenHdfDataset(fileId, FLOW_GROUP_NAME, HEALTHY_CONDUCTANCE_DATASET);
 	CollectiveHdfBlockImport<double, 2, 1>(conductanceDatasetId, start, numRows, conductanceArray);
 	H5Dclose(conductanceDatasetId);

	

 	H5Fclose(fileId);
}