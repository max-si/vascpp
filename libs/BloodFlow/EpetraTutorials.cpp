
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
#include "ml_epetra_preconditioner.h"

#include "Trilinos_Util_CrsMatrixGallery.h"

using namespace Teuchos;
using namespace Trilinos_Util;

#include "EpetraTutorials.h"
#include "PressureListCreator.h"

// Do something with the given communicator.  In this case, we just
// print Epetra's version to the given output stream, on Process 0.
void exampleRoutine (const Epetra_Comm& comm, std::ostream& out)
{
  	if (comm.MyPID () == 0) {
    	// On (MPI) Process 0, print out the Epetra software version.
    	out << Epetra_Version () << std::endl << std::endl;
  	}

  // The type of global indices.  You could just set this to int,
  // but we want the example to work for Epetra64 as well.
#ifdef EPETRA_NO_32BIT_GLOBAL_INDICES
  // Epetra was compiled only with 64-bit global index support, so use
  // 64-bit global indices.
  typedef long long global_ordinal_type;
#else
  // Epetra was compiled with 32-bit global index support.  If
  // EPETRA_NO_64BIT_GLOBAL_INDICES is defined, it does not also
  // support 64-bit indices.
  typedef int global_ordinal_type;
#endif // EPETRA_NO_32BIT_GLOBAL_INDICES

	// 1. Create some Epetra_Map objects
	if (comm.MyPID () == 0) {out << "Creating Epetra_Map objects." << std::endl;}

	// The total (global, i.e., over all MPI processes) number of entries in the Map.
  	const global_ordinal_type numGlobalEntries = comm.NumProc () * 5;

	// Index base of the Map.  We choose zero-based (C-style) indexing.
  	const global_ordinal_type indexBase = 0;

	// Construct a Map that puts the same number of equations on each MPI process contiguosly.
  	Epetra_Map contigMap (numGlobalEntries, indexBase, comm);

	// contigMap is contiguous by construction.
  	if (! contigMap.LinearMap ()) {
    	throw std::logic_error ("The supposedly contiguous Map isn't contiguous.");
  	}

	// create new map that puts equations in round-robin style (1D cyclic)
	const int numGblIndsPerProc = 5;
	global_ordinal_type* gblIndList = new global_ordinal_type [numGblIndsPerProc];

	const int numProcs = comm.NumProc ();
	const int myRank = comm.MyPID ();
	for (int k = 0; k < numGblIndsPerProc; ++k) {
		gblIndList[k] = myRank + k*numProcs;
	}

	Epetra_Map cyclicMap (numGlobalEntries, numGblIndsPerProc, gblIndList, indexBase, comm);

	// deallocate list after construction
	if (gblIndList != NULL) {
    	delete [] gblIndList;
    	gblIndList = NULL;
  	}

	// If there's more than one MPI process in the communicator,
	// then cyclicMap is definitely NOT contiguous.
	if (comm.NumProc () > 1 && cyclicMap.LinearMap ()) {
		throw std::logic_error ("The cyclic Map claims to be contiguous.");
	}

	if (comm.NumProc () > 1 && contigMap.SameAs (cyclicMap)) {
    	throw std::logic_error ("contigMap should not be the same as cyclicMap.");
  	}

	// 2. Create some Epetra_Vector objects
	if (comm.MyPID () == 0) {out << "Epetra Maps created. Creating Epetra Vectors." << std::endl;}

  	// Create an Epetra_Vector with the contiguous Map we created above.
	Epetra_Vector x (contigMap);

	// The copy constructor performs a deep copy. x and y have the same Map.
  	Epetra_Vector y (x);

	// Create a Vector with the 1-D cyclic Map. 'false' leaves data uninitialized
	Epetra_Vector z (cyclicMap, false);

	// Set the entries of z to (pseudo)random numbers.
	(void) z.Random ();

	// Set the entries of x to all ones.
  	(void) x.PutScalar (1.0);

	// Define some constants for use below.
  	const double alpha = 3.14159;
  	const double beta = 2.71828;
  	const double gamma = -10.0;

	// x = beta*x + alpha*z   --> this is legal since x and z have compatible maps
	(void) x.Update (alpha, z, beta);

	// Set all entries of y to 42.0
	(void) y.PutScalar (42.0);

	// y = gamma*y + alpha*x + beta*z
  	y.Update (alpha, x, beta, z, gamma); 

	// Compute the 2-norm of y.
	double theNorm = 0.0;
  	(void) y.Norm2 (&theNorm);

	// Print the norm of y on Proc 0.
  	if (comm.MyPID () == 0) {out << "Norm of y: " << theNorm << std::endl;}

	// 3. Read the entries of the vector
	if (comm.MyPID () == 0) {out << "Vectors created. Reading vector x's entries." << std::endl;}

	{
		const int localLength = y.MyLength ();

		// Count the local number of entries less than 0.5.
		// Use local indices to access the entries of x_data.
		int localCount = 0;
		for (int localIndex = 0; localIndex < localLength; ++localIndex) {
			if (y[localIndex] < 0.5) {
				++localCount;
			}
		}

		int globalCount = 0;
		(void) comm.SumAll (&localCount, &globalCount, 1);

		// Find the total number of entries less than 0.5,
		// over all processes in the Vector's communicator.
		if (comm.MyPID () == 0) {out << "y has " << globalCount << " entr"
			<< (globalCount != 1 ? "ies" : "y")
			<< " less than 0.5." << std::endl;}
  	}

	// 4. Modify vector entries
	if (comm.MyPID () == 0) {out << "Read vector entries. Modifyiing vector." << std::endl;}


	{
		// Use local indices to access the entries of x_data.
		const int localLength = y.MyLength ();
		for (int localIndex = 0; localIndex < localLength; ++localIndex) {
			// Add the value of the local index to every entry of x.
			y[localIndex] += static_cast<double> (localIndex);
		}

		// Print the norm of x.
		theNorm = 0.0;
		(void) y.Norm2 (&theNorm);
		if (comm.MyPID () == 0) {out << "Norm of y (modified random numbers): " << theNorm << std::endl;}
	}
}


double powerMethod(const Epetra_Operator&A, const int niters, const double tolerance) {
	using std::cout;
	using std::endl;

	// an operator doesn't have a comm, but its domain map does
	const Epetra_Comm& comm = A.OperatorDomainMap().Comm();
	const int myRank = comm.MyPID();

	// create 3 vectors for iterating the power method
	// power method computes z = A*q, q should be in the domain of A and z should be in the range
	Epetra_Vector q (A.OperatorDomainMap());
	Epetra_Vector z (A.OperatorRangeMap());
	Epetra_Vector resid (A.OperatorRangeMap());

	// local error code
	int lclerr = 0;

	// fill the iteration vector z with random numbers to start
	lclerr = z.Random();

	// lambda: current approx. of eigenvalue of maximum magnitude
	// normz: the 2-norm of current iteration vector z
	// residual: the 2-norm of the current residual vector resid
	double lambda = 0.0;
	double normz = 0.0;
	double residual = 0.0;

	const double zero = 0.0;
	const double one = 1.0;

	// reporting frequency
	const int reportFrequency = 10;

	// do power method, until method has converged or max iter count reached
	for (int iter = 0; iter < niters; ++iter) {
		// If you feel confident that your code is correct, you may omit
		// everything having to do with lclerr, and just write the following:
		//
		// z.Norm2 (&normz);         // Compute the 2-norm of z
		// q.Scale (one / normz, z); // q := z / normz
		// A.Apply (q, z);           // z := A * q
		// q.Dot (z, &lambda);       // Approx. max eigenvalue

		lclerr = (lclerr == 0) ? z.Norm2 (&normz) : lclerr;
		lclerr = (lclerr == 0) ? q.Scale (one / normz, z) : lclerr;
		lclerr = (lclerr == 0) ? A.Apply (q, z) : lclerr;
		lclerr = (lclerr == 0) ? q.Dot (z, &lambda) : lclerr;

		// Compute and report the residual norm every reportFrequency
		// iterations, or if we've reached the maximum iteration count.
		if (iter % reportFrequency == 0 || iter + 1 == niters) {
			// If you feel confident that your code is correct, you may omit
			// everything having to do with lclerr, and just write the
			// following:
			//
			// resid.Update (one, z, -lambda, q, zero); // z := A*q - lambda*q
			// resid.Norm2 (&residual); // 2-norm of the residual vector

			lclerr = (lclerr == 0) ? resid.Update (one, z, -lambda, q, zero) : lclerr;
			lclerr = (lclerr == 0) ? resid.Norm2 (&residual) : lclerr;

			if (myRank == 0) {
				cout << "Iteration " << iter << ":" << endl
					<< "- lambda = " << lambda << endl
					<< "- ||A*q - lambda*q||_2 = " << residual << endl;
			}
		}
		if (residual < tolerance) {
			if (myRank == 0) {
				cout << "Converged after " << iter << " iterations" << endl;
			}
			break;
		} else if (iter + 1 == niters) {
			if (myRank == 0) {
				cout << "Failed to converge after " << niters << " iterations" << endl;
			}
			break;
		}
  	}

	// if any process failed to insert at least one entry, throw
	int gblerr = 0;
	(void) comm.MaxAll (&lclerr, &gblerr, 1);
	if (gblerr != 0) {
		throw std::runtime_error("The power method failed in some way.");
	}
	
	return lambda;
}


void EPetraCrsMatrix() {
	using std::cout;
  	using std::endl;

	MPI_Comm yourComm = MPI_COMM_WORLD;
	Epetra_MpiComm comm (yourComm);
	const int myRank = comm.MyPID ();
  	const int numProcs = comm.NumProc ();

	if (myRank == 0) {
    	cout << Epetra_Version () << endl << endl << "Total number of processes: " << numProcs << endl;
  	}

	// type of global indices
#ifdef EPETRA_NO_32BIT_GLOBAL_INDICES
	//64-bit global indices
	typedef long long global_ordinal_type;
	std::cout << "Using 64-bit Indices" << std::endl;
#else
	// no 64-bit indices
	typedef int global_ordinal_type;
#endif

	// num rows and columns in matrix
	const global_ordinal_type numGlobalElements = 50;

	// construct a map that puts approx. same num equations per proc
	const global_ordinal_type indexBase = 0;
	Epetra_Map map (numGlobalElements, indexBase, comm);

	// get list of global indices this process owns
	const int numMyElements = map.NumMyElements();

	global_ordinal_type* myGlobalElements = NULL;

#ifdef EPETRA_NO_32BIT_GLOBAL_INDICES
	myGlobalElements = map.MyGlobalElements64();
#else
	myGlobalElements = map.MyGlobalElements();
#endif

	if (numMyElements > 0 && myGlobalElements == NULL) {
		throw std::logic_error ("Failed to get the list of global indices");
	}

	if (myRank == 0) {
		cout << endl << "Creating the sparse matrix." << endl;
	}

	// create Epetra sparse matrix whose rows have distribution given by the map
	// with max num entries per row = 3
	Epetra_CrsMatrix A (Copy, map, 3);

	// local error code for use below
	int lclerr = 0;

	// Fill the sparse matrix one row at a time. InsertGlobalValues add entries
	// to the sparse matrix, using global column indices.
	// It changes both the graph structure and the values
	double tempVals[3];
	global_ordinal_type tempGblInds[3];
	for (int i = 0; i < numMyElements; ++i) {
		//std::cout << "mpiRank = " << myRank << ", myGlocalElements[i] = " << myGlobalElements[i] << std::endl;
		// A(0, 0:1) = [2, -1]
		if (myGlobalElements[i] == 0) {
			tempVals[0] = 2.0;
			tempVals[1] = -1.0;
			tempGblInds[0] = myGlobalElements[i];
			tempGblInds[1] = myGlobalElements[i] + 1;
			if (lclerr == 0) {
				lclerr = A.InsertGlobalValues (myGlobalElements[i], 2, tempVals, tempGblInds);
			}
			if (lclerr!=0) {
				break;
			}
		}
		// A(N-1, N-2:N-1) = [-1, 2]
		else if (myGlobalElements[i] == numGlobalElements-1) {
			tempVals[0] = -1.0;
			tempVals[1] = 2.0;
			tempGblInds[0] = myGlobalElements[i] - 1;
			tempGblInds[1] = myGlobalElements[i];
			if (lclerr == 0) {
				lclerr = A.InsertGlobalValues (myGlobalElements[i], 2, tempVals, tempGblInds);
			}
			if (lclerr!=0) {
				break;
			}
		}
		// A(i, i-i:i+1) = [-1, 2, -1]
		else {
			tempVals[0] = -1.0;
			tempVals[1] = 2.0;
			tempVals[2] = -1.0;
			tempGblInds[0] = myGlobalElements[i] - 1;
			tempGblInds[1] = myGlobalElements[i];
			tempGblInds[2] = myGlobalElements[i] + 1;
			if (lclerr == 0) {
				lclerr = A.InsertGlobalValues (myGlobalElements[i], 3, tempVals, tempGblInds);
			}
			if (lclerr!=0) {
				break;
			}
		}
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

	// output sparse matrix to txt file
	char file[16] = "test_matrix.txt";
	int mat = EpetraExt::RowMatrixToMatrixMarketFile(file, A);

	// Number of iterations
	// const int niters = 500;
	// // desired (abs) residual tolerance
	// const double tolerance = 1.0e-2;

	// // run the power method and report the results
	// double lambda = powerMethod(A, niters, tolerance);
	// if (myRank == 0) {
	// 	cout << endl << "Estimated max eigenvalue: " << lambda << endl;
	// }

	//
	// now change values in sparse matrix and run again
	//

	//
	// increase diagonal dominance
	//
	// if(myRank == 0) {
	// 	cout << endl << "Increasing magnitude of A(0,0), solving again." << endl;
	// }

	// if (A.RowMap().MyGID(0)) {
	// 	// get copy of row with global index 0 and modify diag entries of that row
	// 	const global_ordinal_type gidOfFirstRow = 0;
	// 	const int lidOfFirstRow = A.RowMap().LID(gidOfFirstRow);
	// 	int numEntriesInRow = A.NumMyEntries(lidOfFirstRow);
	// 	double* rowvals = new double[numEntriesInRow];
	// 	global_ordinal_type* rowinds = new global_ordinal_type [numEntriesInRow];

	// 	// get copy of entries and column indices of global row 0.
	// 	if (lclerr == 0) {
	// 		lclerr = A.ExtractGlobalRowCopy (gidOfFirstRow, numEntriesInRow, numEntriesInRow, rowvals, rowinds);
	// 	}
	// 	if (lclerr == 0) {
	// 		for (int i = 0; i < numEntriesInRow; ++i) {
	// 			if (rowinds[i] == gidOfFirstRow) {
	// 				// we found the diagonal entry. modify it
	// 				rowvals[i] *= 10.0;
	// 			}
	// 		}
	// 		// modify values, but not structure of matrix
	// 		// since we have already called FillComplete(), we can't modify its graph structure
	// 		if (lclerr == 0) {
	// 			lclerr = A.ReplaceGlobalValues (gidOfFirstRow, numEntriesInRow, rowvals, rowinds);
	// 		}
	// 	}

	// 	if (rowvals != NULL) {
	// 		delete[] rowvals;
	// 	}
	// 	if (rowinds != NULL) {
	// 		delete[] rowinds;
	// 	}
	// }

	// // if the owning process(es) of global row 0 failed to replace the one entry, throw
	// gblerr = 0;
	// (void) comm.MaxAll(&lclerr, &gblerr, 1);
	// if (gblerr != 0) {
	// 	throw std::runtime_error ("One of the owning process(es) of global row 0 failed to replace an entry.");
	// }

	// // run power method again
	// lambda = powerMethod(A, niters, tolerance);
	// if (myRank == 0) {
	// 	cout << endl << "Estimated max eigenvalue: " << lambda << endl;
	// }

	// this tells trilinos test framework that the test passed
	if (myRank == 0) {
		cout << "End Result: TEST PASSED" << endl;
	}
}


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
	// const global_ordinal_type numGlobalElements = 4 * numVessels - 1;
	// const global_ordinal_type numGlobalElements = numVessels;
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

	// output sparse matrix to txt file
	// char file[16] = "test_matrix.txt";
	// int mat = EpetraExt::RowMatrixToMatrixMarketFile(file, A);

	// create vectors x and b
	Epetra_Vector x (map);
	Epetra_Vector b (map);

	//b.Random(); // fill b with random values

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
  	//ML_Epetra::MultiLevelPreconditioner * MLPrec = new ML_Epetra::MultiLevelPreconditioner(A, true);

	// Tell AztecOO to use this preconditioner.
  	//solver.SetPrecOperator(MLPrec);

	// Tell AztecOO to use CG to solve the problem.
  	//solver.SetAztecOption(AZ_solver, AZ_cg);

	// Tell AztecOO to output status information every iteration 
  	// (hence the 1, which is the output frequency in terms of 
  	// number of iterations).
  	//solver.SetAztecOption(AZ_output, 1);

	// Maximum number of iterations to try.
  	//int Niters = 150; 
  	// Convergence tolerance.
  	//double tol = 1e-10;

	// Solve the linear problem.
	std::cout << std::endl << "Solving Linear Problem" << std::endl;
  	//solver.Iterate (Niters, tol);

	// Print out some information about the preconditioner
  	//if (myRank == 0) { std::cout << MLPrec->GetOutputList(); }

	// We're done with the preconditioner now, so we can deallocate it.
  	//delete MLPrec;

	//! CHANGE PRECONDITIONER
	solver.SetAztecOption(AZ_precond, AZ_Jacobi);
	solver.Iterate(100, 1.0E-8);
	std::cout << "Solver performed " << solver.NumIters() << " iterations." << std::endl
		<< "Norm of true residual = " << solver.TrueResidual() << std::endl;

	// output pressure vector to txt file
	// char pressure_file[20] = "pressure_vector.txt";
	// int vect = EpetraExt::VectorToMatrixMarketFile(pressure_file, x);


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
	// char pois_file[20] = "poiseuille_mat.txt";
	// int pois = EpetraExt::RowMatrixToMatrixMarketFile(pois_file, poiseuille);

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