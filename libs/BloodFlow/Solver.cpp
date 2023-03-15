#include <cmath>
#include <fstream>
#include <functional>
#include <iostream>
#include <limits>
#include <numeric>
#include <utility>
#include <vector>
#include <stdexcept>
#include <sstream>

#include <mpi.h>

#include "BloodFlowVessel.h"
#include "NetworkDescription.h"
#include "Triplet.h"

// Trilinos Epetra headers (not all needed?)
#include <Epetra_config.h>
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

using std::cout;
using std::endl;
using std::vector;

//* structure used to contain pointers to the hyper data structures
struct FormatConversionResult
{
    void ClearData()
    {
        colIndexes.clear();
        vals.clear();
        numColumnsInEachRow.clear();
        rowIndexes.clear();
    }

    void Print()
    {
        cout << "\tStart Row: " << startRow << "\n";
        cout << "\tlocal num Rows: " << localNumRows << "\n";
        cout << "\ttotal num Rows: " << totalNumRows << "\n";
        cout << "\tRow Indexes: ";
        for (auto& x : rowIndexes)
        {
            cout << x << ", ";
        }
        cout << "\n";

        cout << "\tNum Cols Per row: ";
        for (auto& x : numColumnsInEachRow)
        {
            cout << x << ", ";
        }
        cout << "\n";

        cout << "\tCol Indexes: ";
        for (auto& x : colIndexes)
        {
            cout << x << ", ";
        }
        cout << "\n";

        cout << "\tValues: ";
        for (auto& x : vals)
        {
            cout << x << ", ";
        }
        cout << "\n";
        cout << endl;
    }
    std::vector<long long> colIndexes;
    std::vector<double> vals;
    std::vector<long long> numColumnsInEachRow;
    std::vector<long long> rowIndexes;
    long long startRow = -1;
    long long localNumRows = -1;
    long long totalNumRows = -1;
    long long minCol = -1;
    long long maxCol = -1;
};

FormatConversionResult ConvertTripletListToHypreFormat(long long startRow, long long localNumRows, long long totalNumRows, std::vector<Triplet>& tripletList) {
    FormatConversionResult output;
    output.startRow = startRow;
    output.localNumRows = localNumRows;
    output.totalNumRows = totalNumRows;
    long long nnz = tripletList.size();

    output.colIndexes.resize(nnz);
    output.vals.resize(nnz);
    output.numColumnsInEachRow.resize(localNumRows);
    long long minCol = std::numeric_limits<long long>::max(), maxCol = 0;

    for (long long i = 0; i < tripletList.size(); ++i)
    {
		//std::cout << "i = " << i << std::endl;
        ++output.numColumnsInEachRow[tripletList[i].first];
        output.colIndexes[i] = tripletList[i].second;
        minCol = (tripletList[i].second < minCol) ? tripletList[i].second : minCol;
        maxCol = (tripletList[i].second > maxCol) ? tripletList[i].second : maxCol;
        output.vals[i] = tripletList[i].third;
    }
    output.minCol = minCol;
    output.maxCol = maxCol;
    tripletList.clear();
    output.rowIndexes.resize(localNumRows);
    std::iota(output.rowIndexes.begin(), output.rowIndexes.end(), output.startRow);

    return std::move(output);
}


void SolvePressures(const NetworkDescription& network, std::vector<Triplet>& tripletList, std::vector<double>& rhsOutputVector) {
	MPI_Comm yourComm = MPI_COMM_WORLD;
	Epetra_MpiComm comm (yourComm);
	const int myRank = comm.MyPID ();
  	const int numProcs = comm.NumProc ();

	//! may need to recompile Epetra with 64 bit (long long) support
	typedef int global_ordinal_type;

	// global num rows and columns in matrix
	global_ordinal_type totalNumRows = 0;

	for (int i = 0; i < numProcs; i++) {
		const std::pair<long long, long long>& rowBoundPair = network.getRowBoundPair(i);
		totalNumRows += rowBoundPair.second - rowBoundPair.first + 1;
	}

	// local num rows in matrix
	const global_ordinal_type localNumRows = network.getLocalRankEndRow() - network.getLocalRankStartRow() + 1;

	// construct a map that puts approx. same num equations per proc
	const global_ordinal_type indexBase = 0;
	Epetra_Map map (totalNumRows, localNumRows, indexBase, comm);

	// get list of global indices this process owns
	const int numMyElements = map.NumMyElements();
	global_ordinal_type* myGlobalElements = NULL;
	myGlobalElements = map.MyGlobalElements();

	// create Epetra sparse matrix whose rows have distribution given by the map with max num entries per row = 6
	Epetra_CrsMatrix A (Copy, map, 6);

	// local error code for use below
	int lclerr = 0;

	// Fill the sparse matrix one row at a time. InsertGlobalValues add entries
	// to the sparse matrix, using global column indices.
	// It changes both the graph structure and the values
	global_ordinal_type globalRow;
	double tempVals[1];
	global_ordinal_type tempGblInds[1];
	for(int i = 1; i < tripletList.size(); ++i) {
		globalRow = tripletList[i].first;
		tempGblInds[0] = tripletList[i].second;
		tempVals[0] = tripletList[i].third;
		
		if (lclerr == 0) {
				lclerr = A.InsertGlobalValues(globalRow, 1, tempVals, tempGblInds);
		}
		if (lclerr!=0) {
			std::cout << "local error" << std::endl;
			break;
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
	// char file[16] = "test_matrix.txt";
	// int mat = EpetraExt::RowMatrixToMatrixMarketFile(file, A);

	// create vectors x and b
	Epetra_Vector x (map);
	Epetra_Vector b (map);

	// create pressure solution vector
	//if(!myRank) {std::cout << std::endl << "Creating Pressure Solution Vector." << std::endl;}
	for (int i = 0; i < rhsOutputVector.size(); i++) {
		//int 	ReplaceGlobalValues (int NumEntries, const double *Values, const int *Indices)
		tempVals[0] = rhsOutputVector[i];
		tempGblInds[0] = i;
		if (lclerr == 0) {
			lclerr = b.ReplaceGlobalValues(1, tempVals, tempGblInds);
		}
		if (lclerr!=0) {
			break;
		}
	}

	// output pressure vector to txt file
	// char pvec_file[16] = "pvec_vector.txt";
	// int pvec = EpetraExt::VectorToMatrixMarketFile(pvec_file, b);

	// create linear problem
	//if(!myRank) {std::cout << std::endl << "Initializing Linear Problem and Setting Preconditioner" << std::endl;}
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
	//if (myRank == 0) {std::cout << std::endl << "Solving Linear Problem" << std::endl;}
  	solver.Iterate (Niters, tol);

	// Print out some information about the preconditioner
  	if (myRank == 0) { std::cout << MLPrec->GetOutputList(); }

	// deallocate preconditioner
  	delete MLPrec;

	if (myRank == 0) {std::cout << "Solver performed " << solver.NumIters() << " iterations." << std::endl
		<< "Norm of true residual = " << solver.TrueResidual() << std::endl;}

	// output pressure vector to txt file
	// char pressure_file[20] = "par_pressures.txt";
	// int vect = EpetraExt::VectorToMatrixMarketFile(pressure_file, x);

	//! May need to change this... ExtractCopy truncates pressures?
	rhsOutputVector.clear();
	rhsOutputVector.resize(localNumRows);
	x.ExtractCopy(rhsOutputVector.data());
}