/* lin_sys_solver.cpp */

#include <iomanip>
#include <iostream>
#include <fstream>
using namespace std;

// dgeev_ is a symbol in the LAPACK library files
extern "C" {
  extern int dgesv_(int*,int*,double*,int*,int*,double*,int*,int*);
}

// returns index of a matrix stored in row-major format
// i = row #, j = col # (0-indexed)
// m = # of rows, n = # of cols
int ind2D(int i, int j, int m, int n) {
	return j*m+i;
}

// print m x n matrix
void print_matrix(double* mtx, int m, int n) {
	cout << m << "x" << n << " matrix: " << endl;
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			cout << setw(4) << mtx[ind2D(i,j,m,n)] << " ";
		}
		cout << endl;
	}
}

int main(int argc, char** argv) {

	// check for filename
	if (argc < 3) {
		cout << "Usage: " << argv[0] << " [matrix filename] [vector filename]" << endl;
		return EXIT_FAILURE;
	}

	// define matrix dimensions: m x n 
	int m,n;
	double *A;

	// read in txt file containing a matrix
	// store matrix in col-major format
	ifstream fin1(argv[1]);
	if (!fin1.is_open()) {
		cout << "Failed to open " << argv[1] << endl;
		return -1;
	}

	fin1 >> m >> n;	// m = rows, n = cols
	A = new double[m*n];
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			fin1 >> A[ind2D(i,j,m,n)];
		}
	}
	if (fin1.fail() || fin1.eof()){
  	  cerr << "Error while reading " << argv[1] << endl;
  	  return EXIT_FAILURE;
  	}
  	fin1.close();

  	// print matrix A
  	cout << "Matrix A: " << endl;
  	print_matrix(A, m, n);


  	// define vector 
  	int p;
  	double *b;

  	// read in txt file containing a vector
	// first value is length of vector
	ifstream fin2(argv[2]);
	if (!fin2.is_open()) {
		cerr << "Failed to open " << argv[2] << endl;
		return EXIT_FAILURE;
	}
	fin2 >> p;
	if (n != p) {
		cerr << "Matrix and vector dimensions do not match!" << endl;
		return EXIT_FAILURE;
	}
	b = new double[p];
	for (int i = 0; i < p; i++) {
		fin2 >> b[i];
	}
	if (fin2.fail() || fin2.eof()) {
		cerr << "Error while reading " << argv[2] << endl;
		return EXIT_FAILURE;
	}
	fin2.close();

	// print vector b
	cout << "Vector b: " << endl;
	print_matrix(b, p, 1);


	// allocate data for dgesv
	int N = m;				// # of eqns
	int NRHS = 1;			// # of cols in vector b
	//double *A = A 		// matrix A
	int LDA = m;			// leading dimension of A
	int *IPIV = new int[p];	// pivot indices defining permutation matrix
	double *B = b;			// right hand side vector b -> exits as solution
	int LDB = p;			// leading dimension of b
	int INFO;				// info about exit: 0 = successful exit, <0 = error, >0 = singular

	dgesv_(&N,&NRHS,A,&LDA,IPIV,B,&LDB,&INFO);


	// print solution
	cout << "Solution: " << endl;
	print_matrix(b, p, 1);

	delete [] A;
	delete [] b;
	delete [] IPIV;

	return EXIT_SUCCESS;
}
