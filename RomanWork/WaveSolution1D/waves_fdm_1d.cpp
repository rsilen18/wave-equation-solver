/* waves_fdm_1d.cpp */

#define _USE_MATH_DEFINES
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
using namespace std;

// dgeev_ is a symbol in the LAPACK library files
extern "C" {
  extern int dgesv_(int*,int*,double*,int*,int*,double*,int*,int*);
}

// Problem definition class
class Def {
public:
	double a;	// left bound
	double b;	// right bound
	double c;	// wave number
	int N;		// # steps in x-direction
	double t_f;	// final time

	// initial displacement
	double f(double x) const { return sin(x); }
	// initial velocity
	double g(double x) const { return -sin(x); }
	// left Dirichlet BC
	double left(double x) const { return 0; }
	// right Neumann BC
	double right(double x) const { return 0; }
};

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

// Fill in problem definition
void fill_def(Def &def, int N, double t_f, int icase) {
	switch (icase) {
		case 1: 
			def.a = 0;
			def.b = 3*M_PI/2;
			def.c = 1;
			def.N = N;
			def.t_f = t_f;
			break;
		default:
			def.a = 0;
			def.b = 3*M_PI/2;
			def.c = 1;
			def.N = N;
			def.t_f = t_f;
	}
}

// allocate memory for x and initialize
void init_x(double **x, const Def &def, double dx, int order) {
	*x = (double*) calloc(def.N+order+1, sizeof(double));
	// set bounds
	(*x)[order/2] = def.a;
	(*x)[def.N+order/2] = def.b;
	// left ghost points
	for (int i = 0; i < order/2; i++) {
		(*x)[i] = def.a-dx*(order/2-i);
	}
	// right ghost points
	for (int i = def.N+order/2+1; i < def.N+order+1; i++) {
		(*x)[i] = def.b+dx*(i-(def.N+3*order/2));
	}
	// in-between points
	for (int i = order/2+1; i < def.N+order/2; i++) {
		(*x)[i] = def.a+dx*(i-order/2);
	}
}

// allocate memory for unm1, un, unp1
void alloc_time_steps(double** unm1, double** un, double** unp1, const Def &def, int order) {
	*unm1 = (double*) calloc(def.N+order+1, sizeof(double));
	*un = (double*) calloc(def.N+order+1, sizeof(double));
	*unp1 = (double*) calloc(def.N+order+1, sizeof(double));
}

// initial conditions
void ICs(double* unm1, double* x, const Def &def, int order) {
	for (int i = 0; i < def.N+order+1; i++) {
		unm1[i] = def.f(x[i]);
	}
}

void first_time_step(double* un, double* unm1, double* x, const Def &def, 
	double sigma, int order, double dt) {
	for (int i = order/2; i <= def.N+order/2; i++) {
		un[i] = unm1[i]
			   +dt*def.g(x[i])
			   +pow(sigma,2)/2*(
			   		   unm1[i-1]
			   		-2*unm1[i]
			   		+  unm1[i+1]
			   	);
		if (order >= 4) {
			un[i] += pow(sigma,2)/2*(
						-1/12*(
							   unm1[i+2]
							-4*unm1[i+1]
							+6*unm1[i]
							-4*unm1[i-1]
							+  unm1[i-2]
						)
					)
					+dt*pow(sigma,2)/6*(
						   def.g(x[i+1])
						-2*def.g(x[i])
						+  def.g(x[i-1])
						-1/12*(
							   def.g(x[i+2])
							-4*def.g(x[i+1])
							+6*def.g(x[i])
							-4*def.g(x[i-1])
							+  def.g(x[i-2])
						)
					)
					+pow(sigma,4)/24*(
						   unm1[i+2]
						-4*unm1[i+1]
						+6*unm1[i]
						-4*unm1[i-1]
						+  unm1[i-2]
					);
		}
	}
}

/** fill in boundary conditions
	code for Dirichlet left BC, Neumann right BC
	n = timestep
*/
void BCs(double* un, double* x, Def &def, double sigma, int order, 
	double dx, double dt, int n) {
	int ja = order/2;
	int jb = def.N+order/2;
	if (order == 2) {
		un[ja-1] = 2*un[ja]-un[ja+1];
		un[jb+1] = un[jb-1]+2*dx*def.right(n*dt);
	} else if (order == 4) {
		
	}
}

/**
	time steps for n >= 2
*/
void main_time_step(double* unp1, double* un, double* unm1, double* x, const Def &def, 
	double sigma, int order, double dt) {
	int ja = order/2;
	int jb = def.N+order/2;
	for (int i = ja; i <= jb; i++) {
		unp1[i] = 2*un[i]
				  - unm1[i]
				  + pow(sigma,2)*(
				  		un[i+1]
				  	 -2*un[i]
				  	 +  un[i-1]
				  	);
	}
}

/**
	copy contents of src into dest
*/
void copy_array(double* dest, double* src, Def &def, int order) {
	for (int i = 0; i < def.N+order+1; i++) {
		dest[i] = src[i];
	}
}

/**
	print out vector of values into CSV format:
	{x_value, u_value}
*/
void export_results(double* x, double* u, Def &def, int order) {
	cout << "x,u" << endl;
	int ja = order/2;
	int jb = def.N+order/2;
	for (int i = ja; i <= jb; i++) {
		cout << x[i] << "," << u[i] << endl;
	}
}

int main(int argc, char** argv) {

	/* ------ read in command line args --------- */

	/*	
	argv[1] = order
	argv[2] = sigma
	argv[3] = N
	argv[4] = t_f
	argv[5] (optional) = icase (default is 1)
	*/
	Def def = Def();
	if (argc == 6) {
		fill_def(def, atoi(argv[3]), stod(argv[4]), atoi(argv[5]));
	} else if (argc == 5) {
		fill_def(def, atoi(argv[3]), stod(argv[4]), 1);
	} else {
		fprintf(stderr, "Usage: ./waves_fdm_1d.out [order] [sigma] [N] [t_f] [icase (optional)]\n");
		return EXIT_FAILURE;
	}
	int order = atoi(argv[1]);
	double sigma = stod(argv[2]);

	/* -------- setup ------------------ */

	double dx = (def.b-def.a)/def.N;
	double nt = def.t_f/(sigma*dx/def.c);
	nt = ceil(nt);
	double dt = def.t_f/nt;
	sigma = def.c*dt/dx;

	double *x;
	init_x(&x, def, dx, order);

	double *unm1, *un, *unp1;
	alloc_time_steps(&unm1, &un, &unp1, def, order);

	/* ------- main time stepping code ------- */

	ICs(unm1, x, def, order);
	first_time_step(un, unm1, x, def, sigma, order, dt);
	BCs(un, x, def, sigma, order, dx, dt, 1);
	int n = 2;
	while (n*dt <= def.t_f) {
		main_time_step(unp1, un, unm1, x, def, sigma, order, dt);
		BCs(unp1, x, def, sigma, order, dx, dt, n);
		copy_array(unm1, un, def, order);
		copy_array(un, unp1, def, order);
		n++;
	}

	export_results(x, unp1, def, order);

	delete[] x;
	delete[] unm1;
	delete[] un;
	delete[] unp1;
	
	return EXIT_SUCCESS;
}