/* waves_fdm_2d.cpp */

#define _USE_MATH_DEFINES
#include <cmath>
#include <fcntl.h>
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
	double a_x;	// left bound
	double b_x;	// right bound
	double a_y;	// lower bound on y
	double b_y;	// upper bound on y
	double c;	// wave number
	int N;		// # steps in x/y-dimension
	double t_f;	// final time

	// initial displacement
	double f(double x, double y) const { return x*y*(2-x)*(3-y); }
	// initial velocity
	double g(double x, double y) const { return sin(2*M_PI*x); }
	// left Dirichlet BC
	double left(double x, double y) const { return 0; }
	// right Dirichlet BC
	double right(double x, double y) const { return 0; }
	// bottom Dirichlet BC
	double bottom(double x, double y) const { return 0; }
	// top Dirichlet BC
	double top(double x, double y) const { return 0; }
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
			def.a_x = 0;
			def.b_x = 2;
			def.a_y = 0;
			def.b_y = 3;
			def.c = 6;
			def.N = N;
			def.t_f = t_f;
			break;
		default:
			def.a_x = 0;
			def.b_x = 2;
			def.a_y = 0;
			def.b_y = 3;
			def.c = 6;
			def.N = N;
			def.t_f = t_f;
	}
}

// allocate memory for x and initialize
void init_x(double **x, const Def &def, double dx, int order) {
	*x = (double*) calloc(def.N+order+1, sizeof(double));
	// set bounds
	(*x)[order/2] = def.a_x;
	(*x)[def.N+order/2] = def.b_x;
	// left ghost points
	for (int i = 0; i < order/2; i++) {
		(*x)[i] = def.a_x-dx*(order/2-i);
	}
	// right ghost points
	for (int i = def.N+order/2+1; i < def.N+order+1; i++) {
		(*x)[i] = def.b_x+dx*(i-(def.N+order/2));
	}
	// in-between points
	for (int i = order/2+1; i < def.N+order/2; i++) {
		(*x)[i] = def.a_x+dx*(i-order/2);
	}
}

// allocate memory for x and initialize
void init_y(double **y, const Def &def, double dy, int order) {
	*y = (double*) calloc(def.N+order+1, sizeof(double));
	// set bounds
	(*y)[order/2] = def.a_y;
	(*y)[def.N+order/2] = def.b_y;
	// left ghost points
	for (int i = 0; i < order/2; i++) {
		(*y)[i] = def.a_y-dy*(order/2-i);
	}
	// right ghost points
	for (int i = def.N+order/2+1; i < def.N+order+1; i++) {
		(*y)[i] = def.b_y+dy*(i-(def.N+order/2));
	}
	// in-between points
	for (int i = order/2+1; i < def.N+order/2; i++) {
		(*y)[i] = def.a_y+dy*(i-order/2);
	}
}

// allocate memory for unm1, un, unp1
void alloc_time_steps(double** unm1, double** un, double** unp1, const Def &def, int order) {
	*unm1 = (double*) calloc(pow(def.N+order+1,2), sizeof(double));
	*un = (double*) calloc(pow(def.N+order+1,2), sizeof(double));
	*unp1 = (double*) calloc(pow(def.N+order+1,2), sizeof(double));
}

// initial conditions
void ICs(double* unm1, double* x, double* y, const Def &def, int order) {
	for (int i = 0; i < def.N+order+1; i++) {
		for (int j = 0; j < def.N+order+1; j++) {
			unm1[ind2D(i,j,def.N+order+1,def.N+order+1)] = def.f(x[i],y[j]);
		}
	}
}

void first_time_step(double* un, double* unm1, double* x, double* y, const Def &def, 
	double sigma_x, double sigma_y, int order, double dt) {
	int dim = def.N+order+1;	// dimension of un/unm1
	for (int i = order/2; i <= def.N+order/2; i++) {
		for (int j = order/2; j <= def.N+order/2; j++) {
			un[ind2D(i,j,dim,dim)] = 
				 (1-pow(sigma_x,2)-pow(sigma_y,2))*unm1[ind2D(i,j,dim,dim)]
				+dt*def.g(x[i],y[j])
				+pow(sigma_x,2)/2*(unm1[ind2D(i+1,j,dim,dim)]+unm1[ind2D(i-1,j,dim,dim)])
				+pow(sigma_y,2)/2*(unm1[ind2D(i,j+1,dim,dim)]+unm1[ind2D(i,j-1,dim,dim)]);
		}
	}
}

void BCs(double* un, double* x, double* y, Def &def, double sigma_x, double sigma_y, 
	int order, double dx, double dy, double dt, int n) {
	// dimensions for x and y
	int dim = def.N+order+1;
	// non-ghost bounds for indices
	int ia = order/2;
	int ib = def.N+order/2;
	int ja = order/2;
	int jb = def.N+order/2;
	if (order == 2) {
		for (int i = 0; i < dim; i++) {
			un[ind2D(i,ja-1,dim,dim)] = 2*un[ind2D(i,ja,dim,dim)]-un[ind2D(i,ja+1,dim,dim)];
			un[ind2D(i,jb+1,dim,dim)] = 2*un[ind2D(i,jb,dim,dim)]-un[ind2D(i,jb-1,dim,dim)];
			un[ind2D(i,ja,dim,dim)] = def.bottom(x[i], y[ja]);
			un[ind2D(i,jb,dim,dim)] = def.top(x[i], y[jb]);
		}
		for (int j = 0; j < dim; j++) {
			un[ind2D(ia-1,j,dim,dim)] = 2*un[ind2D(ia,j,dim,dim)]-un[ind2D(ia+1,j,dim,dim)];
			un[ind2D(ib+1,j,dim,dim)] = 2*un[ind2D(ib,j,dim,dim)]-un[ind2D(ib-1,j,dim,dim)];
			un[ind2D(ia,j,dim,dim)] = def.left(x[ia], y[j]);
			un[ind2D(ib,j,dim,dim)] = def.right(x[ib], y[j]);
		}
	}
}

/**
	time steps for n >= 2
*/
void main_time_step(double* unp1, double* un, double* unm1, double* x, double* y, const Def &def, 
	double sigma_x, double sigma_y, int order, double dt) {
	// dimensions for x and y
	int dim = def.N+order+1;
	for (int i = order/2; i < def.N+order/2; i++) {
		for (int j = order/2; j < def.N+order/2; j++) {
			unp1[ind2D(i,j,dim,dim)] = 2*un[ind2D(i,j,dim,dim)]
										-unm1[ind2D(i,j,dim,dim)]
										+pow(sigma_x,2)*(
											   un[ind2D(i-1,j,dim,dim)]
										  	-2*un[ind2D(i,j,dim,dim)]
										  	  +un[ind2D(i+1,j,dim,dim)]
										)
										+pow(sigma_y,2)*(
											   un[ind2D(i,j-1,dim,dim)]
											-2*un[ind2D(i,j,dim,dim)]
											  +un[ind2D(i,j+1,dim,dim)]
										);
		}
	}
}

/**
	copy contents of src into dest
*/
void copy_array(double* dest, double* src, Def &def, int order) {
	for (int i = 0; i < pow(def.N+order+1,2); i++) {
		dest[i] = src[i];
	}
}

/**
	print out vector of values into CSV format:
	{x_value, y_value, u_value}
*/
void print_results(double* x, double* y, double* u, Def &def, int order, FILE *fp = stdout) {
	fprintf(fp,"x,y,u\n");
	int ia = order/2;
	int ib = def.N+order/2;
	int ja = order/2;
	int jb = def.N+order/2;
	for (int i = ia; i <= ib; i++) {
		for (int j = ja; j <= jb; j++) {
			fprintf(fp,"%.16f,%.16f,%.16f\n",x[i],y[j],u[ind2D(i,j,def.N+order+1,def.N+order+1)]);
		}
	}
}

int main(int argc, char** argv) {

	/* ------ read in command line args --------- */

	/*	
	argv[1] = order
	argv[2] = cfl (<1 for convergence)
	argv[3] = N (resolution in x/y-dimension)
	argv[4] = t_f
	argv[5] = icase 
	argv[6] = (optional) convergence study flag (Y/N): default NO 
	*/
	Def def = Def();
	char convergence_flag;
	if (argc == 7) {
		convergence_flag = argv[6][0];
	} else if (argc == 6) {
		convergence_flag = 'N';
	} else {
		fprintf(stderr, "Usage: ./waves_fdm_2d.out [order] [cfl] [N] [t_f] [icase] [convergence_flag]\n");
		return EXIT_FAILURE;
	}
	fill_def(def, atoi(argv[3]), stod(argv[4]), atoi(argv[5]));
	int order = atoi(argv[1]);
	double cfl = stod(argv[2]);

	/* -------- setup ------------------ */

	double dx = (def.b_x-def.a_x)/def.N;
	double dy = (def.b_y-def.a_y)/def.N;
	double dt = cfl*dx*dy/def.c*sqrt(1/(pow(dx,2)+pow(dy,2)));
	double nt = ceil(def.t_f/dt);
	dt = def.t_f/nt;
	double sigma_x = def.c*dt/dx;
	double sigma_y = def.c*dt/dy;

	double *x;
	init_x(&x, def, dx, order);
	double *y;
	init_y(&y, def, dy, order);

	double *unm1, *un, *unp1;
	alloc_time_steps(&unm1, &un, &unp1, def, order);

	/* ------- main time stepping code ------- */

	ICs(unm1, x, y, def, order);
	first_time_step(un, unm1, x, y, def, sigma_x, sigma_y, order, dt);
	BCs(un, x, y, def, sigma_x, sigma_y, order, dx, dy, dt, 1);
	int n = 2;
	while (n*dt <= def.t_f) {
		main_time_step(unp1, un, unm1, x, y, def, sigma_x, sigma_y, order, dt);
		BCs(unp1, x, y, def, sigma_x, sigma_y, order, dx, dy, dt, n);
		copy_array(unm1, un, def, order);
		copy_array(un, unp1, def, order);
		n++;
	}

	print_results(x, y, unp1, def, order);

	delete[] x;
	delete[] y;
	delete[] unm1;
	delete[] un;
	delete[] unp1;

	// if (convergence_flag == 'Y')
	// 	convergence_study(def, order, sigma_);
	
	return EXIT_SUCCESS;
}