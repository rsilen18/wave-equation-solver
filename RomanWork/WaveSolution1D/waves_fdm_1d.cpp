/* waves_fdm_1d.cpp */

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
		(*x)[i] = def.b+dx*(i-(def.N+order/2));
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

// discrete delta fcn for left Dirichlet BC
double* f_left(double* u, double* un, int ja, Def &def, double dx) {
	double *val = (double*) calloc(2, sizeof(double));
	val[0] = pow(def.c,2)/pow(dx,2)*(
				un[ja+1]
			 -2*un[ja]
			 +  u[0]
			 -1/12*(
			  	   un[ja+2]
			    -4*un[ja+1]
			    +6*un[ja]
			    -4*u[0]
			    +  u[1]
			  )
			);
	val[1] = pow(def.c,4)/pow(dx,4)*(
				un[ja+2]
			 -4*un[ja+1]
			 +6*un[ja]
			 -4*u[0]
			 +  u[1]
			);
	return val;
}

// discrete delta fcn for right Neumann BC
double* f_right(double* u, double* un, int jb, Def &def, double dx, int n, double dt) {
	double *val = (double*) calloc(2, sizeof(double));
	val[0] = 	u[0]
			   -un[jb-1]
			   -1/3*(
			   		u[1]
			   	 -2*u[0]
			   	 +2*un[jb-1]
			   	 -  un[jb-2]
			   	)
			   -2*dx*def.right(n*dt);
    val[1] = 	u[1]
    		 -2*u[0]
    		 +2*un[jb-1]
    		 -  un[jb-2];
    return val;
}

// solve 2x2 linear system using dgesv
void solve_2x2(double *A, double *x, double *b) {
	int N = 2;
	int NRHS = 1;
	int LDA = 2;
	int *IPIV = (int*) calloc(2, sizeof(int));
	int LDB = 2;
	int INFO;

	dgesv_(&N,&NRHS,A,&LDA,IPIV,b,&LDB,&INFO);
	x[0] = b[0];
	x[1] = b[1];
	free(IPIV);
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
		// Dirichlet left BC - discrete delta fcn
		// u[0] = un[ja-1]
		// u[1] = un[ja-2]
		double *u = (double*) calloc(2, sizeof(double));
		u[0] = 0;
		u[1] = 0;
		double *f0 = f_left(u,un,ja,def,dx);
		u[0] = 1;
		double *f1 = f_left(u,un,ja,def,dx);
		u[0] = 0;
		u[1] = 1;
		double *f2 = f_left(u,un,ja,def,dx);
		// set up 'A' matrix
		double *df_du = (double*) calloc(4, sizeof(double));
		df_du[0] = f1[0]-f0[0];
		df_du[1] = f1[1]-f0[1];
		df_du[2] = f2[0]-f0[0];
		df_du[3] = f2[1]-f0[1];
		// negate f0
		double *b = (double*) calloc(2, sizeof(double));
		b[0] = -1*f0[0];
		b[1] = -1*f0[1];
		solve_2x2(df_du,u,b);
		// assign ghost points
		un[ja-1] = u[0];
		un[ja-2] = u[1];
		// free memory
		free(u);
		free(b);
		free(f0);
		free(f1);
		free(f2);
		free(df_du);

		// Neumann right BC
		// double *A = (double*) calloc(4, sizeof(double));
		// A[0] = 2/3;
		// A[1] = -2;
		// A[2] = -1/12;
		// A[3] = 1;
		// b = (double*) calloc(2, sizeof(double));
		// b[0] =   2/3*un[jb-1]
		// 	   -1/12*un[jb-2]
		// 	   +dx*def.right(n*dt);
		// b[1] = -2*un[jb-1]
		// 	   +  un[jb-2]
		// 	   +2*dx/pow(sigma,2)*(
		// 	   		def.right((n+1)*dt)
		// 	   	 -2*def.right(n*dt)
		// 	   	 +  def.right((n-1)*dt)
		// 	   	);
		// u = (double*) calloc(2, sizeof(double));
		// solve_2x2(A,u,b);
		// un[jb+1] = u[0];
		// un[jb+2] = u[1];
		// free(A);
		// free(b);
		// free(u);
		u = (double*) calloc(2, sizeof(double));
		u[0] = 0;
		u[1] = 0;
		f0 = f_right(u,un,jb,def,dx,n,dt);
		u[0] = 1;
		f1 = f_right(u,un,jb,def,dx,n,dt);
		u[0] = 0;
		u[1] = 1;
		f2 = f_right(u,un,jb,def,dx,n,dt);
		// set up 'A' matrix
		df_du = (double*) calloc(4, sizeof(double));
		df_du[0] = f1[0]-f0[0];
		df_du[1] = f1[1]-f0[1];
		df_du[2] = f2[0]-f0[0];
		df_du[3] = f2[1]-f0[1];
		// negate f0
		b = (double*) calloc(2, sizeof(double));
		b[0] = -1*f0[0];
		b[1] = -1*f0[1];
		solve_2x2(df_du,u,b);
		// assign ghost points
		un[jb+1] = u[0];
		un[jb+2] = u[1];
		// free memory
		free(u);
		free(b);
		free(f0);
		free(f1);
		free(f2);
		free(df_du);
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
		if (order >= 4) {
			unp1[i] += (pow(sigma,4)-pow(sigma,2))/12*(
						  	  un[i+2]
						   -4*un[i+1]
						   +6*un[i]
						   -4*un[i-1]
						   +  un[i-2]
						);
		}
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
void print_results(double* x, double* u, Def &def, int order, FILE *fp = stdout) {
	fprintf(fp,"x,u\n");
	int ja = order/2;
	int jb = def.N+order/2;
	for (int i = ja; i <= jb; i++) {
		fprintf(fp,"%.16f,%.16f\n",x[i],u[i]);
	}
}

/**
	Perform convergence study
	Export results to CSVs
*/
void convergence_study(Def& def, int order, double sigma_) {
	char* file_names[5] = {"icase1_1d_10.csv","icase1_1d_50.csv","icase1_1d_100.csv",
		"icase1_1d_500.csv","icase1_1d_1000.csv"};
	int n_values[5] = {10,50,100,500,1000};
	for (int i = 0; i < 5; i++) {
		FILE *fp = fopen(file_names[i],"w");
		// setup
		def.N = n_values[i];
		double dx = (def.b-def.a)/def.N;
		double nt = def.t_f/(sigma_*dx/def.c);
		nt = ceil(nt);
		double dt = def.t_f/nt;
		double sigma = def.c*dt/dx;
	
		double *x;
		init_x(&x, def, dx, order);

		double *unm1, *un, *unp1;
		alloc_time_steps(&unm1, &un, &unp1, def, order);

		// main time stepping code
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

		print_results(x, unp1, def, order, fp);
		fclose(fp);

		delete[] x;
		delete[] unm1;
		delete[] un;
		delete[] unp1;
	}
}

int main(int argc, char** argv) {

	/* ------ read in command line args --------- */

	/*	
	argv[1] = order
	argv[2] = sigma
	argv[3] = N
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
		fprintf(stderr, "Usage: ./waves_fdm_1d.out [order] [sigma] [N] [t_f] [icase] [convergence_flag]\n");
		return EXIT_FAILURE;
	}
	fill_def(def, atoi(argv[3]), stod(argv[4]), atoi(argv[5]));
	int order = atoi(argv[1]);
	double sigma_ = stod(argv[2]);

	/* -------- setup ------------------ */

	double dx = (def.b-def.a)/def.N;
	double nt = def.t_f/(sigma_*dx/def.c);
	nt = ceil(nt);
	double dt = def.t_f/nt;
	double sigma = def.c*dt/dx;

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

	print_results(x, unp1, def, order);

	delete[] x;
	delete[] unm1;
	delete[] un;
	delete[] unp1;

	if (convergence_flag == 'Y')
		convergence_study(def, order, sigma_);
	
	return EXIT_SUCCESS;
}