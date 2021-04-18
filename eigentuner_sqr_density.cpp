// EigenTuner v0.9
//
// This program finds a 1D potential energy function V(x), with property V(x)=V(-x), which results in the first N eigenenergies
// of a particle trapped in that kind of a potential well having given user-inputted values. The base system is a finite square
// well of given depth and width, and the program seeks a perturbation term that shifts the lowest energies to the appropriate 
// values given by the user.
//
// The computation is done with the shooting method and 4th order Runge-Kutta integration.
// Teemu Isojärvi 18.04.2021


#include <math.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <ctime>

using namespace std;

int idx=0;
int npoints=60000;	// Number of points in the x-discretization
int ni;		// Number of iterations used in the tuning
int Enum; 	// Number of eigenvalues to adjust
int output_setting;

double psi;	// Temporary real wave function value to be stored
double psi1st;	// First derivative of the wave function
double psi2nd;	// Second derivative of the wave function
double E;	// Temporary eigenenergy value
double psiprev;	// Another stored wave function value
double length;	// Width of the square well
double V0;	// Depth of the square well
double step;	// x-step size
double F;	// A number between 0 and 1 telling how much to perturb the well on one iteration
double E1lb;	// User-inputted lower bound for ground state energy. This determines the step size when scanning for correct eigenenergies.
double flag=1.0;	// A variable that signals the end of shooting method integration

double Evalues[30];	// Array for storing the eigenenergies, 30 max
double Etarget[30];	// Target values for energy eigenvalue tuning
double Ediff[30];	// Difference between current energy values and the target values
double linsys[30][30];	// A symmetric matrix for the linear system where the perturbation term multipliers are solved
double rhs[30][1];	// The RHS vector in the linear system
double norm[30];	// An array of normalization constants for the solved energy eigenfunctions
double *pe = new double[100000];	// An array for potential energy values, made large enough to contain them for any grid size needed in practice
double (*psivals)[100000] = new double[30][100000]; // A matrix containing the stored wave function values in its rows. The row and column number is intentionally made larger than needed
double *xarray = new double[100000]; // An array for the x-coordinates of n:th discrete grid point


//******************
//******************

double potential(double x)
{
// Potential energy function V(x) for the QPW

	double V;
	double num=x/step;

	V=pe[(int)num];
	return V;
}

//*****************
//*****************

double k1c1; double k1c2; // Define the variables containing the components
double k2c1; double k2c2; // of vectors k1-k4 in the 2-dimensional RK4 method
double k3c1; double k3c2;
double k4c1; double k4c2;

//*************************
//*************************

double k1(double x,double h) // Parameter 'k1' in RK4 integration
{
// The "k1" in RK4
	k1c1=psi1st;
	k1c2=2.0*(potential(x)-E)*psi;
}

double k2(double x,double h) // Parameter 'k2' in RK4 integration
{
// The "k2" in RK4
	k2c1=psi1st+h*k1c1/2.0;
	k2c2=2.0*(potential(x+h/2.0)-E)*(psi+h*k1c2/2.0);
}

double k3(double x,double h) // Parameter 'k3' in RK4 integration
{
// The "k3" in RK4
	k3c1=psi1st+h*k2c1/2.0;
	k3c2=2.0*(potential(x+h/2.0)-E)*(psi+h*k2c2/2.0);
}

double k4(double x,double h) // Parameter 'k4' in RK4 integration
{
// The "k4" in RK4
	k4c1=psi1st+h*k3c1;
	k4c2=2.0*(potential(x+h)-E)*(psi+h*k3c2);
}

//**************************
//**************************

double solve_eigenstates() { // A subprogram for solving the N energies and probability densities on each iteration

	E=-10.0;
	int n;
	double x;

	for(int qq=0; qq<Enum; qq++) norm[qq]=0.0;

	for(int k=0; k<30; k++) {
		for(int l=0; l<npoints; l++) psivals[k][l]=0.0;
		}

	for(int m=0; m<Enum; m++) { // Loop over the first N eigenfunctions and energies

		flag=1.0;
		n=0;

		while(flag>0.5) // Find a coarse estimate for E_m
		{
			E+=0.05*E1lb; // This is the step length when scanning on the energy axis

			if(m%2==0) {
				psi=1.0;
				psi1st=0.0;
			}
			else {
				psi=0.0;
				psi1st=1.0;
			}

			for(x=0.0; (x<2.0*length && flag>0.5); x+=step)
			{
			// 4th order Runge-Kutta integration with appropriate initial condition at x=0 and with the current
			// energy trial value E.

				k1(x,step);
				k2(x,step);
				k3(x,step);
				k4(x,step);
				if(abs(psi)<100.0) {
					psi1st+=(1.0/6.0)*step*(k1c2+2.0*k2c2+2.0*k3c2+k4c2);
					psi+=(1.0/6.0)*step*(k1c1+2.0*k2c1+2.0*k3c1+k4c1);
				}
			

				if(n==1 && ((psiprev>0.0 && psi<0.0) || (psiprev<0.0 && psi>0.0))) {
				// Is the value of the wave function at the endpoint of the domain of different sign than with the previous
				// tested energy value? If yes, then the correct eigenenergy is between this and the last tested value of E.
					flag=0.0;
				}
				psiprev=psi;
				n=1;
			}
		}

		flag=1.0;
		n=0;

		E-=0.05*E1lb;

		while(flag>0.5) // Approach E_m with short E-step
		{
			E+=0.002*E1lb;	// The energy step is much shorter in this second energy scanning loop

		if(m%2==0) {
			psi=1.0;
			psi1st=0.0;
		}
		else {
			psi=0.0;
			psi1st=1.0;
		}

		idx=0;

		for(x=0.0; (x<2.0*length && flag>0.5); x+=step)
		{
			idx++;
			k1(x,step);
			k2(x,step);
			k3(x,step);
			k4(x,step);
			if(abs(psi)<100.0) {
				psi1st+=(1.0/6.0)*step*(k1c2+2.0*k2c2+2.0*k3c2+k4c2);
				psi+=(1.0/6.0)*step*(k1c1+2.0*k2c1+2.0*k3c1+k4c1);
			}
			if(((double)idx-1.0)*step<0.7*length) psivals[m][idx-1]=psi;
		}
		idx=0;
		if(n==1 && ((psiprev>0.0 && psi<0.0) || (psiprev<0.0 && psi>0.0))) { 
			flag=0.0;
		}
		psiprev=psi;
		n=1;
	}

	norm[m]=0.0;
	for(int count=0; (count<npoints && (double)count*step<0.6*length); count++) norm[m]+=2.0*psivals[m][count]*psivals[m][count]*step;
	for(int count=0; (count<npoints && (double)count*step<0.6*length); count++) psivals[m][count]/=sqrt(norm[m]);

	Evalues[m]=E;

	}
}

//*******************************************
//*******************************************

void gauss_jordan()
{
// This routine calculates the solution of a linear system of equations with Gaussian elimination. It is needed
// for calculating the multipliers of the perturbation terms on each iteration.
// The code has been mimicked from the book "Numerical recipes in C" and modified enough to not be a copyright problem.

	int n = Enum;
	int m = 1;
	int column_indices[n];
	int row_indices[n];
	int pivot_indices[n];
	int i,j,k,l,q;
	int column_index;
	int row_index;

	double largest;
	double pivot_element;
	double temp_var1;
	double temp_var2;

	for (j=0;j<n;j++) pivot_indices[j]=0;
	for (i=0;i<n;i++) {
		largest=0.0;
		for (j=0;j<n;j++) {
			if (pivot_indices[j] != 1) {
				for (k=0;k<n;k++) {
					if (pivot_indices[k] == 0) {
						if (abs(linsys[j][k]) >= largest) {
							row_index=j;
							column_index=k;
							largest=abs(linsys[j][k]);	
						}
					}
				}
			}
		}

		pivot_indices[column_index]++;

		if (row_index != column_index) {
			for (l=0;l<n;l++) {
				temp_var2=linsys[row_index][l];
				linsys[row_index][l]=linsys[column_index][l];
				linsys[column_index][l]=temp_var2;
			}

			for (l=0;l<m;l++) {
				temp_var2=rhs[row_index][l];
				rhs[row_index][l]=rhs[column_index][l];
				rhs[column_index][l]=temp_var2;
			}
		}
		row_indices[i]=row_index;
		column_indices[i]=column_index;
		
		pivot_element=linsys[column_index][column_index];
		linsys[column_index][column_index]=1.0;
		for (l=0;l<n;l++) linsys[column_index][l]/=pivot_element;
		for (l=0;l<m;l++) rhs[column_index][l]/=pivot_element;
		for (q=0;q<n;q++) {
			if (q != column_index) {
				temp_var1=linsys[q][column_index];
				linsys[q][column_index]=0.0;
				for (l=0;l<n;l++) linsys[q][l] -= linsys[column_index][l]*temp_var1;
				for (l=0;l<m;l++) rhs[q][l] -= rhs[column_index][l]*temp_var1;
			}
		}
	}

	for (l=n-1;l>=0;l--) {
		if (row_indices[l] != column_indices[l]) {
			for (k=0;k<n;k++) {
				temp_var2=linsys[k][column_indices[l]];
				linsys[k][column_indices[l]]=linsys[k][row_indices[l]];
				linsys[k][row_indices[l]]=temp_var2;
			}
		}
	}

}

//*******************************************
//*******************************************

main()
{
// This is the 'main' subprogram for this application

	double intgrl;

	cout << "This program forms a 1D quantum potential well with the first N eigenvalues of the Hamiltonian having some pre-defined\n";
	cout << "values given by the user. It begins with a finite square well centered at x=0 and having length L and depth V0. Then\n";
	cout << "it repeatedly adds perturbation terms that can be expected to move the eigenvalues to the direction of the correct target\n";
	cout << "values. If the parameters are chosen correctly, the result is a potential energy function and a set of wave function\n";
	cout << "probability densities of a system where the energies E_1 to E_N are close to the user-defined target values.\n\n";

	cout << "As an input, the program takes the length and depth of the initial square well, a lower bound estimate for the ground state energy\n";
	cout << "of the system, and a fraction F which is a number between 0.0 and 1.0 and tells how much the program is to perturb the system\n";
	cout << "during a single iteration. These values should be given in decimal form, such as V0 = 200.0 instead of V0 = 200. In addition to this,\n";
	cout << "the program wants the number N which is the number of lowest eigenenergies to tune, the number of iterations to be made (this, of\n";
	cout << "course, should be a large number whenever F is small), and the values of target energies E_1 to E_N. The integer values N and the\n";
	cout << "number of iterations have to be written in integer form such as 4 instead of 4.0, and the energies in decimal form such as\n";
	cout << "E_2 = 4.0.\n\n";

	cout << "All of the target eigenenergies should be numbers between 0 and V0 and given in increasing order (e.g. not E_1 = 5.0 followed by E_2 = 2.0).\n";
	cout << "It is also not a good idea to set the first eigenenergy E_1 to a value that is much closer to zero than the other energies, e.g.\n";
	cout << "E_1 = 0.01 and E_2 = 6.0. This is not necessary, either, because in practice the spacings between energies have more importance\n";
	cout << "than their absolute values.\n\n";

	cout << "As an option, you are also asked whether the data for the potential energy function V(x) should be saved in an output file after\n";
	cout << "every iteration or only after the last one. The first N wave function probability densities of the resulting potential well are\n";
	cout << "saved in an output file only after the last iteration\n\n";

	cout << "The length of the square well: ";
	cin >> length;
	if(length <= 0.0) {
		cout << "\nError: The well length should be greater than zero.\n\n";
		exit(0);
	}

	cout << "\nThe depth of the square well: ";
	cin >> V0;
	if(V0 <= 0.0) {
		cout << "\nError: The well depth should be greater than zero.\n\n";
		exit(0);
	}

	cout << "\nGive a lower bound estimate for ground state energy: ";
	cin >> E1lb;
	if(E1lb<=0.0) {
		cout << "\nError: The lower bound has to be greater than zero.\n\n";
		exit(0);
	}

	cout << "\nNumber of energy eigenvalues to calculate and tune: ";
	cin >> Enum;
	if(Enum<1) {
		cout << "\nError: The number of considered eigenvalues has to be at least 1.\n\n";
		exit(0);
	}

	cout << "\nFraction F: ";
	cin >> F;
	if(F<=0.0 || F>1.0) {
		cout << "\nError: The value of F has to be between 0 and 1.\n\n";
		exit(0);
	}

	cout << "\nNumber of iterations: ";
	cin >> ni;
	if(ni<1) {
		cout << "\nError: The number of iterations has to be at least 1.\n\n";
		exit(0);
		}
	

	for(int n=0; n<Enum; n++) {
		cout << "\nEnergy " << (n+1) << ": ";
		cin >> Etarget[n];
		if(n>0 && Etarget[n]<Etarget[n-1]) {
			cout << "\nError: Eigenenergies must be given in increasing order.\n\n";
			exit(0);
		}

		if(Etarget[n]<=0.0) {
			cout << "\nError: Only positive eigenenergy values are currently supported.\n\n";
			exit(0);
		}

		if(Etarget[n]>=V0) {
			cout << "\nError: Do not set eigenenergies greater than the well depth.\n\n";
			exit(0);
		}
	}

	cout << "\nProduce V(x) output file on every step (1) or only the last (2)? Enter 1 or 2: ";
	cin >> output_setting;
	if(output_setting!=1 && output_setting!=2) output_setting=1;

	cout << "\n";

	step = 2.0*length/(double)npoints;


	// First, solve the energies and wave functions of the unperturbed potential well

	for(int n = 0; n<Enum; n++) {
		for(int m = 0; m<npoints; m++) {
			psivals[n][m]=0.0;
			if(n==0) xarray[m]=step*(double)m;
			if(n==0) pe[m]=0.0;
			if(n==0 && step*(double)m > length/2.0) pe[m]+=V0;
		}
	}

	solve_eigenstates();	// Find the first eigenenergies of the original system

	cout << "Energies of the unperturbed system:\n\n";
	for(int n = 0; n<Enum; n++) {
		cout << "Energy " << (n+1) << ": " << Evalues[n] << "\n";
	}

	cout << "\n";

	if(Evalues[Enum-1]>V0) {
		cout << "Error: The square well of this depth and length can't hold " << Enum << " bound states.\n\n";
		exit(0);
	}

	for(int iter=0; iter<ni; iter++) {
	// The loop for performing the appropriate number of energy tuning iterations

		for(int n=0; n<Enum; n++) {
			Ediff[n]=Etarget[n]-Evalues[n];	// Calculate the required changes in energy eigenvalues
			rhs[n][0]=Ediff[n];
		}

		for(int n=0; n<Enum; n++) {
			for(int m=0; m<Enum; m++) {
// These loops form the square matrix for the linear system to be solved on each iteration

				intgrl=0.0;
				for(int qq=0; (double)qq*step<0.6*length; qq++) {
					intgrl += 2.0*step*psivals[n][qq]*psivals[n][qq]*psivals[n][qq]*psivals[n][qq]*psivals[m][qq]*psivals[m][qq];
				}
				linsys[m][n]=intgrl;
				
			}
		}

		gauss_jordan();

		for(int p = 0; p<Enum; p++) {	// Adjust the potential energy function for the next iteration
			for(int r = 0; (double)r*step<0.6*length; r++) {
				pe[r] += F*rhs[p][0]*psivals[p][r]*psivals[p][r]*psivals[p][r]*psivals[p][r]; // The perturbation terms are proportional to 4th power moduli of previous wave functions
			}
		}

		solve_eigenstates();

		cout << "Energies after step " << (iter+1) << ":\n\n";
		for(int n = 0; n<Enum; n++) {
			cout << "Energy " << (n+1) << ": " << Evalues[n] << "\n";
		}

		cout << "\n";

		if(output_setting==1 || iter==ni-1) {
			ofstream pe_out("potential_energy_step" + to_string(iter+1) + ".txt");
			if (pe_out.is_open())	// Produce an output text file that contains x-coordinates and V(x) values in two columns
			{
				for(int count = (int)npoints-1; count >= 0; count--) {
					pe_out << -xarray[count] << " " << potential((double)count*step) << "\n";
    				}

				for(int count = 1; count < (int)npoints; count++) {
					pe_out << xarray[count] << " " << potential((double)count*step) << "\n";
    				}

				pe_out << "\n";

				pe_out.close();
			}
		}


	}

	ofstream pdf_out("probability_densities.txt");

	if (pdf_out.is_open())	// Produce an output text file that contains x-coordinates and V(x) values in two columns
		{
			for(int m = (int)npoints-1; m >= 0; m--) {
    				
				pdf_out << -xarray[m];

				for(int n = 0; n < Enum; n++) {
					pdf_out << " " << 20.0*psivals[n][m]*psivals[n][m];
    				}

				pdf_out << "\n";
			}

			for(int m = 1; m < (int)npoints; m++) {

				pdf_out << xarray[m];
    				
				for(int n = 0; n < Enum; n++) {
					pdf_out << " " << 20.0*psivals[n][m]*psivals[n][m];
    				}

				pdf_out << "\n";
			}

				
			pdf_out.close();
		}

	

}



