#include <cmath>
#include <iostream>

#ifndef SPECTRAL_H
#define SPECTRAL_H

//Chebyshev spectral methods functions

double ch_x(double xst, double xend, int n, int i);
//Gives out Chebyshev-Gauss-Lobatto point adjusted to [xst, xend] segment

int ch_x_inv_low(double xst, double xend, int n, double x);
//Inverts ch_x expression and finds the closest to x CGL point from the left

//Matrix routines

double cheb_Dn1_el(int n, int i, int j);
//Gives out an element of 1st derivative Chebyshev differentiation matrix

double cheb_Dn2_el(int n, int i, int j);
//Gives out an element of 2nd derivative Chebyshev differentiation matrix

double cheb_Dn_el(int npoint, unsigned int ord, int i, int j);
//Gives out an element of ord-th derivative Chebyshev differentiation matrix, explicit for ord = 1 and ord = 2, recursively for higer orders

double cheb_Ex_el(int n, int j, double xst, double xend, double x);
//Gives extrapolation coeffient to point x of the collocation polynomial 

template<unsigned int ord>
inline double cheb_Dn_el(int npoint, int i, int j){//Very inefficient, only works reasonably for ords 1 and 2
	if(ord == 0){
		if(i == j) return 1.0;
		else return 0.0;
	}
	else if(ord == 1){
		return cheb_Dn1_el(npoint, i, j);
	}
	else if(ord == 2){
		return cheb_Dn2_el(npoint, i, j);
	}	
	else{
		std::cout << "Chebyshev differentiation order too high" << std::endl;
		throw("Wrong differention order");
	}
}

template<unsigned int ord>
double cheb_Dn(int npoint, int i, int j){
	return cheb_Dn_el<ord>(npoint, npoint-i, npoint-j);
}

template<unsigned int ord, bool side>
double cheb_Dn_double(int n, bool parity, int i, int j){
	if(side){
		if(parity) return cheb_Dn_el<ord>(2*n-1, 2*n-i-1, 2*n-j-1) - cheb_Dn_el<ord>(2*n-1, 2*n-i-1, j);
		else return cheb_Dn_el<ord>(2*n-1, 2*n-i-1, 2*n-j-1) + cheb_Dn_el<ord>(2*n-1, 2*n-i-1, j);
	}
	else{
		if(parity) return cheb_Dn_el<ord>(2*n-1, n-i-1, n-j-1) - cheb_Dn_el<ord>(2*n-1, n-i-1, n+j);
		else return cheb_Dn_el<ord>(2*n-1, n-i-1, n-j-1) + cheb_Dn_el<ord>(2*n-1, n-i-1, n+j);
	}
}

//This part relies on fft transformations

void cheb_diff(double xst, double xend, double *specti, int n, double *specto);
//Finds Chebyshev polynomial coefficients for the derivative *spect given coefficients for the function

void cheb_int(double xst, double xend, double *specti, int n, double *specto);
//Finds Chebyshev polynomial coefficients for the antiderivative *spect given coefficients for the function

double chebev(double xst, double xend, double *spect, int n, double x);
//Evaluates a function at a arbitrary point x between xst and xend given Chebyshev polynomial expansion coefficients

double chebev_d(double xst, double xend, double *spect, int n, double x);
//Evaluates function's 1st derivative at an arbitrary point x between xst and xend given Chebyshev polynomial expansion coefficients

double chebev_int(double xst, double xend, double *spect, int n, double x);
//Evaluates function's integral at an arbitrary point x between xst and xend given Chebyshev polynomial expansion coefficients

double cheb_integrate(double xst, double xend, double *spect, int n);
//Evaluates function's integral at an arbitrary point x between xst and xend given Chebyshev polynomial expansion coefficients

#endif