#include <spectral.h>

#define ODDM(a) ((a)%2 ? -1.0 : 1.0)

//Chebyshev spectral methods functions

inline double chpoint(int n, int i){//Gives out Chebyshev-Gauss-Lobatto point
	return cos(M_PI * i/(1.0 * n));
}

double ch_x(double xst, double xend, int n, int i){//Gives out Chebyshev-Gauss-Lobatto point
	return xst + 0.5 * (xend - xst) * (cos(M_PI * (n - i)/(1.0 * n)) + 1.0);
}

int ch_x_inv_low(double xst, double xend, int n, double x){
	double tmp = (2.0 * x - xst - xend)/(xend - xst);
	if(tmp >= 1.0) tmp = 0.0;
	else if(tmp <= -1.0) tmp = 1.0;
	else tmp = acos(tmp)/M_PI;
	return floor(n*(1.0 - tmp));
}

double cheb_Dn1_el(int n, int i, int j){ //Gives out an element of 1st derivative Chebyshev-Gauss-Lobatto matrix
	if(i == 0){
		if(j == 0) return (2.0 * n * n + 1.0)/6.0;
		else if(j == n) return 0.5 * ODDM(n);
		else return 2.0 * ODDM(j)/(1.0 - chpoint(n, j));
	}
	else if (i == n){
		if(j == 0) return - 0.5 * ODDM(n);
		else if(j == n) return -(2.0 * n * n + 1.0)/6.0;
		else return -2.0 * ODDM(n+j)/(1.0 + chpoint(n, j));
	}
	else{
		if(j == 0) return -0.5 * ODDM(i)/(1.0 - chpoint(n, i));
		else if(j == n) return 0.5 * ODDM(i+n)/(1.0 + chpoint(n, i));
		else if(j == i)	return - 0.5 * chpoint(n, j)/(1.0 - pow(chpoint(n, j), 2.0));
		else return ODDM(i+j)/(chpoint(n, i) - chpoint(n, j));
	}
}

double cheb_Dn2_el(int n, int i, int j){ //Gives out an element of 2nd derivative Chebyshev-Gauss-Lobatto matrix
	double xi, xj;
	if(i == 0){
		if(j == 0) return (1.0 * n * n * n * n - 1.0)/15.0;
		else if(j == n) return ODDM(n) * (n * n - 1.0)/3.0;
		else{
			xj = chpoint(n, j);
			return 2.0/3.0 * ODDM(j) * ((2.0 * n * n + 1.0) * (1.0 - xj) - 6.0)/pow(1.0 - xj, 2.0);}
	}
	else if (i == n){
		if(j == 0) return ODDM(n) * (n * n - 1.0)/3.0;
		else if(j == n) return (1.0 * n * n * n * n - 1.0)/15.0;
		else{
			xj = chpoint(n, j);
			return 2.0/3.0 * ODDM(j+n) * ((2.0 * n * n + 1.0) * (1.0 + xj) - 6.0)/pow(1.0 + xj, 2.0);}
	}
	else{
		xi = chpoint(n, i);
		if(j == 0) return 0.5 * ODDM(i) * (xi*xi + xi - 2.0)/(1.0 - xi*xi)/pow(xi - 1.0, 2.0);
		else if(j == n) return 0.5 * ODDM(i+n) * (xi*xi - xi - 2.0)/(1.0 - xi*xi)/pow(xi + 1.0, 2.0);
		else if(j == i)	return -((n*n - 1.0) * (1.0 - xi*xi) + 3.0)/3.0/pow(1.0 - xi*xi, 2.0);
		else {
			xj = chpoint(n, j);
			return ODDM(i+j) * (xi*xi + xi*xj - 2.0)/(1.0 - xi*xi)/pow(xi - xj, 2.0);}			
	}
}

double cheb_Ex_el(int n, int j, double xst, double xend, double x){
	double xc = (xst + xend - 2*x)/(xend - xst);
	if(xc <= -1.0) xc = M_PI;
	else if(xc >= 1.0) xc = 0.0;
	else xc = acos(xc);
	double denom = xc - M_PI*j/(1.0*n);
	if(denom == 0) return 1.0;

	if(j == 0){
		if(xc == 0.0) return 1.0;
		else return sin(n * xc) / tan(0.5 * xc) / (2.0 * n);
	}

	double sign = -1.0;
	if(j%2) sign = 1.0;

	if(j == n){
		if(xc == M_PI) return 1.0;
		else return sign * sin(n * xc) * tan(0.5 * xc) / (2.0 * n);
	}

	return (sin((n - 0.5)*denom)/sin(0.5*denom) +
				sin((n - 0.5)*(xc + M_PI*j/(1.0*n)))/sin(0.5*(xc + M_PI*j/(1.0*n)))
				- 2.0 * sign * cos(n*xc)) / (2.0 * n);
}

//This part relies on fft transformations

void cheb_diff(double xst, double xend, double *specti, int n, double *specto){
	double tmp;
	double sp1 = specti[n-1];
	specto[n-1] = - 4.0*n*specti[n]/(xend-xst);
	specto[n] = 0.0;
	for(int i = n-2; i >= 0; i--){
		tmp = specti[i];
		specto[i] = specto[i+2] - 4.0*(i+1)*sp1/(xend-xst);
		sp1 = tmp;
	}
	return;
}

void cheb_int(double xst, double xend, double *specti, int n, double *specto){
	double tmp1 = specti[0], tmp2;
	for(int i = 1; i < n; i++){
		tmp2 = specti[i];
		specto[i] = -(xend-xst)*(tmp1 - specti[i+1])/(4.0*i);
		tmp1 = tmp2;
	}
	specto[0] = 0.0;
	specto[n] = -(xend-xst)*tmp1/(4.0*n);
	return;
}

double chebev(double xst, double xend, double *spect, int n, double x){
	double d = 0.0, dd = 0.0, sv, y, y2;
	int j;
	if ((x-xst)*(x-xend) > 0.0){
		printf("x is not in the range xst to xend\n");
		return NAN;
	}

	y2 = 2.0 * (y = (xst + xend - 2.0 * x)/(xend - xst)); //Change of variable.
	
	for (j = n-1; j >= 1; j--){ //Clenshaw’s recurrence.
		sv = d;
		d = y2 * d - dd + spect[j];
		dd = sv;
	}

	return y * d - dd + 0.5 * spect[0]; //Last step is different.
}

double chebev_d(double xst, double xend, double *spect, int n, double x){
	double d = 0.0, dd = 0.0, sv, y, y2;
	int j;
	if ((x-xst)*(x-xend) > 0.0){
		printf("x is not in the range xst to xend\n");
		return NAN;
	}

	y2 = 2.0 * (y = (xst + xend - 2.0 * x)/(xend - xst)); //Change of variable.
	
	for (j = n-2; j >= 1; j--){ //Clenshaw’s recurrence.
		sv = d;
		d = (y2 * (j+2) * d - (j+3) * dd)/(1.0*(j+1)) + spect[j+1];
		dd = sv;
	}

	return (4.0 * y * d - 3.0 * dd + spect[1]) * 2.0 / (xst - xend); //Last step is different.
}

double chebev_int(double xst, double xend, double *spect, int n, double x){
	double d = 0.0, dd = 0.0, sv, y, y2;
	int j;
	if ((x-xst)*(x-xend) > 0.0){
		printf("x is not in the range xst to xend\n");
		return NAN;
	}

	y2 = 2.0 * (y = (xst + xend - 2.0 * x)/(xend - xst)); //Change of variable.
	
	for (j = n-1; j >= 1; j--){ //Clenshaw’s recurrence.
		sv = d;
		d = y2 * d - dd - (spect[j-1] - spect[j+1])/(2.0*j);
		dd = sv;
	}

	return 0.5*(y * d - dd)*(xend-xst); //Last step is different.
}

double cheb_integrate(double xst, double xend, double *spect, int n){
	double sum = 0.5*spect[0];
	for(int i = 2; i <= n; i+=2) sum -= spect[i]/(1.0*i*i - 1.0);
	return (xend - xst)*sum;
}