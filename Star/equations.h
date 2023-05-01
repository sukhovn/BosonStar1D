#include <cmath>
#include <limits>

#ifndef EQUATIONS_H
#define EQUATIONS_H

class Linear_Fractional_Compactification{
public:
	double cc;

	Linear_Fractional_Compactification(): cc(0.0) {};
	Linear_Fractional_Compactification(double cci): cc(cci) {};

	double xc(double x) const{
		if(x*cc == 1.0) return std::numeric_limits<double>::infinity();
		else return x/(1 - cc*x);
	}

	double dxc(double x) const{
		if(x*cc == 1.0) return std::numeric_limits<double>::infinity();
		else return 1.0/std::pow(1 - cc*x, 2.0);
	}

	double ddxc(double x) const{
		if(x*cc == 1.0) return std::numeric_limits<double>::infinity();
		else return 2.0*cc/std::pow(1 - cc*x, 3.0);
	}

	double xc_inv(double x) const{
		return x/(1 + cc*x);
	}
};


class Boson_Star_Equations{
public:
	std::vector<double> poly_coef;

	double kappa;
	double omega;

	Boson_Star_Equations(): poly_coef(4, 0.0), kappa(8.0*M_PI) {};

	double vpot(double phi) const{
		int deg = poly_coef.size()-1;
		double tmp = poly_coef[deg];
		for(int i = deg; i > 0; i--){
			tmp *= phi;
			tmp += poly_coef[i-1];
		}
		return phi*tmp;
	}

	double dvpot(double phi) const{
		int deg = poly_coef.size()-1;
		double tmp = (deg+1)*poly_coef[deg];
		for(int i = deg; i > 0; i--){
			tmp *= phi;
			tmp += i*poly_coef[i-1];
		}
		return tmp;
	}

	double equations(int fn, double x, double *phi, double *dphi, double *ddphi) const{
		double rhs, phi_sq = phi[2]*phi[2];
		switch(fn){
			case 0:
				rhs = exp(-phi[1])*omega*omega*phi_sq + exp(-phi[0])*pow(dphi[2], 2.0) + vpot(phi_sq);
				return dphi[0] - (1.0 - exp(phi[0])*(1.0 - kappa*x*x*rhs))/x;
				break;
			case 1:
				rhs =  exp(-phi[1])*omega*omega*phi_sq + exp(-phi[0])*pow(dphi[2], 2.0) - vpot(phi_sq);
				return dphi[1] - (exp(phi[0])*(kappa*x*x*rhs + 1.0) - 1.0)/x;
				break;
			case 2:
				rhs = exp(phi[0])*(dvpot(phi_sq) - exp(-phi[1])*omega*omega)*phi[2];
				return ddphi[2] + (1.0 + exp(phi[0])*(1.0 - kappa*x*x*vpot(phi_sq)))*dphi[2]/x - rhs;
				break; 
			default:
				std::cout << "Wrong equation number" << std::endl;
				throw("Wrong equation number");
		}
	}

	double boundary_condition_left(int fn, double *phi, double *dphi) const{
		switch(fn){
			case 0:
				return phi[0];
				break;
			case 1:
				return dphi[1];
				break;
			case 2:
				return dphi[2];
				break;
		}
		std::cout << "Wrong equation number" << std::endl;
		throw("Wrong equation number");
	}

	double boundary_condition_right(int fn, double *phi, double *dphi) const{
		switch(fn){
			case 0:
				return phi[0];
				break;
			case 1:
				return phi[1];
				break;
			case 2:
				return phi[2];
				break;	
		}
		std::cout << "Wrong equation number" << std::endl;
		throw("Wrong equation number");
	}

	double residual(int fn, double x, double *phi, double *dphi, double *ddphi) const{
		double phi_sq = phi[2]*phi[2];
		double dense = exp(-phi[1])*omega*omega*phi_sq - exp(-phi[0])*pow(dphi[2], 2.0) - vpot(phi_sq);
		return std::exp(-phi[0])*(ddphi[1] + 0.5*(dphi[1] - dphi[0])*(2.0/x + dphi[1])) - 2.0*kappa*dense;
	}
};

#endif