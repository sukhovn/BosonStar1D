#include <ode_solver.h>
#include <cmath>
#include <iomanip>

#ifndef STAR_EQUATIONS_H
#define STAR_EQUATIONS_H

class Boson_Star_RK: public RKDP_stepper{
public:
	double kappa;

	double mass;
	double sigma;

	double omega;
	double v0;
	double phi0;

	Boson_Star_RK() = delete;
	Boson_Star_RK(double mass_, double sigma_, double omega_, double phi0_): RKDP_stepper(4),
		kappa(1.0), mass(mass_), sigma(sigma_), omega(omega_), v0(0.0), phi0(phi0_) {};

	// Potential and potential derivatives
	double vpot(double phi2) const{
		if(sigma == 0)
			return mass*mass*phi2;
		return mass*mass*phi2*std::pow(1.0 - 2.0*phi2/(sigma*sigma), 2.0);
	}
	
	double dvpot(double phi2) const{
		if(sigma == 0)
			return mass*mass;
		return mass*mass*(1.0 - 2.0*phi2/(sigma*sigma))*(1.0 - 6.0*phi2/(sigma*sigma));
	}

	double ddvpot(double phi2) const{
		if(sigma == 0)
			return 0.0;
		return 8.0*mass*mass*(3.0 * phi2/(sigma*sigma) - 1.0)/(sigma*sigma);
	}

	void derivatives(const double x, const double *y, double *dy) const{
		double rho = std::exp(-y[1])*omega*omega*y[2]*y[2] + std::exp(-y[0])*y[3]*y[3] + vpot(y[2]*y[2]);
		double pr = std::exp(-y[1])*omega*omega*y[2]*y[2] + std::exp(-y[0])*y[3]*y[3] - vpot(y[2]*y[2]);
		
		dy[0] = (1 - std::exp(y[0])*(1 - kappa*x*x*rho))/x;
		dy[1] = (std::exp(y[0])*(1 + kappa*x*x*pr) - 1)/x;
		dy[2] = y[3];
		dy[3] = -(1 + std::exp(y[0])*(1 - 0.5*kappa*x*x*(rho-pr)))*y[3]/x - std::exp(y[0])*(std::exp(-y[1])*omega*omega - dvpot(y[2]*y[2]))*y[2];
	}

	double dv_infinty(const double x, const double *y) const{
		return y[0] + y[1];
	}

	struct Origin_Asymptotic{
		double v0;
		double phi0;

		double u2;
		double v2;
		double phi2;

		double u4;
		double v4;
		double phi4;

		Origin_Asymptotic(const Boson_Star_RK &bstar): v0(bstar.v0), phi0(bstar.phi0){
			double kappa = bstar.kappa;
			double omega = bstar.omega;

			u2 = 1.0/3.0*kappa*(std::exp(-v0)*phi0*phi0*omega*omega + bstar.vpot(phi0*phi0));
			v2 = 1.0/3.0*kappa*(2*std::exp(-v0)*phi0*phi0*omega*omega - bstar.vpot(phi0*phi0));
			phi2 = 1.0/6.0*phi0*(-std::exp(-v0)*omega*omega + bstar.dvpot(phi0*phi0));

			u4 = 0.1*(-u2*u2 + 8*kappa*phi2*phi2 + 2*std::exp(-v0)*kappa*phi0*(u2*phi0 - v2*phi0 + 2*phi2)*omega*omega + 2*u2*kappa*bstar.vpot(phi0*phi0) + 4*kappa*phi0*phi2*bstar.dvpot(phi0*phi0));
			v4 = 0.1*(u2*u2 + 12*kappa*phi2*phi2 + 3*std::exp(-v0)*kappa*phi0*(u2*phi0 - v2*phi0 + 2*phi2)*omega*omega - 2*u2*kappa*bstar.vpot(phi0*phi0) - 4*kappa*phi0*phi2*bstar.dvpot(phi0*phi0));
			phi4 = 0.05*(-std::exp(-v0)*(u2*phi0 - v2*phi0 + phi2)*omega*omega + (u2*phi0 + phi2)*bstar.dvpot(phi0*phi0) + 2*phi2*(u2 - v2 + phi0*phi0*bstar.ddvpot(phi0*phi0)));
		};

		void evaluate(double xa, double *ya){
			ya[0] = u2*xa*xa + u4*xa*xa*xa*xa;
			ya[1] = v0 + v2*xa*xa + v4*xa*xa*xa*xa;
			ya[2] = phi0 + phi2*xa*xa + phi4*xa*xa*xa*xa;
			ya[3] = 2.0*phi2*xa + 4.0*phi4*xa*xa*xa;
		}
	};

	void set_initial_conditions(const double xi, double *yi){
		x = xi;

		Origin_Asymptotic asym(*this);

		asym.evaluate(xi, y);
		derivatives(xi, y, dy);
	}

	friend class Boson_Star_Shooter;
};

#endif