#include <star_equations.h>
#include <cstring>

#ifndef STAR_SHOOTER_H
#define STAR_SHOOTER_H

#define SHOOT_TRIES 50

class Boson_Star_Shooter{
	double xst;
	double xend;

	int npoints;
	double *xstorage;
	double *ystorage;

public:
	Boson_Star_RK bstar;
	
	double *xarr;
	double *yarr;

	Boson_Star_Shooter(double mass, double sigma, double phi0, double xst_, double xend_):
		bstar(mass, sigma, 1.0, phi0), xst(xst_), xend(xend_), xarr(nullptr), yarr(nullptr), xstorage(nullptr), ystorage(nullptr) {};
	Boson_Star_Shooter(const Boson_Star_Shooter &bshoot) = delete;
	Boson_Star_Shooter& operator=(const Boson_Star_Shooter& that) = delete;

	~Boson_Star_Shooter(){
		if(xstorage != nullptr)
			delete [] xstorage;
		if(ystorage != nullptr)
			delete [] ystorage;
	}

	std::pair<double, double> bracket_omega(double omega_step){
		double omega_st, omega_end;
		bstar.omega = 0.0;

		double *xint = new double[2];
		xint[0] = xst; xint[1] = xend;
		double *yint = new double[8];
		for(int i = 0; i < 8; ++i) yint[i] = 0.0;

		while(1){
			try{
		 		bstar.integrate_grid(1, xint, yint);
			}
			catch(const char *error){
				if(strcmp(error, "stepsize underflow in RKDP_stepper") != 0){
					std::cout << "Bracketing exited with unknown exception" << std::endl;
					throw(error);
				}
			}
			if(yint[6] == bstar.y[2])
				break;
			else if(bstar.y[2] < 0.0){
				std::cout << "Bracketing failed to find nonsingular interval" << std::endl;
				delete [] xint; delete [] yint;
				throw("bracketing fail");
			}
			bstar.omega += omega_step;
		}

		try{
			if(yint[6] < 0.0){
				omega_end = bstar.omega;
				while(1){
					bstar.omega -= 0.1*omega_step;
					bstar.integrate_grid(1, xint, yint);
					if(bstar.y[2] > 0.0){
						omega_st = bstar.omega;
						break;
					}
				}
			}
			if(yint[6] > 0.0){
				omega_st = bstar.omega;
				while(1){
					bstar.omega += 0.1*omega_step;
					bstar.integrate_grid(1, xint, yint);
					if(bstar.y[2] < 0.0){
						omega_end = bstar.omega;
						break;
					}
				}
			}
		}	
		catch(const char *error){
			delete [] xint, delete [] yint;
			if(strcmp(error, "stepsize underflow in RKDP_stepper") == 0){
				std::cout << "Bracketing failed to find one of the interval boundaries" << std::endl;
				throw("bracketing fail");
			}
			else{
				std::cout << "Bracketing exited with unknown exception" << std::endl;
				throw(error);
			}
		}

		delete [] xint;
		delete [] yint;

		return std::make_pair(omega_st, omega_end);
	}

	double shoot_interval(double omega_st, double omega_end){
		double *xint = new double [2];
		xint[0] = xst; xint[1] = xend;
		double *yint = new double[8];

		// Shooting for omega_st
		bstar.omega = omega_st;
		bstar.integrate_grid(1, xint, yint);
		double phi_st = yint[6];

		// Shooting for omega_st
		bstar.omega = omega_end;
		bstar.integrate_grid(1, xint, yint);
		double phi_end = yint[6];

		if(phi_st*phi_end > 0.0){
			std::cout << "Improper omega interval chosen" << std::endl;
			std::cout << "phi' difference for omega " << omega_st << " is " << phi_st << std::endl;
			std::cout << "phi' difference for omega " << omega_end << " is " << phi_end << std::endl;
			std::cout << "Differences should have opposite signs" << std::endl;
			delete [] xint; delete [] yint;
			throw("improper omega interval");
		}

		double omega_try, phi_try;
		int at;
		for(at = 0; at < SHOOT_TRIES; ++at){
			omega_try = 0.5*(omega_st + omega_end);
			bstar.omega = omega_try;
			bstar.integrate_grid(1, xint, yint);
			phi_try = yint[6];
			if(phi_try*phi_st > 0.0){
				phi_st = phi_try;
				omega_st = omega_try;
			}
			else{
				phi_end = phi_try;
				omega_end = omega_try;
			}

			if(std::fabs(phi_try) < 1.0e-10)
				break;
		}

		double dv0 = bstar.dv_infinty(xend, yint+4);

		delete [] xint;
		delete [] yint;

		if(at == SHOOT_TRIES){
			std::cout << std::setprecision(16);
			std::cout << "Final omega interval [" << omega_st << ", " << omega_end << "]" << std::endl;
			std::cout << "Final phi' differences [" << phi_st << ", " << phi_end << "]" << std::endl;
			std::cout << "Shooting didn't converge" << std::endl;
			throw("shooting didn't converge");
		}

		return dv0;
	}

	void shoot_for_omega(double omega_step = 0.01, bool comments = true){
		if(comments){
			std::cout << "Finding the solution by shooting\n";
			std::cout << "The value of phi0 is: " << bstar.phi0 << std::endl;
		}
		std::pair<double, double> bracket = bracket_omega(omega_step);
		if(comments) std::cout << "The interval bracketed [" << bracket.first << "," << bracket.second << ']' << std::endl;
		double dv0 = shoot_interval(bracket.first, bracket.second);
		bstar.v0 = -dv0;
		bstar.omega *= std::exp(-0.5*dv0);
		if(comments){
			std::cout << "The value of omega is: " << bstar.omega << std::endl;
			std::cout << "The value of v0 is: " << bstar.v0 << std::endl;
		}
	}

	template<class Fill_Func>
	void fill_solution(int npoints_, const Fill_Func &func){
		npoints = npoints_;
		xstorage = new double[npoints+1];
		ystorage = new double[4*(npoints+1)];
		xarr = xstorage+1;
		yarr = ystorage+4;

		xstorage[0] = 0.0;
		for(int i = 0; i < npoints; ++i){
			xarr[i] = func(i);
			if(xarr[i] > xend){
				std::cout << "Attempting to find a solution beyond xend" << std::endl;
				throw("fill_solution error");
			}
		}

		for(int i = 0; i < 4*(npoints+1); ++i)
			ystorage[i] = 0.0;

		int ist = 0, iend;
		double *xint = xarr;
		double *yint = yarr;
		double xst_tmp, xend_tmp;

		while(xarr[ist] < xst)
			++ist;
		--ist;

		xint += ist;
		yint += 4*ist;
		xst_tmp = xint[0];
		xint[0] = xst;

		Boson_Star_RK::Origin_Asymptotic asym_orig(bstar);

		// Integrating
 		bstar.integrate_grid(npoints-ist-1, xint, yint);

 		// Filling values before xst with asymptotics
 		xint[0] = xst_tmp;
 		yint = yarr;
 		for(int i = 0; i <= ist; ++i){
 			asym_orig.evaluate(xarr[i], yint);
			yint += 4;
 		}
	}
};

#endif