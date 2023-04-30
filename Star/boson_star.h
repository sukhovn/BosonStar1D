#include <tuple>
#include <container.h>
#include <parameter_file.h>
#include <min_root.h>
#include <collocation_solver.h>
#include <equations.h>
#include <fstream>

#ifndef BOSON_STAR_H
#define BOSON_STAR_H

#define RAD_FRACTION 0.99

class Boson_Star: public Container, public Boson_Star_Equations{
private:
	double phi0;
	double msmass;
	double charge;
	double rad;
	double compactness;
	bool bound_null_geodesics;

	double residual_error;

	Parameter_File pfile;
	void read_parameters(const std::string &folder_name);
	void record_parameters(const std::string &folder_name);

public:
	double coordinate_komar_radius;
	double ps_length;

	//Standard constructor
	Boson_Star(const std::vector<int> &nxi, const std::vector<double> &boundi);

	void load(const std::string &folder_name){
		read_parameters(folder_name);
		Container::load(folder_name);
		return;
	}

	//Load from file constructor
	Boson_Star(const std::string &folder_name);

	//Copy constructor
	// Boson_Star(const Boson_Star &star);

	// Boson_Star &operator=(const Boson_Star &star);

	//Saves star solution
	void save(const std::string &folder_name){
		std::cout << "Saving the solution\n";
		Container::save(folder_name);
		record_parameters(folder_name);
	}

	friend class Boson_Star_Solver;
	friend class Solution_Saver;
	friend class Extrapolation;
};

class Boson_Star_Solver: private Linear_Fractional_Compactification, public Collocation_Solver<Linear_Fractional_Compactification, Boson_Star_Equations>{
private:
	Boson_Star *star;

public:
	Boson_Star_Solver(Boson_Star *stari):
		Linear_Fractional_Compactification(),
		Collocation_Solver<Linear_Fractional_Compactification, Boson_Star_Equations>(*this, stari, stari),
		star(stari)
	{
		const std::vector<double> &boundaries = star->boundaries();
		cc = 1.0/(boundaries[boundaries.size()-1] - boundaries[boundaries.size()-2]);
	};

	void load_spline(const std::string &folder_name);

	// void change_star(Boson_Star *stari){
	// 	star = stari;
	// 	eq_system = stari;
	// 	change_container(*stari);
	// }

	void find_star_parameters(void){
		Evaluator profile(*star);

		//Evaluating mass and charge
		double domain_start = star->boundaries()[0];
		double comp_bound = star->boundaries()[star->n_domains()-1];
		double domain_end = star->boundaries()[star->n_domains()];
		Evaluator mass_func(*star,2,[&profile, this, &comp_bound, &domain_end](int fn, double x){
			if(x == domain_end)
				return 0.0;

			//Evaluating fields
			double u = profile(x, 0);
			double v = profile(x, 1);
			double phi = profile(x, 2);
			double dphi = profile.derivative(x, 2);

			//Evaluating coordinate and rescaling factor;
			double dx = 1.0;
			if(x >= comp_bound){
				dx = dxc(x - comp_bound);
				x = comp_bound + xc(x - comp_bound);
			}

			if(fn == 0){
				double dens = std::exp(-v)*std::pow(star->omega*phi, 2.0) + std::exp(-u)*std::pow(dphi/dx, 2.0) + star->vpot(phi*phi);
				return 0.5*star->kappa*dens*x*x*dx;
			}
			if(fn == 1){
				double dens = std::exp(0.5*(u-v))*star->omega*phi*phi;
				return std::sqrt(2.0*star->kappa)*dens*x*x*dx;
			}
			std::cout << "Wrong density number" << std::endl;
			throw("integral calculation error");
		});
		mass_func.integrate(0.0);

		star->msmass = mass_func(domain_end, 0);
		star->charge = mass_func(domain_end, 1);

		//Evaluating radius
		Coordinate_Compactified<Linear_Fractional_Compactification> coord(*this, star);
		double xst_rad = domain_start, xend_rad = domain_end;
		while(xend_rad - xst_rad > 1.0e-6){
			double tmp = mass_func(0.5*(xst_rad+xend_rad), 0);
			if(tmp > RAD_FRACTION*star->msmass)
				xend_rad = 0.5*(xst_rad + xend_rad);
			else
				xst_rad = 0.5*(xst_rad + xend_rad);
		}
		star->rad = coord.coord(0.5*(xst_rad + xend_rad));

		//Evaluating compactness
		Evaluator mass_radius(*star,1,[&profile, this, &domain_end](int fn, double x){
			if(x == domain_end)
				return 0.0;

			//Evaluating fields
			double u = profile(x, 0);

			return 0.5*(1.0 - std::exp(-u));
		});

		//Finding maximal compactness
		auto lambda_comp = [&mass_radius](double x){return -mass_radius(x, 0);};
		Minimization max_compactness(lambda_comp);
		max_compactness.set_bracket(domain_start, 0.5*(domain_start+domain_end), domain_end);
		max_compactness.brent();
		star->compactness = -max_compactness.fmin;

		std::cout << "Misner-Sharp mass is: " << star->msmass << std::endl;
		std::cout << "Star charge is: " << star->charge << std::endl;
		std::cout << RAD_FRACTION * 100.0 << "% mass reached at the radius " << star->rad << std::endl;
		std::cout << "The solution has maximal compactness: " << star->compactness << std::endl;

		//Checking if bound null geodesics are present
		Evaluator effective_compactness(*star,1,[&profile, this, &domain_end, &comp_bound](int fn, double x){
			if(x == domain_end)
				return 2.0;

			//Evaluating fields
			double dv = profile.derivative(x, 1);

			//Evaluating coordinate and rescaling factor;
			double dx = 1.0;
			if(x >= comp_bound){
				dx = dxc(x - comp_bound);
				x = comp_bound + xc(x - comp_bound);
			}

			return 2.0 - dv*x/dx;
		});

		auto lambda_geod = [&effective_compactness](double x){return effective_compactness(x);};
		Minimization geodesic_condition(lambda_geod);
		geodesic_condition.set_bracket(domain_start, 0.5*(domain_start+domain_end), domain_end);
		geodesic_condition.brent();
		if(geodesic_condition.fmin < 0.0){
			std::cout << "The star solution contains bound null geodesics\n";
			// std::cout << geodesic_condition.fmin << std::endl;
			star->bound_null_geodesics = true;
		}
		else{
			std::cout << "The star solution does not contain bound null geodesics\n";
			// std::cout << geodesic_condition.fmin << std::endl;
			star->bound_null_geodesics = false;
		}
	}

	int solve(void){
		std::cout << "Solving star equations" << std::endl;
		int err = Collocation_Solver<Linear_Fractional_Compactification, Boson_Star_Equations>::solve();
		if(err) return err;

		//Calculating maximal residual values
		auto residual_data = find_function_maximum(&Boson_Star_Equations::residual);
		star->residual_error = std::get<0>(residual_data)[0];
		std::cout << "The residual maximum is: " << star->residual_error << std::endl;

		find_star_parameters();
		return 0;
	}

	friend class Extrapolation;
};

#endif