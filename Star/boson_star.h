#include <tuple>
#include <container.h>
#include <parameter_file.h>
#include <min_root.h>
#include <collocation_solver.h>
#include <equations.h>
#include <fstream>
#include <unordered_map>

#ifndef BOSON_STAR_H
#define BOSON_STAR_H

#define RAD_FRACTION 0.99

class Boson_Star: public Container, public Boson_Star_Equations{
private:
	std::unordered_map<std::string, double> params;

	double residual_error;

	Parameter_File pfile;
	void read_parameters(const std::string &folder_name);
	void record_parameters(const std::string &folder_name);

public:
	double coordinate_radius;
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

	//Saves star solution
	void save(const std::string &folder_name){
		std::cout << "Saving the solution\n";
		Container::save(folder_name);
		record_parameters(folder_name);
	}

	void regrid_solution(const std::vector<int> &_nx, const std::vector<double> &_bound){
		double comp_bound = boundaries()[n_domains()-1];
		double domain_end = boundaries()[n_domains()];
		Linear_Fractional_Compactification compact_old(1.0/(domain_end - comp_bound));
		Coordinate_Compactified<Linear_Fractional_Compactification> coord_old(compact_old, this);

		Evaluator profile(*this);
		adjust_partitioning(_nx, _bound);

		comp_bound = boundaries()[n_domains()-1];
		domain_end = boundaries()[n_domains()];
		Linear_Fractional_Compactification compact_new(1.0/(domain_end - comp_bound));
		Coordinate_Compactified<Linear_Fractional_Compactification> coord_new(compact_new, this);

		auto lambda = [&profile, &coord_old, &coord_new, &domain_end](int fn, double x){
			if(x == domain_end)
				return 0.0;
			return profile(coord_old.comp_coord(coord_new.coord(x)), fn);
		};

		fill_function(lambda);
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

	void change_star(Boson_Star *star_){
		star = star_;
		change_container(star_, star_);
		const std::vector<double> &boundaries = star->boundaries();
		cc = 1.0/(boundaries[boundaries.size()-1] - boundaries[boundaries.size()-2]);
	}

	void find_star_parameters(void);

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

	template<class Parameter_Calculator>
	int parameter_change(Parameter_Calculator &calc, double par_new){
		double par = calc.calculate_parameter(star);
		double par_past = par;
		double mult_coef, direction;
		int ord = 1;

		if(par_new > par){
			direction = 1.0;
			mult_coef = 1.1;
		}
		else{
			direction = -1.0;
			mult_coef = 0.9;
		}

		std::cout << "Starting at multiplication coefficient " << mult_coef << std::endl;

		while(1){
			std::cout << std::endl << "Solving for " + calc.name + ": " << par << std::endl;
			if(solve() != 0){
				if(par == par_past){
					std::cout << "Star incorrectly initialized" << std::endl;
					return -1;
				}
				mult_coef = std::sqrt(mult_coef);
				par = par_past * mult_coef;
				calc.update_parameter(star, par);
				std::cout << std::endl << "Changing multiplication coefficient to " << mult_coef << std::endl;
				ord++;
			}
			else{
				par_past = par;
				par *= mult_coef;
				calc.update_parameter(star, par);
			}
			if(ord > 10){
				std::cout << std::endl << "Minimal increment in " + calc.name + " reached, terminating" << std::endl;
				return 1;
			}
			if((par_new - par) * direction < 0.0){
				par = par_new;
				calc.update_parameter(star, par);
				std::cout << std::endl << "Solving for " + calc.name + " = " << par << std::endl;
				if(solve() == 0)
					break;
				else{
					while((par_new - par_past * mult_coef) * direction < 0.0)
						mult_coef = std::sqrt(mult_coef);
					par = par_past * mult_coef;
					calc.update_parameter(star, par);
				}
			}
		}

		std::cout << std::endl << "We reached " + calc.name + " value of " << par << std::endl;
		return 0;
	}

	friend class Extrapolation;
};

#endif