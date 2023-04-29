#include <tuple>
#include <container.h>
#include <parameter_file.h>
#include <min_root.h>
#include <collocation_solver.h>
#include <equations.h>

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

	int solve(void){
		std::cout << "Solving star equations" << std::endl;
		int err = Collocation_Solver<Linear_Fractional_Compactification, Boson_Star_Equations>::solve();
		if(err) return err;

		//Calculating maximal residual values
		auto residual_data = find_function_maximum(&Boson_Star_Equations::residual);
		star->residual_error = std::get<0>(residual_data)[0];
		std::cout << "The residual maximum is: " << star->residual_error << std::endl;

		return 0;
	}

	friend class Extrapolation;
};

#endif