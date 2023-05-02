#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <string.h>
#include <time.h>
#include <boson_star.h>
#include <supplement.h>


int main(int argc, char const *argv[]){
	std::string inp_name = "../Setup_init";
	std::string domains_file;

	bool change_sigma = false;
	double sigma;

	bool save_out = false;
	std::string out_folder;

	bool load_spectral = false;

	for(int i = 1; i < argc; i++){
		if(strcmp(argv[i], "-init") == 0){
			inp_name = argv[i+1];
			i++;
		}
		else if(strcmp(argv[i], "-in") == 0){
			inp_name = argv[i+1];
			load_spectral = true;
			i++;
		}
		else if(strcmp(argv[i], "-domains") == 0){
			domains_file = argv[i+1];
			i++;
		}
		else if(strcmp(argv[i], "-sigma") == 0){
			change_sigma = true;
			sigma = atof(argv[i+1]);
			i++;
		}
		else if(strcmp(argv[i], "-out") == 0){
			save_out = true;
			out_folder = argv[i+1];
			i++;
		}
		else{
			std::cout << "Wrong key option received" << std::endl;
			return -1;
		}
	}

	std::vector<int> nx;
	std::vector<double> bound;

	if(!domains_file.empty())
		read_raw_domains(domains_file, nx, bound);
	else
		read_raw_domains(inp_name + "/domains.dat", nx, bound);

	Boson_Star star(nx, bound);
	star.set_doubling_begin({false, false, false});

	if(load_spectral){
		star.load(inp_name);
		if(!domains_file.empty()){
			std::cout << "Regridding starting solution to new domain configuration:" << std::endl;
			for(int i = 0; i < nx.size(); ++i)
				std::cout << nx[i] << ", [" << bound[i] << ',' << bound[i+1] << ']' << std::endl;
			std::cout << std::endl;
			star.regrid_solution(nx, bound);
		}
	}

	Boson_Star_Solver solver(&star);

	if(!load_spectral)
		solver.load_spline(inp_name);

	struct{
		std::string name = "sigma";

		double calculate_parameter(const Boson_Star *star){
			return std::sqrt(-4.0*star->poly_coef[0]/star->poly_coef[1]);
		}

		void update_parameter(Boson_Star *star, double sigma){
			double mass = std::sqrt(star->poly_coef[0]);
			star->poly_coef[1] = -4.0*mass*mass/(sigma*sigma);
			star->poly_coef[2] = 4.0*mass*mass/std::pow(sigma, 4.0);
		}
	} sigma_calculator;

	if(change_sigma)
		solver.parameter_change(sigma_calculator, sigma);
	else
		solver.solve();

	if(save_out){
		star.save(out_folder);
		solver.save_csv(out_folder + "/fields.csv", 0.0, 20.0, 1001);
		solver.save_function_to_file(out_folder + "/Residual", &Boson_Star_Equations::residual, 1, true);
	}

	return 0;
}