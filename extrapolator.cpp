#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <string.h>
#include <time.h>
#include <boson_star.h>
#include <extrapolation.h>

static bool isNumber(const std::string& str)
{
	bool had_point = false;
	auto c = str.begin();
    if(*c == '-') ++c;
    for(; c != str.end(); ++c){
        if(*c == '.'){
        	if(had_point) return false;
        	else had_point = true;
        }
        else if (std::isdigit(*c) == 0) return false;
    }
    return true;
}

void read_raw_domains(const std::string &file_name, std::vector<int> &nx, std::vector<double> &bound){
	std::ifstream file(file_name);
    if(!file.is_open()){
		std::cout << "No domain file found" << std::endl;
		throw("Read fail");
	}
    
	std::string line, token; 
	std::string::size_type sztp;

	nx.clear();
	bound.clear();

	for(int i = 0; i < 2; ++i){
		std::getline(file, line);
		std::istringstream iss(line);

		while(iss.peek() != EOF){
			std::getline(iss, token, ' ');
			if(isNumber(token)){
				if(i) bound.push_back(std::stoi(token, &sztp));
				else nx.push_back(std::stoi(token, &sztp));
			}
			else{
				std::cout << "Wrong key loaded while reading domain structure" << std::endl;
				throw("Wrong domain structure");
			}
		}
	}

	if(nx.size() < 1 || nx.size() + 1 != bound.size()){
		std::cout << "Wrong domain structure" << std::endl;
		throw("Wrong domain structure");
	}
}

int main(int argc, char const *argv[]){
	double sigma = 1.0;
	double increment = 0.01;
	double ps_end = 0.5;
	std::string out_folder_name = "Tracing";

	for(int i = 1; i < argc; i++){
		if(strcmp(argv[i], "-sigma") == 0){
			sigma = atof(argv[i+1]);
			i++;
		}
		if(strcmp(argv[i], "-ps_end") == 0){
			ps_end = atof(argv[i+1]);
			i++;
		}
		if(strcmp(argv[i], "-increment") == 0){
			increment = atof(argv[i+1]);
			i++;
		}
	}

	std::string inp_name = "Setup_init";
	std::vector<int> nx;
	std::vector<double> bound;

	read_raw_domains(inp_name + "/domains.dat", nx, bound);

	Boson_Star star(nx, bound);
	star.set_doubling_begin({false, false, false});

	Boson_Star_Solver solver(&star);
	solver.load_spline(inp_name);

	struct{
		std::string name = "sigma";

		double calculate_parameter(const Boson_Star *star){
			if(star->poly_coef[1] == 0.0)
				return std::numeric_limits<double>::infinity();
			return std::sqrt(-4.0*star->poly_coef[0]/star->poly_coef[1]);
		}

		void update_parameter(Boson_Star *star, double sigma){
			double mass = std::sqrt(star->poly_coef[0]);
			star->poly_coef[1] = -4.0*mass*mass/(sigma*sigma);
			star->poly_coef[2] = 4.0*mass*mass/std::pow(sigma, 4.0);
		}
	} sigma_calculator;

	double sigma_solution = sigma_calculator.calculate_parameter(&star);
	
	if(std::isinf(sigma_solution)){
		std::cout << "Sigma parameter is switched off, sigma adjustment not possible" << std::endl;
		solver.solve();
	}
	else if(std::fabs(sigma_solution - sigma) > 1.0e-10)
		solver.parameter_change(sigma_calculator, sigma);
	else
		solver.solve();

	#include <conditions.h>

	Extrapolation extra(star, 5, out_folder_name);
	extra.conditions.push_back(&check_maximum_mass);
	extra.conditions.push_back(&check_bound_null_geodesics);
	// extra.conditions.push_back(&change_domain_size);
	extra.pslenght_tracer(increment, ps_end);

	return 0;
}