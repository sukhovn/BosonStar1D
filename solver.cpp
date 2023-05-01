#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <string.h>
#include <time.h>
#include <boson_star.h>

void read_raw_domains(const std::string &file_name, std::vector<int> &nx, std::vector<double> &bound);

int main(int argc, char const *argv[]){
	std::string inp_name = "Setup_init";
	std::string out_folder = "Test_solution";
	std::string domains_file;

	for(int i = 1; i < argc; i++){
		if(strcmp(argv[i], "-in") == 0){
			inp_name = argv[i+1];
			i++;
		}
		else if(strcmp(argv[i], "-domains") == 0){
			domains_file = argv[i+1];
			i++;
		}
		else{
			std::cout << "Wrong key option received" << std::endl;
			return -1;
		}
	}

	std::vector<int> nx = {80, 20};
	std::vector<double> bound = {0.0, 20.0, 50.0};

	read_raw_domains("Setup_init/domains.dat", nx, bound);

	if(!domains_file.empty())
		read_raw_domains("../" + domains_file, nx, bound);

	Boson_Star star(nx, bound);
	star.set_doubling_begin({false, false, false});

	Boson_Star_Solver solver(&star);
	solver.load_spline(inp_name);
	solver.solve();

	std::string folder_name(out_folder);
	star.save(folder_name);
	solver.save_csv(folder_name + "/fields.csv", 0.0, 15.0, 1001);
	solver.save_function_to_file(folder_name + "/Residual", &Boson_Star_Equations::residual, 1, true);
	std::cout << std::endl;

	// nx = {70, 40};
	// bound = {0.0, 5.0, 10.0};
	nx = {15, 70, 20};
	bound = {0.0, 2.0, 15.0, 25.0};
	// nx = {25, 70, 25};
	// bound = {0.0, 1.0, 10.0, 20.0};

	Boson_Star star2(star);
	star2.regrid_solution(nx, bound);

	solver.change_star(&star2);
	solver.solve();
	star2.save("Test_solution2");
	solver.save_csv("Test_solution2/fields.csv", 0.0, 15.0, 1001);
	solver.save_function_to_file("Test_solution2/Residual", &Boson_Star_Equations::residual, 1, true);

	return 0;
}

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