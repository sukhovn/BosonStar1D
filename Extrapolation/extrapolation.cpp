#include <iostream>
#include <iomanip>
#include <fstream>
#include <new>
#include <math.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include <interpolation.h>
#include <extrapolation.h>

#define ZERO_CHECK(a) (a < 1.0e-5 ? 1.0e-5 : a)
#define MAX_FAILS 5

//Here we define parameter recorder

//This routine initializes a structure for parameter recorder, it creates or opens a folder with files with phi, w, msmass and rad
//If append != 0 it appends solutions to already existing files
Solution_Saver::Solution_Saver(const std::string &file_name, bool append) : is_on(false), records(0),
	save_parameters({"omega", "phi0", "mass", "charge", "radius", "compactness"})
{	
	//If append == 0 we are overwriting old files, if not we are appending new data to the old files

	std::string head_line = "ps_length";
	for(auto par : save_parameters)
		head_line = head_line + ',' + par;

	//Creating or opening output files
	std::ifstream f_test(file_name);
	if(f_test.good() == false) append = false;

	if(append)
		out.open(file_name, std::ofstream::app);
	else{
		out.open(file_name);
		out << head_line << std::endl;
	}

	out << std::setprecision(15) << std::fixed;
}

//This routine stores solution for parameter recorder
Solution_Saver& Solution_Saver::operator<< (const Boson_Star& star){
	if(is_on){
		std::cout << "Recording star solution\n";
		out << star.ps_length;
		for(std::string par : save_parameters)
			out << ',' << star.params.at(par);
		out << std::endl;
		records++;}
	return *this;
}


Extrapolation::Extrapolation(const Boson_Star &sol, int num, const std::string &_save_folder):
	Boson_Star_Solver(new Boson_Star(sol)), sol_num(num), is_initialized(false),
	ps_length_components({"omega", "mass", "phi0"}),
	save_folder(_save_folder)
{
	sol_current = star;
	sol_array = new Boson_Star*[num];

	stepper_increment = 0.01;
}

int Extrapolation::initial_fill(const double dpar){
	if(is_initialized){
		std::cout << "The extrapolation has already been initialized, aborting" << std::endl;
		return 1;
	}

	double mult_coef;

	if(dpar > 0.0){
		mult_coef = 1.0 + stepper_increment;
	}
	else{
		mult_coef = 1.0 - stepper_increment;
	}

	double omega_past = sol_current->omega, par_last = -1.0e5;
	int ord = 1, sol_record = 0;

	for(int i = 0; i < 70; ++i)
		std::cout << '*';
	std::cout << std::endl;
	std::cout << "Initializing extrapolation" << std::endl;;

	while(1){
		for(int i = 0; i < 70; ++i)
			std::cout << '-';
		std::cout << std::endl;
		std::cout << "Solving for omega = " << sol_current->omega << "\n";
		if(solve()){
			for(int i = 0; i < 70; ++i)
				std::cout << '-';
			std::cout << std::endl;
			mult_coef = std::sqrt(mult_coef);
			std::cout << "Changing multiplication coefficient to " << mult_coef << "\n";
			ord++;
		}
		else{
			omega_past = sol_current->omega;
			if(fabs(sol_current->ps_length - par_last) > fabs(dpar)){
				if(sol_record == 0)
					sol_current->ps_length = 0.0;
				else{
					double tmp, dlen = 0.0;
					for(auto par : ps_length_components){
						tmp = (sol_current->params.at(par) - sol_array[sol_record-1]->params.at(par))/ZERO_CHECK(sol_array[sol_record-1]->params.at(par));
						dlen += tmp*tmp;
					}
					sol_current->ps_length = sol_array[sol_record-1]->ps_length + std::sqrt(dlen);
				}

				std::cout << "Recording the solution " << sol_record+1 << " out of " << sol_num << " with pseudolengh " << sol_current->ps_length << " for extrapolation\n";
				sol_array[sol_record] = new Boson_Star(*sol_current);

				for(int i = 0; i < 70; ++i)
					std::cout << '-';
				std::cout << std::endl;

				sol_record++;

				if(sol_record > sol_num-1)
					break;
			}
		}
		if(ord > 20){
			printf("Minimal increment in omega reached, terminating\n");
			return 1;
			for(int i = 0; i < sol_record; i++) delete sol_array[i];
		}
		sol_current->omega = omega_past * mult_coef;
	}

	std::cout << "Extrapolation initialized\n";
	for(int i = 0; i < 70; ++i)
		std::cout << '*';
	std::cout << std::endl;

	is_initialized = true;

	return 0;
}

int Extrapolation::extrapolate_solution(double *param_array, const double param){
	if(!is_initialized){
		std::cout << "The extrapolation has not been initialized, aborting\n";
		return 1;
	}
	std::cout << "Attempting extrapolation to pseudolength " << param << "\n";

	double tmp, dyi;
	Boson_Star *sol_pt;
	std::vector<double> tmparr(sol_num);

	std::cout << "Extrapolating solution\n";

	for(int i = 0; i < sol_num; i++)
		tmparr[i] = sol_array[i]->omega;
	if(polint(param_array, tmparr.data(), sol_num, param, &tmp, &dyi) != 0){
		std::cout << "Could not extrapolate omega, extrapolation failed\n";
		return 1;
	}

	sol_current->omega = tmp;
	std::cout << "Extrapolated omega value is " << tmp << "\n";

	for(int i = 0; i < sol_current->size(); i++){
		for(int ii = 0; ii < sol_num; ii++)
			tmparr[ii] = (*sol_array[ii])[i];
		if(polint(param_array, tmparr.data(), sol_num, param, &tmp, &dyi) != 0){
			std::cout << "Could not extrapolate solution, extrapolation failed\n";
			return 1;
		}
		(*sol_current)[i] = tmp;
	}

	std::cout << "Extrapolation is done, getting a solution from extrapolation\n";
	if(solve() != 0){
		std::cout << "Could not get a solution from extrapolation\n";
		sol_current->Boson_Star::operator=(*sol_array[sol_num-1]);
		return 1;
	}

	double dlen = 0.0;
	for(auto par : ps_length_components){
		tmp = (sol_current->params.at(par) - sol_array[sol_num-1]->params.at(par))/ZERO_CHECK(sol_array[sol_num-1]->params.at(par));
		dlen += tmp*tmp;
	}
	sol_current->ps_length = sol_array[sol_num-1]->ps_length + std::sqrt(dlen);
	std::cout << "The pseudolength is " << sol_current->ps_length << "\n";

	//Interchanging solutions
	delete sol_array[0];

	for(int i = 0; i < sol_num-1; i++){
		sol_array[i] = sol_array[i+1];
	}

	sol_array[sol_num-1] = new Boson_Star(*sol_current);

	return 0;
}

int Extrapolation::regrid_solutions(const std::vector<int> &_nx, const std::vector<double> &_bound){
	for(int i = 0; i < 70; ++i)
		std::cout << '*';
	std::cout << std::endl;

	std::cout << "New domains configuration:" << std::endl;
	for(int i = 0; i < _nx.size(); ++i)
		std::cout << _nx[i] << ", [" << _bound[i] << ',' << _bound[i+1] << ']' << std::endl;

	Boson_Star **sol_backup = new Boson_Star*[sol_num];
	for(int i = 0; i < sol_num; ++i){
		std::cout << "\nRegridding solution: " << i << std::endl;
		sol_backup[i] = new Boson_Star(*sol_array[i]);
		sol_array[i]->regrid_solution(_nx, _bound);
		change_star(sol_array[i]);
		if(solve() != 0){
			std::cout << "Regridding failed" << std::endl;
			for(int j = 0; j <= i; ++j){
				delete sol_array[j];
				sol_array[j] = sol_backup[j];
			}
			change_star(sol_current);
			delete [] sol_backup;
			return 1;
		}
	}
	std::cout << "\nRegridding was successful" << std::endl;
	for(int i = 0; i < 70; ++i)
		std::cout << '*';
	std::cout << std::endl;

	delete sol_current;
	sol_current = new Boson_Star(*sol_array[sol_num-1]);
	change_star(sol_current);
	return 0;
}

int Extrapolation::param_change(double param_new, Solution_Saver &saver){
	bool sample_save = true;
	double param_last = sol_current->ps_length;
	double mult_coef, mult_coef_init, direction, param_try;
	std::vector<double> param_array(sol_num);
	int fails = 1, easy_steps = 0;
	bool max_comp_change = false;

	if(param_new > param_last){
		direction = 1.0;
		mult_coef = mult_coef_init = 1.0 + stepper_increment;
	}
	else{
		direction = -1.0;
		mult_coef = mult_coef_init = 1.0 - stepper_increment;
	}

	while(1){
		param_try = param_last * mult_coef;
		for(int i = 0; i < sol_num; i++)
			param_array[i] = sol_array[i]->ps_length;

		for(int i = 0; i < 70; ++i)
			std::cout << '-';
		std::cout << std::endl;
		if(extrapolate_solution(param_array.data(), param_try) != 0){
			for(int i = 0; i < 70; ++i)
				std::cout << '-';
			std::cout << std::endl;

			for(int i = 0; i < 70; ++i)
				std::cout << '+';
			std::cout << std::endl;
			mult_coef = std::sqrt(mult_coef);
			std::cout << "Changing multiplication coefficient to " << mult_coef << "\n";
			fails++;
			for(int i = 0; i < 70; ++i)
				std::cout << '+';
			std::cout << std::endl;
		}
		else{
			saver << *sol_current;
			for(int i = 0; i < 70; ++i)
				std::cout << '-';
			std::cout << std::endl;

			if(iterations > 2){
				mult_coef = std::sqrt(mult_coef);
				std::cout << "Too many iterations, changing multiplication coefficient to " << mult_coef << "\n";
			}
			param_last = sol_current->ps_length;

			if(iterations == 1) easy_steps++;
			else easy_steps = 0;

			for(auto cond = conditions.begin(); cond != conditions.end(); ++cond)
				(*cond)->check_condition(*this);
		}

		//No further extrapolation termination
		if(fails > MAX_FAILS){
			for(int i = 0; i < 70; ++i)
				std::cout << '+';
			std::cout << std::endl;
			std::cout << "Minimal increment is reached, terminating\n";
			for(int i = 0; i < 70; ++i)
				std::cout << '+';
			std::cout << std::endl;
			return 1;
		}

		//If there were too many easy steps we increase the stepsize
		if(easy_steps > 10){
			if(mult_coef < mult_coef_init - 1.0e-6){
				mult_coef = std::pow(mult_coef, 2.0);
				std::cout << "Increasing interpolation step, changing multiplication coefficient to " << mult_coef << std::endl;
			}
			easy_steps = 0;
			fails = 0;
		}

		//Termination condition
		if((param_new - param_last * mult_coef) * direction < 0.0){
			if((param_new - param_last) * direction < 0.0){
				break;
			}
			param_try = param_new;
			for(int i = 0; i < sol_num; i++)
				param_array[i] = sol_array[i]->ps_length;

			for(int i = 0; i < 70; ++i)
				std::cout << '-';
			std::cout << std::endl;
			if(extrapolate_solution(param_array.data(), param_try) == 0){
				saver << *sol_current;
				for(int i = 0; i < 70; ++i)
					std::cout << '-';
				std::cout << std::endl;
				break;
			}
			else{
				for(int i = 0; i < 70; ++i)
					std::cout << '-';
				std::cout << std::endl;

				fails++;
				while((param_new - (param_last * mult_coef)) * direction < 0.0)
					mult_coef = std::sqrt(mult_coef);
			}
		}
	}

	std::cout << "The extrapolation is done, omega is: " << sol_current->omega << "\n";
	return 0;
}

int Extrapolation::pslenght_tracer(double increment, double pslength_end){
	std::string file_name;
	stepper_increment = increment;

	//Creating output folder
	struct stat sttus = {0};
	if (stat(save_folder.c_str(), &sttus) == -1) {
		mkdir(save_folder.c_str(), 0700);
	}

	file_name = save_folder + "/tracing.csv";
	Solution_Saver saver(file_name, false);
	saver.on();
	
	std::cout << "Tracing pslength\n";
	if(param_change(pslength_end, saver)){
		std::cout << "\nTracing interrupted" << std::endl;
		for(auto cond = conditions.begin(); cond != conditions.end(); ++cond)
			(*cond)->final_action(*this, false);
		return 1;
	}

	std::cout << "\nTracing is finished" << std::endl;

	for(auto cond = conditions.begin(); cond != conditions.end(); ++cond)
		(*cond)->final_action(*this, true);

	save_current_solution(save_folder + "/Final", 0);

	return 0;
}