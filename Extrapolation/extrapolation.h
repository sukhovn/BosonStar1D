#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdlib.h>
#include <math.h>

#include <boson_star.h>

#ifndef EXTRAPOLATION_H
#define EXTRAPOLATION_H

class Solution_Saver{
private:
	std::ofstream out;

	bool is_on; //If this parameter is zero, no new records are added, if not, records are being added
	int records;

	std::vector<std::string> save_parameters;

public:
	//This routine initializes a structure for parameter recorder, it creates or opens a folder with files with phi, w, msmass and rad
	//If append != 0 it appends solutions to already existing files
	Solution_Saver(const std::string &file_name, bool append);

	//Turning on and off recording
	void on(void){
		std::cout << "Recording star solutions parameters\n";
		is_on = true;};
	void off(void){
		std::cout << "Stopping recording star solutions parameters\n";
		is_on = false;};

	Solution_Saver& operator<< (const Boson_Star& sol);
};

class Extrapolation_Condition;

class Extrapolation: public Boson_Star_Solver{
	Solution_Saver *saver;

	bool is_initialized;

	std::vector<std::string> ps_length_components;
	double mult_coef;
	double mult_coef_max;

public:
	int sol_num;
	Boson_Star *sol_current;
	Boson_Star **sol_array; //Recorded solutions

	std::vector<Extrapolation_Condition*> conditions;
	double stepper_increment; //Initial stepper increment

	const std::string save_folder;

public:
	Extrapolation(const Boson_Star &sol, int num, const std::string &_save_folder);

	~Extrapolation(){
		if(is_initialized){
			for(int i = 0; i < sol_num; i++)
				delete sol_array[i];
			delete saver;
		}
		delete [] sol_array;
		delete sol_current;
	}

	void save_current_solution(const std::string &folder_name){
		sol_current->save(folder_name);

		double radius = 0.0;
		while(radius < 100.0 && radius < sol_current->params.at("radius"))
			radius += 5.0;

		save_csv(folder_name + "/fields.csv", 0.0, radius, 1001);
		save_function_to_file(folder_name + "/Residual", &Boson_Star_Equations::residual, 1, true);
	}

	void save_solution(const std::string &folder_name, int num){
		if(!is_initialized){
			std::cout << "The extrapolation has not been initialized, solution cannot be saved\n" << std::endl;
			return;
		}
		if(num >= sol_num, num < 0){
			std::cout << "Wrong solution number provided for saving" << std::endl;
			return;
		}

		sol_array[sol_num-num-1]->save(folder_name);

		double radius = 0.0;
		while(radius < 100.0 && radius < sol_array[sol_num-num-1]->params.at("radius"))
			radius += 5.0;

		change_star(sol_array[sol_num-num-1]);
		save_csv(folder_name + "/fields.csv", 0.0, radius, 1001);
		save_function_to_file(folder_name + "/Residual", &Boson_Star_Equations::residual, 1, true);
		change_star(sol_current);
	}

	//Fills initial sol_num solutions for further extrapolation
	int initial_fill(const double dpar);

	//Loads previous extrapolator state
	void load_tracing_state(const std::string &folder_name);

	//Regrids all stored solution to the new grid
	int regrid_solutions(const std::vector<int> &_nx, const std::vector<double> &_bound);

	//Extrapolates the solution to ps_length value param
	int extrapolate_solution(double *param_array, const double param, bool swap_solutions);

	//Gradually increases ps_length value via extrapolation, saving star solutions parametes into saver along the way
	int param_change(double param_new);

	//A wrap up for tracing ps_length
	int pslenght_tracer(double pslength_end);

	//Saves the current extrapolator state
	void save_tracing_state(const std::string &folder_name);
};

struct Extrapolation_Condition{
	virtual void check_condition(Extrapolation &extra) = 0;
	virtual void final_action(Extrapolation &extra, bool normal_exit) = 0;
};

#endif