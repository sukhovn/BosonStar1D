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

class Extrapolation: public Boson_Star_Solver{
	public:
		int sol_num;
		Boson_Star *sol_current;
		Boson_Star **sol_array; //Recorded solutions
		bool is_initialized;

		// double max_int_size;
		double stepper_increment; //Initial stepper increment
		std::vector<std::string> ps_length_components;

		//Tracking maximal compactness and maximal stable compactness
		double max_compactness;
		bool bound_null_geodesics;

	public:
		Extrapolation(const Boson_Star &sol, int num);
		~Extrapolation(){
			if(is_initialized)
				for(int i = 0; i < sol_num; i++) delete sol_array[i];
			delete [] sol_array;
			delete sol_current;
		}

		//Fills initial sol_num solutions for further extrapolation
		int initial_fill(const double dpar);

		//Extrapolates the solution to ps_length value param
		int extrapolate_solution(double *param_array, const double param);

		//Gradually increases ps_length value via extrapolation, saving star solutions parametes into saver along the way
		int param_change(double param_new, Solution_Saver &saver);

		//A wrap up for tracing ps_length
		int pslenght_tracer(double increment, double pslength_end, const std::string &folder_path);
};

#endif