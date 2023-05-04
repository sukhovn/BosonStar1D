struct: public Extrapolation_Condition{
	double previous_local_maximum_observed = 0.0;
	double global_maximum = 0.0;

	void check_previous_max_save(const std::string &max_path){
		struct stat sttus = {0};
		std::string::size_type sztp;
		std::string line, token;
		if(stat(max_path.c_str(), &sttus) == -1) return;
		std::ifstream inp_file(max_path + "/parameters.dat");
		if(!inp_file.is_open()){
			std::cout << "No parameters file found in previous maximal mass solution" << std::endl;
			throw("reading error");
		}
		while(std::getline(inp_file, line)){
			std::istringstream iss(line);
			std::getline(iss, token, ':');
			if(token.compare("The Misner-Sharp mass is") == 0){
				std::getline(iss, token, ':');
				global_maximum = std::stod(token, &sztp);
				break;
			}
		}
	}

	void check_condition(Extrapolation &extra){
		if(previous_local_maximum_observed > extra.sol_current->params.at("mass")){
			if(previous_local_maximum_observed > global_maximum){
				global_maximum = previous_local_maximum_observed;
				std::cout << "Local mass maximum reached at " << previous_local_maximum_observed << std::endl;
				extra.save_solution(extra.save_folder + "/Maximal_mass",1);
			}
		}
		else
			previous_local_maximum_observed = extra.sol_current->params.at("mass");
	}

	void final_action(Extrapolation &extra, bool normal_exit){
		if(previous_local_maximum_observed > global_maximum){
			global_maximum = previous_local_maximum_observed;
			std::cout << "A maximum in mass was reached at the end: " << global_maximum << std::endl;
			extra.save_current_solution(extra.save_folder + "/Maximal_mass");
		}
		else{
			std::cout << "A maximum in mass was reached at: " << global_maximum << std::endl;
		}
	}

} check_maximum_mass;

struct: public Extrapolation_Condition{
	bool bound_null_geodesics = false;

	void check_condition(Extrapolation &extra){
		if(bound_null_geodesics)
			return;

		if(extra.sol_current->params.at("bound_null_geodesics") == 1.0)
			bound_null_geodesics = true;
	}

	void final_action(Extrapolation &extra, bool normal_exit){
		if(bound_null_geodesics){
			std::cout << "Stars with bound null geodesics detected" << std::endl;
		}
		else{
			std::cout << "No stars with bound null geodesics detected" << std::endl;
		}
	}

} check_bound_null_geodesics;

struct: public Extrapolation_Condition{
	double residual_max;
	double residual_max_ps_length;

	void check_condition(Extrapolation &extra){
		if(extra.sol_current->residual_error > residual_max){
			residual_max = extra.sol_current->residual_error;
			residual_max_ps_length = extra.sol_current->ps_length;
		}
	}

	void final_action(Extrapolation &extra, bool normal_exit){
		std::cout << "Maximal residual: " << residual_max << std::endl;
		std::cout << "Maximal residual reached at pseudolength: " << residual_max_ps_length << std::endl;
	}

} maximal_residual;