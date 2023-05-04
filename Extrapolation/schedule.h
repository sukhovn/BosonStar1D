#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
#include <utility>
#include <algorithm>
#include <map>
#include <extrapolation.h>

class Tracing_Schedule: public Extrapolation_Condition{
public:
	std::map<double, std::pair<std::vector<int>, std::vector<double>>>::iterator curr;
	std::map<double, std::pair<std::vector<int>, std::vector<double>>> schedule;

	Tracing_Schedule(): curr(schedule.end()) {};

	void load_schedule(const std::string &file_name){
		double curr_pslength, tmp;
		std::vector<int> curr_nx;
		std::vector<double> curr_bound;
		int read_state = 0;

		std::ifstream file(file_name);
		if(!file){
			throw("file not found");
		}
		std::string line;
		std::string token;
		std::string::size_type sztp;

		while(std::getline(file, line)){
			if(line.size() == 0)
				continue;
			std::istringstream iss(line);
			std::getline(iss, token, ':');
			if(token.compare("Pseudolength") == 0){
				if(read_state != 0){
					std::cout << "Corrupt schedule file, expected domain entries instead of line" << std::endl;
					std::cout << line << std::endl;
					throw("corrupt file");
				}
				std::getline(iss, token, ':');
				curr_pslength = std::stod(token, &sztp);
				read_state = 1;
			}

			iss.seekg(0, iss.beg);
			std::getline(iss, token, ',');
			try{
				tmp = std::stoi(token, &sztp);
			}
			catch(std::invalid_argument const& ex){
				continue;
			}

			iss.seekg(0, iss.beg);
			if(read_state == 1){
				while(std::getline(iss, token, ',')){
					if(!std::all_of(token.begin(), token.end(), [](char c){return ::isdigit(c) || c == ' ';})){
						std::cout << "Corrupt schedule file, expected all integers in line" << std::endl;
						std::cout << line << std::endl;
						throw("corrupt file");
					}
					curr_nx.push_back(std::stoi(token, &sztp));
				}
				read_state = 2;
			}
			else if(read_state == 2){
				while(std::getline(iss, token, ','))
					curr_bound.push_back(std::stod(token, &sztp));
				
				//Recording schedule entry
				bool correct_domains = (curr_nx.size()+1 == curr_bound.size());
				for(int i = 0; i < curr_nx.size(); ++i){
					if(curr_bound[i] >= curr_bound[i+1]){
						correct_domains = false;
						break;
					}
				}
				if(!correct_domains){
					std::cout << "Incorrect domain configuration" << std::endl;
					
					for(int i = 0; i < curr_nx.size(); ++i){
						if(i) std::cout << ',';
						std::cout << curr_nx[i];
					}
					std::cout << std::endl;

					for(int i = 0; i < curr_bound.size(); ++i){
						if(i) std::cout << ',';
						std::cout << curr_bound[i];
					}
					std::cout << std::endl;
					
					throw("corrupt file");
				}

				schedule[curr_pslength] = std::make_pair(curr_nx, curr_bound);
				curr_nx.clear();
				curr_bound.clear();
				read_state = 0;
			}
			else{
				std::cout << "Corrupt schedule file, expected next pseudolength entry instead of line" << std::endl;
				std::cout << line << std::endl;
				throw("corrupt file");
			}
		}
	}

	void setup_scheduling(double ps_length){
		curr = schedule.begin();
		while(curr != schedule.end()){
			if(curr->first > ps_length)
				break;
			++curr;
		}
	}

	void check_condition(Extrapolation &extra){
		if(curr == schedule.end()) return;
		if(extra.sol_current->ps_length > curr->first){
			extra.regrid_solutions(curr->second.first, curr->second.second);
			++curr;
		}
	};

	void final_action(Extrapolation &extra, bool normal_exit){return;}
	friend std::ostream &operator<<(std::ostream &os, const Tracing_Schedule &sch);
};

std::ostream &operator<<(std::ostream &os, const Tracing_Schedule &sch){
	bool first;
	for(auto el = sch.schedule.begin(); el != sch.schedule.end(); ++el){
		if(el != sch.schedule.begin()) os << '\n';
		os << "Pseudolength: " << el->first << std::endl;
		os << "Domains:\n";

		first = true;
		for(int nx : el->second.first){
			if(!first) os << ',';
			os << nx;
			if(first) first = false;
		}
		os << std::endl;

		first = true;
		for(double bound : el->second.second){
			if(!first) os << ',';
			os << bound;
			if(first) first = false;
		}
		os << std::endl;
	}
	return os;
}