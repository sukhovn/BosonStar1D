#include <iostream>
#include <vector>

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
			std::getline(iss, token, ',');
			if(isNumber(token)){
				if(i) bound.push_back(std::stod(token, &sztp));
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