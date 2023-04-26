#include <container.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include <fftw3.h>

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

int count_func_files(const std::string &folder_name){
	std::string path;
	int nfunc_load = 0;
	struct stat sttus = {0};
	while(1){
		path = folder_name + "/f" + std::to_string(nfunc_load) + ".dat";
		if (stat(path.c_str(), &sttus) == -1) break;
		else nfunc_load++;
	}
	return nfunc_load;
}

void Container::read_domain_structure(const std::string &file_name, int fn, bool fill){
	std::vector<int> nx_inp;
	std::vector<double> bound_inp;
	std::ifstream file(file_name);
	std::string line, token, par; 
	std::getline(file, line);
	std::istringstream iss(line);
	std::string::size_type sztp;

	bool first = false;
	while(iss.peek() != EOF){
		std::getline(iss, token, ',');
		if(isNumber(token)){
			if(nx_inp.size() == bound_inp.size())
				bound_inp.push_back(std::stod(token, &sztp));
			else
				nx_inp.push_back(std::stoi(token, &sztp));
		}
		else{
			if(token.size() != 1){
				std::cout << "Wrong key loaded while reading domain structure" << std::endl;
				throw("Wrong domain structure");
			}
			par.push_back(token[0]);
		}
	}

	if(par.size() != 2){
		std::cout << "Wrong number of keys loaded while reading domain structure" << std::endl;
		throw("Wrong domain structure");
	}

	if(fill){
		if(par[0] == 'n') parity_st.resize(0);
		else parity_st.resize(nfunc);
		if(par[1] == 'n') parity_end.resize(0);
		else parity_end.resize(nfunc);
		bound = bound_inp;
		nx = nx_inp;
	}
	else{
		if(par[0] == 'n' && !parity_st.empty() || par[1] == 'n' && !parity_end.empty() || 
			par[0] != 'n' && parity_st.empty()  || par[1] != 'n' && parity_end.empty() || 
			nx != nx_inp || bound != bound_inp)
		{
			std::cout << "Different data files have different domain structure" << std::endl;
			throw("Domain structure mismatch");
		}
	}

	if(!parity_st.empty()){
		if(par[0] == 'e') parity_st[fn] = false;
		else if(par[0] == 'o') parity_st[fn] = true;
		else{
			std::cout << "Wrong key loaded while reading domain structure" << std::endl;
			throw("Wrong domain structure");
		}
	}

	if(!parity_end.empty()){
		if(par[1] == 'e') parity_end[fn] = false;
		else if(par[1] == 'o') parity_end[fn] = true;
		else{
			std::cout << "Wrong key loaded while reading domain structure" << std::endl;
			throw("Wrong domain structure");
		}
	}
}

void Container::save(const std::string &folder_name) const{
	struct stat sttus = {0};
	if (stat(folder_name.c_str(), &sttus) == -1) {
		mkdir(folder_name.c_str(), 0700);
	}

	std::ofstream out;
	std::string file_name;

	for(int fn = 0; fn < nfunc; fn++){
		file_name = folder_name + "/f" + std::to_string(fn) + ".dat";
		out.open(file_name);
		out << std::setprecision(15) << std::fixed;

		for(int dm = 0; dm < ndomains; dm++)
			out << bound[dm] << "," << nx[dm] << ",";
		out << bound[ndomains] << ",";

		if(parity_st.empty())
			out << "n";
		else{
			if(parity_st[fn])
				out << "o";
			else
				out << "e";
		}	
		out << ",";
		if(parity_end.empty())
			out << "n";
		else{
			if(parity_end[fn])
				out << "o";
			else
				out << "e";
		}
			
		out << std::endl;

		out << data[fn]; 
		for(int i = fn+nfunc; i < sz; i+=nfunc)
			out << "," << data[i];
		out.close();
	}
}

void Container::load(const std::string &folder_name){
	nfunc = count_func_files(folder_name);
	parity_st.resize(nfunc);
	parity_end.resize(nfunc);

	std::string path, line, token;
	std::ifstream file;
	std::string::size_type sztp;

	for(int fn = 0; fn < nfunc; ++fn){
		path = folder_name + "/f" + std::to_string(fn) + ".dat";
		read_domain_structure(path, fn, fn == 0);
		if(fn == 0)
			reallocate_partitioning(nx);

		file.open(path);
		file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
		std::getline(file, line);
		std::istringstream iss(line);
		for(int i = fn; i < sz; i+=nfunc){
			std::getline(iss, token, ',');
			data[i] = std::stod(token, &sztp);
		}
		file.close();
	}
}

void Evaluator::array_fft(void){
	double *ptr;
	
	int fft_size;
	double fft_coef;
	fftw_plan planCOS;
	fftw_r2r_kind kind = FFTW_REDFT00;

	for(int dm = 0; dm < ndomains; ++dm){
		ptr = data + index_st[dm];
		fft_coef = nx[dm]-1;

		planCOS = fftw_plan_many_r2r(1, &(nx[dm]), nfunc, ptr, NULL, nfunc, 1, ptr, NULL, nfunc, 1, &kind, FFTW_ESTIMATE);

		for(int i = 0; i < nfunc*nx[dm]; ++i) ptr[i] /= fft_coef;

		fftw_execute_r2r(planCOS, ptr, ptr);

		fftw_destroy_plan(planCOS);
	}
}

Evaluator::Evaluator(const Container &cont): 
	nfunc(cont.nfunc), ndomains(cont.ndomains), nx(cont.ndomains), bound(cont.bound), doubling_init(!cont.parity_st.empty()), sz(0)
{	
	for(int dm = 0; dm < ndomains; ++dm){
		index_st.push_back(sz);
		if((dm == 0 && cont.parity_st.size()) || (dm == ndomains-1 && cont.parity_end.size()))
			sz += nfunc*(nx[dm] = 2*cont.nx[dm]);
		else sz += nfunc*(nx[dm] = cont.nx[dm]);
	}

	if(cont.parity_st.size()) bound[0] = 2*bound[0] - bound[1];
	if(cont.parity_end.size()) bound[ndomains] = 2*bound[ndomains] - bound[ndomains-1];

	data = new double[sz];

	int indx = 0;
	for(int dm = 0; dm < ndomains; ++dm){
		if(dm == 0 && cont.parity_st.size()){
			for(int i = nfunc*(cont.nx[0]-1); i >= 0; i -= 2*nfunc){
				for(int fn = 0; fn < nfunc; ++fn){
					if(cont.parity_st[fn]) data[indx++] = -cont.data[i++];
					else data[indx++] = cont.data[i++];
				}
			}
		}
		
		for(int i = 0; i < nfunc*cont.nx[dm]; ++i)
			data[indx++] = cont.data[cont.index_st[dm]+i];
		
		if(dm == ndomains-1 && cont.parity_end.size()){
			for(int i = cont.sz-nfunc; i >= cont.index_st[dm]; i -= 2*nfunc){
				for(int fn = 0; fn < nfunc; ++fn){
					if(cont.parity_end[fn]) data[indx++] = -cont.data[i++];
					else data[indx++] = cont.data[i++];
				}
			}
		}
	}

	array_fft();
}

double Evaluator::operator()(double coord, int fn) const{
	if(coord < bound[0] || coord > bound[ndomains]){
		std::cout << "Coordinate is out of bounds" << std::endl;
		throw("Coordinate is out of bounds");
	}
	int dm = 0;
	while(coord > bound[dm+1]) ++dm;

	double d = 0.0, dd = 0.0, sv, y, y2;
	int ind;

	y2 = 2.0 * (y = (bound[dm] + bound[dm+1] - 2.0 * coord)/(bound[dm+1] - bound[dm])); //Change of variable.
	
	for(ind = index_st[dm]+nfunc*(nx[dm]-1)+fn; ind >= index_st[dm]+nfunc; ind -= nfunc){ //Clenshaw’s recurrence.
		sv = d;
		d = y2 * d - dd + data[ind];
		dd = sv;
	}

	return y * d - dd + 0.5 * data[ind]; //Last step is different.
}

double Evaluator::derivative(double coord, int fn) const{
	if(coord < bound[0] || coord > bound[ndomains]){
		std::cout << "Coordinate is out of bounds" << std::endl;
		throw("Coordinate is out of bounds");
	}
	int dm = 0;
	while(coord > bound[dm+1]) ++dm;

	double d = 0.0, dd = 0.0, sv, y, y2;
	int ind = index_st[dm]+nfunc*nx[dm]+fn;

	y2 = 2.0 * (y = (bound[dm] + bound[dm+1] - 2.0 * coord)/(bound[dm+1] - bound[dm])); //Change of variable.
	
	for(int j = nx[dm]-1; j > 1; j--){ //Clenshaw’s recurrence.
		sv = d;
		d = (y2 * (j+1) * d - (j+2) * dd)/(1.0*j) + data[ind -= nfunc];
		dd = sv;
	}

	return (4.0 * y * d - 3.0 * dd + data[ind -= nfunc]) * 2.0 / (bound[dm] - bound[dm+1]); //Last step is different.
}

void Evaluator::differentiate(void){
	double *sp1 = new double[nfunc];

	double length, tmp;
	int indx = sz-1;
	for(int dm = ndomains-1; dm >= 0; --dm){
		indx -= nfunc;
		length = bound[dm+1] - bound[dm];
		for(int fn = 0; fn < nfunc; ++fn){
			sp1[fn] = data[indx];
			data[indx] = - 4.0*(nx[dm]-1)*data[indx+nfunc]/length;
			data[indx+nfunc] = 0.0;
			--indx;
		}
		for(int i = nx[dm]-3; i >= 0; --i){
			for(int fn = 0; fn < nfunc; ++fn){
				tmp = data[indx];
				data[indx] = data[indx+2*nfunc] - 4.0*(i+1)*sp1[fn]/length;
				sp1[fn] = tmp;
				--indx;
			}
		}
	}

	delete [] sp1;
}

void Evaluator::integrate(double zero_val){
	double *tmp = new double[nfunc];
	double *next = new double[nfunc];
	double *prev = new double[nfunc];

	for(int fn = 0; fn < nfunc; ++fn)
		prev[fn] = zero_val;

	int indx;
	double length, tmp2;
	for(int dm = 0; dm < ndomains; ++dm){
		indx = index_st[dm];
		length = bound[dm+1] - bound[dm];
		
		//Doing forward integration
		for(int fn = 0; fn < nfunc; ++fn)
			tmp[fn] = data[indx++];

		for(int i = 1; i < nx[dm]-1; ++i){
			for(int fn = 0; fn < nfunc; ++fn){
				tmp2 = data[indx];
				data[indx] = -length*(tmp[fn] - data[indx+nfunc])/(4.0*i);
				tmp[fn] = tmp2;
				++indx;
			}
		}

		for(int fn = 0; fn < nfunc; ++fn){
			data[indx++] = -length*tmp[fn]/(4.0*(nx[dm]-1));
			tmp[fn] = 0.0;
			next[fn] = 0.0;
		}

		//Computing first and last values
		indx = index_st[dm]+nfunc-1;
		for(int i = 1; i < nx[dm]; ++i){
			for(int fn = 0; fn < nfunc; ++fn){
				next[fn] -= (2.0*(i%2) - 1.0)*data[++indx];
				if(doubling_init && (dm == 0)){
					if(i%2 == 0){
						tmp[fn] -= ((i%4) - 1.0)*data[indx];
					}
				}
				else
					tmp[fn] += data[indx];
			}
		}

		for(int fn = 0; fn < nfunc; ++fn){
			next[fn] += 0.5*(data[index_st[dm]+fn] = 2.0*(prev[fn] - tmp[fn]));
			prev[fn] = next[fn];
		}
	}
	
	delete [] next;
	delete [] tmp;
}