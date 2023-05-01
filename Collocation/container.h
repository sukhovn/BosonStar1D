#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <iomanip>
#include <limits>
#include <cctype>

#include <spectral.h>

#ifndef CONTAINER_H
#define CONTAINER_H

int count_func_files(const std::string &folder_name);
class Collocation_Matrix;


class Container{
private:
	int nfunc;

	int ndomains;
	std::vector<int> nx;
	std::vector<int> index_st;
	std::vector<double> bound;

	std::vector<bool> parity_st;//false - even, true - odd
	std::vector<bool> parity_end;

	int sz;
	double *data;

	void reallocate_partitioning(const std::vector<int> &nxi){
		if(nxi.size() < 1){
			throw("Zero domains provided");
			std::cout << "Zero domains provided" << std::endl;
		}

		ndomains = nxi.size();
		nx = nxi;
		index_st.resize(0);
		int sz_tmp = 0;

		for(int nxd : nx){
			index_st.push_back(sz_tmp);
			sz_tmp += nfunc*(nxd-1);
		}
		sz_tmp += nfunc;

		if(sz != sz_tmp){
			if(data != NULL) delete [] data;
			data = new double[sz_tmp];
			sz = sz_tmp;
			for(int i = 0; i < sz; ++i)
				data[i] = 0.0;
		}
	}

public:
	Container(int nfunci, const std::vector<int> &nxi, const std::vector<double> &boundi):
		nfunc(nfunci), bound(boundi), sz(0), data(NULL)
	{
		if(nxi.size() + 1 != bound.size()){
			std::cout << "Wrong domain configuration" << std::endl;
			throw("Wrong domain configuration");
		}

		reallocate_partitioning(nxi);
	}

	void set_doubling_begin(const std::vector<bool> &parity){
		if(!parity_end.empty() && ndomains == 1){
			std::cout << "Cannot set doubling on both sized of the same domain" << std::endl;
			throw("Cannot set doubling on both sized of the same domain");
		}
		parity_st = parity;
	}

	void set_doubling_end(const std::vector<bool> &parity){
		if(!parity_st.empty() && ndomains == 1){
			std::cout << "Cannot set doubling on both sized of the same domain" << std::endl;
			throw("Cannot set doubling on both sized of the same domain");
		}
		parity_end = parity;
	}

	Container(const std::string &folder_name): nfunc(count_func_files(folder_name)), sz(0), data(NULL){
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

	void read_domain_structure(const std::string &file_name, int fn, bool fill = true);

	Container(const Container &cont):
		nfunc(cont.nfunc), ndomains(cont.ndomains), nx(cont.nx), index_st(cont.index_st), bound(cont.bound), sz(cont.sz),
		parity_st(cont.parity_st), parity_end(cont.parity_end)
	{
		data = new double[sz];
		for(int i = 0; i < sz; ++i)
			data[i] = cont.data[i];
	}

	Container &operator=(const Container &cont){
		if(cont.nx != nx || nfunc != cont.nfunc){
			nfunc = cont.nfunc;
			reallocate_partitioning(cont.nx);
		}
		bound = cont.bound;
		for(int i = 0; i < sz; ++i)
			data[i] = cont.data[i];

		parity_st = cont.parity_st;
		parity_end = cont.parity_end;

		return *this;
	}

	~Container(){
		delete [] data;
	}

	double norm() const{
		double tmp = 0.0;
		for(int i = 0; i < sz; ++i){
			tmp += data[i]*data[i];
		}
		return std::sqrt(tmp);
	}

	int n_func() const{return nfunc;}
	int n_domains() const{return ndomains;}
	int n_x(int i) const{return nx[i];}
	int size() const{return sz;}

	const std::vector<int> &n_x() const{return nx;}
	const std::vector<double> &boundaries() const {return bound;}
	const double &operator[](int i) const {return data[i];}
	double &operator[](int i){return data[i];}
	double last_domain_length() const {return bound[ndomains] - bound[ndomains-1];}

	double coord(unsigned int dm, unsigned int i) const{
		if(!parity_st.empty() && dm == 0) return ch_x(2.0*bound[0] - bound[1], bound[1], 2*nx[0]-1, nx[0]+i);
		if(!parity_end.empty() && dm == ndomains-1) return ch_x(bound[dm], 2*bound[dm+1]-bound[dm], 2*nx[dm]-1, i);
		return ch_x(bound[dm], bound[dm+1], nx[dm]-1, i);
	}

	template<class Callable>
	void fill_function(const Callable &func){
		int indx = 0;
		double x;
		for(int dm = 0; dm < ndomains; ++dm){
			for(int i = (dm ? 1 : 0); i < nx[dm]; ++i){
				x = coord(dm, i);
				for(int fn = 0; fn < nfunc; fn++){
					data[indx++] = func(fn, x);
				}
			}
		}
	}

	template<class Callable>
	double compare_function(const Callable &func) const{
		int indx = 0;
		double x, tmp, max = 0.0;
		for(int dm = 0; dm < ndomains; ++dm){
			for(int i = (dm ? 1 : 0); i < nx[dm]; ++i){
				x = coord(dm, i);
				for(int fn = 0; fn < nfunc; fn++){
					tmp = std::fabs(data[indx++] - func(fn, x));
					// data[indx++] -= func(fn, x);
					if(tmp > max) max = tmp;
				}
			}
		}
		return max;
	}

	void subtract_delta(const Collocation_Matrix &colm, double coeff = 1.0);

	//Changes partitioning while zeroing solution
	void adjust_partitioning(const std::vector<int> &nxi, const std::vector<double> &boundi){
		bound = boundi;
		reallocate_partitioning(nxi);
	}

	void save(const std::string &folder_name) const;
	void load(const std::string &folder_name);

	friend class Evaluator;
	friend class Collocation_Matrix;

	friend class Point;
	template<class>
	friend class Coordinate_Compactified;
	template<class>
	friend class Point_Compactified;
	template<class>
	friend class Collocation_Solution;

	template<class, class>
	friend class Collocation_Solver;
};


class Evaluator{
private:
	int nfunc;

	int ndomains;
	std::vector<int> nx;
	std::vector<int> index_st;
	std::vector<double> bound;

	bool doubling_init;

	int sz;
	double *data;

	void array_fft(void);

	Evaluator(int nfunci, const std::vector<int> &nxi, const std::vector<double> &boundi):
		nfunc(nfunci), ndomains(nxi.size()), nx(nxi), bound(boundi), sz(0), data(NULL), doubling_init(false)
	{
		if(nxi.size() + 1 != bound.size()){
			std::cout << "Wrong domain configuration" << std::endl;
			throw("Wrong domain configuration");
		}

		if(nxi.size() < 1){
			throw("Zero domains provided");
			std::cout << "Zero domains provided" << std::endl;
		}

		for(int dm = 0; dm < ndomains; ++dm){
			index_st.push_back(sz);
			sz += nfunc*nx[dm];
		}

		data = new double[sz];
	}

public:
	Evaluator(const Container &cont);

	template<class Callable>
	Evaluator(const Container &cont, int nfunc, const Callable &func): Evaluator(nfunc, cont.nx, cont.bound){
		int indx = -1;
		double x;
		for(int dm = 0; dm < ndomains; ++dm){
			if(dm){
				for(int fn = 0; fn < nfunc; ++fn){
					++indx;
					data[indx] = data[indx-nfunc];
				}
			}
			for(int i = (dm ? 1 : 0); i < nx[dm]; ++i){
				x = ch_x(bound[dm], bound[dm+1], nx[dm]-1, i);
				for(int fn = 0; fn < nfunc; ++fn){
					data[++indx] = func(fn, x);
				}
			}
		}

		array_fft();
	}

	Evaluator(Evaluator &&evali): nfunc(evali.nfunc), nx(evali.nx), bound(evali.bound),
		ndomains(evali.ndomains), index_st(evali.index_st), sz(evali.sz), data(evali.data) {};

	~Evaluator(){
		delete [] data;
	}

	double operator()(double coord, int fn = 0) const;
	double derivative(double coord, int fn = 0) const;

	void test(const Container &cont) const{
		int indx = -1;
		double max = 0.0, tmp, coord;
		for(int dm = 0; dm < cont.ndomains; ++dm){
			for(int i = (dm ? 1 : 0); i < cont.nx[dm]; ++i){
				coord = cont.coord(dm, i);
				for(int fn = 0; fn < cont.nfunc; ++fn){
					tmp = std::fabs(cont.data[++indx] - (*this)(coord, fn));
					if(tmp > max)
						max = tmp;
				}
			}
		}
		std::cout << "Maximal discrepancy is " << max << std::endl;
	}

	void differentiate(void);
	void integrate(double zero_val = 0.0);
};

#endif