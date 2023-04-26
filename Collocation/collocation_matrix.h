#include <container.h>

// #define EIGEN_SUPERLU_SUPPORT
// #define EIGEN_USE_BLAS
// #define EIGEN_USE_LAPACKE

// #define LAPACK_COMPLEX_CUSTOM
// #define lapack_complex_float std::complex<float>
// #define lapack_complex_double std::complex<double>

#include <eigen3/Eigen/Dense>

#ifndef COLLOCATION_MATRIX_H
#define COLLOCATION_MATRIX_H

class Collocation_Matrix{
private:
	int nfunc;

	int ndomains;
	std::vector<int> nx;
	std::vector<int> index_st;
	std::vector<double> bound;

	bool doubling_st;
	std::vector<bool> parity_st;//false - even, true - odd
	bool doubling_end;
	std::vector<bool> parity_end;

	int sz;

public:
	Eigen::MatrixXd matr;
	Eigen::VectorXd rhs;
	Eigen::VectorXd delta;

	Collocation_Matrix(const Container *cont):
		nfunc(cont->nfunc), ndomains(cont->ndomains), nx(cont->nx), index_st(cont->index_st), bound(cont->bound), sz(cont->sz),
		parity_st(cont->parity_st), parity_end(cont->parity_end)
	{
		matr.resize(sz, sz);
		rhs.resize(sz);
		delta.resize(sz);
		matr.setZero();
		rhs.setZero();
		delta.setZero();
	}

	int size() const{return sz;}

	//This functions zeroes the matrix
	void setZero_jacobian(void){
		matr.setZero();
	}

	void setZero_rhs(void){
		rhs.setZero();
	}

	//This function zeroes the matrix row
	void setZero_row(unsigned int dm, unsigned int ind, unsigned int fn){
		matr.row(index_st[dm] + nfunc*ind + fn).setZero();
	}

	//This function acceses right hand side array
	double &right_hand_side(unsigned int dm, unsigned int ind, unsigned int fn){
		return rhs(index_st[dm] + nfunc*ind + fn);
	}

	//This function adds a differention operator to the (row_dm, row_ind) row, multiplying it by coef
	//der1 and der2 indicate derivative dimensions, der1 = -1 and/or der2 = -1 are used to indicate no differention
	//This is the main function used to construct the Jacobian
	template<unsigned int ord>
	void row_add_Dn_factor(unsigned int row_dm, unsigned int row_ind, unsigned int row_fn,
		unsigned int dm, unsigned int ind, unsigned int fn, double coef)
	{
		int row_index = index_st[row_dm] + nfunc*row_ind + row_fn;

		if(ord == 0){
			matr(row_index, index_st[dm] + nfunc*ind + fn) += coef;
			return;
		}

		double factor;
		int indx = index_st[dm] + fn;

		if(dm == 0 && !parity_st.empty()){
			factor = coef * std::pow(1.0/(bound[dm+1] - bound[dm]), ord);
			for(int i = 0; i < nx[dm]; ++i){
				matr(row_index, indx) += factor * cheb_Dn_double<ord,false>(nx[0], parity_st[fn], ind, i);
				indx += nfunc;
			}
		}
		else if(dm == ndomains-1 && !parity_end.empty()){
			factor = coef * std::pow(1.0/(bound[dm+1] - bound[dm]), ord);
			for(int i = 0; i < nx[dm]; ++i){
				matr(row_index, indx) += factor * cheb_Dn_double<ord,true>(nx[dm], parity_end[fn], ind, i);
				indx += nfunc;
			}
		}
		else{
			factor = coef * std::pow(2.0/(bound[dm+1] - bound[dm]), ord);
			for(int i = 0; i < nx[dm]; ++i){
				matr(row_index, indx) += factor * cheb_Dn<ord>(nx[dm]-1, ind, i);
				indx += nfunc;
			}
		}
		return;
	}


	//Computes the next step given the Jacobian and the equation difference
	void compute_delta(void){
		delta = matr.partialPivLu().solve(rhs);
	}

	// Changes domains configuration
	void regrid(const Container *cont){
		nfunc = cont->nfunc; ndomains = cont->ndomains; nx = cont->nx;
		index_st = cont->index_st; bound = cont->bound; sz = cont->sz;
		parity_st = cont->parity_st; parity_end = cont->parity_end;
		matr.resize(sz, sz);
		rhs.resize(sz);
		delta.resize(sz);
		matr.setZero();
		rhs.setZero();
		delta.setZero();
	}

	friend class Container;
};

#endif