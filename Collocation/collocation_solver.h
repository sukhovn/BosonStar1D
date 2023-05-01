#include <point_class.h>
#include <collocation_solution.h>
#include <tuple>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#ifndef COLLOCATION_SOLVER_H
#define COLLOCATION_SOLVER_H

//This functions find maximum and minimum values out of a and b.
#ifndef FMAX
#define FMAX(a, b) (a > b ? a : b)
#endif

//Difference with respect to which the equations are numerically differentiated
#define H_DERIV 1.0e-7

//Collocation iteration parameters
#define ALPHA 1.0e-4 //Lambda check parameter
#define TOLX 1.0e-10 //Minimum step parameter

#define STPMX 100.0 //Maximal Newton step

#define COLLOC_THRSHLD 1.0e-13
#define MAX_ITER 100 //Limit for maximal number of iterations

template<class Compactification, class Equations_Class>
class Collocation_Solver: public Collocation_Solution<Compactification>{
private:
	Collocation_Matrix colm;
	double eqn_error;
	int iterations;
	const Equations_Class *eq_system;

	using Collocation_Solution<Compactification>::cont;
	using Collocation_Solution<Compactification>::pnt;

	typedef double (Equations_Class::*EqnPointer)(int, double, double *, double *, double *) const;
	typedef double (Equations_Class::*BcPointer)(int, double *, double *) const;

	double deq_d_phi(const EqnPointer eqn_func, int fn_i, int fn_j){
			double tmp;
			pnt.phi[fn_j] += H_DERIV;
			tmp = (eq_system->*eqn_func)(fn_i, pnt.x, pnt.phi, pnt.dphi, pnt.ddphi);
			pnt.phi[fn_j] -= 2.0*H_DERIV;
			tmp -= (eq_system->*eqn_func)(fn_i, pnt.x, pnt.phi, pnt.dphi, pnt.ddphi);
			pnt.phi[fn_j] += H_DERIV;
			return tmp/(2.0*H_DERIV);
	}
	
	double deq_d_dphi(const EqnPointer eqn_func, int fn_i, int fn_j){
			double tmp;
			pnt.dphi[fn_j] += H_DERIV;
			tmp = (eq_system->*eqn_func)(fn_i, pnt.x, pnt.phi, pnt.dphi, pnt.ddphi);
			pnt.dphi[fn_j] -= 2.0*H_DERIV;
			tmp -= (eq_system->*eqn_func)(fn_i, pnt.x, pnt.phi, pnt.dphi, pnt.ddphi);
			pnt.dphi[fn_j] += H_DERIV;
			return tmp/(2.0*H_DERIV);
	}

	double deq_d_ddphi(const EqnPointer eqn_func, int fn_i, int fn_j){
			double tmp;
			pnt.ddphi[fn_j] += H_DERIV;
			tmp = (eq_system->*eqn_func)(fn_i, pnt.x, pnt.phi, pnt.dphi, pnt.ddphi);
			pnt.ddphi[fn_j] -= 2.0*H_DERIV;
			tmp -= (eq_system->*eqn_func)(fn_i, pnt.x, pnt.phi, pnt.dphi, pnt.ddphi);
			pnt.ddphi[fn_j] += H_DERIV;
			return tmp/(2.0*H_DERIV);
	}

	double dbc_d_phi(const BcPointer bc_func, int fn_i, int fn_j){
			double tmp;
			pnt.phi[fn_j] += H_DERIV;
			tmp = (eq_system->*bc_func)(fn_i, pnt.phi, pnt.dphi);
			pnt.phi[fn_j] -= 2.0*H_DERIV;
			tmp -= (eq_system->*bc_func)(fn_i, pnt.phi, pnt.dphi);
			pnt.phi[fn_j] += H_DERIV;
			return tmp/(2.0*H_DERIV);
	}

	double dbc_d_dphi(const BcPointer bc_func, int fn_i, int fn_j){
			double tmp;
			pnt.dphi[fn_j] += H_DERIV;
			tmp = (eq_system->*bc_func)(fn_i, pnt.phi, pnt.dphi);
			pnt.dphi[fn_j] -= 2.0*H_DERIV;
			tmp -= (eq_system->*bc_func)(fn_i, pnt.phi, pnt.dphi);
			pnt.dphi[fn_j] += H_DERIV;
			return tmp/(2.0*H_DERIV);
	}

	template<bool FILL_RHS, bool FILL_JACOBIAN>
	void fill_left_boundary_condition(){
		double tmp;
		pnt.get_point(0, 0);
		for(int fn = 0; fn < cont->nfunc; ++fn){
			if(FILL_RHS) colm.right_hand_side(0, 0, fn) = eq_system->boundary_condition_left(fn, pnt.phi, pnt.dphi);

			if(FILL_JACOBIAN){
				for(int fnd = 0; fnd < cont->nfunc; fnd++){
					tmp = dbc_d_phi(&Equations_Class::boundary_condition_left, fn, fnd);
					if(tmp != 0) colm.row_add_Dn_factor<0>(0, 0, fn, 0, 0, fnd, tmp);
					tmp = dbc_d_dphi(&Equations_Class::boundary_condition_left, fn, fnd);
					if(tmp != 0) colm.row_add_Dn_factor<1>(0, 0, fn, 0, 0, fnd, tmp);
				}
			}
		}
	}

	template<bool FILL_RHS, bool FILL_JACOBIAN>
	void fill_right_boundary_condition(){
		double tmp;
		int dm = cont->ndomains-1;
		int ind = cont->nx[dm]-1;
		pnt.get_point(dm, ind);
		for(int fn = 0; fn < cont->nfunc; ++fn){
			if(FILL_RHS) colm.right_hand_side(dm, ind, fn) = eq_system->boundary_condition_right(fn, pnt.phi, pnt.dphi);
			if(FILL_JACOBIAN){
				for(int fnd = 0; fnd < cont->nfunc; fnd++){
					tmp = dbc_d_phi(&Equations_Class::boundary_condition_right, fn, fnd);
					if(tmp != 0) colm.row_add_Dn_factor<0>(dm, ind, fn, dm, ind, fnd, tmp);
					tmp = dbc_d_dphi(&Equations_Class::boundary_condition_right, fn, fnd);
					if(tmp != 0) colm.row_add_Dn_factor<1>(dm, ind, fn, dm, ind, fnd, tmp);
				}
			}
		}
	}

	template<bool FILL_RHS, bool FILL_JACOBIAN>
	void fill_interdomain_boundary(int dm){
		pnt.get_point(dm, cont->nx[dm]-1);
		for(int fn = 0; fn < cont->nfunc; ++fn){
			if(FILL_RHS) colm.right_hand_side(dm+1, 0, fn) = -pnt.dphi[fn];
			if(FILL_JACOBIAN) colm.row_add_Dn_factor<1>(dm+1, 0, fn, dm, cont->nx[dm]-1, fn, -1.0/pnt.dx);
		}
		pnt.get_point(dm+1, 0);
		for(int fn = 0; fn < cont->nfunc; ++fn){
			if(FILL_RHS) colm.right_hand_side(dm+1, 0, fn) += pnt.dphi[fn];
			if(FILL_JACOBIAN) colm.row_add_Dn_factor<1>(dm+1, 0, fn, dm+1, 0, fn, 1.0/pnt.dx);
		}
	}

	template<bool FILL_RHS, bool FILL_JACOBIAN>
	void fill_bulk(int dm){
		int ist = !cont->parity_st.empty() && dm == 0 ? 0 : 1;
		int iend = !cont->parity_end.empty() && dm == cont->ndomains-1 ? cont->nx[dm] : cont->nx[dm]-1;
		double tmp, deq, ddeq;

		for(int i = ist; i < iend; ++i){
			pnt.get_point(dm, i);

			for(int fn = 0; fn < cont->nfunc; ++fn){
				if(FILL_RHS) colm.right_hand_side(dm, i, fn) = eq_system->equations(fn, pnt.x, pnt.phi, pnt.dphi, pnt.ddphi);

				if(FILL_JACOBIAN){
					//Adding field derivatives to the Jacobian
					for(int fnd = 0; fnd < cont->nfunc; ++fnd){
						tmp = deq_d_phi(&Equations_Class::equations, fn, fnd);
						colm.row_add_Dn_factor<0>(dm, i, fn, dm, i, fnd, tmp);
					}

					//Adding field derivatives' derivatives to the Jacobian
					for(int fnd = 0; fnd < cont->nfunc; ++fnd){
						deq = deq_d_dphi(&Equations_Class::equations, fn, fnd);
						ddeq = deq_d_ddphi(&Equations_Class::equations, fn, fnd);

						if(deq != 0.0 || (ddeq != 0.0 && pnt.ddx != 0)){
							tmp = deq/pnt.dx - ddeq*pnt.ddx/std::pow(pnt.dx, 3.0);
							colm.row_add_Dn_factor<1>(dm, i, fn, dm, i, fnd, tmp);
						}
						if(ddeq != 0) colm.row_add_Dn_factor<2>(dm, i, fn, dm, i, fnd, ddeq/(pnt.dx*pnt.dx));
					}
				}
			}
		}
	}

public:
	Collocation_Solver(const Compactification &compact, const Equations_Class *eq_sys, Container *conti):
		Collocation_Solution<Compactification> (compact, conti), eq_system(eq_sys), colm(conti) {};

	void change_container(Container *cont_change,	const Equations_Class *eq_system_change){
		Collocation_Solution<Compactification>::change_container(cont_change);
		colm.regrid(cont_change);
		eq_system = eq_system_change;
	}

	template<bool FILL_RHS, bool FILL_JACOBIAN>
	void collocation_fill(){
		int ist, iend;
		if(FILL_RHS) colm.setZero_rhs();
		if(FILL_JACOBIAN) colm.setZero_jacobian();

		if(cont->parity_st.empty()) fill_left_boundary_condition<FILL_RHS, FILL_JACOBIAN>();
		for(int dm = 0; dm < cont->ndomains; ++dm){
			fill_bulk<FILL_RHS, FILL_JACOBIAN>(dm);
			if(dm != cont->ndomains-1) fill_interdomain_boundary<FILL_RHS, FILL_JACOBIAN>(dm); 
		}
		if(cont->parity_end.empty()) fill_right_boundary_condition<FILL_RHS, FILL_JACOBIAN>();
	}

	std::tuple<std::vector<double>, std::vector<int>, std::vector<int>> 
		find_function_maximum(EqnPointer func, int nfunc = 1, bool interior_check_only = true){
		int ist, iend;
		double tmp;
		std::vector<int> dm_max(nfunc, 0);
		std::vector<int> i_max(nfunc, 0);
		std::vector<double> max_diff(nfunc, 0.0);
		for(int dm = 0; dm < cont->ndomains; ++dm){
            ist = dm == 0 && (!interior_check_only || !cont->parity_st.empty()) ? 0 : 1;
			iend = interior_check_only && dm == cont->ndomains-1 && cont->parity_end.empty() ? 
							cont->nx[dm]-1 : cont->nx[dm];
			for(int i = ist; i < iend; ++i){
				pnt.get_point(dm, i);
				for(int fn = 0; fn < cont->nfunc; fn++){
					tmp = std::fabs((eq_system->*func)(fn, pnt.x, pnt.phi, pnt.dphi, pnt.ddphi));
					if(tmp > max_diff[fn]){
						max_diff[fn] = tmp;
						dm_max[fn] = dm;
						i_max[fn] = i;
					}
				}
			}
		}
		return std::make_tuple(max_diff, dm_max, i_max);
	}

	void save_function_to_file(const std::string &folder_name, EqnPointer func, int nfunc = 1, bool interior_check_only = true){
		struct stat sttus = {0};
		if (stat(folder_name.c_str(), &sttus) == -1) {
			mkdir(folder_name.c_str(), 0700);
		}

        std::string file_name;
		std::ofstream **out = new std::ofstream*[nfunc];
		for(int fn = 0; fn < nfunc; ++fn){
			file_name = folder_name + "/eval" + std::to_string(fn) + ".dat";
			out[fn] = new std::ofstream(file_name);
		}

		for(int fn = 0; fn < nfunc; ++fn){
			for(int dm = 0; dm < cont->ndomains; dm++)
				*out[fn] << cont->bound[dm] << "," << cont->nx[dm] << ",";
			*out[fn] << cont->bound[cont->ndomains] << ",";

			if(cont->parity_st.empty())
				*out[fn] << "n";
			else
				*out[fn] << "a";
			*out[fn] << ",";

			if(cont->parity_end.empty())
				*out[fn] << "n";
			else
				*out[fn] << "a";

			*out[fn] << std::endl;
		}

		int ist, iend;
		double tmp;

		for(int fn = 0; fn < nfunc; ++fn){
			pnt.get_point(0, interior_check_only && cont->parity_st.empty() ? 1 : 0);
			*out[fn] << (eq_system->*func)(fn, pnt.x, pnt.phi, pnt.dphi, pnt.ddphi);
		}

		for(int dm = 0; dm < cont->ndomains; ++dm){
            ist = dm == 0 && interior_check_only && cont->parity_st.empty() ? 2 : 1;
			iend = interior_check_only && dm == cont->ndomains-1 && cont->parity_end.empty() ? cont->nx[dm]-1 : cont->nx[dm];

			for(int i = ist; i < iend; ++i){
				pnt.get_point(dm, i);
				for(int fn = 0; fn < nfunc; fn++)
					*out[fn] << "," << (eq_system->*func)(fn, pnt.x, pnt.phi, pnt.dphi, pnt.ddphi);
			}
		}

		for(int fn = 0; fn < nfunc; ++fn)
			delete out[fn];
		delete [] out;
	}

private:
	int lnsearch(double stepmax){
		double temp = colm.delta.norm(); //Limiting maximal step 
		if(temp > stepmax)
			colm.delta *= stepmax/temp;

		double test = 0.0; //Adjusting minimal lambda
		for(int i = 0; i < cont->sz; ++i){
			temp = std::fabs(colm.delta(i))/FMAX(std::fabs(cont->data[i]), 1.0);
			if(temp > test) test = temp;
		}

		double alamin = TOLX/test;
		double alam1, alam2, alamc = 1.0;
		double f1, f2, fc;
		double f = eqn_error;
		double a, b, disc;

		double diff_to_origin = 1.0;
		cont->subtract_delta(colm);

		collocation_fill<true, false>();
		fc = 0.5*colm.rhs.squaredNorm();

		while(1){
			if(fc <= (1.0 - 2.0*ALPHA) * f){ //Checking if the step with current lambda is adequatly small
				eqn_error = fc;
				return 0;
			}
	//		std::cout << "Continuing Lambda search\n";

			alam2 = alam1;
			f2 = f1;
			alam1 = alamc;
			f1 = fc;

			if(alamc == 1.0){ //Initial lambda step
				alamc = f/(f + fc);
			}
			else{ //Secondary lambda step
				a = ((f1 + (2.0 * alam1 - 1.0) * f)/(alam1 * alam1) - (f2 + (2.0 * alam2 - 1.0) * f)/(alam2 * alam2))/(alam1-alam2);
				b = (alam1 * (f2 + (2.0 * alam2 - 1.0) * f)/(alam2 * alam2) - alam2 * (f1 + (2.0 * alam1 - 1.0) * f)/(alam1 * alam1))/(alam1-alam2);
				if(a == 0.0)
					alamc = f/b;
				else{
					disc = b * b + 6.0 * a * f;
					if(disc < 0.0)
						alamc = 0.5 * alam1;
					else if(b <= 0.0)
						alamc = (-b + std::sqrt(disc))/(3.0 * a);
					else
						alamc = 2.0 * f/(b + std::sqrt(disc));
				}
				if(alamc > 0.5 * alam1)
					alamc = 0.5 * alam1;
			}
			alamc = FMAX(alamc, 0.1 * alam1);

			if(alamc < alamin){
				std::cout << "Reached minimal lambda " << alamc << '\n';
				cont->subtract_delta(colm, -diff_to_origin);
				collocation_fill<true, false>();
				return 1;
			}

			cont->subtract_delta(colm, alamc-diff_to_origin);
			diff_to_origin = alamc;

			collocation_fill<true, false>();
			fc = 0.5*colm.rhs.squaredNorm();
		}
	}

	int iteration(double stepmax){
		int code;//._.

		collocation_fill<false, true>();

		colm.compute_delta();

		// cont->subtract_delta(colm);
		// collocation_fill<true, false>();
		// eqn_error = 0.5*colm.rhs.squaredNorm();

		// std::cout << "Searching for the next step\n";

		code = lnsearch(stepmax);

		iterations++;

		return 0;//code;
	}

public:
	int solve(void){
		int code;
		double err, stepmax;
		std::ios oldState(nullptr);
		oldState.copyfmt(std::cout);

		collocation_fill<true, false>();
		eqn_error = 0.5*colm.rhs.squaredNorm();

		stepmax = STPMX*FMAX(std::sqrt(cont->norm()), (double) cont->size());
		iterations = 0;

		for(int i = 0; i < MAX_ITER; i++){
			code = iteration(stepmax);
			err = std::sqrt(2.0*eqn_error)/(1.0*colm.size());
			std::cout << "Step is done, the error is: " << std::scientific << err << std::endl;
			if(err < COLLOC_THRSHLD){
				std::cout << "Calculation is done, the error is: " << std::scientific << err << std::endl;
				std::cout << "Number of iterations is: " << iterations << "\n" << std::endl;
				std::cout.copyfmt(oldState);
				return 0;
			}
			if(code != 0){
				std::cout << "Reached the minimal newton step, the error is: " << std::scientific << err << std::endl;
				std::cout << "Number of iterations is: " << iterations << "\n" << std::endl;
				std::cout.copyfmt(oldState);
				return 1;
			}
		}

		std::cout << "Solution did not converge, no data is recorded\n" << std::endl;
		std::cout.copyfmt(oldState);
		return 1;
	}
};
#endif