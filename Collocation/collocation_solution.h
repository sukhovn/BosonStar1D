#include <container.h>
#include <collocation_matrix.h>

#ifndef COLLOCATION_SOLUTION_H
#define COLLOCATION_SOLUTION_H

template<class Compactification>
class Collocation_Solution{
private:
	Point_Compactified<Compactification> pnt;
	Container *cont;

public:
	const Compactification &compact;

	Collocation_Solution(const Compactification &compact_, Container *conti):
		compact(compact_), cont(conti), pnt(compact_, conti) {};

	void change_container(Container *cont_change){
		pnt.change_container(cont_change);
		cont = cont_change;
	}

	template<class Callable>
	void fill_seed(Callable func, double value_at_infinity = 0.0){
		int indx = 0;
		double x_tmp;

		for(int dm = 0; dm < cont->ndomains; ++dm){
			for(int i = (dm ? 1 : 0); i <  cont->nx[dm]; ++i){
				x_tmp = pnt.get_coordinate(dm, i);

				if(std::isinf(x_tmp)){
					for(int fn = 0; fn < cont->nfunc; ++fn) cont->data[indx++] = value_at_infinity;
					continue;
				}

				for(int fn = 0; fn <  cont->nfunc; fn++){
					 cont->data[indx++] = func(fn, x_tmp);
				}
			}
		}
	}

	void save_csv(const std::string &file_name, double xmin, double xmax, int npoints){
		if(xmin < cont->bound[0]){
			std::cout << "Interpolation reaching beyond domain range" << std::endl;
			throw("xmin or xmax out of bounds");
		}

		std::ofstream out(file_name);

		out << 'x';
		for(int fn = 0; fn < cont->nfunc; ++fn)
			out << ",f" << fn;
		out << std::endl;

		double x, x_c;
		double dx = (xmax - xmin)/(1.0*(npoints-1));

		Evaluator evol(*cont);
		Coordinate_Compactified<Compactification> coord(compact, cont);

		out << std::setprecision(16);
		for(int i = 0; i < npoints; ++i){
			x_c = coord.comp_coord(x);

			out << x;
			for(int fn = 0; fn < cont->nfunc; ++fn)
				out << ',' << evol(x_c, fn);
			out << std::endl;

			x += dx;
		}
	}

	template<class, class>
	friend class Collocation_Solver;
};

#endif