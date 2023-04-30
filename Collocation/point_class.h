#include <container.h>
#include <limits>

#ifndef POINT_CLASS_H
#define POINT_CLASS_H

class Point{
public:
	const Container *cont;

	unsigned int dm;
	unsigned int ind;
	double x;

	double *phi;

	double *dphi;
	double *ddphi;

	Point(const Container *conti): cont(conti){
		phi = new double[cont->nfunc];
		dphi = new double[cont->nfunc];
		ddphi = new double[cont->nfunc];
	}
	~Point(){
		delete [] phi; delete [] dphi; delete [] ddphi;
	}

	void get_point_fields_only(unsigned int dmi, unsigned int indi){
		if(indi >= cont->nx[dmi]){
			std::cout << "Point " << indi << " is out of range " << cont->nx[dmi]-1 << " in domain " << dmi << "\n";
			throw("Point i is out of range");
		}

		dm = dmi;
		ind = indi;
		x = cont->coord(dm, ind);

		for(int fn = 0; fn < cont->nfunc; ++fn){
			phi[fn] = cont->data[cont->index_st[dm] + cont->nfunc*ind + fn];
		}
		for(int fn = 0; fn < cont->nfunc; ++fn) dphi[fn] = 0.0;
		for(int fn = 0; fn < cont->nfunc; ++fn) ddphi[fn] = 0.0;
	}

	void get_derivatives(){
		int indx = 0;
		double *pnt = cont->data + cont->index_st[dm];
		double rescale;

		// Case when the domain is doubled at the beginning
		if(dm == 0 && !cont->parity_st.empty()){
			for(int i = 0; i < cont->nx[0]; ++i){
				for(int fn = 0; fn < cont->nfunc; ++fn){
					dphi[fn] += cheb_Dn_double<1,false>(cont->nx[0], cont->parity_st[fn], ind, i) * pnt[indx];
					ddphi[fn] += cheb_Dn_double<2,false>(cont->nx[0], cont->parity_st[fn], ind, i) * pnt[indx++];
				}
			}
			rescale = 1.0/(cont->bound[1] - cont->bound[0]);
		}
		// Case when the domain is doubled at the end
		else if(dm == cont->ndomains-1 && !cont->parity_end.empty()){
			for(int i = 0; i < cont->nx[dm]; ++i){
				for(int fn = 0; fn < cont->nfunc; ++fn){
					dphi[fn] += cheb_Dn_double<1,true>(cont->nx[dm], cont->parity_end[fn], ind, i) * pnt[indx];
					ddphi[fn] += cheb_Dn_double<2,true>(cont->nx[dm], cont->parity_end[fn], ind, i) * pnt[indx++];
				}
			}
			rescale = 1.0/(cont->bound[dm+1] - cont->bound[dm]);
		}
		// Default no doubling case
		else{
			for(int i = 0; i < cont->nx[dm]; ++i){
				for(int fn = 0; fn < cont->nfunc; ++fn){
					dphi[fn] += cheb_Dn<1>(cont->nx[dm]-1, ind, i) * pnt[indx];
					ddphi[fn] += cheb_Dn<2>(cont->nx[dm]-1, ind, i) * pnt[indx++];
				}
			}
			rescale = 2.0/(cont->bound[dm+1] - cont->bound[dm]);
		}
		for(int fn = 0; fn < cont->nfunc; ++fn){
			dphi[fn] *= rescale;
			ddphi[fn] *= rescale*rescale;
		}
	}

	void get_point(unsigned int dmi, unsigned int indi){
		get_point_fields_only(dmi, indi);
		get_derivatives();
	}
};


template<class Compactification>
class Coordinate_Compactified{
	double comp_end;
	double xend;
public:
	const Compactification &compact;

	Coordinate_Compactified(const Compactification &comp_, const Container *cont):
		compact(comp_), comp_end(cont->bound[cont->ndomains-1]), xend(cont->bound[cont->ndomains]) {};

	double coord(double x){
		if(x >= xend)
			return std::numeric_limits<double>::infinity();
		else if(x < comp_end)
			return x;
		else
			return comp_end + compact.xc(x - comp_end);
	}

	double comp_coord(double x){
		if(x < comp_end)
			return x;
		else
			return comp_end + compact.xc_inv(x - comp_end);
	}
};


template<class Compactification>
class Point_Compactified: public Point{
public:
	const Compactification &compact;

	double x_c;
	double dx;
	double ddx;

	Point_Compactified(const Compactification &compi, const Container *conti): Point(conti), compact(compi) {};

	double get_coordinate(unsigned int dmi, unsigned int indi){
		if(dmi == cont->ndomains-1){
			return compact.xc(cont->coord(dmi, indi) - cont->bound[dmi]) + cont->bound[dmi];
		}
		else
			return cont->coord(dmi, indi);
	}

	void get_point_fields_only(unsigned int dmi, unsigned int indi){
		Point::get_point_fields_only(dmi, indi);

		if(dmi == cont->ndomains-1){
			x_c = x;
			x = compact.xc(x_c - cont->bound[dmi]) + cont->bound[dmi];
		}
		else{
			x_c = x;
		}
	}

	bool get_point(unsigned int dmi, unsigned int indi){
		Point::get_point_fields_only(dmi, indi);

		if(dmi == cont->ndomains-1){
			x_c = x;
			x = compact.xc(x_c - cont->bound[dmi]) + cont->bound[dmi];
			if(std::isinf(x)) return true;
			dx = compact.dxc(x_c - cont->bound[dmi]);
			ddx = compact.ddxc(x_c - cont->bound[dmi]);
		}
		else{
			x_c = x;
			dx = 1.0;
			ddx = 0.0;
		}

		Point::get_derivatives();

		if(dmi == cont->ndomains-1){
			for(int fn = 0; fn < cont->nfunc; ++fn){
				dphi[fn] /= dx;
				ddphi[fn] = (ddphi[fn] - ddx*dphi[fn])/(dx*dx);
			}
		}

		return false;
	}
};

#endif