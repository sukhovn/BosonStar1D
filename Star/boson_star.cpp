#include <boson_star.h>
#include <par_solitonic.h>
#include <interpolation.h>
#include <limits>

void Boson_Star::read_parameters(const std::string &folder_name){
	pfile.read_file(folder_name + "/parameters.dat");

	int shift;
	poly_coef.resize(4);

	if(pfile.name() == "Solitonic"){
		shift = 2;
		double mass = pfile[0];
		double sigma = pfile[1];

		if(sigma > 0.0){
			poly_coef[0] = mass*mass;
			poly_coef[1] = -4.0*mass*mass/(sigma*sigma);
			poly_coef[2] = 4.0*mass*mass/std::pow(sigma, 4.0);
		}
		else{
			poly_coef[0] = mass*mass;
		}
	}
	else{
		shift = 4;
		for(int i = 0; i < 5; i++)
			poly_coef[i] = pfile[i];
	}

	kappa = pfile[shift+0];
	params["omega"] = omega = pfile[shift+1];
	
	params["phi0"] = pfile[shift+2];
	params["mass"] = pfile[shift+3];
	params["charge"] = pfile[shift+4];
	params["radius"] = pfile[shift+5];
	params["compactness"] = pfile[shift+6];
	params["bound_null_geodesics"] = pfile[shift+7];
}

void Boson_Star::record_parameters(const std::string &folder_name){
	int shift;

	if(pfile.name() == "Solitonic"){
		shift = 2;
		if(std::fabs(poly_coef[1]) < 1.0e-14 && std::fabs(poly_coef[2]) < 1.0e-14){
			pfile[0] = std::sqrt(poly_coef[0]);
			pfile[1] = 0.0;
		}
		else{
			bool coef_3_cond = (std::fabs(poly_coef[2] - 0.25*poly_coef[1]*poly_coef[1]/poly_coef[0])/poly_coef[2] < 1.0e-14); 
			if(coef_3_cond){
				pfile[0] = std::sqrt(poly_coef[0]);
				pfile[1] = std::sqrt(-4.0*poly_coef[0]/poly_coef[1]);
			}
			else{
				pfile[0] = nan("");
				pfile[1] = nan("");
			}
		}
	}
	else{
		shift = 4;
		for(int i = 0; i < 4; i++)
			pfile[i] = poly_coef[i];
	}

	pfile[shift+0] = kappa;
	pfile[shift+1] = params["omega"];

	pfile[shift+2] = params["phi0"];
	pfile[shift+3] = params["mass"];
	pfile[shift+4] = params["charge"];
	pfile[shift+5] = params["radius"];
	pfile[shift+6] = params["compactness"];
	pfile[shift+7] = params["bound_null_geodesics"];

	pfile.print_file(folder_name + "/parameters.dat");
}

Boson_Star::Boson_Star(const std::vector<int> &nxi, const std::vector<double> &boundi):
	Boson_Star_Equations(), Container(3, nxi, boundi),
	params({{"omega", 1.0}, {"phi0", 0.0}, {"mass", 0.0}, {"charge", 0.0}, {"radius", 0.0}, {"compactness", 0.0}, {"bound_null_geodesics", 0.0}}),
	pfile((const char*) __Star_par_solitonic_dat, __Star_par_solitonic_dat_len),
	coordinate_radius(0.0), ps_length(0.0) {};

Boson_Star::Boson_Star(const std::string &folder_name):
	Boson_Star_Equations(), Container(folder_name),
	params({{"omega", 1.0}, {"phi0", 0.0}, {"mass", 0.0}, {"charge", 0.0}, {"radius", 0.0}, {"compactness", 0.0}, {"bound_null_geodesics", 0.0}}),
	pfile((const char*) __Star_par_solitonic_dat, __Star_par_solitonic_dat_len),
	coordinate_radius(0.0), ps_length(0.0)
{
	load(folder_name);
};

void Boson_Star_Solver::load_spline(const std::string &folder_name){
	star->read_parameters(folder_name);

	std::fstream file(folder_name + "/fields.csv");
	if(!file){
		std::cout << "Fields file not found" << std::endl;
		throw("file missing");
	}

	std::string line;
	std::getline(file, line);

	//Counting number of field in the file
	int nfields = -1;
	std::string token;
	std::istringstream iss(line);
	while(iss.peek() != EOF){
		std::getline(iss, token, ',');
		if(nfields == -1){
			if(token.compare("x")){
				std::cout << "Fields file head is corrupted" << std::endl;
				throw("file corrupted");
			}
		}
		else if(token.compare("y" + std::to_string(nfields))){
			std::cout << "Fields file head is corrupted" << std::endl;
			throw("file corrupted");
		}
		++nfields;
	}

	if(nfields < star->n_func()){
		std::cout << "Number of fields in the fields file doesn't match number of fields in BS solver" << std::endl;
		throw("file mismatch");
	}

	//Counting number of points in the file
	int pos_data = file.tellg();
	int npoints = 0;
	while(file.peek() != EOF){
		std::getline(file, line);
		if(!line.compare(""))
			break;
		++npoints;
	}
	file.clear();
	file.seekg(pos_data, file.beg);

	//Reading file
	double *xcoord = new double[npoints];
	double **fields = new double*[nfields];
	for(int i = 0; i < nfields; ++i)
		fields[i] = new double[npoints];
	double *fields_spline = new double[npoints];

	std::string::size_type sztp;
	for(int i = 0; i < npoints; ++i){
		std::getline(file, line);
		std::istringstream iss(line);
		for(int j = -1; j < nfields; ++j){
			if(iss.peek() == EOF){
				std::cout << "Fields file corrupted" << std::endl;
				std::cout << "Missing the entry [" << i << ',' << j << "]" << std::endl;
				throw("file corrupted");
			}
			std::getline(iss, token, ',');
			if(j == -1){
				xcoord[i] = std::stod(token, &sztp);
			}
			else{
				fields[j][i] = std::stod(token, &sztp);
			}
		}
	}

	double xi, tmp;
	double inter_bound = xcoord[npoints-1];
	double m_est = xcoord[npoints-1]*(1.0 - std::exp(-fields[0][npoints-1]));
	Point_Compactified<Linear_Fractional_Compactification> point(*this, star);

	for(int fn = 0; fn < star->n_func(); ++fn){
		spline(xcoord, fields[fn], npoints, 0.0, 0.0, fields_spline);
		int indx = fn;
		for(int dm = 0; dm < star->n_domains(); ++dm){
			for(int i = (dm ? 1 : 0); i < star->n_x(dm); ++i){
				xi = point.get_coordinate(dm, i);
				if(xi < inter_bound){
					if(splint(xcoord, fields[fn], fields_spline, npoints, xi, &tmp)){
						delete [] fields_spline;
						for(int i = 0; i < nfields; ++i) delete [] fields[i];
						delete [] fields;
						delete [] xcoord;
						throw("Spline interpolation failed, initialization stopped\n");
					}
				}
				else{
					if(std::isinf(xi)){
						tmp = 0.0;
					}
					else{
						if(fn == 0)
							tmp = -std::log(1.0 - m_est/xi);
						else if(fn == 1)
							tmp = std::log(1.0 - m_est/xi);
						else
							tmp = 0.0;
					}
				}
				(*star)[indx] = tmp;
				indx += star->n_func();
			}
		}
	}

	delete [] fields_spline;
	for(int i = 0; i < nfields; ++i) delete [] fields[i];
	delete [] fields;
	delete [] xcoord;
}

void Boson_Star_Solver::find_star_parameters(void){
	star->params["omega"] = star->omega;
	Evaluator profile(*star);

	//Finding the central value of the scalar field
	double domain_start = star->boundaries()[0];
	star->params["phi0"] = profile(domain_start, 2);

	//Evaluating mass and charge
	double comp_bound = star->boundaries()[star->n_domains()-1];
	double domain_end = star->boundaries()[star->n_domains()];
	Evaluator mass_func(*star,2,[&profile, this, &comp_bound, &domain_end](int fn, double x){
		if(x == domain_end)
			return 0.0;

		//Evaluating fields
		double u = profile(x, 0);
		double v = profile(x, 1);
		double phi = profile(x, 2);
		double dphi = profile.derivative(x, 2);

		//Evaluating coordinate and rescaling factor;
		double dx = 1.0;
		if(x >= comp_bound){
			dx = dxc(x - comp_bound);
			x = comp_bound + xc(x - comp_bound);
		}

		if(fn == 0){
			double dens = std::exp(-v)*std::pow(star->omega*phi, 2.0) + std::exp(-u)*std::pow(dphi/dx, 2.0) + star->vpot(phi*phi);
			return 0.5*star->kappa*dens*x*x*dx;
		}
		if(fn == 1){
			double dens = std::exp(0.5*(u-v))*star->omega*phi*phi;
			return std::sqrt(2.0*star->kappa)*dens*x*x*dx;
		}
		std::cout << "Wrong density number" << std::endl;
		throw("integral calculation error");
	});
	mass_func.integrate(0.0);

	star->params["mass"] = mass_func(domain_end, 0);
	star->params["charge"] = mass_func(domain_end, 1);

	//Evaluating radius
	Coordinate_Compactified<Linear_Fractional_Compactification> coord(*this, star);
	double xst_rad = domain_start, xend_rad = domain_end;
	while(xend_rad - xst_rad > 1.0e-6){
		double tmp = mass_func(0.5*(xst_rad+xend_rad), 0);
		if(tmp > RAD_FRACTION*star->params["mass"])
			xend_rad = 0.5*(xst_rad + xend_rad);
		else
			xst_rad = 0.5*(xst_rad + xend_rad);
	}
	star->coordinate_radius = 0.5*(xst_rad + xend_rad);
	star->params["radius"] = coord.coord(star->coordinate_radius);

	//Evaluating compactness
	Evaluator mass_radius(*star,1,[&profile, this, &domain_end](int fn, double x){
		if(x == domain_end)
			return 0.0;

		//Evaluating fields
		double u = profile(x, 0);

		return 0.5*(1.0 - std::exp(-u));
	});

	//Finding maximal compactness
	auto lambda_comp = [&mass_radius](double x){return -mass_radius(x, 0);};
	Minimization max_compactness(lambda_comp);
	max_compactness.set_bracket(domain_start, 0.5*(domain_start+domain_end), domain_end);
	max_compactness.brent();
	star->params["compactness"] = -max_compactness.fmin;

	std::cout << "Central value of the scalar field is: " << star->params["phi0"] << std::endl;
	std::cout << "Misner-Sharp mass is: " << star->params["mass"] << std::endl;
	std::cout << "Star charge is: " << star->params["charge"] << std::endl;
	std::cout << RAD_FRACTION * 100.0 << "% mass reached at the radius " << star->params["radius"] << std::endl;
	std::cout << "The solution has maximal compactness: " << star->params["compactness"] << std::endl;

	//Checking if bound null geodesics are present
	Evaluator effective_compactness(*star,1,[&profile, this, &domain_end, &comp_bound](int fn, double x){
		if(x == domain_end)
			return 2.0;

		//Evaluating fields
		double dv = profile.derivative(x, 1);

		//Evaluating coordinate and rescaling factor;
		double dx = 1.0;
		if(x >= comp_bound){
			dx = dxc(x - comp_bound);
			x = comp_bound + xc(x - comp_bound);
		}

		return 2.0 - dv*x/dx;
	});

	auto lambda_geod = [&effective_compactness](double x){return effective_compactness(x);};
	Minimization geodesic_condition(lambda_geod);
	geodesic_condition.set_bracket(domain_start, 0.5*(domain_start+domain_end), domain_end);
	geodesic_condition.brent();
	if(geodesic_condition.fmin < 0.0){
		std::cout << "The star solution contains bound null geodesics\n";
		// std::cout << geodesic_condition.fmin << std::endl;
		star->params["bound_null_geodesics"] = 1.0;
	}
	else{
		std::cout << "The star solution does not contain bound null geodesics\n";
		// std::cout << geodesic_condition.fmin << std::endl;
		star->params["bound_null_geodesics"] = 0.0;
	}
}