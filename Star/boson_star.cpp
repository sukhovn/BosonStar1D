#include <boson_star.h>
#include <par_solitonic.h>
#include <interpolation.h>
#include <limits>

Boson_Star &Boson_Star::operator<<(const Parameter_File &reader){
	int shift;

	poly_coef.resize(4);

	if(reader.name() == "Solitonic"){
		shift = 2;
		double mass = reader[0];
		double sigma = reader[1];

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
			poly_coef[i] = reader[i];
	}

	const std::vector<double> &boundaries = this->boundaries();
	kappa = reader[shift+0];
	if(std::fabs(reader[shift+1] - 1.0/boundaries[boundaries.size()-1]) > 1.0e-14){
		std::cout << "Warning! Incorrect compactification constant in the parameter file" << std::endl;
	}
	omega = reader[shift+2];
	
	phi0 = reader[shift+3];
	msmass = reader[shift+4];
	charge = reader[shift+5];
	rad = reader[shift+6];
	compactness = reader[shift+7];
	bound_null_geodesics = reader[shift+8];

	return *this;
}

void Boson_Star::operator>>(Parameter_File &reader) const{
	int shift;

	if(reader.name() == "Solitonic"){
		shift = 2;
		if(std::fabs(poly_coef[1]) < 1.0e-14 && std::fabs(poly_coef[2]) < 1.0e-14){
			reader[0] = std::sqrt(poly_coef[0]);
			reader[1] = 0.0;
		}
		else{
			bool coef_3_cond = (std::fabs(poly_coef[2] - 0.25*poly_coef[1]*poly_coef[1]/poly_coef[0])/poly_coef[2] < 1.0e-14); 
			if(coef_3_cond){
				reader[0] = std::sqrt(poly_coef[0]);
				reader[1] = std::sqrt(-4.0*poly_coef[0]/poly_coef[1]);
			}
			else{
				reader[0] = nan("");
				reader[1] = nan("");
			}
		}
	}
	else{
		shift = 4;
		for(int i = 0; i < 4; i++)
			reader[i] = poly_coef[i];
	}

	const std::vector<double> &boundaries = this->boundaries();
	reader[shift+0] = kappa;
	reader[shift+1] = 1.0/boundaries[boundaries.size()-1];
	reader[shift+2] = omega;
	reader[shift+3] = phi0;
	reader[shift+4] = msmass;
	reader[shift+5] = charge;
	reader[shift+6] = rad;
	reader[shift+7] = compactness;
	reader[shift+8] = (double) bound_null_geodesics;
}

Boson_Star::Boson_Star(const std::vector<int> &nxi, const std::vector<double> &boundi):
	Boson_Star_Equations(), Container(3, nxi, boundi),
	pfile((const char*) __Star_par_solitonic_dat, __Star_par_solitonic_dat_len),
	coordinate_komar_radius(0.0), ps_length(0.0) {};

Boson_Star::Boson_Star(const std::string &folder_name):
	Boson_Star_Equations(), Container(folder_name),
	pfile((const char*) __Star_par_solitonic_dat, __Star_par_solitonic_dat_len),
	coordinate_komar_radius(0.0), ps_length(0.0)
{
	load(folder_name);
};

void Boson_Star::load_spline(const std::string &folder_name){
	int size_x = 0;
	double xi, tmp;
	std::string path, line, token;
	std::ifstream file;
	std::string::size_type sztp;

	pfile.read_file(folder_name + "/parameters.dat");
	*this << pfile;

	file.open(folder_name + "/x.dat");
	if(!file.is_open()){
		std::cout << "No coordinate file found" << std::endl;
		throw("Read fail");
	}
	std::getline(file, line);
	std::istringstream iss(line);

	//Counting points
	while(iss.peek() != EOF){
		std::getline(iss, token, ' ');
		try{std::stod(token, &sztp);}
		catch(const std::invalid_argument& ia){break;}
		++size_x;
	}
	iss.clear();
	iss.seekg(0);

	double *x_array = new double[size_x];
	double *f_array = new double[size_x];
	double *f_sp_array = new double[size_x];

	for(int i = 0; i < size_x; ++i){
		std::getline(iss, token, ' ');
		x_array[i] = std::stod(token, &sztp);
	}

	file.close();

	for(int fn = 0; fn < n_func(); ++fn){
		path = folder_name + "/f" + std::to_string(fn) + ".dat";

		file.open(path);
		if(!file.is_open()){
			std::cout << "No function file found" << std::endl;
			throw("Read fail");
		}

		std::getline(file, line);
		iss.clear();
		iss.str(line);

		for(int i = 0; i < size_x; ++i){
			if(iss.peek() == EOF){
				std::cout << "Not enough data in file " << fn << std::endl;
				throw("Read fail");
			}
			std::getline(iss, token, ' ');
			try{f_array[i] = std::stod(token, &sztp);}
			catch(const std::invalid_argument& ia){
				std::cout << "Improper number read in file " << fn << std::endl;
				throw("Read fail");
			}
		}

		spline(x_array, f_array, size_x, 0.0, 0.0, f_sp_array);
		int indx = fn;
		for(int dm = 0; dm < n_domains(); ++dm){
			for(int i = (dm ? 1 : 0); i < n_x(dm); ++i){
				xi = coord(dm, i);
				if(splint(x_array, f_array, f_sp_array, size_x, xi, &tmp)){
					delete [] f_array; delete [] f_sp_array; delete [] x_array;
					throw("Spline interpolation failed, initialization stopped\n");
				}
				(*this)[indx] = tmp;
				indx += n_func();
			}
		}

		file.close();
	}

	delete [] x_array;
	delete [] f_array;
	delete [] f_sp_array;
}