#include <star_shooter.h>
#include <fstream>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

void save_grid(const char *file_name, const int nfunc, const int npoints, const double *xarr, const double *yarr){
	std::ofstream out;
	out.open(file_name);

	out << "x";
	for(int i = 0; i < nfunc; i++) out << "," << "y" << i;
	out << std::endl;
	
	out << std::setprecision(16);

	const double *ycurr = yarr;
	for(int i = 0; i < npoints; i++){
		out << xarr[i];
		for(int i = 0; i < nfunc; i++) out << "," << ycurr[i];
		out << std::endl;
		ycurr += nfunc;
	}

	out.close();
	return;
}

int main(int argc, char const *argv[]){
	std::string folder_out;
	double mass = 1.0, sigma = 0.0;
	double xst = 0.01, xend = 10.0, phi0 = 0.4;

	for(int i = 1; i < argc; i++){
		if(strcmp(argv[i], "-mass") == 0){
			mass = atof(argv[i+1]);
			i++;
		}
		else if(strcmp(argv[i], "-sigma") == 0){
			sigma = atof(argv[i+1]);
			i++;
		}
		else if(strcmp(argv[i], "-phi0") == 0){
			phi0 = atof(argv[i+1]);
			i++;
		}
		else if(strcmp(argv[i], "-xend") == 0){
			xend = atof(argv[i+1]);
			i++;
		}
		else if(strcmp(argv[i], "-out") == 0){
			folder_out = argv[i+1];
			i++;
		}
		else{
			std::cout << "Wrong key option received" << std::endl;
			return -1;
		}
	}

	double *xint = new double[2];
	xint[0] = xst; xint[1] = xend;
	double *yint = new double[8];
	for(int i = 0; i < 8; ++i) yint[i] = 0.0;

	Boson_Star_Shooter shooter(mass, sigma, phi0, xst, xend);
	std::cout << "Lagrangian mass is: " << mass << std::endl;
	std::cout << "Sigma is: " << sigma << std::endl;
	std::cout << "End of integration domain is: " << xend << std::endl;
	shooter.shoot_for_omega();

	if(!folder_out.empty()){
		//Creating output folder
		struct stat sttus = {0};
		if(stat(folder_out.c_str(), &sttus) == -1){
	    	mkdir(folder_out.c_str(), 0700);
		}

		int npoints = 1001;
		double dx = xend/(npoints - 1);
		shooter.fill_solution(npoints, [&dx](int i){return dx*i;});

		std::string file_out = folder_out + "/fields.csv";
		save_grid(file_out.c_str(), 4, npoints, shooter.xarr, shooter.yarr);

		std::ofstream pfile(folder_out + "/parameters.dat");
		pfile << "The mass value is: " << mass << std::endl;
		pfile << "The sigma value is: " << sigma << std::endl;
		pfile << "The gravitational constant kappa is: " << shooter.bstar.kappa << std::endl;
		pfile << "The start of the coordinate domain is: " << 0.0 << std::endl;
		pfile << "The end of the coordinate domain is: " << xend << std::endl;
		pfile << "The number of spatial points is: " << npoints << std::endl;
		pfile << "The central scalar field value is: " << phi0 << std::endl;
		pfile << std::setprecision(16);
		pfile << "The frequency w value is: " << shooter.bstar.omega << std::endl;
		pfile.close();
	}

	return 0;
}