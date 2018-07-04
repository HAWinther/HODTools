#ifndef _HOD_HEADER
#define _HOD_HEADER
#include <iostream>
#include <random>
#include <vector>
#include <math.h>
#include <algorithm>  
#include <numeric>  
#include "random_methods.h"
#include "rockstar_halo.h"
#include "galaxy.h"

class HodModel{
  private:

    // HOD parameters
    double logMmin;
    double M1;
    double M0;
    double sigma;
    double alpha;
    double Mhalo_min;

    // Fiducial parameters
    double logMmin_fiducial   = 13.09;
    double M0_fiducial        = pow(10.0,14.0);
    double M1_fiducial        = pow(10.0,13.077);
    double sigma_fiducial     = 0.596;
    double alpha_fiducial     = 1.0127;
    double Mhalo_min_fiducial = 1e10;

    // The number of parameters to fit for (if 5 we don't include Mhalo_min in simplex search)
    int num_hod_param = 5;

  public:

    HodModel(){
      logMmin   = logMmin_fiducial;
      M1        = M1_fiducial;
      M0        = M0_fiducial;
      alpha     = alpha_fiducial;
      sigma     = sigma_fiducial;
      Mhalo_min = Mhalo_min_fiducial;
    }

    HodModel(std::vector<double> &param){
      set_hod_param(param);
    }

    // This method converts the parameters we fit for in the simplex search to
    // the true parameters
    void simplex_param_conversion(std::vector<double> &param){
      logMmin   = logMmin_fiducial + param[0];
      M1        = pow(10.0, log10(M1_fiducial) + param[1]);
      M0        = pow(10.0, log10(M0_fiducial) + param[2]);
      sigma     = sigma_fiducial + param[3];
      alpha     = alpha_fiducial + param[4];
      Mhalo_min = Mhalo_min_fiducial;
      if(num_hod_param == 6) Mhalo_min += param[5];
    }

    void set_hod_param(std::vector<double> &param){
      logMmin   = param[0];
      M1        = param[1];
      M0        = param[2];
      sigma     = param[3];
      alpha     = param[4];
      Mhalo_min = param[5];
    }

    std::vector<double> get_hod_param(){
      std::vector<double> hod_param;
      hod_param.push_back(logMmin);
      hod_param.push_back(M1);
      hod_param.push_back(M0);
      hod_param.push_back(sigma);
      hod_param.push_back(alpha);
      hod_param.push_back(Mhalo_min);
      return hod_param;
    }

    double get_logMmin()  { return logMmin;   };
    double get_M1()       { return M1;        };
    double get_M0()       { return M0;        };
    double get_sigma()    { return sigma;     };
    double get_alpha()    { return alpha;     };
    double get_Mhalo_min(){ return Mhalo_min; };
    int get_num_hod_param(){ return num_hod_param; } 

    void set_logMmin(double p)  { logMmin   = p; };
    void set_M1(double p)       { M1        = p; };
    void set_M0(double p)       { M0        = p; };
    void set_sigma(double p)    { sigma     = p; };
    void set_alpha(double p)    { alpha     = p; };
    void set_Mhalo_min(double p){ Mhalo_min = p; };

    double NcentralGalaxy(double M){
      return 0.5 * (1.0 + std::erf((std::log10(M)-logMmin)/sigma));
    }

    double NsatelitteGalaxy(double M, double Ncentral){
      if(M < M0) return 0.0;
      return Ncentral * pow( (M-M0)/M1, alpha);
    }

    void print(){
      std::cout << logMmin << " vs " << logMmin_fiducial << "\n";
      std::cout << M1 << " vs " << M1_fiducial << "\n";
      std::cout << M0 << " vs " << M0_fiducial << "\n";
      std::cout << sigma << " vs " << sigma_fiducial << "\n";
      std::cout << alpha << " vs " << alpha_fiducial <<"\n";
      std::cout << Mhalo_min << " vs " << Mhalo_min_fiducial << "\n";
    }
};

double ExpectedNumberDensity(std::vector<Halo> &halos, HodModel &hod, double box);
void generateMock(std::vector<Halo> &halos, std::vector<Galaxy> &mock, HodModel &hod, double box);
void generateMock_from_halofile(std::string filename_rockstar_halos, std::string outputname, HodModel &hod, double box);
void output_mock(std::string filename, std::vector<Galaxy> &mock);
int value_locate(std::vector<double> vec, double val);

#endif
