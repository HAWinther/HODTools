#include <iostream>
#include <random>
#include <fstream>
#include <vector>
#include <math.h>
#include <algorithm>  
#include <numeric>  
#include "src/rockstar_halo.h"
#include "src/cuter_library.h"
#include "src/galaxy.h"
#include "src/random_methods.h"
#include "src/simplex.h"
#include "src/hod.h"

//===========================================================
// Parameters and global variables
//===========================================================
// r-limits for evaluation of chi-square for correlation function
const double rmin_corr_chisquared_comp = 0.5;
const double rmax_corr_chisquared_comp = 100.0;

// Parameters for computation of correlation function
const double rmax_global = 100.0;
const int    nbins_global = 20;

// Parameters for simplex-search or mcmc search
double chi2_convergence_criterion = 0.01;
int nstep_max = 5000;

// Global variables: boxsize of simulations, the halo-catalog and zeta(r) for fiducial model
double box_global;
std::vector<Halo> halos_global;
BinnedCorrelationFunction *bcf_fiducial;
//===========================================================

//===========================================================
// Create a mock given HOD parameters. Compute zeta(r)
// Evaluate this relative to the reference and return the
// chi-squared value. This is used in simplex or MCMC search
//===========================================================
double create_mock_compute_2pcf_return_chi2(std::vector<double> &param){
 
#ifdef TESTING
  //===============================================
  // For testing that the simplex/mcmc search works
  //===============================================
  double chi2test = 0.0;
  for(int i = 0; i < param.size(); i++){
    chi2test += pow2(param[i] / (0.1*(i+1)) - 1.0) / pow2(0.01);
  }
  return chi2test / double(param.size());

  //===============================================
#endif

  const int nbins   = nbins_global;
  const double rmax = rmax_global;
  const double box  = box_global;

  // Set up HOD model with fiducial parameters and convert simplex-params to HOD params
  HodModel hod;
  hod.simplex_param_conversion(param);
 
  // Generate mock
  std::vector<Galaxy> mock;
  generator.seed(STANDARD_SEED);
#ifdef VELOCITY
  generator2.seed(STANDARD_SEED);
#endif
  generateMock(halos_global, mock, hod, box);

  // Compute chi(r) of mock
  BinnedCorrelationFunction *bcf = CUTER_correlation_function_periodic_from_galaxies(&mock[0], mock.size(), nbins, rmax, box);

  // Compute chi2
  double chi2 = 0.0;
  int np = 0;
  for(int i = 0; i < nbins; i++){
    if(bcf->r[i] > rmin_corr_chisquared_comp && bcf->r[i] < rmax_corr_chisquared_comp){
      chi2 += pow((bcf->corr_func[i] - bcf_fiducial->corr_func[i])/bcf_fiducial->corr_func[i],2);
      np++;
    }
  }
  free_binned_correlation_function(bcf);

  // For MCMC we need to divide this by pow2(0.01) or something to make sure we don't accept every step
  return chi2;
}

int main(int argc, char** argv){
#ifdef USE_MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
#endif
  
  // The boxsize
  box_global = 1024.0;
  
  // Read in fiducial halo-catalog
  std::string filename_halo_fiducial = "halo_cat.dat";
  std::vector<Halo> halos_fiducial;
  readRockstarHalos(filename_halo_fiducial, halos_fiducial);

  // Set up HOD model (standard constructor is fiducial parameters)
  HodModel hod_fiducial;

  // Compute number density of tracers in fiducial model
  std::cout << "Expected number density in fiducial model: " << ExpectedNumberDensity(halos_fiducial, hod_fiducial, box_global) << std::endl;

  // Compute mock using fiducial parameters
  generator.seed(STANDARD_SEED);
#ifdef VELOCITY
  generator2.seed(STANDARD_SEED);
#endif
  std::vector<Galaxy> mock_fiducial;
  generateMock(halos_fiducial, mock_fiducial, hod_fiducial, box_global);
  
  // Compute correlation function for the fiducial mock
  bcf_fiducial = CUTER_correlation_function_periodic_from_galaxies(&mock_fiducial[0], mock_fiducial.size(), nbins_global, rmax_global, box_global);

  // Read in halo-catalog for which we are to match to the fiducial one
  std::string filename_halo = "halo_cat.dat";
  readRockstarHalos(filename_halo, halos_global);

  // Perform simplex-search
  std::vector<double> start(hod_fiducial.get_num_hod_param(), 0.0); // Initial guess
  std::vector<double> dx   (hod_fiducial.get_num_hod_param(), 0.1); // Step-size
  Evaluation_function eval_func = create_mock_compute_2pcf_return_chi2;
  simplex_search(start, dx, eval_func, chi2_convergence_criterion, nstep_max);
  //mcmc_search(start, dx, eval_func, chi2_convergence_criterion, nstep_max);

#ifdef USE_MPI
  MPI_Finalize();
#endif
}

