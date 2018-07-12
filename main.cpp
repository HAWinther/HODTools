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
const double rmin_global =  0.1;
const double rmax_global = 80.0;
const int    nbins_global = 60;

// Parameters for simplex-search or mcmc search
double chi2_convergence_criterion = 0.001;
int nstep_max = 1000;

// Global variables: boxsize of simulations, the halo-catalog and zeta(r) for fiducial model
double box_global, box_fiducial;
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
  // Best-fit is: (0.1, 0.2, 0.3, ...)
  //===============================================
  double chi2test = 0.0;
  for(int i = 0; i < param.size(); i++){
    chi2test += pow2(param[i] / (0.1*(i+1)) - 1.0) / pow2(0.01);
  }
  return chi2test / double(param.size());

  //===============================================
#endif

  const int nbins   = nbins_global;
  const double rmin = rmin_global;
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
  BinnedCorrelationFunction *bcf = CUTER_correlation_function_periodic_from_galaxies(&mock[0], mock.size(), nbins, rmin, rmax, box);
  
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
  return chi2 / pow2(0.01) / double(np);
}

int main(int argc, char** argv){
  CUTER_set_bintype_log();
  CUTER_set_verbose( int(false) );

#ifdef USE_MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
#endif

  // The fiducial halo-catalog
  std::string filename_halo_fiducial = "halo_cat.dat";
  box_fiducial = 1024.0;

  // The halo-catalog we are to match the HOD to
  std::string filename_halo          = "halo_cat.dat";
  box_global = 1024.0;

  // The outputname
  std::string filename_output_fiducial = "mock_fiducial.txt";
  std::string filename_output          = "mock.txt";

  int do_matching = 0;

  if(argc < 8){
    std::cout << "\n====================================================\n";
    std::cout << "Run as ./HODMatch param[]\n";
    std::cout << "====================================================\n";
    std::cout << "Param [1]: filename_halocat_fiducial\n";
    std::cout << "Param [2]: boxsize_fiducial (Mpc/h)\n";
    std::cout << "Param [3]: filename_halocat\n";
    std::cout << "Param [4]: boxsize (Mpc/h)\n";
    std::cout << "Param [5]: filename_mock_output_fiducial\n";
    std::cout << "Param [6]: filename_mock_output\n";
    std::cout << "Param [7]: do_matching (0,1,2)\n";
    std::cout << "If do_matching = 0 we compute the HOD for the fiducial halo-catalog (the other input is ignored) and output it\n";
    std::cout << "If do_matching = 1 we try to do the matching using simplex search. If do_matching = 2 we do a MCMC search\n";
    std::cout << "NB: one probably *has* to edit main.cpp to set the initial guess/step-size for the parameters to get it to converge\n";
    std::cout << "====================================================\n\n";
    exit(1);
  } else {
    filename_halo_fiducial   = std::string(argv[1]);
    box_fiducial             = atof(argv[2]);
    filename_halo            = std::string(argv[3]);
    box_global               = atof(argv[4]);
    filename_output_fiducial = std::string(argv[5]);
    filename_output          = std::string(argv[6]);
    do_matching              = atoi(argv[7]);

    std::cout << "\n====================================================\n";
    std::cout << "HODMatching Parameters:                               \n";
    std::cout << "====================================================\n";
    if(do_matching == 0){
      std::cout << "Filename Fiducial HaloCat: [ " << filename_halo_fiducial   << " ]\n";
      std::cout << "Boxsize Fiducial         : [ " << box_fiducial             << " ] Mpc/h\n";
      std::cout << "Filename Output          : [ " << filename_output_fiducial << " ]\n";
    } else {
      std::cout << "Filename Fiducial HaloCat: [ " << filename_halo_fiducial   << " ]\n";
      std::cout << "Filename HaloCat         : [ " << filename_halo            << " ]\n";
      std::cout << "Boxsize Fiducial         : [ " << box_fiducial             << " ] Mpc/h\n";
      std::cout << "Boxsize                  : [ " << box_global               << " ] Mpc/h\n";
      std::cout << "Filename Fiducial Output : [ " << filename_output_fiducial << " ]\n";
      std::cout << "Filename Output          : [ " << filename_output          << " ]\n";
      if(do_matching == 1) std::cout << "Will match using simplex search\n";
      else std::cout << "Will match using mcmc search\n";
    }
    std::cout << "\n";
  }

  // Read in fiducial halo-catalog
  std::vector<Halo> halos_fiducial;
  readRockstarHalos(filename_halo_fiducial, halos_fiducial);

  // Set up HOD model (standard constructor is fiducial parameters)
  HodModel hod_fiducial;

  // Compute number density of tracers in fiducial model
  std::cout << "Expected number density in fiducial model: " << ExpectedNumberDensity(halos_fiducial, hod_fiducial, box_fiducial) << " (h/Mpc)^3\n";

  // Compute mock using fiducial parameters
  generator.seed(STANDARD_SEED);
#ifdef VELOCITY
  generator2.seed(STANDARD_SEED);
#endif
  std::vector<Galaxy> mock_fiducial;
  generateMock(halos_fiducial, mock_fiducial, hod_fiducial, box_fiducial);

  // Output fiducial mock
  output_mock(filename_output_fiducial, mock_fiducial); 
  if(do_matching == 0) return 1;

  // Compute correlation function for the fiducial mock needed for th matching
  bcf_fiducial = CUTER_correlation_function_periodic_from_galaxies(&mock_fiducial[0], mock_fiducial.size(), nbins_global, rmin_global, rmax_global, box_fiducial);

  // Read in halo-catalog for which we are to match to the fiducial one
  readRockstarHalos(filename_halo, halos_global);

  // Do matching
  std::vector<double> param, step_size;
  if(do_matching == 1){
    // Perform simplex-search
    Evaluation_function eval_func = create_mock_compute_2pcf_return_chi2;
    param     = std::vector<double>(hod_fiducial.get_num_hod_param(), 0.001); // Initial guess
    step_size = std::vector<double>(hod_fiducial.get_num_hod_param(), 0.1); // Step-size
    simplex_search(param, step_size, eval_func, chi2_convergence_criterion, nstep_max);
  } else {
    // Perform mcmc-search
    Evaluation_function eval_func = create_mock_compute_2pcf_return_chi2;
    param     = std::vector<double> (hod_fiducial.get_num_hod_param(), 0.001);   // Initial guess
    step_size = std::vector<double> (hod_fiducial.get_num_hod_param(), 0.001); // Step-size
    mcmc_search(param, step_size, eval_func, chi2_convergence_criterion, nstep_max);
  }

  // Convert simplex/mcmc params to true parameters
  HodModel hod;
  hod.simplex_param_conversion(param);
 
  // Output best fit parameters
  hod.print();

  // Generate mock
  std::vector<Galaxy> mock;
  generator.seed(STANDARD_SEED);
#ifdef VELOCITY
  generator2.seed(STANDARD_SEED);
#endif
  generateMock(halos_global, mock, hod, box_global);

  // Output mock
  output_mock(filename_output, mock); 

#ifdef USE_MPI
  MPI_Finalize();
#endif
}

