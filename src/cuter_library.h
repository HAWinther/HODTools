#ifndef _CUTER_HEADER
#define _CUTER_HEADER
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#if defined(USE_OMP)
#include <omp.h>
#elif defined(USE_MPI)
#include <mpi.h>
#endif
#include "gsl_spline_wrapper.h"
#include "galaxy.h"
#define MIN(x,y) (x > y ? y : x)
#define MAX(x,y) (x > y ? x : y)
#define pow2(x) ((x)*(x))
#define pow3(x) ((x)*(x)*(x))
#define SPEED_OF_LIGHT_IN_KM_PER_SEC 2.99792458e5
#define FILEFORMAT_RA_DEC_Z_POSITIONS 1
#define FILEFORMAT_PHYSICAL_POSITIONS 2

//====================================================
// 
// Very simple code to compute the radial monopole correlation
// function for data from a galaxy survey / mocks. Using
// grids to speed it up.
// Same (or a bit faster) speed as CUTE with OpenMP
// Made to be used as a library for easy calls from other 
// C or C++ code
//
// For the non-periodic option the input files are assumed 
// to have format: [RA, DEC, z, (weight)]
// with RA/DEC in degree if file_format = 1 and
// [X, Y, Z, (weight)] with X,Y,Z in Mpc/h if file_format = 2
// Rmax in Mpc/h
//
// Defines:
// WEIGHTS    : Use weights in galaxy mock
// USE_OMP    : Parallelize using OpenMP
// USE_MPI    : Parallelize using MPI
// BRUTEFORCE : Brute-force pair-counts
// PERIODIC   : Periodic box, no random catalog
//
// Written by Hans Winther (2018)
//
//====================================================

extern int mpi_rank;
extern int mpi_size;
extern int cuter_library_verbose;
extern int cuter_library_logbin;

//====================================================
// This is a single grid-cell
//====================================================
typedef struct Cell{
  int np;           // Number of galaxies in the cell
  Galaxy *galaxy;   // Array of galaxies
} Cell;

//====================================================
// A grid containing a list of cells and general info
//====================================================
typedef struct Grid{
  Cell *cells;                // List of cells
  double cell_size;           // Size of cells in Mpc/h
  int ngrid;                  // Number of gridcells per dimension
  int ngalaxies;              // Number of galaxies in cell
  int max_ix, max_iy, max_iz; // No galaxies in cells (i,j,k) where i>max_ix or j>max_iy etc. 
  int allocated;              // Is cells allocated or not?
} Grid;

//====================================================
// Container for a binning
//====================================================
typedef struct PairCountBinning{
  int nbins;             // Number of linear bins between r=0 and r=RMAX
  double rmin;           // The RMIN we bin down to
  double rmax;           // The RMAX we bin up to
  double norm;           // Normalization of paircount: (weighted) total number of pairs
  double *paircount;     // The paircounts
  int allocated;         // Is paircount allocated or not?
} PairCountBinning;

//====================================================
// Container for a binning
//====================================================
typedef struct BinnedCorrelationFunction{
  int nbins;             // Number of linear bins between r=0 and r=RMAX
  double rmin;           // The RMIN we bin down to
  double rmax;           // The RMAX we bin up to
  double *r;             // The r of the bin
  double *DD;            // The paircounts
  double *DR;            // The paircounts
  double *RR;            // The paircounts
  double *corr_func;     // Correlation function
  double *err_corr;      // Poisson error for correlation function
  int allocated;         // Is binning allocated or not?
} BinnedCorrelationFunction;

//====================================================
// Global r(z) spline
//====================================================
extern GSL_Spline *global_spline_rofz;

//====================================================
// Internal methods
//====================================================

GalaxyCatalog *read_galaxies_from_file(char *filename, int npart, int file_format);
GalaxyCatalog *create_galaxy_catalog_from_galaxies(Galaxy *galaxies, int ngalaxies);
GalaxyCatalog *create_galaxy_catalog_from_galaxies_copy(Galaxy *galaxies, int ngalaxies);

Grid *create_grid(int ngalaxies, double rmax, double box);
PairCountBinning *create_binning(int nbins, double rmin, double rmax);
GSL_Spline *create_rofz_spline(double OmegaM);

double r_of_z(double z);

int  ode_rofz(double z, const double r[], double drdz[], void *params);
int  count_lines_in_file(char *filename);

void add_galaxies_to_cells(Grid *grid, GalaxyCatalog *cat);
void brute_force_pair_counting(GalaxyCatalog *cat, PairCountBinning *pc, double box);
void brute_force_cross_pair_counting(GalaxyCatalog *cat, GalaxyCatalog *cat2, PairCountBinning *pc, double box);
void grid_pair_counting(Grid *grid, PairCountBinning *pc, double box);
void grid_cross_pair_counting(Grid *grid, Grid *grid2, PairCountBinning *pc, double box);
void compute_correlation_function_lz(PairCountBinning *DD, PairCountBinning *DR, PairCountBinning *RR, char *filename);
void compute_correlation_function(PairCountBinning *DD, char *filename, double box);
void outputGalaxies(GalaxyCatalog *cat, char *filename);
void compute_boxsize_shift_positions(GalaxyCatalog *cat, GalaxyCatalog *cat2, double *box);

void free_cat(GalaxyCatalog *cat);
void free_grid(Grid *grid);
void free_binning(PairCountBinning *pc);
void free_binned_correlation_function(BinnedCorrelationFunction *bcf);

BinnedCorrelationFunction* create_binning_correlation_function(
    int nbins, double rmin, double rmax, double *DD, double *DR, double *RR, double *corr, double *err_corr);

double r_binning(double ibin, int nbins, double rmin, double rmax);

//=========================================================
// Methods for periodic box
//=========================================================
BinnedCorrelationFunction *CUTER_correlation_function_periodic(
    char *filename_galaxies, int nbins, double rmin, double rmax, double box);
BinnedCorrelationFunction *CUTER_correlation_function_periodic_from_catalog(
    GalaxyCatalog *galaxy_cat, int nbins, double rmin, double rmax, double box);
BinnedCorrelationFunction *CUTER_correlation_function_periodic_from_galaxies(
    Galaxy *galaxies, int ngalaxies, int nbins, double rmin, double rmax, double box);

//=========================================================
// Methods for survey data (LZ estimator)
// NB: CUTER_correlation_function_from_galaxies modifies galaxy positions!
// Use the _copy method to avoid changing them
//=========================================================
BinnedCorrelationFunction *CUTER_correlation_function(
    char *filename_galaxies, char* filename_random, int file_format, int nbins, double rmin, double rmax, double OmegaM);
BinnedCorrelationFunction *CUTER_correlation_function_from_catalog(
    GalaxyCatalog *galaxy_cat, GalaxyCatalog *random_cat, int nbins, double rmin, double rmax, double OmegaM);
BinnedCorrelationFunction *CUTER_correlation_function_from_galaxies(
    Galaxy *galaxies, Galaxy *random, int ngalaxies, int nrandom, int nbins, double rmin, double rmax, double OmegaM);
BinnedCorrelationFunction *CUTER_correlation_function_from_galaxies_copy(
    Galaxy *galaxies, Galaxy *random, int ngalaxies, int nrandom, int nbins, double rmin, double rmax, double OmegaM);

//=========================================================
// Set binning and verbose options externally
//=========================================================
void CUTER_set_bintype_log();
void CUTER_set_bintype_lin();
void CUTER_set_verbose(int verbose);

#endif
