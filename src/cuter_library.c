#include "cuter_library.h"

int mpi_rank = 0;
int mpi_size = 1;
int cuter_library_verbose = 0;
int cuter_library_logbin  = 0;
GSL_Spline *global_spline_rofz;

//====================================================
// The position of the bins 
// (Left edge of bin i  : ibin = i      )
// (Center of bin i     : ibin = i + 0.5)
// (Right edge of bin i : ibin = i + 1.0)
//====================================================
double r_binning(double ibin, int nbins, double rmin, double rmax){
  double r;
  if(cuter_library_logbin){
    // Log-bins
    r = exp( log(rmin) + log(rmax/rmin) / (double) nbins * (ibin) );
  } else {
    // Linear binning
    r = rmin + (rmax - rmin)/(double)(nbins) * (ibin);
  }
  return r;
}

void CUTER_set_bintype_log(){
  cuter_library_logbin = 1;
}

void CUTER_set_bintype_lin(){
  cuter_library_logbin = 0;
}

void CUTER_set_verbose(int verbose){
  cuter_library_verbose = verbose;
}

//====================================================
// Spline of r(z) = Int_0^z dz/H(z)
//====================================================
double r_of_z(double z){
  return Lookup_GSL_Spline(global_spline_rofz, z);
}

//====================================================
// Brute-force pair counting
//====================================================
void brute_force_pair_counting(GalaxyCatalog *cat, PairCountBinning *pc, double box){
  // Fetch data from galaxy catalog
  Galaxy *allgalaxies = cat->galaxies;
  int ngalaxies       = cat->ngalaxies;
  
  // Fetch data from binning
  int nbins   = pc->nbins;
  double rmin = pc->rmin;
  double rmax = pc->rmax;
  double *XX  = pc->paircount;
 
  // Other variables
  double pairs_dist = 0.0, pairs_dist2 = 0.0;
  double nbins_over_rmax = nbins / rmax, rmax2 = rmax * rmax, rmin2 = rmin * rmin;
  double nbins_over_logrmaxrmin = nbins / log(rmax/rmin);
  int i;
  
  // Initialize OpenMP
  int nthreads = 1, id = 0;
#if defined(USE_OMP) 
#pragma omp parallel
  {
    if(omp_get_thread_num() == 0) nthreads = omp_get_num_threads();
  }
#elif defined(USE_MPI)
  nthreads = mpi_size;
#endif

  // Initialize binning
  double **XX_threads = (double **) malloc(sizeof(double *)*nthreads);
  for(id = 0; id < nthreads; id++){
    XX_threads[id] = (double *) malloc(sizeof(double)*nbins);
    for(i = 0; i < nbins; i++){
      XX_threads[id][i] = 0.0;
      XX[i] = 0.0;
    }
  }

  if(mpi_rank == 0 && cuter_library_verbose){
    printf("\n====================================\n");
    printf("Brute-force pair counting:\n");
    printf("====================================\n");
    printf("Using n = %i threads\n", nthreads);
  }

  // Double loop over all pairs of galaxies
  clock_t start = clock(), diff;
  int istart = 0, iend = ngalaxies;
#if defined(USE_OMP)
#pragma omp parallel for private(id) reduction(+:pairs_dist2) reduction(+:pairs_dist) schedule(dynamic)
#elif defined(USE_MPI)
  int i_per_task = ngalaxies / nthreads;
  istart = i_per_task * mpi_rank;
  iend   = MIN( i_per_task * (mpi_rank + 1), ngalaxies);
#endif
  for(i = istart; i < iend; i++){
#if defined(USE_OMP)
    id = omp_get_thread_num();
#elif defined(USE_MPI)
    id = mpi_rank;
#else
    id = 0;
#endif
    Galaxy *p1 = &allgalaxies[i];
    int j;
    for(j = i+1; j < ngalaxies; j++){
      Galaxy *p2 = &allgalaxies[j];

      // Distance between galaxies
      double dx = fabs(p1->x[0] - p2->x[0]);
      double dy = fabs(p1->x[1] - p2->x[1]);
      double dz = fabs(p1->x[2] - p2->x[2]);
#ifdef PERIODIC
      if(dx > box/2.0) dx -= box;
      if(dy > box/2.0) dy -= box;
      if(dz > box/2.0) dz -= box;
#endif
      double dist2 = pow2(dx)+pow2(dy)+pow2(dz);

      // Add to bin
      if(dist2 < rmax2 && dist2 > rmin2){
        double dist = sqrt(dist2);

        // The index in the binning
        int ibin;
        if(cuter_library_logbin){
          ibin = (int) (log(dist/rmin) * nbins_over_logrmaxrmin);
        } else {
          ibin = (int) ((dist - rmin) * nbins_over_rmax);
        }

#ifdef WEIGHTS
        XX_threads[id][ibin] += p1->w * p2->w;
#else
        XX_threads[id][ibin] += 1.0;
#endif
        // Total number of pairs we have computed distances for
        pairs_dist += 1.0;
      }

      // Total number of pairs we have computed square distances for
      pairs_dist2 += 1.0;
    }
  }

#ifdef USE_MPI
  // Gather data from all CPUs
  double *recv = (double *) malloc(sizeof(double)*nbins*nthreads);
  MPI_Allgather(XX_threads[mpi_rank], nbins, MPI_DOUBLE, recv, nbins, MPI_DOUBLE, MPI_COMM_WORLD);
  for(id = 0; id < nthreads; id++)
    for(i = 0; i < nbins; i++)
      XX_threads[id][i] = recv[id * nbins + i];
  free(recv);
#endif

  // Sum up results from different threads
  for(id = 0; id < nthreads; id++){
    for(i = 0; i < nbins; i++){
      XX[i] += XX_threads[id][i];
    }
    free(XX_threads[id]);
  }
  free(XX_threads);
#ifdef PERIODIC
  for(i = 0; i < nbins; i++) XX[i] *= 2;
#endif

  // Show timing and how many dist^2 and square roots we had to compute
  diff = clock() - start;
  if(mpi_rank == 0 && cuter_library_verbose){
    printf("We have computed dist^2 for [%.0lf] pairs out of [%.0lf]. That is %lf%%\n", pairs_dist2, 
        (double)ngalaxies * (double)ngalaxies/2.0, pairs_dist2/((double)ngalaxies*(double) ngalaxies/2.0) * 100.0);
    printf("We have computed sqrt(dist^2) for [%.0lf] pairs out of [%.0lf]. That is %lf%%\n", pairs_dist, 
        (double)ngalaxies * (double)ngalaxies/2.0, pairs_dist/((double)ngalaxies*(double) ngalaxies/2.0) * 100.0);
    printf("This took: [%8.3lf] sec\n", (double)(1000 * diff / CLOCKS_PER_SEC)/1000.0);
  }
}

//====================================================
// Brute-force cross pair counting
//====================================================
void brute_force_cross_pair_counting(GalaxyCatalog *cat, GalaxyCatalog *cat2, PairCountBinning *pc, double box){
  // Fetch data from galaxy catalog
  Galaxy *allgalaxies = cat->galaxies;
  int ngalaxies       = cat->ngalaxies;

  // Fetch data from galaxy catalog2
  Galaxy *allgalaxies2 = cat2->galaxies;
  int ngalaxies2       = cat2->ngalaxies;

  // Fetch data from binning
  int nbins   = pc->nbins;
  double rmin = pc->rmin;
  double rmax = pc->rmax;
  double *XY  = pc->paircount;

  // Other variables
  double pairs_dist = 0.0, pairs_dist2 = 0.0;
  double nbins_over_rmax = nbins / rmax, rmax2 = rmax * rmax, rmin2 = rmin * rmin;
  double nbins_over_logrmaxrmin = nbins / log(rmax/rmin);
  int i;

  // Initialize OpenMP
  int nthreads = 1, id;
#if defined(USE_OMP) 
#pragma omp parallel
  {
    if(omp_get_thread_num() == 0) nthreads = omp_get_num_threads();
  }
#elif defined(USE_MPI)
  nthreads = mpi_size;
#endif

  // Initialize binning
  double **XY_threads = (double **) malloc(sizeof(double *)*nthreads);
  for(id = 0; id < nthreads; id++){
    XY_threads[id] = (double *) malloc(sizeof(double)*nbins);
    for(i = 0; i < nbins; i++){
      XY_threads[id][i] = 0.0;
      XY[i] = 0.0;
    }
  }

  if(mpi_rank == 0 && cuter_library_verbose){
    printf("\n====================================\n");
    printf("Brute-force cross pair counting:\n");
    printf("====================================\n");
    printf("Using n = %i threads\n", nthreads);
  }

  // Double loop over all pairs of galaxies
  clock_t start = clock(), diff;
  int istart = 0, iend = ngalaxies;
#if defined(USE_OMP)
#pragma omp parallel for private(id) reduction(+:pairs_dist2) reduction(+:pairs_dist) schedule(dynamic)
#elif defined(USE_MPI)
  int i_per_task = ngalaxies / nthreads;
  istart = i_per_task * mpi_rank;
  iend   = MIN( i_per_task * (mpi_rank + 1), ngalaxies);
#endif
  for(i = istart; i < iend; i++){
#if defined(USE_OMP)
    id = omp_get_thread_num();
#elif defined(USE_MPI)
    id = mpi_rank;
#else
    id = 0;
#endif
    Galaxy *p1 = &allgalaxies[i];
    int j;
    for(j = 0; j < ngalaxies2; j++){
      Galaxy *p2 = &allgalaxies2[j];

      // Distance between galaxies
      double dx = fabs(p1->x[0] - p2->x[0]);
      double dy = fabs(p1->x[1] - p2->x[1]);
      double dz = fabs(p1->x[2] - p2->x[2]);
#ifdef PERIODIC
      if(dx > box/2.0) dx -= box;
      if(dy > box/2.0) dy -= box;
      if(dz > box/2.0) dz -= box;
#endif
      double dist2 = pow2(dx)+pow2(dy)+pow2(dz);

      // Add to bin
      if(dist2 < rmax2 && dist2 > rmin2){
        double dist = sqrt(dist2);

        // The index in the binning
        int ibin;
        if(cuter_library_logbin){
          ibin = (int) (log(dist/rmin) * nbins_over_logrmaxrmin);
        } else {
          ibin = (int) ((dist-rmin) * nbins_over_rmax);
        }

#ifdef WEIGHTS
        XY_threads[id][ibin] += p1->w * p2->w;
#else
        XY_threads[id][ibin] += 1.0;
#endif
        // Total number of pairs we have computed distances for
        pairs_dist += 1.0;
      }

      // Total number of pairs we have computed square distances for
      pairs_dist2 += 1.0;
    }
  }

#ifdef USE_MPI
  // Gather data from all CPUs
  double *recv = (double *) malloc(sizeof(double)*nbins*nthreads);
  MPI_Allgather(XY_threads[mpi_rank], nbins, MPI_DOUBLE, recv, nbins, MPI_DOUBLE, MPI_COMM_WORLD);
  for(id = 0; id < nthreads; id++)
    for(i = 0; i < nbins; i++)
      XY_threads[id][i] = recv[id * nbins + i];
  free(recv);
#endif

  // Sum up results from different threads
  for(id = 0; id < nthreads; id++){
    for(i = 0; i < nbins; i++){
      XY[i] += XY_threads[id][i];
    }
    free(XY_threads[id]);
  }
  free(XY_threads);

  // Show timing and how many dist^2 and square roots we had to compute
  diff = clock() - start;
  if(mpi_rank == 0 && cuter_library_verbose){
    printf("We have computed dist^2 for [%.0lf] pairs out of [%.0lf]. That is %lf%%\n", pairs_dist2, 
        (double)ngalaxies * (double)ngalaxies2, pairs_dist2/((double)ngalaxies*(double) ngalaxies2) * 100.0);
    printf("We have computed sqrt(dist^2) for [%.0lf] pairs out of [%.0lf]. That is %lf%%\n", pairs_dist, 
        (double)ngalaxies * (double)ngalaxies2, pairs_dist/((double)ngalaxies*(double) ngalaxies2) * 100.0);
    printf("This took: [%8.3lf] sec\n", (double)(1000 * diff / CLOCKS_PER_SEC)/1000.0);
  }
}

//====================================================
// Pair counting using grid to speed it up
//====================================================
void grid_pair_counting(Grid *grid, PairCountBinning *pc, double box){
  // Fetch data from grid
  Cell *cells   = grid->cells;
  int ngrid     = grid->ngrid;
  int ngalaxies = grid->ngalaxies;
  int max_ix    = grid->max_ix;
  int max_iy    = grid->max_iy;
  int max_iz    = grid->max_iz;
  double cell_size = grid->cell_size;

  // Fetch data from binning
  int nbins   = pc->nbins;
  double rmin = pc->rmin;
  double rmax = pc->rmax;
  double *XX  = pc->paircount;

  // Other variables
  double pairs_dist2 = 0.0, pairs_dist = 0.0;
  double nbins_over_rmax = nbins / rmax, rmax2 = rmax * rmax, rmin2 = rmin * rmin;
  double nbins_over_logrmaxrmin = nbins / log(rmax/rmin);
  int i;

  // Initialize OpenMP
  int nthreads = 1, id = 0;
#if defined(USE_OMP) 
#pragma omp parallel
  {
    if(omp_get_thread_num() == 0) nthreads = omp_get_num_threads();
  }
#elif defined(USE_MPI)
  nthreads = mpi_size;
#endif

  // Initialize binning
  double **XX_threads = (double **) malloc(sizeof(double *)*nthreads);
  for(id = 0; id < nthreads; id++){
    XX_threads[id] = (double *) malloc(sizeof(double)*nbins);
    for(i = 0; i < nbins; i++){
      XX_threads[id][i] = 0.0;
      XX[i] = 0.0;
    }
  }

  if(mpi_rank == 0 && cuter_library_verbose){
    printf("\n====================================\n");
    printf("Pair counting using grid:\n");
    printf("====================================\n");
    printf("Using n = %i threads\n", nthreads);
  }

  // How many cells in each direction we must search
  int delta_ncells = (int)(ceil(rmax / cell_size)) + 1;
  if(mpi_rank == 0 && cuter_library_verbose)
    printf("We will go left and right: [%i] cells. The corresponding distance: +/-[%lf] Mpc/h\n", 
        delta_ncells,  delta_ncells * cell_size);

  //==========================================================
  // Loop over all the cells
  // The loops only go up to max_ix etc. since we wan to skip 
  // looping over cells that we know are empty
  //==========================================================

  int ix0, num_processed = 0;
  clock_t start = clock(), diff;
  int istart = 0, iend = max_ix + 1;
#if defined(USE_OMP) 
#pragma omp parallel for private(id) reduction(+:pairs_dist2) reduction(+:pairs_dist) schedule(dynamic)
#elif defined(USE_MPI)
  int i_per_task = (max_ix + 1) / nthreads;
  istart = i_per_task * mpi_rank;
  iend   = i_per_task * (mpi_rank + 1);
  if(mpi_rank == mpi_size - 1) iend = max_ix + 1;
#endif
  for(ix0 = istart; ix0 < iend; ix0++){
#if defined(USE_OMP) 
    id = omp_get_thread_num();
#elif defined(USE_MPI)
    id = mpi_rank;
#else
    id = 0;
#endif
    int iy0;
    for(iy0 = 0; iy0 <= max_iy; iy0++){
      int iz0;
      for(iz0 = 0; iz0 <= max_iz; iz0++){

        // Index of current cell
        int index = (ix0*ngrid + iy0)*ngrid + iz0;

        // Pointer to current cell
        Cell *curcell = &cells[index];

        // Number of galaxies in current cell
        int np_cell = curcell->np;

        // Loop over all galaxies in current cell
        int ipart_cell;
        for(ipart_cell = 0; ipart_cell < np_cell; ipart_cell++){

          // Current particle
          Galaxy *curpart_cell = &curcell->galaxy[ipart_cell];

          // We now want to loop over nearby cells by looking at cube of cells around current cell
#ifdef PERIODIC
          int ix_left  = -delta_ncells, ix_right = delta_ncells;
          int iy_left  = -delta_ncells, iy_right = delta_ncells;
          int iz_left  = -delta_ncells, iz_right = delta_ncells;
#else
          int ix_right = ix0 + delta_ncells <= max_ix  ? ix0 + delta_ncells : max_ix;
          int iy_right = iy0 + delta_ncells <= max_iy  ? iy0 + delta_ncells : max_iy;
          int iz_right = iz0 + delta_ncells <= max_iz  ? iz0 + delta_ncells : max_iz;
          int ix_left  = ix0 - delta_ncells >= 0       ? ix0 - delta_ncells : 0;
          int iy_left  = iy0 - delta_ncells >= 0       ? iy0 - delta_ncells : 0;
          int iz_left  = iz0 - delta_ncells >= 0       ? iz0 - delta_ncells : 0;
#endif

          // Loop over neightbor cells
          int delta_ix, delta_iy, delta_iz;
          for(delta_ix = ix_left; delta_ix <= ix_right; delta_ix++){
#ifdef PERIODIC
            int ix = ix0 + delta_ix;
            while(ix >= ngrid) ix -= ngrid;
            while(ix < 0) ix += ngrid;
#else
            int ix = delta_ix;
#endif
#ifndef PERIODIC
            // Avoid double counting so we skip cells that have been correlated with this one before
            if(ix < ix0) continue;
#endif
            for(delta_iy = iy_left; delta_iy <= iy_right; delta_iy++){
#ifdef PERIODIC
              int iy = iy0 + delta_iy;
              while(iy >= ngrid) iy -= ngrid;
              while(iy < 0) iy += ngrid;
#else
              int iy = delta_iy;
#endif
#ifndef PERIODIC
              // Avoid double counting so we skip cells that have been correlated with this one before
              if(ix == ix0 && iy < iy0) continue;
#endif
              for(delta_iz = iz_left; delta_iz <= iz_right; delta_iz++){
#ifdef PERIODIC
                int iz = iz0 + delta_iz;
                while(iz >= ngrid) iz -= ngrid;
                while(iz < 0) iz += ngrid;
#else
                int iz = delta_iz;
#endif

#ifndef PERIODIC
                // Avoid double counting so we skip cells that have been correlated with this one before
                if(ix == ix0 && iy == iy0 && iz < iz0) continue;
#endif
                // Index of neighboring cell
                int index_neighbor_cell = (ix*ngrid + iy)*ngrid + iz;

                // Pointer to neighboring cell
                Cell *neighborcell = &cells[index_neighbor_cell];

                // Number of galaxies in neighboring cell
                int npart_neighbor_cell = neighborcell->np;

                // Careful: if the nbor cell is the same as the current cell then 
                // we will overcount if we do all particles so only correlate with partices we haven't touched yet 
                int istart_nbor_cell = 0;
#ifndef PERIODIC
                if(curcell == neighborcell) istart_nbor_cell = ipart_cell + 1;
#endif
                // Loop over galaxies in neighbor cells
                int ipart_neighbor_cell;
                for(ipart_neighbor_cell = istart_nbor_cell; ipart_neighbor_cell < npart_neighbor_cell; ipart_neighbor_cell++){

                  // Galaxy in neighboring cell
                  Galaxy *curpart_neighbor_cell = &neighborcell->galaxy[ipart_neighbor_cell];

                  // ==================================================================
                  // We now count up the pair [curpart_cell] x [curpart_neighbor_cell]
                  // ==================================================================

                  // The distance between the two galaxies
                  double dx = fabs(curpart_cell->x[0] - curpart_neighbor_cell->x[0]);
                  double dy = fabs(curpart_cell->x[1] - curpart_neighbor_cell->x[1]);
                  double dz = fabs(curpart_cell->x[2] - curpart_neighbor_cell->x[2]);
#ifdef PERIODIC
                  if(dx > box/2.0) dx -= box;
                  if(dy > box/2.0) dy -= box;
                  if(dz > box/2.0) dz -= box;
#endif
                  double dist2 = pow2(dx)+pow2(dy)+pow2(dz);

                  // Add to bin
                  if(dist2 < rmax2 && dist2 > rmin2 && dist2 > 0.0){
                    double dist = sqrt(dist2);

                    // The index in the binning
                    int ibin;
                    if(cuter_library_logbin){
                      ibin = (int) (log(dist/rmin) * nbins_over_logrmaxrmin);
                    } else {
                      ibin = (int) ((dist-rmin) * nbins_over_rmax);
                    }

#ifdef WEIGHTS
                    XX_threads[id][ibin] += curpart_cell->w * curpart_neighbor_cell->w;
#else
                    XX_threads[id][ibin] += 1.0;
#endif

                    // Total number of pairs we have computed distances for
                    pairs_dist += 1.0;
                  }

                  // Total number of pairs we have computed square distances for
                  pairs_dist2 += 1.0;
                }
              }
            }
          }
        }
      }
    }

    // Show progress...
#if defined(USE_OMP)
#pragma omp critical
    {
      if(cuter_library_verbose) printf("Processed (%4i / %4i)   (ThreadID: %3i)\n", num_processed, max_ix, id);
      num_processed++;
    }
#endif
  }

#ifdef USE_MPI
  // Gather results from all CPUs
  double *recv = (double *) malloc(sizeof(double)*nbins*nthreads);
  MPI_Allgather(XX_threads[mpi_rank], nbins, MPI_DOUBLE, recv, nbins, MPI_DOUBLE, MPI_COMM_WORLD);
  for(id = 0; id < nthreads; id++)
    for(i = 0; i < nbins; i++)
      XX_threads[id][i] = recv[id * nbins + i];
  free(recv);
#endif

  // Sum up results from different threads
  for(id = 0; id < nthreads; id++){
    for(i = 0; i < nbins; i++){
      XX[i] += XX_threads[id][i];
    }
    free(XX_threads[id]);
  }
  free(XX_threads);

  // Show timing and how many dist^2 and square roots we had to compute
  diff = clock() - start;
  if(mpi_rank == 0 && cuter_library_verbose){
    printf("We have computed dist^2 for [%.0lf] pairs out of [%.0lf]. That is %lf%%\n", pairs_dist2, 
        (double)ngalaxies * (double)ngalaxies/2.0, pairs_dist2/((double)ngalaxies*(double) ngalaxies/2.0) * 100.0);
    printf("We have computed sqrt(dist^2) for [%.0lf] pairs out of [%.0lf]. That is %lf%%\n", pairs_dist, 
        (double)ngalaxies * (double)ngalaxies/2.0, pairs_dist/((double)ngalaxies*(double) ngalaxies/2.0) * 100.0);
    printf("This took: [%8.3lf] sec\n", (double)(1000 * diff / CLOCKS_PER_SEC)/1000.0);
  }
}

//====================================================
// Read galaxies from file. Format: RA DEC z
//====================================================
GalaxyCatalog *read_galaxies_from_file(char *filename, int ngalaxies, int file_format){
  if(ngalaxies == 0) return NULL;
  int i;

  if(mpi_rank == 0 && cuter_library_verbose){
    printf("\n====================================\n"); 
    printf("Reading from file [%s]\n", filename); 
    printf("====================================\n"); 
    printf("Galaxy file has %d galaxies\n",ngalaxies);
  }

  // Allocate particle array
  GalaxyCatalog *cat = (GalaxyCatalog *) malloc(sizeof(GalaxyCatalog));
  cat->galaxies = (Galaxy *) malloc(sizeof(Galaxy)*ngalaxies);
  Galaxy *allgalaxies = cat->galaxies;
  cat->ngalaxies = ngalaxies;
  cat->allocated = 1;

  double sum_w = 0.0, sum_w2 = 0.0;

  // Read the data
  FILE *fp = fopen(filename, "r");
  if(fp == NULL){
    printf("Error in read_galaxies_from_file. File [ %s ] does not exist\n", filename);
    exit(1);
  }

  for(i = 0; i < ngalaxies; i++){
    char line[1024];
    double Pos[3];

    // Read galaxy
    fgets(line, sizeof(line),fp);
   
    double w = 1.0;
    if(file_format == FILEFORMAT_RA_DEC_Z_POSITIONS){
      double RA, DEC, z;
#ifdef WEIGHTS
      sscanf(line, "%lf %lf %lf %lf", &RA, &DEC, &z, &w);
#else
      sscanf(line, "%lf %lf %lf", &RA, &DEC, &z);
#endif

      // Convert to positions in Mpc/h
      double costheta = cos(2.0 * M_PI / 360.0 * (90.0 - DEC));
      double sintheta = sqrt(1.0 - costheta*costheta);
      double phi      = RA * 2.0 * M_PI / 360.0;
      double r        = r_of_z(z);

      // Set positions
      Pos[0] = r * sintheta * cos(phi);
      Pos[1] = r * sintheta * sin(phi);
      Pos[2] = r * costheta;

    } else if (file_format == FILEFORMAT_PHYSICAL_POSITIONS) {
#ifdef WEIGHTS
      sscanf(line, "%lf %lf %lf %lf", &Pos[0], &Pos[1], &Pos[2], &w);
#else
      sscanf(line, "%lf %lf %lf", &Pos[0], &Pos[1], &Pos[2]);
#endif
    } else {
      printf("Error: unknown fileformat [%i]\n", file_format);
      exit(1);
    }

    // Store galaxies
    allgalaxies[i].x[0] = Pos[0];
    allgalaxies[i].x[1] = Pos[1];
    allgalaxies[i].x[2] = Pos[2];
#ifdef WEIGHTS
    allgalaxies[i].w = w;
#endif  
    sum_w  += w;
    sum_w2 += w*w;
  }
  fclose(fp);

  // The mean weight and RMS
  if(mpi_rank == 0 && cuter_library_verbose)
    printf("Mean weight: %lf  RMS: %lf\n", sum_w/(double)ngalaxies, sqrt(sum_w2/(double)ngalaxies));
  cat->sum_w = sum_w;
  cat->sum_w2 = sum_w2;

  return cat;
}

//====================================================
// Compute the boxsize we need to encompas all galaxies
// in both catalogs and shift the positions so that 
// they are inside [0, BOX]^3
//====================================================
void compute_boxsize_shift_positions(GalaxyCatalog *cat, GalaxyCatalog *cat2, double *box){
  int i;

  // Compute max and min position
  double max_x = -1e100, min_x = 1e100;
  double max_y = -1e100, min_y = 1e100;
  double max_z = -1e100, min_z = 1e100;

  if(cat != NULL){
    int ngalaxies  = cat->ngalaxies;
    Galaxy *galaxies = cat->galaxies;
    for(i = 0; i < ngalaxies; i++){
      double *Pos = &galaxies[i].x[0];
      if(Pos[0] > max_x) max_x = Pos[0];
      if(Pos[1] > max_y) max_y = Pos[1];
      if(Pos[2] > max_z) max_z = Pos[2];
      if(Pos[0] < min_x) min_x = Pos[0];
      if(Pos[1] < min_y) min_y = Pos[1];
      if(Pos[2] < min_z) min_z = Pos[2];
    }
  }

  if(cat2 != NULL){
    int ngalaxies2 = cat2->ngalaxies;
    Galaxy *galaxies2 = cat2->galaxies;
    for(i = 0; i < ngalaxies2; i++){
      double *Pos = &galaxies2[i].x[0];
      if(Pos[0] > max_x) max_x = Pos[0];
      if(Pos[1] > max_y) max_y = Pos[1];
      if(Pos[2] > max_z) max_z = Pos[2];
      if(Pos[0] < min_x) min_x = Pos[0];
      if(Pos[1] < min_y) min_y = Pos[1];
      if(Pos[2] < min_z) min_z = Pos[2];
    }
  }

  if(cat != NULL){
    int ngalaxies = cat->ngalaxies;
    Galaxy *galaxies = cat->galaxies;
    // Shift positions
    for(i = 0; i < ngalaxies; i++){
      double *Pos = &galaxies[i].x[0];
      Pos[0] -= min_x;
      Pos[1] -= min_y;
      Pos[2] -= min_z;
    }
  }

  if(cat2 != NULL){
    int ngalaxies2 = cat2->ngalaxies;
    Galaxy *galaxies2 = cat2->galaxies;
    for(i = 0; i < ngalaxies2; i++){
      double *Pos = &galaxies2[i].x[0];
      Pos[0] -= min_x;
      Pos[1] -= min_y;
      Pos[2] -= min_z;
    }
  }

  // New max_x,y,z values
  max_x = max_x-min_x;
  max_y = max_y-min_y;
  max_z = max_z-min_z;

  // The min/max positions (separations)
  if(mpi_rank == 0 && cuter_library_verbose){
    printf("\n====================================\n");
    printf("Shifting particles and computing boxsize:\n");
    printf("====================================\n");
    printf("Min/Max X position: 0.0 -> %5.2lf Mpc/h\n", max_x);
    printf("Min/Max Y position: 0.0 -> %5.2lf Mpc/h\n", max_y);
    printf("Min/Max Z position: 0.0 -> %5.2lf Mpc/h\n", max_z);
  }

  // The largest displacement in either direction
  *box = 1.01 * MAX(max_x, MAX(max_y, max_z));
  if(mpi_rank == 0 && cuter_library_verbose)
    printf("The boxsize we will use is %5.2lf Mpc/h\n", *box);
}

//====================================================
// Function to create a new grid
// Estimates the gridsize from rmax, box and nparticles
// and initalizes all th cells
//====================================================
Grid *create_grid(int ngalaxies, double rmax, double box){
  if(ngalaxies == 0) return NULL;

  //====================================================
  // Estimate optimal size of grid
  //====================================================
  int ngrid1 = (int)(8.0*box/rmax);                  // 8 cells to get to rmax
  int ngrid2 = (int)(pow(0.5*ngalaxies,0.3333333));  // 2 galaxies per cells on average
  int ngrid  = (int) ceil( MIN( ngrid1, ngrid2 ) );
  if(ngrid < 10) ngrid = 10;                         // Make sure we have atleast 1000 cells in total

  //====================================================
  // Make our grid. This is an array of N^3 cells
  //====================================================
  int NcellTotal  = pow3(ngrid);
  Grid *grid      = (Grid *) malloc(sizeof(Grid));
  grid->cells     = (Cell *) malloc(sizeof(Cell) * NcellTotal);
  Cell *cells     = grid->cells;
  grid->ngrid     = ngrid;
  grid->ngalaxies = ngalaxies;
  grid->cell_size = box/(double) ngrid;
  grid->allocated = 1;

  if(mpi_rank == 0 && cuter_library_verbose){
    printf("\n====================================\n"); 
    printf("Creating new grid\n"); 
    printf("====================================\n"); 
    printf("ngrid = %i [min(%i,%i)] CellSize: %lf Mpc/h\n", ngrid, ngrid1, ngrid2, box/(double) ngrid);
    printf("Total number of cells: [%i] Galaxies that will be added to grid: [%i]\n", NcellTotal, ngalaxies);
  }

  //====================================================
  // Initialize the number of galaxies in all cells
  //====================================================
  int i;
  for(i = 0; i < NcellTotal; i++){
    cells[i].np = 0;
  }

  return grid;
}

//====================================================
// Count the number of lines in a file
//====================================================
int count_lines_in_file(char *filename){
  int n = 0;
  FILE *fp = fopen(filename, "r");
  if(fp == NULL){
    printf("Error in count_lines_in_file. File [ %s ] does not exist. Return 0\n", filename);
    return 0;
  }
  while(!feof(fp)){
    char ch = fgetc(fp);
    if(ch == '\n') n++;
  }
  return n;
}

//====================================================
// Add a list of galaxies to a grid. We first
// count how many are in each cell, allocate and then add
//====================================================
void add_galaxies_to_cells(Grid *grid, GalaxyCatalog *cat){
  if(cat == NULL) return; 

  // Fetch data from grid
  Cell *cells = grid->cells;
  int ngrid   = grid->ngrid;
  double cell_size = grid->cell_size;

  // Fetch data from galaxy catalog
  Galaxy *allgalaxies = cat->galaxies;
  int ngalaxies       = cat->ngalaxies;

  // Other variables
  int NcellTotal = pow3(ngrid), i;
  int max_ix = 0, max_iy = 0, max_iz = 0;

  //====================================================
  // Count how many particles there are in each cell
  //====================================================
  for(i = 0; i < ngalaxies; i++){
    // Position of current galaxy
    double *Pos = &allgalaxies[i].x[0];

    // Determine the cell-coordinate the galaxies belongs to
    int ix = (int)(Pos[0]/cell_size);
    int iy = (int)(Pos[1]/cell_size);
    int iz = (int)(Pos[2]/cell_size);

    // Check for errors
    if(ix >= ngrid || iy >= ngrid || iz >= ngrid){
      printf("Error in add_galaxies_to_grid (%i,%i,%i)  ngrid = %i\n", ix, iy, iz, ngrid);
      exit(1);
    }

    // Find the maximum values (so we can skip looping through them later)
    if(ix > max_ix) max_ix = ix;
    if(iy > max_iy) max_iy = iy;
    if(iz > max_iz) max_iz = iz;

    // Index of cell particle belongs to in grid
    int index = (ix*ngrid + iy)*ngrid + iz;

    // Increase the number of galaxies we have in the cell
    cells[index].np += 1; 
  }

  // Max grid coordinates where we have galaxies
  max_ix = MIN(max_ix,ngrid-1);
  max_iy = MIN(max_iy,ngrid-1);
  max_iz = MIN(max_iz,ngrid-1);
#ifdef PERIODIC
  max_ix = max_iy = max_iz = ngrid-1;
#endif
  grid->max_ix = max_ix;
  grid->max_iy = max_iy;
  grid->max_iz = max_iz;

  //====================================================
  // Now that we know how many galaxies are in each cell
  // Allocate the particle array in each cell
  //====================================================
  int nempty = 0;
  for(i = 0; i < NcellTotal; i++){
    grid->cells[i].galaxy = (Galaxy *) malloc(sizeof(Galaxy) * cells[i].np);
    if(cells[i].np == 0) nempty += 1;
  }

  if(mpi_rank == 0 && cuter_library_verbose){
    printf("\n====================================\n");
    printf("Adding galaxies to grid\n");
    printf("====================================\n");
    printf("There are [%i] empty cells (%lf %%) in the grid\n", nempty, (double)nempty / (double) NcellTotal * 100.0);
  }

  //====================================================
  // Go through galaxies and add them to the grid cells
  // where they belong
  //====================================================

  // Book-keeping array when adding particles (tells how many we have added in each cell so far)
  int *ncurrent = (int *) malloc(sizeof(int)*NcellTotal);
  for(i = 0; i < NcellTotal; i++) ncurrent[i] = 0;

  for(i = 0; i < ngalaxies; i++){

    // Current particle position
    double *Pos = &allgalaxies[i].x[0];

    // Determine the cell-coordinate the galaxies belongs to
    int ix = (int)(Pos[0]/cell_size);
    int iy = (int)(Pos[1]/cell_size);
    int iz = (int)(Pos[2]/cell_size);

    // Index of cell particle belongs to in grid
    int index = (ix*ngrid + iy)*ngrid + iz;

    // Assign particle to the array in the given cell
    int partnumber = ncurrent[index];
    cells[index].galaxy[partnumber].x[0] = Pos[0];
    cells[index].galaxy[partnumber].x[1] = Pos[1];
    cells[index].galaxy[partnumber].x[2] = Pos[2];
#ifdef WEIGHTS
    cells[index].galaxy[partnumber].w = allgalaxies[i].w;
#endif

    // Increase book-keeping variable 
    ncurrent[index]++;
  }
  free(ncurrent);
}

//====================================================
// Cross pair counts using grid to speed it up
// Cross seems to be faster if we loop over the coarsest 
// grid first so call in order (galaxy_grid, random_grid)
//====================================================
void grid_cross_pair_counting(Grid *grid, Grid *grid2, PairCountBinning *pc, double box){
  // Fetch data from the grid
  Cell *cells    = grid->cells;
  int ngrid      = grid->ngrid;
  int ngalaxies  = grid->ngalaxies;
  int max_ix     = grid->max_ix;
  int max_iy     = grid->max_iy;
  int max_iz     = grid->max_iz;

  // Fetch data from the grid2
  Cell *cells2   = grid2->cells;
  int ngrid2     = grid2->ngrid;
  int ngalaxies2 = grid2->ngalaxies;
  int max_ix2    = grid2->max_ix;
  int max_iy2    = grid2->max_iy;
  int max_iz2    = grid2->max_iz;
  double cell_size2 = grid2->cell_size;

  // Fetch data from the binning
  int nbins      = pc->nbins;
  double rmin    = pc->rmin;
  double rmax    = pc->rmax;
  double *XY     = pc->paircount;

  // Other variables
  double pairs_dist2 = 0.0, pairs_dist = 0.0;
  double nbins_over_rmax = nbins / rmax, rmax2 = rmax * rmax, rmin2 = rmin * rmin;
  double nbins_over_logrmaxrmin = nbins / log(rmax/rmin);
  int i;

  // Initialize OpenMP
  int nthreads = 1, id = 0;
#if defined(USE_OMP) 
#pragma omp parallel
  {
    if(omp_get_thread_num() == 0) nthreads = omp_get_num_threads();
  }
#elif defined(USE_MPI)
  nthreads = mpi_size;
#endif

  // Initialize binning
  double **XY_threads = (double **) malloc(sizeof(double *)*nthreads);
  for(id = 0; id < nthreads; id++){
    XY_threads[id] = (double *) malloc(sizeof(double)*nbins);
    for(i = 0; i < nbins; i++){
      XY_threads[id][i] = 0.0;
      XY[i] = 0.0;
    }
  }

  if(mpi_rank == 0 && cuter_library_verbose){
    printf("\n====================================\n");
    printf("Cross pair counts using grid:\n");
    printf("====================================\n");
    printf("Using n = %i threads\n", nthreads);
  }

  // How many cells in each direction we must search in the second grid
  int delta_ncells2 = (int)(ceil(rmax / cell_size2)) + 2;
  if(mpi_rank == 0 && cuter_library_verbose)
    printf("We will go left and right: [%i] cells. The corresponding distance: +/-[%lf] Mpc/h\n", 
        delta_ncells2,  delta_ncells2 * cell_size2);

  //==========================================================
  // Loop over all the cells in grid1
  // The loops only go up to max_ix etc. since we wan to skip 
  // looping over cells that we know are empty
  //==========================================================

  int ix0, num_processed = 0;
  clock_t start = clock(), diff;
  int istart = 0, iend = max_ix + 1;
#if defined(USE_OMP) 
#pragma omp parallel for private(id) reduction(+:pairs_dist2) reduction(+:pairs_dist) schedule(dynamic)
#elif defined(USE_MPI)
  int i_per_task = (max_ix + 1) / nthreads;
  istart = i_per_task * mpi_rank;
  iend   = i_per_task * (mpi_rank + 1);
  if(mpi_rank == mpi_size - 1) iend = max_ix + 1;
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  for(ix0 = istart; ix0 < iend; ix0++){
#if defined(USE_OMP) 
    id = omp_get_thread_num();
#elif defined(USE_MPI)
    id = mpi_rank;
#else
    id = 0;
#endif
    int iy0;
    for(iy0 = 0; iy0 <= max_iy; iy0++){
      int iz0;
      for(iz0 = 0; iz0 <= max_iz; iz0++){

        // Index of current cell
        int index = (ix0*ngrid + iy0)*ngrid + iz0;

        // Pointer to current cell
        Cell *curcell = &cells[index];

        // Number of galaxies in current cell
        int np_cell = curcell->np;

        // Loop over all galaxies in current cell
        int ipart_cell;
        for(ipart_cell = 0; ipart_cell < np_cell; ipart_cell++){

          // Current particle
          Galaxy *curpart_cell = &curcell->galaxy[ipart_cell];

          //========================================
          // Now we find the index of the second grid 
          // this grid corresponds to
          //========================================
          int ix_grid2 = (int)(ix0 * ngrid2/(double)ngrid);
          int iy_grid2 = (int)(iy0 * ngrid2/(double)ngrid);
          int iz_grid2 = (int)(iz0 * ngrid2/(double)ngrid);

          // We now want to loop over nearby cells by looking at cube of cells around current cell
#ifdef PERIODIC
          int ix2_left  = -delta_ncells2, ix2_right = delta_ncells2;
          int iy2_left  = -delta_ncells2, iy2_right = delta_ncells2;
          int iz2_left  = -delta_ncells2, iz2_right = delta_ncells2;
#else
          int ix2_right = ix_grid2 + delta_ncells2 <= max_ix2  ? ix_grid2 + delta_ncells2 : max_ix2;
          int iy2_right = iy_grid2 + delta_ncells2 <= max_iy2  ? iy_grid2 + delta_ncells2 : max_iy2;
          int iz2_right = iz_grid2 + delta_ncells2 <= max_iz2  ? iz_grid2 + delta_ncells2 : max_iz2;
          int ix2_left  = ix_grid2 - delta_ncells2 >= 0        ? ix_grid2 - delta_ncells2 : 0;
          int iy2_left  = iy_grid2 - delta_ncells2 >= 0        ? iy_grid2 - delta_ncells2 : 0;
          int iz2_left  = iz_grid2 - delta_ncells2 >= 0        ? iz_grid2 - delta_ncells2 : 0;
#endif

          // Loop over neightbor cells
          int delta_ix2, delta_iy2, delta_iz2;
          for(delta_ix2 = ix2_left; delta_ix2 <= ix2_right; delta_ix2++){
#ifdef PERIODIC
            int ix2 = ix0 + delta_ix2;
            ix2 = ix2 < 0       ? ngrid2 - ix2 : ix2;
            ix2 = ix2 >= ngrid2 ? ix2 - ngrid2 : ix2;
#else
            int ix2 = delta_ix2;
#endif
            for(delta_iy2 = iy2_left; delta_iy2 <= iy2_right; delta_iy2++){
#ifdef PERIODIC
              int iy2 = iy0 + delta_iy2;
              iy2 = iy2 < 0       ? ngrid2 - iy2 : iy2;
              iy2 = iy2 >= ngrid2 ? iy2 - ngrid2 : iy2;
#else
              int iy2 = delta_iy2;
#endif
              for(delta_iz2 = iz2_left; delta_iz2 <= iz2_right; delta_iz2++){
#ifdef PERIODIC
                int iz2 = iz0 + delta_iz2;
                iz2 = iz2 < 0       ? ngrid2 - iz2 : iz2;
                iz2 = iz2 >= ngrid2 ? iz2 - ngrid2 : iz2;
#else
                int iz2 = delta_iz2;
#endif

                // Index of neighboring cell
                int index_neighbor_cell = (ix2*ngrid2 + iy2)*ngrid2 + iz2;

                // Pointer to neighboring cell
                Cell *neighborcell = &cells2[index_neighbor_cell];

                // Number of galaxies in neighboring cell
                int npart_neighbor_cell = neighborcell->np;

                // Loop over galaxies in neighbor cells
                int ipart_neighbor_cell;
                for(ipart_neighbor_cell = 0; ipart_neighbor_cell < npart_neighbor_cell; ipart_neighbor_cell++){

                  // Galaxy in neighboring cell
                  Galaxy *curpart_neighbor_cell = &neighborcell->galaxy[ipart_neighbor_cell];

                  // ==================================================================
                  // We now count up the pair [curpart_cell] x [curpart_neighbor_cell]
                  // ==================================================================
                  
                  // The distance between the two galaxies
                  double dx = fabs(curpart_cell->x[0] - curpart_neighbor_cell->x[0]);
                  double dy = fabs(curpart_cell->x[1] - curpart_neighbor_cell->x[1]);
                  double dz = fabs(curpart_cell->x[2] - curpart_neighbor_cell->x[2]);
#ifdef PERIODIC
                  if(dx > box/2.0) dx -= box;
                  if(dy > box/2.0) dy -= box;
                  if(dz > box/2.0) dz -= box;
#endif
                  double dist2 = pow2(dx)+pow2(dy)+pow2(dz);

                  // Add to bin
                  if(dist2 < rmax2 && dist2 > rmin2){
                    double dist = sqrt(dist2);

                    // The index in the binning
                    int ibin;
                    if(cuter_library_logbin){
                      ibin = (int) (log(dist/rmin) * nbins_over_logrmaxrmin);
                    } else {
                      ibin = (int) ((dist-rmin) * nbins_over_rmax);
                    }

#ifdef WEIGHTS
                    XY_threads[id][ibin] += curpart_cell->w * curpart_neighbor_cell->w;
#else
                    XY_threads[id][ibin] += 1.0;
#endif

                    // Total number of pairs we have computed distances for
                    pairs_dist += 1.0;
                  }

                  // Total number of pairs we have computed square distances for
                  pairs_dist2 += 1.0;
                }
              }
            }
          }
        }
      }
    }

    // Show progress...
#if defined(USE_OMP)
#pragma omp critical
    {
      if(cuter_library_verbose) printf("Processed (%4i / %4i)   (ThreadID: %3i)\n", num_processed, max_ix, id);
      num_processed++;
    }
#endif
  }

#ifdef USE_MPI
  // Gather results from all CPUs
  double *recv = (double *) malloc(sizeof(double)*nbins*nthreads);
  MPI_Allgather(XY_threads[mpi_rank], nbins, MPI_DOUBLE, recv, nbins, MPI_DOUBLE, MPI_COMM_WORLD);
  for(id = 0; id < nthreads; id++)
    for(i = 0; i < nbins; i++)
      XY_threads[id][i] = recv[id * nbins + i];
  free(recv);
#endif

  // Sum up results from different threads
  for(id = 0; id < nthreads; id++){
    for(i = 0; i < nbins; i++){
      XY[i] += XY_threads[id][i];
    }
    free(XY_threads[id]);
  }
  free(XY_threads);

  // Show timing and how many dist^2 and square roots we had to compute
  diff = clock() - start;
  if(mpi_rank == 0 && cuter_library_verbose){
    printf("We have computed dist^2 for [%.0lf] pairs out of [%.0lf]. That is %lf%%\n", pairs_dist2, 
        (double)ngalaxies * (double)ngalaxies2, pairs_dist2/((double)ngalaxies*(double) ngalaxies2) * 100.0);
    printf("We have computed sqrt(dist^2) for [%.0lf] pairs out of [%.0lf]. That is %lf%%\n", pairs_dist, 
        (double)ngalaxies * (double)ngalaxies2, pairs_dist/((double)ngalaxies*(double) ngalaxies2) * 100.0);
    printf("This took: [%8.3lf] sec\n", (double)(1000 * diff / CLOCKS_PER_SEC)/1000.0);
  }
}

//====================================================
// Compute correlation function directly
//====================================================
void compute_correlation_function(PairCountBinning *DD, char *filename, double box){
  // Fetch data from binning
  int nbins      = DD->nbins;
  double rmin    = DD->rmin;
  double rmax    = DD->rmax;
  double norm_DD = DD->norm;
  int i;

  if(mpi_rank == 0 && cuter_library_verbose){
    printf("\n====================================\n");
    printf("Correlation function:\n");
    printf("====================================\n");
    printf("Outputfile [%s] has format [r  xi  err_xi  DD]\n", filename);
  }

  FILE *fp;
  if(mpi_rank == 0){
    fp = fopen(filename, "w");
    if(fp == NULL){
      printf("Error in compute_correlation_function. Cannot open file [ %s ]. Will print to screen\n", filename);
    }
  }

  double rho_avg = norm_DD / pow3(box);

  for(i = 0; i < nbins; i++){
    // Left edge of bin
    double r = r_binning(i+0.0, nbins, rmin, rmax);

    // Right edge of bin
    double r1 = r_binning(i+1.0, nbins, rmin, rmax);

    // Center of bin
    double rbin = r_binning(i+0.5, nbins, rmin, rmax);

    // Compute correlation function with Poisson errors
    double corr = 0.0, err_corr = 0.0;
    if(DD->paircount[i] != 0.0){
      double one_over_sqrtDD = 1.0/sqrt(DD->paircount[i]);
      double vol = 4.0 * M_PI / 3.0 * (pow3(r1) - pow3(r));
      double rho = DD->paircount[i] / (norm_DD * vol);

      corr = rho/rho_avg - 1.0;
      err_corr = (1.0 + corr) * (one_over_sqrtDD);
    }

    // Output to file
    if(mpi_rank == 0){
#ifdef WEIGHTS
      if(fp != NULL) 
        fprintf(fp,"%le %le %le %le\n",  rbin, corr, err_corr, DD->paircount[i]);
      else
        printf("%le %le %le %le\n",  rbin, corr, err_corr, DD->paircount[i]);
#else
      if(fp != NULL) 
        fprintf(fp,"%le %le %le %d\n",  rbin, corr, err_corr, (int)DD->paircount[i]);
      else  
        printf("%le %le %le %le\n",  rbin, corr, err_corr, DD->paircount[i]);
#endif
    }
  }
  if(mpi_rank == 0)
    fclose(fp);
}

//====================================================
// Compute correlation function using LS
//====================================================
void compute_correlation_function_lz(PairCountBinning *DD, PairCountBinning *DR, PairCountBinning *RR, char *filename){
  // Fetch data from binning
  int nbins      = DD->nbins;
  double rmin    = DD->rmin;
  double rmax    = DD->rmax;
  double norm_DD = DD->norm;
  double norm_DR = DR->norm;
  double norm_RR = RR->norm;
  int i;

  if(mpi_rank == 0 && cuter_library_verbose){
    printf("\n====================================\n");
    printf("Correlation function using LS estimator:\n");
    printf("====================================\n");
    printf("Outputfile [%s] has format [r  xi  err_xi  DD  DR  RR]\n", filename);
  }

  FILE *fp;
  if(mpi_rank == 0){
    fp = fopen(filename, "w");
    if(fp == NULL){
      printf("Error in compute_correlation_function. Cannot open file [ %s ]. Will print to screen\n", filename);
    }
  }

  for(i = 0; i < nbins; i++){
    // Center of bin
    double r = r_binning(i+0.0, nbins, rmin, rmax);

    // Compute correlation function with Poisson errors
    double corr = 0.0, err_corr = 0.0;
    if(DD->paircount[i] != 0.0){

      double one_over_sqrtDD = 1.0/sqrt(DD->paircount[i]);
      double one_over_sqrtDR = 1.0/sqrt(DR->paircount[i]);
      double one_over_sqrtRR = 1.0/sqrt(RR->paircount[i]);

      double normed_DD = DD->paircount[i] / norm_DD;
      double normed_DR = DR->paircount[i] / norm_DR;
      double normed_RR = RR->paircount[i] / norm_RR;

      corr = (normed_DD - 2.0*normed_DR + normed_RR) / normed_RR;
      err_corr = (1.0 + corr) * (one_over_sqrtDD + one_over_sqrtDR + one_over_sqrtRR);
    }

    // Output to file
    if(mpi_rank == 0){
#ifdef WEIGHTS
      if(fp != NULL) 
        fprintf(fp,"%le  %le  %le  %le  %le  %le\n", r, corr, err_corr, DD->paircount[i], DR->paircount[i], RR->paircount[i]);
      else  
        printf("%le  %le  %le  %le  %le  %le\n", r, corr, err_corr, DD->paircount[i], DR->paircount[i], RR->paircount[i]);
#else
      if(fp != NULL) 
        fprintf(fp,"%le  %le  %le  %d   %d   %d\n",  r, corr, err_corr, (int)DD->paircount[i], (int)DR->paircount[i], (int)RR->paircount[i]);
      else  
        printf("%le  %le  %le  %d   %d   %d\n",  r, corr, err_corr, (int)DD->paircount[i], (int)DR->paircount[i], (int)RR->paircount[i]);
#endif
    }
  }
  if(mpi_rank == 0)
    fclose(fp);
}

//====================================================
// Allocate memory for a binning
//====================================================
PairCountBinning *create_binning(int nbins, double rmin, double rmax){
  PairCountBinning *pc = (PairCountBinning *) malloc(sizeof(PairCountBinning));
  pc->paircount = (double *) malloc(sizeof(double) * nbins);
  pc->nbins = nbins;
  pc->rmin  = rmin;
  pc->rmax  = rmax;
  pc->norm  = 1.0;
  pc->allocated = 1;
  return pc;
}

//====================================================
// Allocate memory for a binning
//====================================================
BinnedCorrelationFunction* create_binning_correlation_function(int nbins, double rmin, double rmax, double *DD, double *DR, double *RR, double *corr_func, double *err_corr){
  BinnedCorrelationFunction *tmp = (BinnedCorrelationFunction *) malloc(sizeof(BinnedCorrelationFunction));
  tmp->nbins = nbins;
  tmp->rmin  = rmin;
  tmp->rmax  = rmax;
  tmp->r     = (double *) malloc(sizeof(double)*nbins);
  tmp->DD    = (double *) malloc(sizeof(double)*nbins);
  tmp->DR    = (double *) malloc(sizeof(double)*nbins);
  tmp->RR    = (double *) malloc(sizeof(double)*nbins);
  tmp->corr_func = (double *) malloc(sizeof(double)*nbins);
  tmp->err_corr  = (double *) malloc(sizeof(double)*nbins);
  tmp->allocated = 1;
  int i;
  for(i = 0; i < nbins; i++){
    tmp->r[i] = r_binning(i+0.5, nbins, rmin, rmax);
    if(DD != NULL)   tmp->DD[i] = DD[i];
    if(DR != NULL)   tmp->DR[i] = DR[i];
    if(RR != NULL)   tmp->RR[i] = RR[i];
    if(corr_func != NULL) tmp->corr_func[i] = corr_func[i];
    if(err_corr != NULL)  tmp->err_corr[i]  = err_corr[i];
  }
  return tmp;
}

//====================================================
// Free memory associated with a binning
//====================================================
void free_binning(PairCountBinning *pc){
  if(pc == NULL) return;
  if(pc->allocated){
    free(pc->paircount);
  }
  free(pc);
}

//====================================================
// Free memory associated with a binning
//====================================================
void free_binned_correlation_function(BinnedCorrelationFunction *bcf){
  if(bcf == NULL) return;
  if(bcf->allocated){
    free(bcf->r);
    free(bcf->DD);
    free(bcf->DR);
    free(bcf->RR);
    free(bcf->corr_func);
    free(bcf->err_corr);
  }
  free(bcf);
}

//====================================================
// Free the memory associated with a galaxy catalog
//====================================================
void free_cat(GalaxyCatalog *cat){
  if(cat == NULL) return;
  if(cat->allocated){
    free(cat->galaxies);
  }
  free(cat);
}

//====================================================
// Free the memory associated with a grid
//====================================================
void free_grid(Grid *grid){
  if(grid == NULL) return;
  if(grid->allocated){
    int ngrid = grid->ngrid;
    int NcellTotal = pow3(ngrid), i;
    for(i = 0; i < NcellTotal; i++){
      free(grid->cells[i].galaxy);
    }
    free(grid->cells);
  }
  free(grid);
}

//====================================================
// Output galaxy catalog in physical coordinates
//====================================================
void outputGalaxies(GalaxyCatalog *cat, char *filename){
  if(mpi_rank != 0) return;

  int ngalaxies = cat->ngalaxies;
  Galaxy *galaxies = cat->galaxies;
  FILE *fp = fopen(filename, "w");

  int i;
  for(i = 0; i < ngalaxies; i++){
    Galaxy *curgalaxy = &galaxies[i];
#ifdef WEIGHTS
    fprintf(fp, "%lf  %lf  %lf  %lf\n", curgalaxy->x[0], curgalaxy->x[1], curgalaxy->x[2], curgalaxy->w);
#else
    fprintf(fp, "%lf  %lf  %lf\n", curgalaxy->x[0], curgalaxy->x[1], curgalaxy->x[2]);
#endif
  }
  fclose(fp);
}

//====================================================
// The ODE dr/dz = 1/H(z) in units of Mpc/h for LCDM
// Hubble_Length = c/H0
//====================================================
int ode_rofz(double z, const double r[], double drdz[], void *params){
  const double Hubble_Length_in_Mpch = SPEED_OF_LIGHT_IN_KM_PER_SEC / 100.0;
  double OmegaM = *(double *) params;

  drdz[0] =  Hubble_Length_in_Mpch / sqrt(OmegaM*pow3(1.0+z) + 1.0 - OmegaM);
  return GSL_SUCCESS;
}

//====================================================
// Integrate the ODE for r(z) and return a spline of the result
//====================================================
GSL_Spline* create_rofz_spline(double OmegaM){
  const int npts      = 1000;
  const double zmax   = 5.0;
  const double deltaz = zmax/(double) (npts-1);

  // Set up ODE system
  gsl_odeiv2_system sys_rofz = {ode_rofz, NULL, 1, &OmegaM};
  gsl_odeiv2_driver *ode = gsl_odeiv2_driver_alloc_y_new (&sys_rofz, gsl_odeiv2_step_rk2, 1e-10, 1e-10, 0.0);

  double *z_arr = (double *) malloc(sizeof(double) * npts);
  double *rofz_arr = (double *) malloc(sizeof(double) * npts);

  // Initial conditions r(z=0) = 0
  double ode_z = z_arr[0] = 0.0;
  double rofz_now[1] = {0.0};
  rofz_arr[0] = 0.0;

  // Integration over redshift
  int i;
  for(i = 1; i < npts; i++){
    double znow = i * deltaz;

    // Integrate one step
    int status = gsl_odeiv2_driver_apply(ode, &ode_z, znow, rofz_now);
    if(status != GSL_SUCCESS){
      printf("Error in integrating z = %f  r = %f\n", znow, rofz_now[0]);
      exit(1);
    }

    // Store values
    z_arr[i] = znow;
    rofz_arr[i] = rofz_now[0];
  }

  // Spline up the results
  GSL_Spline *rofz_spline = Create_GSL_Spline(z_arr, rofz_arr, npts);

  // Free memory
  free(z_arr);
  free(rofz_arr);

  return rofz_spline;
}

//====================================================
// Take a list of galaxies and copy this into a galaxy catalog
//====================================================
GalaxyCatalog *create_galaxy_catalog_from_galaxies_copy(Galaxy *galaxies, int ngalaxies){
  // Allocate particle array
  GalaxyCatalog *cat = (GalaxyCatalog *) malloc(sizeof(GalaxyCatalog));
  cat->galaxies = (Galaxy *) malloc(sizeof(Galaxy)*ngalaxies);
  Galaxy *allgalaxies = cat->galaxies;
  cat->ngalaxies = ngalaxies;
  cat->allocated = 1;

  int i;
  double sum_w = 0.0, sum_w2 = 0.0, w = 1.0;
  for(i = 0; i < ngalaxies; i++){
    allgalaxies[i].x[0] = galaxies[i].x[0];
    allgalaxies[i].x[1] = galaxies[i].x[1];
    allgalaxies[i].x[2] = galaxies[i].x[2];
#ifdef WEIGHTS
    w = galaxies[i].w;
    allgalaxies[i].w = w;
#endif  
    sum_w  += w;
    sum_w2 += w*w;
  }

  // The mean weight and RMS
  if(mpi_rank == 0 && cuter_library_verbose)
    printf("Create galaxy catalog from galaxies. Mean weight: %lf  RMS: %lf\n", sum_w/(double)ngalaxies, sqrt(sum_w2/(double)ngalaxies));
  cat->sum_w = sum_w;
  cat->sum_w2 = sum_w2;
  return cat;
}
  
//====================================================
// Take a list of galaxies and assign this to a galaxy catalog
//====================================================
GalaxyCatalog *create_galaxy_catalog_from_galaxies(Galaxy *galaxies, int ngalaxies){
  // Allocate particle array
  GalaxyCatalog *cat = (GalaxyCatalog *) malloc(sizeof(GalaxyCatalog));
  cat->galaxies = galaxies;
  cat->ngalaxies = ngalaxies;
  cat->allocated = 0;

  int i;
  double sum_w = 0.0, sum_w2 = 0.0, w = 1.0;
  for(i = 0; i < ngalaxies; i++){
#ifdef WEIGHTS
    w = galaxies[i].w;
#endif  
    sum_w  += w;
    sum_w2 += w*w;
  }

  // The mean weight and RMS
  if(mpi_rank == 0 && cuter_library_verbose)
    printf("Create galaxy catalog from galaxies. Mean weight: %lf  RMS: %lf\n", sum_w/(double)ngalaxies, sqrt(sum_w2/(double)ngalaxies));
  cat->sum_w = sum_w;
  cat->sum_w2 = sum_w2;
  return cat;
}

//====================================================
// Compute correlation function from a given galaxy catalog
//====================================================
BinnedCorrelationFunction *CUTER_correlation_function_periodic_from_catalog(GalaxyCatalog *galaxy_cat, int nbins, double rmin, double rmax, double box){
#ifdef USE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
#endif
  
  if(cuter_library_logbin && rmin <= 0.0) {
    printf("Error: CUTER called with cuter_library_logbin = 1 and rmin = %lf <= 0.0\n", rmin);
    return NULL;
  } else if (rmin < 0.0){
    printf("Error: CUTER called with rmin = %lf < 0.0\n", rmin);
    return NULL;
  }

  int ngalaxies = galaxy_cat->ngalaxies;
  Grid *galaxy_grid = create_grid(ngalaxies, rmax, box);
  add_galaxies_to_cells(galaxy_grid, galaxy_cat);
  PairCountBinning *DD = create_binning(nbins, rmin, rmax);
  DD->norm = galaxy_cat->sum_w;
#ifdef BRUTEFORCE
  brute_force_pair_counting(galaxy_cat, DD, box);
#else
  grid_pair_counting(galaxy_grid, DD, box);
#endif
  int i;
  double *corr_func = (double *) malloc(sizeof(double)*nbins);
  double *err_corr  = (double *) malloc(sizeof(double)*nbins);
  double norm_DD = DD->norm, rho_avg = norm_DD / pow3(box);
  for(i = 0; i < nbins; i++){
    double rleft  = r_binning(i+0.0, nbins, rmin, rmax);
    double rbin   = r_binning(i+0.5, nbins, rmin, rmax);
    double rright = r_binning(i+1.0, nbins, rmin, rmax);
    corr_func[i] = err_corr[i] = 0.0;
    if(DD->paircount[i] != 0.0){
      double one_over_sqrtDD = 1.0/sqrt(DD->paircount[i]);
      double vol = 4.0 * M_PI / 3.0 * (pow3(rright) - pow3(rleft));
      double rho = DD->paircount[i] / (norm_DD * vol);
      corr_func[i] = rho/rho_avg - 1.0;
      err_corr[i] = (1.0 + corr_func[i]) * (one_over_sqrtDD);
    }
  }

  BinnedCorrelationFunction *corr_binning = create_binning_correlation_function(
      nbins, rmin, rmax, DD->paircount, NULL, NULL, corr_func, err_corr);

  free_grid(galaxy_grid);
  free(corr_func);
  free(err_corr);
  free_binning(DD);

  return corr_binning;
}

//====================================================
// Compute correlation function from an array of galaxies
//====================================================
BinnedCorrelationFunction *CUTER_correlation_function_periodic_from_galaxies(Galaxy *galaxies, int ngalaxies, int nbins, double rmin, double rmax, double box){
  GalaxyCatalog *galaxy_cat = create_galaxy_catalog_from_galaxies(galaxies, ngalaxies);

  BinnedCorrelationFunction *corr_binning = CUTER_correlation_function_periodic_from_catalog(galaxy_cat, nbins, rmin, rmax, box);
  
  free(galaxy_cat);
  return corr_binning;
}

//====================================================
// Compute correction function from data in a file
//====================================================
BinnedCorrelationFunction *CUTER_correlation_function_periodic(char *filename_galaxies, int nbins, double rmin, double rmax, double box){
  int file_format = FILEFORMAT_PHYSICAL_POSITIONS;
  int ngalaxies = count_lines_in_file(filename_galaxies);
  
  GalaxyCatalog *galaxy_cat = read_galaxies_from_file(filename_galaxies, ngalaxies, file_format);
  
  BinnedCorrelationFunction *corr_binning = CUTER_correlation_function_periodic_from_catalog(galaxy_cat, nbins, rmin, rmax, box);
  
  free_cat(galaxy_cat);
  return corr_binning;
}

//====================================================
// Compute correction function from a galaxy catalog
//====================================================
BinnedCorrelationFunction *CUTER_correlation_function_from_catalog(GalaxyCatalog *galaxy_cat, GalaxyCatalog *random_cat, int nbins, double rmin, double rmax, double OmegaM){
#ifdef USE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
#endif
  
  if(cuter_library_logbin && rmin <= 0.0) {
    printf("Error: CUTER called with cuter_library_logbin = 1 and rmin = %lf <= 0.0\n", rmin);
    return NULL;
  } else if (rmin < 0.0){
    printf("Error: CUTER called with rmin = %lf < 0.0\n", rmin);
    return NULL;
  }

  int ngalaxies = galaxy_cat->ngalaxies;
  int nrandom   = random_cat->ngalaxies;

  double box;
  compute_boxsize_shift_positions(galaxy_cat, random_cat, &box);
  
  Grid *galaxy_grid = create_grid(ngalaxies, rmax, box);
  Grid *random_grid = create_grid(nrandom, rmax, box);
  add_galaxies_to_cells(galaxy_grid, galaxy_cat);
  add_galaxies_to_cells(random_grid, random_cat);
  
  PairCountBinning *DD = create_binning(nbins, rmin, rmax);
  PairCountBinning *DR = create_binning(nbins, rmin, rmax);
  PairCountBinning *RR = create_binning(nbins, rmin, rmax);
  
  DD->norm = 0.5*(pow2(galaxy_cat->sum_w) - galaxy_cat->sum_w2);
  RR->norm = 0.5*(pow2(random_cat->sum_w) - random_cat->sum_w2);
  DR->norm = galaxy_cat->sum_w * random_cat->sum_w;
#ifdef BRUTEFORCE
  brute_force_pair_counting(galaxy_cat, DD, box);
  brute_force_cross_pair_counting(random_cat, galaxy_cat, DR, box);
  brute_force_pair_counting(random_cat, RR, box);
#else
  grid_pair_counting(galaxy_grid, DD, box);
  grid_cross_pair_counting(galaxy_grid, random_grid, DR, box);
  grid_pair_counting(random_grid, RR, box);
#endif
  int i;
  double *corr_func = (double *) malloc(sizeof(double)*nbins);
  double *err_corr  = (double *) malloc(sizeof(double)*nbins);
  double norm_DD = DD->norm;
  double norm_DR = DR->norm;
  double norm_RR = RR->norm;
  for(i = 0; i < nbins; i++){
    double r = r_binning(i+0.5, nbins, rmin, rmax);
    corr_func[i] = err_corr[i] = 0.0;
    if(DD->paircount[i] != 0.0){
      double one_over_sqrtDD = 1.0/sqrt(DD->paircount[i]);
      double one_over_sqrtDR = 1.0/sqrt(DR->paircount[i]);
      double one_over_sqrtRR = 1.0/sqrt(RR->paircount[i]);
      double normed_DD = DD->paircount[i] / norm_DD;
      double normed_DR = DR->paircount[i] / norm_DR;
      double normed_RR = RR->paircount[i] / norm_RR;
      corr_func[i] = (normed_DD - 2.0*normed_DR + normed_RR) / normed_RR;
      err_corr[i] = (1.0 + corr_func[i]) * (one_over_sqrtDD + one_over_sqrtDR + one_over_sqrtRR);
    }
  }

  BinnedCorrelationFunction *corr_binning = create_binning_correlation_function(
      nbins, rmin, rmax, DD->paircount, DR->paircount, RR->paircount, corr_func, err_corr);

  free_grid(galaxy_grid);
  free_grid(random_grid);
  free_binning(DD);
  free_binning(DR);
  free_binning(RR);
  free(corr_func);
  free(err_corr);
  return corr_binning;
}

//====================================================
// Compute correction function from data in files
//====================================================
BinnedCorrelationFunction *CUTER_correlation_function(char *filename_galaxies, char* filename_random, int file_format, int nbins, double rmin, double rmax, double OmegaM){
  global_spline_rofz = create_rofz_spline(OmegaM); 
  
  int ngalaxies = count_lines_in_file(filename_galaxies);
  int nrandom   = count_lines_in_file(filename_random);
  
  GalaxyCatalog *galaxy_cat = read_galaxies_from_file(filename_galaxies, ngalaxies, file_format);
  GalaxyCatalog *random_cat = read_galaxies_from_file(filename_random, nrandom, file_format);

  BinnedCorrelationFunction *corr_binning = CUTER_correlation_function_from_catalog(galaxy_cat, random_cat, nbins, rmin, rmax, OmegaM);
  
  free_cat(galaxy_cat);
  free_cat(random_cat);
  Free_GSL_Spline(global_spline_rofz);
  return corr_binning;
}

//====================================================
// Compute correction function from an array of galaxies
// NB: this modifies the Galaxy positions by shifting them!
//====================================================
BinnedCorrelationFunction *CUTER_correlation_function_from_galaxies(Galaxy *galaxies, Galaxy *random, int ngalaxies, int nrandom, int nbins, double rmin, double rmax, double OmegaM){
  global_spline_rofz = create_rofz_spline(OmegaM); 
  
  GalaxyCatalog *galaxy_cat = create_galaxy_catalog_from_galaxies(galaxies, ngalaxies);
  GalaxyCatalog *random_cat = create_galaxy_catalog_from_galaxies(random,   nrandom);
  
  BinnedCorrelationFunction *corr_binning = CUTER_correlation_function_from_catalog(galaxy_cat, random_cat, nbins, rmin, rmax, OmegaM);
  
  free(galaxy_cat);
  free(random_cat);
  return corr_binning;
}

//====================================================
// Compute correction function from an array of galaxies
//====================================================
BinnedCorrelationFunction *CUTER_correlation_function_from_galaxies_copy(Galaxy *galaxies, Galaxy *random, int ngalaxies, int nrandom, int nbins, double rmin, double rmax, double OmegaM){
  global_spline_rofz = create_rofz_spline(OmegaM); 
 
  GalaxyCatalog *galaxy_cat = create_galaxy_catalog_from_galaxies_copy(galaxies, ngalaxies);
  GalaxyCatalog *random_cat = create_galaxy_catalog_from_galaxies_copy(random,   nrandom);
  
  BinnedCorrelationFunction *corr_binning = CUTER_correlation_function_from_catalog(galaxy_cat, random_cat, nbins, rmin, rmax, OmegaM);
  
  free_cat(galaxy_cat);
  free_cat(random_cat);
  return corr_binning;
}

