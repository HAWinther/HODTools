#include "rockstar_halo.h"

void readRockstarHalos(std::string filename, std::vector<Halo> &halos){
  clock_t start = clock();
 
  std::cout << "\n====================================================\n";
  std::cout << "Reading Rockstar file: " << filename << std::endl;
  std::cout << "====================================================\n";

  std::ifstream infile(filename.c_str());

  if(!infile.good()){
    std::cout << "Error: file[ " << filename << " ] does not exist\n"; 
    exit(1);
  }

  // Count number of lines in file = number of halos + number of header lines
  int nhalos = std::count(std::istreambuf_iterator<char>(infile), 
    std::istreambuf_iterator<char>(), '\n') - NHEADER;
  
  std::cout << "Filename: " << filename << " nhalos: " << nhalos << std::endl;

  // Allocate memory
  halos = std::vector<Halo>(nhalos);

  // Rewind file
  infile.clear();
  infile.seekg(0);

  // Read headerlines
  for(int i = 0; i < NHEADER; i++){
    std::string headerline; 
    std::getline(infile, headerline);
#ifdef _DEBUG
    std::cout << headerline << std::endl;
    if(i == NHEADER - 1) std::cout << "\n";
#endif
  }

  // Read halo by halo
  double tmp[NCOLS];
  for(int i = 0; i < nhalos; i++){
    for(int j = 0; j < NCOLS; j++) 
      infile >> tmp[j];
    
    Halo *curhalo = &halos[i];

    curhalo->ID   = int(tmp[COL_ID]);
    curhalo->M    = tmp[COL_MASS];
    curhalo->vrms = tmp[COL_VRMS];
    curhalo->Rvir = tmp[COL_RVIR]/1000.0;
    curhalo->Rs   = tmp[COL_RS]/1000.0;
    curhalo->x[0] = tmp[COL_X];
    curhalo->x[1] = tmp[COL_Y];
    curhalo->x[2] = tmp[COL_Z];
    curhalo->v[0] = tmp[COL_VX];
    curhalo->v[1] = tmp[COL_VY];
    curhalo->v[2] = tmp[COL_VZ];
    curhalo->PID  = tmp[COL_PID];

#ifdef _DEBUG
    if(i < 3 || i > nhalos-4)
      std::cout << "Mass: " << curhalo->M << " x: " << curhalo->x[0] << std::endl;
#endif
  }

  // Sort halos by mass from largest to smallest
  std::sort(halos.begin(), halos.end(), compHaloByMass());

  clock_t end = clock();
  double time = (double) (end-start) / CLOCKS_PER_SEC;
  std::cout << "readRockstarHalos done. Time: " << time << " sec" << std::endl; 
  std::cout << "====================================================\n\n";
}
