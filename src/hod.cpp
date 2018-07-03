#include "hod.h"

double ExpectedNumberDensity(std::vector<Halo> &halos, HodModel &hod, double box){ 
  double xmin[3] = {1e100,1e100,1e100}, xmax[3] = {-1e100,-1e100,-1e100};

  int nhalos = halos.size();
  double ngalaxy = 0.0;
  double ncen = 0.0, nsat = 0.0;
  for(int i = 0; i < nhalos; i++){
    double M = halos[i].M;
    if(M < hod.get_Mhalo_min() || halos[i].PID > 0.0) continue;
    for(int j = 0; j < 3; j++){
      if(halos[i].x[j] < xmin[j]) xmin[j] = halos[i].x[j];
      if(halos[i].x[j] > xmax[j]) xmax[j] = halos[i].x[j];
    }

    //double Ncentral_EV     = NcentralGalaxy(M, logMmin, sigma);
    double Ncentral_EV     = hod.NcentralGalaxy(M);
    double Nsatelittes_EV  = hod.NsatelitteGalaxy(M, Ncentral_EV);

    ngalaxy += Ncentral_EV + Nsatelittes_EV;
    ncen += Ncentral_EV;
    nsat += Nsatelittes_EV;
  }

  double volume = pow(box, 3);
  return ngalaxy / volume;
}

void generateMock(std::vector<Halo> &halos, std::vector<Galaxy> &mock, HodModel &hod, double box){
  clock_t start = clock();
  const int nbins = 100;
  const bool is_central_galaxy = true;

  // Loop over halos 
  int nhalos = halos.size();
  int ngalaxy = 0;
  int ncentral = 0;
  int nsatelittes = 0;
  for(int i = 0; i < nhalos; i++){
    // Current halo
    Halo *curhalo = &halos[i];

    // Only process halos above treshold and those that have no parent
    double M = curhalo->M;
    if(M < hod.get_Mhalo_min() || halos[i].PID > 0.0) continue;

    double Ncentral_EV = hod.NcentralGalaxy(M);
    double Nsat_EV     = hod.NsatelitteGalaxy(M, Ncentral_EV);

    double r1 = generateUniform();
    int Nsat  = generatePoisson(Nsat_EV);

    // Central galaxy
    if(r1 < Ncentral_EV){

      // Create galaxy
#ifdef VELOCITY
      Galaxy g(curhalo->x, curhalo->v, is_central_galaxy);
#else
      Galaxy g(curhalo->x, is_central_galaxy);
#endif
      // Add to mock
      mock.push_back(g);

      ++ngalaxy;
      ++ncentral;
    } else {

      if(Nsat > 0){

        // Create galaxy
#ifdef VELOCITY
        Galaxy g(curhalo->x, curhalo->v, is_central_galaxy);
#else
        Galaxy g(curhalo->x, is_central_galaxy);
#endif
        // Add to mock
        mock.push_back(g);

        --Nsat;
        ++ngalaxy;
        ++ncentral;
      }
    }

    // Populate with satelittes
    if(Nsat > 0) {
      std::vector<double> u = generateUniformVector(Nsat);
      std::vector<double> v = generateUniformVector(Nsat);
      std::vector<double> w = generateUniformVector(Nsat);
#ifdef VELOCITY
      std::vector<double> xx = generateNormalVector2(Nsat, 0.0, 1.0);
      std::vector<double> yy = generateNormalVector2(Nsat, 0.0, 1.0);
      std::vector<double> zz = generateNormalVector2(Nsat, 0.0, 1.0);
#endif

      // Consentration of halo
      double cons = std::min(3.0, curhalo->Rvir / curhalo->Rs);
    
      // Make radial bins with mass-fraction
      std::vector<double> rbins(nbins);
      std::vector<double> massfracbins(nbins);
      double massnorm = log(1.0 + cons) - cons / (1.0 + cons);
      for(int j = 0; j < nbins; j++){
        rbins[j] = j/double(nbins);
        massfracbins[j] = (log(1.0 + cons*rbins[j]) - cons*rbins[j]/(1.0 + cons*rbins[j])) / massnorm;
      }

      // Make Nsat satelittes
      for(int j = 0; j < Nsat; j++){
        int loc = value_locate(massfracbins, u[j]);

        // Angle-position of satelitte within halo
        double rr = (rbins[loc] + 0.5/double(nbins)) * curhalo->Rvir;
        double tt = std::acos(v[j]*2.0 - 1.0);
        double pp = w[j]*2.0*M_PI;

        // Physical position of satelitte within halo
        double dx[3];
        dx[0] = rr*sin(tt)*cos(pp);
        dx[1] = rr*sin(tt)*sin(pp);
        dx[2] = rr*cos(tt);
    
        // Create galaxy
        Galaxy g(curhalo->x, curhalo->v, !is_central_galaxy);
       
        // Add displacement
        for(int axis = 0; axis < 3; axis++){
          g.x[axis] += dx[axis];
          if(g.x[axis] > box) g.x[axis] -= box;
          if(g.x[axis] < 0.0) g.x[axis] += box;
        }

#ifdef VELOCITY
        // Add velocity displacement
        double dv[3];
        dv[0] = xx[j] * curhalo->vrms/sqrt(3.0);
        dv[1] = yy[j] * curhalo->vrms/sqrt(3.0);
        dv[2] = zz[j] * curhalo->vrms/sqrt(3.0);
        for(int axis = 0; axis < 3; axis++){
          g.v[axis] += dv[axis];
        }
#endif

        // Add to mock
        mock.push_back(g);

        ++nsatelittes;
        ++ngalaxy;
      }
    }
  }

  clock_t end = clock();
  double time = (double) (end-start) / CLOCKS_PER_SEC;
#ifdef _DEBUG
  std::cout << "Ncentral: " << ncentral << " Nsat: " << nsatelittes << " Total: " << ngalaxy << std::endl; 
  std::cout << "generateMock done. Time: " << time << " sec" << std::endl; 
#endif
}

//===========================================================
// Read a halo-catalog, create a mock and output it
//===========================================================
void generateMock_from_halofile(std::string filename_rockstar_halos, std::string outputname, HodModel &hod, double box){
  std::vector<Halo> halos;
  readRockstarHalos(filename_rockstar_halos, halos);
  std::vector<Galaxy> mock;
  generateMock(halos, mock, hod, box);
  output_mock(outputname, mock);
}

//===========================================================
// Output a mock
//===========================================================
void output_mock(std::string filename, std::vector<Galaxy> &mock){
  std::ofstream fp(filename.c_str());
  for(auto g: mock){
    fp << g.x[0] << " " << g.x[1] << " " << g.x[2] << " ";
#ifdef VELOCITY
    fp << g.v[0] << " " << g.v[1] << " " << g.v[2];
#endif
    fp << "\n";
  }
}

int value_locate(std::vector<double> vec, double val){
  int nvec = vec.size();
  if(val < vec[0]) return -1;
  for(int i = 1; i < nvec; i++)
    if(val < vec[i]) return i-1;
  return nvec-1;
}

