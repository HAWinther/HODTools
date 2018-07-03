#ifndef _ROCKSTAR_HALO
#define _ROCKSTAR_HALO
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>  

// Fileformat
#define NHEADER  16
#define NCOLS    37
#define COL_ID   0
#define COL_MASS 21
#define COL_VRMS 4
#define COL_RVIR 5
#define COL_RS   6
#define COL_X    8
#define COL_Y    9
#define COL_Z    10
#define COL_VX   11
#define COL_VY   12
#define COL_VZ   13
#define COL_PID  36

class Halo {
  public:

    //===============================================
    double M;     // Mass
    double x[3];  // Pos
    double v[3];  // Vel
    int    ID;    // ID
    double PID;   // Parent ID (-1 for no parent)
    double Rvir;  // Virial radius
    double Rs;    // Scale radius
    double vrms;  // Velocity dispersion 
    //===============================================

    Halo(){}

    Halo(double _M, double *_x, double *_v, int _ID, double _PID, double _Rvir, double _Rs){
      M    = _M;
      Rvir = _Rvir;
      Rs   = _Rs;
      ID   = _ID;
      PID  = _PID;
      x[0] = _x[0];
      x[1] = _x[1];
      x[2] = _x[2];
      v[0] = _v[0];
      v[1] = _v[1];
      v[2] = _v[2];
    }
};

struct compHaloByMass {
  inline bool operator() (const Halo& lhs, const Halo& rhs){
    return (lhs.M > rhs.M);
  }
};

void readRockstarHalos(std::string filename, std::vector<Halo> &halos);

#endif
