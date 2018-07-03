#ifndef _GALAXY_HEADER
#define _GALAXY_HEADER

//====================================================
// A single galaxy
//====================================================
typedef struct Galaxy{
  double x[3];
#ifdef VELOCITY
  double v[3];
#endif
#ifdef WEIGHTS
  double w;
#endif
  bool central;

  Galaxy(double *_x, double *_v, bool _central){
    x[0] = _x[0];
    x[1] = _x[1];
    x[2] = _x[2];
#ifdef VELOCITY
    v[0] = _v[0];
    v[1] = _v[1];
    v[2] = _v[2];
#endif
    central = _central;
  }
  
  Galaxy(double *_x, bool _central){
    x[0] = _x[0];
    x[1] = _x[1];
    x[2] = _x[2];
    central = _central;
  }

} Galaxy;

//====================================================
// A galaxy catalog containing a list of galaxies and info
//====================================================
typedef struct GalaxyCatalog{
  int ngalaxies;    // Number of galaxies
  Galaxy *galaxies; // List of galaxies
  double sum_w;     // Sum of weight
  double sum_w2;    // Sum of weights^2
  int allocated;    // Is galaxies allocated or not?
} GalaxyCatalog;

#endif
