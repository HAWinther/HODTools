#ifndef _RANDOM_HEADER
#define _RANDOM_HEADER
#include <random>
#include <vector>

#define STANDARD_SEED  12345

extern std::default_random_engine generator;
extern std::default_random_engine generator2;

std::vector<double> generateUniformVector(int n);
std::vector<double> generateNormalVector(int n, double mu, double sigma);
std::vector<double> generateNormalVector2(int n, double mu, double sigma);
int generatePoisson(double lambda);
double generateUniform();
double generateNormal(double mu, double sigma);

#endif
