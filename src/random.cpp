#ifndef _RANDOM_HEADER
#define _RANDOM_HEADER
#include <random>
#include <vector>

#define SEED     2345

static std::default_random_engine generator(SEED);

std::vector<double> generateUniformVector(int n);
std::vector<double> generateNormalVector(int n, double mu, double sigma);
int generatePoisson(double lambda);
double generateUniform();
double generateNormal(double mu, double sigma)l

#endif
