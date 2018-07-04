#include "random_methods.h"

std::default_random_engine generator(STANDARD_SEED);
std::default_random_engine generator2(STANDARD_SEED);

std::vector<double> generateUniformVector(int n){
  std::vector<double> v(n);
  std::uniform_real_distribution<double> distribution(0.0, 1.0);
  for(int i = 0; i < n; i++){
    v[i] = distribution(generator);
  }
  return v;
}

std::vector<double> generateNormalVector(int n, double mu, double sigma){
  std::vector<double> v(n);
  std::normal_distribution<double> distribution(mu, sigma);
  for(int i = 0; i < n; i++){
    v[i] = distribution(generator);
  }
  return v;
}

// A second one independent of the generator above to ensure consistent results
// if velocity is defined
std::vector<double> generateNormalVector2(int n, double mu, double sigma){
  std::vector<double> v(n);
  std::normal_distribution<double> distribution(mu, sigma);
  for(int i = 0; i < n; i++){
    v[i] = distribution(generator2);
  }
  return v;
}

int generatePoisson(double lambda){
  std::poisson_distribution<int> distribution(lambda);
  return distribution(generator);
}

double generateUniform(){
  std::uniform_real_distribution<double> distribution(0.0, 1.0);
  return distribution(generator);
}

double generateNormal(double mu, double sigma){
  std::normal_distribution<double> distribution(mu,sigma);
  return distribution(generator);
}
