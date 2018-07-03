#ifndef _SIMPLEX_HEADER
#define _SIMPLEX_HEADER
#include <iostream>
#include <vector>
#include <math.h>
#include <algorithm>  
#include <iomanip>
#include <numeric>  
#include "random_methods.h"

typedef double (*Evaluation_function)(std::vector<double> &param);

template <typename Container>
struct compare_indirect_index{
  const Container& container;
  compare_indirect_index( const Container& container ): container( container ) { }
  bool operator () ( size_t lindex, size_t rindex ) const{
    return container[ lindex ] < container[ rindex ];
  }
};

std::vector<size_t> get_sorted_indices(std::vector<double> &numbers);
void simplex_search(std::vector<double> &start, std::vector<double> dx, Evaluation_function &eval_func, double chi2_converge_limit, int istep_max);
void mcmc_search(std::vector<double> &start, std::vector<double> &step_size, Evaluation_function &eval_func, double chi2_converge_limit, int istep_max);

#endif
