#include "simplex.h"

std::vector<size_t> get_sorted_indices(std::vector<double> &numbers){
  std::vector <size_t> indices( numbers.size(), 0 );
  std::iota( indices.begin(), indices.end(), 0 );  // found in <numeric>
  std::sort( indices.begin(), indices.end(), compare_indirect_index <decltype(numbers)> ( numbers ) );
  return indices;
}

void simplex_search(std::vector<double> &start, std::vector<double> dx, Evaluation_function &eval_func, double chi2_converge_limit, int istep_max){
  int npar = start.size();
  double simplex[npar][npar+1];
  double chi2_r, chi2_e, chi2_c, alp = 1.0, gam = 2.0, del = 0.5, bet = 0.5;
  std::vector<double> chi2(npar+1);
  std::vector<double> x_r(npar), x_c(npar), x_e(npar), centroid(npar);

  std::cout << "\n====================================================\n";
  std::cout << "Starting simplex search: \n";
  std::cout << "====================================================\n";
  
  for(int i = 0; i < npar; i++){
    for(int j = 0; j < npar; j++){
      simplex[i][j] = 0.0;
    }
  }

  // Set starting value
  for(int i = 0; i < npar; i++){
    simplex[i][0] = start[i];
  }

  for(int i = 0; i < npar; i++){
    for(int j = 0; j < npar; j++) simplex[j][i+1] = simplex[j][0];
    simplex[i][i+1] = simplex[i][i+1] + dx[i];
  }

  for(int i = 0; i < npar + 1; i++){
    std::vector<double> param_0(npar);
    for(int j = 0; j < npar; j++) param_0[j] = simplex[j][i];
    chi2[i] = eval_func(param_0);
  }
  std::cout << "0.0 Param: ["; for(int j = 0; j < npar; j++) { std::cout << " " << std::setw(12) << simplex[j][0];} std::cout << "] Chi2: " << std::setw(12) << chi2[0] << "\n";

  int icount = 1;
  while(icount < istep_max){
    ++icount;
    std::vector<size_t> sorted_index = get_sorted_indices(chi2);

    int ind_h1 = sorted_index[npar];
    int ind_h2 = sorted_index[npar-1];
    int ind_l1 = sorted_index[0];
  
    std::vector<double> centroid(npar,0.0);
    for(int i = 0; i < npar + 1; i++){
      if(i == ind_h1) continue;
      for(int j = 0; j < npar; j++) centroid[j] = centroid[j] + simplex[j][i]/double(npar);
    }
      
    for(int j = 0; j < npar; j++) x_r[j] = centroid[j] + alp * (centroid[j] - simplex[j][ind_h1]);
    
    std::vector<double> param_1(x_r);
    chi2_r = eval_func(param_1);
    
    if(chi2_r < chi2[ind_l1]){
      for(int j = 0; j < npar; j++) x_e[j] = centroid[j] + gam * (x_r[j] - centroid[j]);
      
      std::vector<double> param_2(x_e);
      chi2_e = eval_func(param_2);
      
      if(chi2_e < chi2_r){
        for(int j = 0; j < npar; j++) simplex[j][ind_h1] = x_e[j];
        chi2[ind_h1] = chi2_e;
      } else {
        for(int j = 0; j < npar; j++) simplex[j][ind_h1] = x_r[j];
        chi2[ind_h1] = chi2_r;
      }
      std::cout << "1.0 Param: ["; for(int j = 0; j < npar; j++) { std::cout << " " << std::setw(12) << simplex[j][ind_h1];} std::cout << "] Chi2: " << std::setw(12) << chi2[ind_h1] << "\n";
    } else if(chi2_r < chi2[ind_h2]) {
      for(int j = 0; j < npar; j++) simplex[j][ind_h1] = x_r[j];
      chi2[ind_h1] = chi2_r;
      std::cout << "2.0 Param: ["; for(int j = 0; j < npar; j++) { std::cout << " " << std::setw(12) << simplex[j][ind_h1];} std::cout << "] Chi2: " << std::setw(12) << chi2[ind_h1] << "\n";
    } else if(chi2_r > chi2[ind_h2] && chi2_r < chi2[ind_h1]){
      for(int j = 0; j < npar; j++) x_c[j] = centroid[j] + bet * (x_r[j] - centroid[j]);
      
      std::vector<double> param_3(x_c);
      chi2_c = eval_func(param_3);
      
      if(chi2_c < chi2_r){
        for(int j = 0; j < npar; j++) simplex[j][ind_h1] = x_c[j];
        chi2[ind_h1] = chi2_c;
        std::cout << "3.0 Param: ["; for(int j = 0; j < npar; j++) { std::cout << " " << std::setw(12) << simplex[j][ind_h1];} std::cout << "] Chi2: " << std::setw(12) << chi2[ind_h1] << "\n";
      } else {
        for(int i = 0; i < npar + 1; i++){
          if(i == ind_l1) continue;
          for(int j = 0; j < npar; j++) simplex[j][i] = simplex[j][ind_l1] + del * (simplex[j][i] - simplex[j][ind_l1]);
          
          std::vector<double> param_4(npar);
          for(int j = 0; j < npar; j++) param_4[j] = simplex[j][i];
          chi2[i] = eval_func(param_4);
          
          std::cout << "3.5 Param: ["; for(int j = 0; j < npar; j++) { std::cout << " " << std::setw(12) << simplex[j][i];} std::cout << "] Chi2: " << std::setw(12) << chi2[i] << "\n";
        }
      }
    } else if(chi2_r > chi2[ind_h1]){
      for(int j = 0; j < npar; j++) x_c[j] = centroid[j] + bet * (simplex[j][ind_h1] - centroid[j]);
          
      std::vector<double> param_5(x_c);
      chi2_c = eval_func(param_5);
      
      if(chi2_c < chi2[ind_h1]){
        for(int j = 0; j < npar; j++) simplex[j][ind_h1] = x_c[j];
        chi2[ind_h1] = chi2_c;
        std::cout << "4.0 Param: ["; for(int j = 0; j < npar; j++) { std::cout << " " << std::setw(12) << simplex[j][ind_h1];} std::cout << "] Chi2: " << std::setw(12) << chi2[ind_h1] << "\n";
      } else {
        for(int i = 0; i < npar + 1; i++){
          if(i == ind_l1) continue;
          for(int j = 0; j < npar; j++) simplex[j][i] = simplex[j][ind_l1] + del * (simplex[j][i] - simplex[j][ind_l1]);
      
          std::vector<double> param_6(npar);
          for(int j = 0; j < npar; j++) param_6[j] = simplex[j][i];
          chi2[i] = eval_func(param_6);
          
          std::cout << "4.5 Param: ["; for(int j = 0; j < npar; j++) { std::cout << " " << std::setw(12) << simplex[j][i];} std::cout << "] Chi2: " << std::setw(12) << chi2[i] << "\n";
        }
      }
    }
    // Break statement:
    bool converged = false;
    for(int j = 0; j < npar; j++)
      if(chi2[j] < chi2_converge_limit) converged = true;
    if(converged) break;
  }
  std::cout << "We finished the search after icount = " << icount << "\n";

  // Copy over the best-fit parameters
  for(int i = 0; i < npar; i++) start[i] = simplex[i][0];
}

void mcmc_search(std::vector<double> &start, std::vector<double> &step_size, Evaluation_function &eval_func, double chi2_converge_limit, int istep_max){
  int nparam = start.size();

  std::vector<double> current_param = start, min_param = start;
  double chi2 = eval_func(current_param), chi2_min = chi2;
  int nstep = 0, naccept = 0;
  while(nstep < istep_max){
    
    // Make a new proposal
    std::vector<double> proposal_param(nparam);
    for(int j = 0; j < nparam; j++){
      proposal_param[j] = current_param[j] + step_size[j] * generateNormal(0.0, 1.0);
    }

    // Compute chi2
    double chi2_new = eval_func(proposal_param);
    
    // Metropolis step
    double ratio = std::exp(-0.5 * (chi2_new - chi2));
    if( ratio > generateUniform() ){
      current_param = proposal_param;
      chi2 = chi2_new;
      if(chi2 < chi2_min){
        chi2_min = chi2;
        min_param = current_param;
      }
      naccept++;
#ifdef _DEBUG
      std::cout << "Ratio: " << ratio << " ACCEPT!" << std::endl;
#endif
    }
    nstep++;

    std::cout << "Acceptance ratio: " << naccept / double(nstep) << " Chi2: " << chi2 << "\nParam: [";
    for(int j = 0; j < nparam; j++) std::cout << std::setw(12) << min_param[j] << " ";
    std::cout << "]\n";
  
    if(chi2 < chi2_converge_limit) break;
  }
  std::cout << "We finished the search after nstep = " << nstep << " MinChi2: " << chi2 << "\nParam:    [";
  for(int j = 0; j < nparam; j++) std::cout << std::setw(12) << min_param[j] << " ";
  std::cout << "]\n";

#ifdef TESTING
  std::cout << "Expected: [";
  for(int j = 0; j < nparam; j++) std::cout << std::setw(12) << 0.1*(j+1) << " ";
  std::cout << "]\n";
#endif

  // Copy over the parameters
  start = min_param;
}
