

//================================================================================
double jalla(std::vector<double> &param){
  return param[0];
}
void test(Evaluation_function &f){
  std::vector<double> test {5.0};
  double x = f(test);
  std::cout << x << "\n";
}
void test2(){
  Evaluation_function f = jalla;
  test(f);
}
//================================================================================
