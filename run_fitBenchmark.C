#include "fitBenchmark.C"

void run_fitBenchmark(int n = 1000) {

   useNumConv = false;
   useLogLFit = true;
   initParamMode = InitParamClose;
   
   int seed = GenerateSeed();
   std::cout << "using seed " << seed << std::endl;
   //ROOT METHODS
   fitBenchmark(n,"Fumili","",seed,0,10);
   fitBenchmark(n,"Fumili2","",seed,0,10);
   fitBenchmark(n,"Minuit2","",seed,0,10);
   fitBenchmark(n,"Minuit2","BFGS",seed,0,10);
   fitBenchmark(n,"Minuit","",seed,0,10);
   strategy=1;
   fitBenchmark(n,"Minuit2","str1",seed,0,10);
   fitBenchmark(n,"Minuit","str1",seed,0,10);
   strategy=0;
   if (!useLogLFit) {
   fitBenchmark(n,"GSLMultiFit","trust_lm",seed,0,10);
   fitBenchmark(n,"GSLMultiFit","trust_lmaccel",seed,0,10);
   fitBenchmark(n,"GSLMultiFit","trust_dogleg",seed,0,10);
   fitBenchmark(n,"GSLMultiFit","trust_ddogleg",seed,0,10);
   fitBenchmark(n,"GSLMultiFit","subspace2D",seed,0,10);
   }

   //SCIPY METHODS
   //trust-krylov, trust-exact, "trust-ncg"  (Jacobian is required for trust region)
   //"Newton-CG"  (Jacobian is required for Newton-CG method)
   // "dogleg" (Jacobian is required for dogleg minimization)
   //std::string methods[]={"Nelder-Mead","L-BFGS-B","Powell","CG","BFGS","TNC","COBYLA","SLSQP","trust-constr","Newton-CG", "dogleg", "trust-ncg","trust-exact","trust-krylov"};
   std::string methods[]={"Nelder-Mead","L-BFGS-B","Powell","CG","BFGS","TNC","COBYLA","SLSQP","trust-constr","Newton-CG", "dogleg","trust-exact","trust-krylov"};
   
   //std::string methods[]={"L-BFGS-B","Powell","CG","BFGS","TNC","COBYLA","SLSQP","trust-constr",};
   for(const std::string &method : methods)
   {
   fitBenchmark(n,"Scipy",method.c_str(),seed,0,10);
   }

   //"ma86" "ma77" fails
   // std::string methods_ipopt[]={"mumps","ma27","ma97","pardisomkl","spral"};
   // for(const std::string &method : methods_ipopt)
   // {
   // fitBenchmark(n,"Ipopt",method.c_str(),seed,0,10);
   // }


}
