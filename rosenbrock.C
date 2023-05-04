#include "Math/ScipyMinimizer.h"
#include "Math/MultiNumGradFunction.h"
#include <Fit/ParameterSettings.h>
#include "Math/Functor.h"
#include <string>
#include "Math/MinimizerOptions.h"
#include "TStopwatch.h"


double RosenBrock(const double *xx )
{
  const Double_t x = xx[0];
  const Double_t y = xx[1];
  const Double_t tmp1 = y-x*x;
  const Double_t tmp2 = 1-x;
  return 100*tmp1*tmp1+tmp2*tmp2;
}

double RosenBrockGrad(const double *x, unsigned int ipar)
{
   if (ipar == 0)
      return -2 * (1 - x[0]) + 200 * (x[1] - x[0] * x[0]) * (-2 * x[0]);
   else
      return 200 * (x[1] - x[0] * x[0]);
}

bool RosenBrockHessian(const std::vector<double> &xx, double *hess)
{
   const double x = xx[0];
   const double y = xx[1];
   
   hess[0] = 1200*x*x - 400*y + 2;
   hess[1] = -400*x;
   hess[2] = -400*x;
   hess[3] = 200;
   
   return true;
}

void fit(const char* minName, const char* algoName, int &success,double &time, double *xs, double &f_value, int &ncalls,int &status)
{
   ROOT::Math::Minimizer* minimizer =
      ROOT::Math::Factory::CreateMinimizer(minName, algoName);
   if (!minimizer) {
      std::cerr << "Error: cannot create minimizer \"" << minName
                << "\". Maybe the required library was not built?" << std::endl;
      return;
   }
      minimizer->SetMaxFunctionCalls(1000000);
      minimizer->SetMaxIterations(100000);
      minimizer->SetTolerance(1e-5); //increased, requested by Lorenzo
      //minimizer->SetExtraOption("gtol",1e-3);
      ROOT::Math::GradFunctor f(&RosenBrock,&RosenBrockGrad,2); 
      double step[2] = {0.01,0.01};
      double variable[2] = { -1.2,1.0};
   
      minimizer->SetFunction(f);
      //minimizer->SetHessianFunction(RosenBrockHessian);
   
      // Set the free variables to be minimized!
      minimizer->SetVariable(0,"x",variable[0], step[0]);
      minimizer->SetVariable(1,"y",variable[1], step[1]);
      minimizer->SetPrintLevel(1);
      TStopwatch t;

      t.Reset();
      t.Start();
      success = minimizer->Minimize(); 
      t.Stop();
      TString name = Form("%s_%s",minName,algoName);
      name.ReplaceAll("-","_");

      xs = (double *)minimizer->X();
      status = minimizer->Status();
      ncalls = minimizer->NCalls();
      f_value = RosenBrock(xs);
      time = t.RealTime();
      cout << "Results Method : " << minName<<" / "<<algoName;
      cout << " \t Minimum: f(" << xs[0] << "," << xs[1] << ")= " 
         << f_value;
      cout << " \t Real Time (sec) = " << time << "\t ncalls = "<<ncalls<< endl;
      cout << endl << "===============" << endl;
}

// methods that requires hessian to work "dogleg", "trust-ncg","trust-exact","trust-krylov"
using namespace std;
int rosenbrock()
{ 
   TFile fileOut("Rosenbrock.root","RECREATE");
   std::string methods[]={"Nelder-Mead","L-BFGS-B","Powell","CG","BFGS","TNC","COBYLA","SLSQP","trust-constr","Newton-CG", "dogleg", "trust-ncg","trust-exact","trust-krylov"};
   
//   ROOT::Math::MinimizerOptions::SetDefaultMaxIterations(10000);
//ROOT::Math::MinimizerOptions::SetDefaultTolerance   	
   
   TTree t1("t1","Performance Results");
   TString minName;
   TString algoName;
   int success;
   int status;
   double time;
   double xs[2];
   double f_value;
   int ncalls;
   t1.Branch("minimizer",&minName);
   t1.Branch("algo",&algoName);
   t1.Branch("success",&success,"success/I");
   t1.Branch("status",&status,"status/I");
   t1.Branch("time",&time,"time/F");
   t1.Branch("x0",&xs[0],"x0/F");
   t1.Branch("x1",&xs[1],"x1/F");
   t1.Branch("minimum",&f_value,"minimum/F");
   t1.Branch("ncalls",&ncalls,"ncalls/I");
   minName="Scipy";
   for(const std::string &text : methods)
   {
      algoName = text;
      fit(minName.Data(),algoName.Data(),success,time,xs,f_value,ncalls,status);
      t1.Fill();
   }

   std::string mmethods[]={"Migrad","Simplex","Combined","MigradImproved","Scan","Seek"};
   minName="Minuit";
   for(const std::string &text : mmethods)
   {
      algoName = text;
      fit(minName.Data(),algoName.Data(),success,time,xs,f_value,ncalls,status);
      t1.Fill();
   }

   minName="Minuit2";
   for(const std::string &text : mmethods)
   {
      algoName = text;
      fit(minName.Data(),algoName.Data(),success,time,xs,f_value,ncalls,status);
      t1.Fill();
   }
   minName="Minuit2";
   algoName ="Fumili2"; 
   fit(minName.Data(),algoName.Data(),success,time,xs,f_value,ncalls,status);
   t1.Fill();
   minName="Fumili";
   algoName =""; 
   fit(minName.Data(),algoName.Data(),success,time,xs,f_value,ncalls,status);
   t1.Fill();
   
   std::string gslmethods[]={"ConjugateFR", "ConjugatePR", "BFGS", "BFGS2", "SteepestDescent"};
     minName = "GSLMultiMin";
     for(const std::string &text : gslmethods)
    {
      algoName = text;
      fit(minName.Data(),algoName.Data(),success,time,xs,f_value,ncalls,status);
      t1.Fill();
    }
    std::string gslloglmethods[]={"trust_lm","trust_lmaccel","trust_dogleg","trust_ddogleg","subspace2D"};
     minName = "GSLMultiFit";
     for(const std::string &text : gslloglmethods)
    {
      algoName = text;
      fit(minName.Data(),algoName.Data(),success,time,xs,f_value,ncalls,status);
      t1.Fill();
    }

   minName="GSLSimAn";
   algoName =""; 
   fit(minName.Data(),algoName.Data(),success,time,xs,f_value,ncalls,status);
   t1.Fill();

   t1.Write();
   fileOut.Close();
   return 0;
}

int main()
{
  return rosenbrock();
}