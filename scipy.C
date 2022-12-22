#include "Math/ScipyMinimizer.h"
#include "Math/MultiNumGradFunction.h"
#include "Math/Functor.h"
#include <string>
#include "Math/MinimizerOptions.h"

class Wood4GradientFunction : public ROOT::Math::IGradientFunctionMultiDim {
public:
   double DoEval(const double *par) const
   {
      const Double_t w = par[0];
      const Double_t x = par[1];
      const Double_t y = par[2];
      const Double_t z = par[3];

     const Double_t w1 = w-1;
   const Double_t x1 = x-1;
   const Double_t y1 = y-1;
   const Double_t z1 = z-1;
   const Double_t tmp1 = x-w*w;
   const Double_t tmp2 = z-y*y;

  return 100*tmp1*tmp1+w1*w1+90*tmp2*tmp2+y1*y1+10.1*(x1*x1+z1*z1)+19.8*x1*z1;
   }
   unsigned int NDim() const { return 4; }
   ROOT::Math::IGradientFunctionMultiDim *Clone() const { return new Wood4GradientFunction(); }
   double DoDerivative(const double *par, unsigned int ipar) const
   {
      const Double_t w = par[0];
      const Double_t x = par[1];
      const Double_t y = par[2];
      const Double_t z = par[3];

      if (ipar == 0)
         return 400*(y-w*w)*w+2*(w-1);
      if (ipar == 1)
         return 20.2*(x-1)+19.8*(z-1);
      if (ipar == 2)
         return 200*(y-w*w)+360*(z-y*y)*y+2*(1-y);
      if (ipar == 3)
         return 180*(z-y*y)+20.2*(z-1)+19.8*(x-1);
   }
};


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

// methods that requires hessian to work "dogleg", "trust-ncg","trust-exact","trust-krylov"
using namespace std;
int scipy()
{ 
   
   std::string methods[]={"Nelder-Mead","L-BFGS-B","Powell","CG","BFGS","TNC","COBYLA","SLSQP","trust-constr","Newton-CG", "dogleg", "trust-ncg","trust-exact","trust-krylov"};
   //Wood4GradientFunction wgf;
   for(const std::string &text : methods)
   {
   ROOT::Math::Experimental::ScipyMinimizer minimizer(text.c_str());
   minimizer.SetMaxFunctionCalls(1000000);
   minimizer.SetMaxIterations(100000);
   minimizer.SetTolerance(1e-3);
   minimizer.SetExtraOption("gtol",1e-3);
   ROOT::Math::GradFunctor f(&RosenBrock,&RosenBrockGrad,2); 
   double step[2] = {0.01,0.01};
   double variable[2] = { -1.2,1.0};
 
   minimizer.SetFunction(f);
   minimizer.SetHessianFunction(RosenBrockHessian);
   //minimizer.SetFunction(rgf);
   //minimizer.SetFunction(wgf);
 
   // Set the free variables to be minimized!
   minimizer.SetVariable(0,"x",variable[0], step[0]);
   minimizer.SetVariable(1,"y",variable[1], step[1]);
 
   minimizer.Minimize(); 
 
   const double *xs = minimizer.X();
   cout << "Minimum: f(" << xs[0] << "," << xs[1] << "): " 
        << RosenBrock(xs) << endl;
   cout << endl << "===============" << endl;
   }
   return 0;
}

int main()
{
  return scipy();
}