//#include "Math/ScipyMinimizer.h"
#include "Math/MultiNumGradFunction.h"
#include <Fit/ParameterSettings.h>
#include "Math/Functor.h"
#include <string>
#include "Math/MinimizerOptions.h"
#include "TStopwatch.h"

/// @brief  Test of Fitting convolution of Landau with Gaussian function
///

const int nbins = 1000;
double xmin = -20.;
double xmax = 80;
int verbose = 0;

TH1 * genLangaus() {

   auto h1 = new TH1D("h1","h1",nbins,xmin,xmax);

   int n = 10000000;

   for (int i = 0; i < n; i++) {
      double y = gRandom->Landau();

      y += gRandom->Gaus(0,1);

      h1->Fill(y);

   }

   h1->Draw();
   std::cout << h1->GetXaxis()->GetBinCenter( h1->GetMaximumBin()) << std::endl;
   return h1;

}

void fitLangaus(TH1 * h1, bool useBounds, bool fft = false,TString image="image.png") {

   TF1Convolution fc("landau","gaus", xmin, xmax, fft);
   if (fft) fc.SetNofPointsFFT(10000);

   TF1 f1("f1",fc,xmin,xmax,5);

   double initParams[] = {100., 1., 2., 0., 2. };

   f1.SetParameters(initParams);
   f1.FixParameter(3,0);
   f1.SetNpx(10000);
   if (useBounds) {
      f1.SetParLimits(1,-1,10);
      f1.SetParLimits(2,0,10);
      f1.SetParLimits(4,0,10);
   }
   TString optFit = "L G";
   if (verbose == 1) optFit += " V";
   if (verbose == 2) optFit += " VV";
   h1->Fit(&f1,optFit);
   TCanvas *c1 = new TCanvas();
   h1->Draw();
   c1->Print(image);
}

void testLanGaus(const char * name="Minuit2", const char * name2="", bool useBounds = true, int seed = 1111) {

   gRandom->SetSeed(seed);

   auto h1 = genLangaus();

   std::cout << "do fitting ...\n";
   TStopwatch w;
   w.Start();
   ROOT::Math::MinimizerOptions::SetDefaultMinimizer(name, name2);
   ROOT::Math::MinimizerOptions::SetDefaultStrategy(0); //hessian is 2 https://root.cern.ch/root/htmldoc/guides/users-guide/FittingHistograms.html
   ROOT::Math::MinimizerOptions::SetDefaultMaxIterations(10000);
   ROOT::Math::MinimizerOptions::SetDefaultTolerance(1e-3);
   ROOT::Math::IntegratorOneDimOptions::SetDefaultIntegrator("Gauss");
   TString filename=name;
   filename+="_";
   filename+=name2;
   filename+=".png";
   
   fitLangaus(h1, useBounds, false,filename);
   w.Print();
}


// methods that requires hessian to work "dogleg", "trust-ncg","trust-exact","trust-krylov"
using namespace std;
int ipopt_langauss()
{ 
   //std::string methods[]={"Nelder-Mead","L-BFGS-B","Powell","CG","BFGS","TNC","COBYLA","SLSQP","trust-constr","Newton-CG", "dogleg", "trust-ncg","trust-exact","trust-krylov"};
   //std::string methods[]={"Nelder-Mead","L-BFGS-B","Powell","BFGS","TNC","COBYLA","SLSQP","trust-constr","Newton-CG", "dogleg", "trust-ncg","trust-exact","trust-krylov"};
 //,"L-BFGS-B","Powell","CG","BFGS","TNC","COBYLA","SLSQP","trust-constr","Newton-CG", "dogleg", "trust-ncg","trust-exact","trust-krylov"};
   //testLanGaus();
// ma27 (load the Harwell routine MA27 from library at runtime)
//         ma57 (load the Harwell routine MA57 from library at runtime)
//         ma77 (load the Harwell routine HSL_MA77 from library at runtime)
//         ma86 (load the Harwell routine MA86 from library at runtime)
//         ma97 (load the Harwell routine MA97 from library at runtime)
//         pardiso (load the Pardiso package from pardiso-project.org from user-provided library at runtime)
//         wsmp (use the Wsmp package)
//         mumps (use the Mumps package)



//    checking for MC19... yes
// checking for MA27... yes
// checking for MA28... yes
// checking for MA57... yes
// checking for MA77... yes
// checking for MA86... yes
// checking for MA97... yes
// checking for MC68... yes
   //testLanGaus("Ipopt","spral");
   testLanGaus("Ipopt","pardisomkl");
   
   //testLanGaus("Ipopt","ma27");
   //std::string methods[]={"mumps","ma57","ma77","ma86","ma97"};
   // std::string methods[]={"mumps","ma27","ma57","ma77","ma86","ma97"};
   // //TStopwatch t;
   //  for(const std::string &text : methods)
   //  {
   //    std::cout<<"=== Ipopt  ===  "<<text.c_str()<<" ==="<<std::endl;
   //     testLanGaus("Ipopt",text.c_str());

   //  }
   return 0;
}

int main()
{
  return ipopt_langauss();
}