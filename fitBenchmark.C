/// \file
/// \ingroup tutorial_fit
/// \notebook
/// Demonstrate performance and usage of Minuit2 and Fumili2 for monodimensional fits.
/// version2 : fit without increasing statistics, always a new histogram
/// \macro_image
/// \macro_output
/// \macro_code
///
/// \author Lorenzo Moneta

#include "TH1.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TStopwatch.h"
#include "TSystem.h"
#include "TRandom3.h"
#include "Math/MinimizerOptions.h"
#include "TPaveLabel.h"
#include "TStyle.h"
#include "TMath.h"
#include "TROOT.h"
#include "TFrame.h"
/*#include "Fit/FitConfig.h"*/


// TF1 *fitFcn;
// TF1 *genFcn;
// TH1 *histo;

int nbins0 = 100;
int nevts0 = 100000;
int nbins = nbins0;
int nevts = nevts0;

int strategy = 0;
int verbose = 0;

int specialPassNumber = -1;

bool useLogLFit = true;

bool useNumConv = false;

enum { NoConv = 0, NumConv = 1, FFTConv = 2 };  // convolution mode

int useConvMode = NoConv;

bool useAD = true;

enum { InitParamFixed = 1, InitParamClose = 2, InitParamFar = 3 };

int initParamMode = InitParamClose;   // select mode for initial parameters
   
bool useSameHisto = false;


bool useNewFile = true;

bool generateHistos = true;


TString fumiliMethod = "tr";


double genPars[] = {10.,1,0.5,1.,0.1, 1, 0.05};
//{1.,1.,1.,1.,0.3,1.}  // 3 poly parameters + 3 BreitWigner (Amplitude, Gamma, Mean} + 1 conv gaussian sigma

double initPars[] = { 1,1,1,100,.3, 0.5, 0.2 };


//far away values
std::vector<double> initParsMin1 = {  3000,   0.7,   0.2,   200,  .05, 0.7, 0.005};
std::vector<double> initParsMax1 = {  5000,   1.3,   0.8,  600 ,  .15 ,1.3,  0.1};

// close valuesd
std::vector<double> initParsMin2 = {  3000,   0.8,   0.4,  200,  .08, 0.9, 0.03};
std::vector<double>initParsMax2 = {  5000,   1.2,   0.6,  500,  .12 ,1.1, 0.08};

std::vector<double> initParsMin = (initParamMode == InitParamFar) ? initParsMin1 : initParsMin2;
std::vector<double> initParsMax = (initParamMode == InitParamFar) ? initParsMax1 : initParsMax2;

TParameter<int> seedValue;

// Quadratic background function
double background(double *x, double *par) {
   return par[0] * std::exp(- (par[1]*x[0] + par[2]*x[0]*x[0] ));
}

// Lorenzian Peak function
double lorentzianPeak(double *x, double *par) {
   return (0.5*par[0]*par[1]/TMath::Pi()) /
   TMath::Max( 1.e-10,(x[0]-par[2])*(x[0]-par[2]) + .25*par[1]*par[1]);
}



// Sum of background and peak function
double fitFunction(double *x, double *par) {
  return background(x,par) + lorentzianPeak(x,&par[3]);
}

std::string bkgFormula = "[0] * exp(-([1]*x + [2]*x*x))";
std::string sigFormula = "0.5*[0]*[1]/TMath::Pi()/TMath::Max(1.e-10,(x-[2])*(x-[2]) + 0.25*[1]*[1])";

std::vector<TH1*> CreateHistograms(TF1 * genFcn, int npass) {
   TRandom & rng = *gRandom;
   TStopwatch timer;
   std::vector<TH1*> histos(npass);


   TH1::AddDirectory(false);
   timer.Start();
   bool ok = true;
   // fill histogram many times
   // every time increase its statistics and re-use previous fitted
   // parameter values as starting point
      genFcn->SetParameters(genPars);
      int nh = (useSameHisto) ? 1 : npass;
      for (int pass=0;pass<nh;pass++) {
         histos[pass] = new TH1D(TString::Format("hist%d",pass),"",nbins,0,3);
         for (int i=0;i<nevts;i++) {
            histos[pass]->Fill(genFcn->GetRandom(&rng));
         }
      }
      if (useSameHisto) {
          for (int pass=1;pass<npass;pass++) {
             histos[pass] = histos[0];
          }
      }

   timer.Print();
   TH1::AddDirectory(true);
   return histos;
}

std::vector<TH1*> ReadHistograms(TFile * f, int npass) {
    std::vector<TH1*> histos;
    histos.reserve(npass);
    for (int i = 0; i < npass; i++) {
       auto h = (TH1*) f->Get(TString::Format("hist%d",i));
       if (!h) {
          Error("ReadHIstogram","Error reading histogram for pass %d",i);
          histos.clear(); return histos;
       }
       histos.push_back(h);
    }
    return histos;
}

void  WriteHistograms(std::vector<TH1*> & histos, TString fileName) {
   auto f = TFile::Open(fileName, "RECREATE");
   int npass = histos.size();
   for (int i = 0; i < npass; i++) {
       histos[i]->Write();
   }
   f->Close();
   delete f;
   return;
}


int GenerateSeed(int s = 0) {
   TRandom3 r(s);
   return r.Integer(TMath::Limits<int>::Max());
}

void RandomizeParams(int npar, double * p) {
   if (verbose) std::cout << "Randomize initial parameter values....." << std::endl;
   for (int i = 0; i < npar; i++) {
      p[i] = gRandom->Uniform(initParsMin[i],initParsMax[i]);
      if (verbose) std::cout << "\t parameter " << i << " set to : " << p[i] << std::endl;
   }
}



std::vector<bool> DoFit(const std::vector<TH1*>& histos, TF1 * fitFcn, const char *   fitter, const char * algoName,  TVirtualPad *pad, bool usegrad = false, TString extraOpt = "") {
   printf("\n*********************************************************************************\n");
   printf("\t %s %s\n",fitter,algoName);
   printf("*********************************************************************************\n");

   // auto hchi2 = new TH1D("hchi2","chi2 distribution",100,0.9,1.5);
   // auto hnfcn = new TH1D("hnfcn","Number of calls",100,1,0);
   // auto htime = new TH1D("htime","Time distribution",100,1,0);

   TString name = fitter;
   TString name2 = algoName;
   TString fitter2 = fitter;
   if (!name2.IsNull()) fitter2 = TString::Format("%s_%s",fitter,algoName);
   if (usegrad) fitter2 += "_G";
   fitter2.ReplaceAll("-","_");
   TString tName = "t_" + fitter2;
   TString ttName = "Fit results for " + name;
   if (!name2.IsNull()) ttName += "/" + name2;
   TNtuple *ntuple = new TNtuple(tName,ttName,"status:time:nfcn:chi2:ndf");
   
   TStopwatch timer;
   //   timer.Start();

   extraOpt += " S ";

  
   ROOT::Math::MinimizerOptions::SetDefaultMinimizer(name, name2);
   if (name == "Fumili" || name == "Fumili2" || name == "GSLMultiFit") 
      ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(100);
   else
      ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(0);  // use default values
   ROOT::Math::MinimizerOptions::SetDefaultStrategy(strategy);
   if ((std::string(name) =="Minuit2" && std::string(name2) == "Fumili") || std::string(name) == "Fumili2") {
      auto & opt = ROOT::Math::MinimizerOptions::Default("Minuit2");
      opt.SetValue("FumiliMethod",fumiliMethod);
   }
   std::cout << "Minimizer options for " << name << std::endl;
   ROOT::Math::MinimizerOptions::PrintDefault(name);

   
   pad->SetGrid();
   pad->SetLogy();
   std::string title = std::string(fitter);
   if (!name2.IsNull()) title += "/" + std::string(algoName);
   title += " fit bench";
   if (usegrad) {
      extraOpt += " G ";
      title += " (G)";
   }
   if (useLogLFit)
      extraOpt += " L ";

   TString fitterType(fitter);


   int npass = histos.size();

   std::vector<bool> ret(npass);
   std::vector<bool> retChi2(npass);
   int nFailChi2Only = 0;

   if (specialPassNumber < 0 || specialPassNumber >= npass )
      specialPassNumber = npass-1;  // use last event for debugging by default
   // copy histogram of special event
   auto histo = new TH1D((const TH1D&) *histos[specialPassNumber]);
   histo->SetName(TString(histo->GetName()) + fitter);
   histo->SetTitle(title.c_str());

   timer.Start();
   bool ok = true;
   // fill histogram many times
   // every time increase its statistics and re-use previous fitted
   // parameter values as starting point
  
   for (int pass=0;pass<npass;pass++) {
      if (pass%100 == 0) printf("pass : %d\n",pass);
      else printf(".");
      // reset also parameter values and error (errors are crucial)
      if (initParamMode != InitParamFixed) RandomizeParams( fitFcn->GetNpar(), initPars);
      fitFcn->SetParameters(initPars);
      for (int i = 0; i < fitFcn->GetNpar(); i++)
         fitFcn->SetParError(i,0.);
      int iret = 0;
      //fitFcn->SetParLimits(0,0,1000);
      fitFcn->SetParLimits(1,0,10);
      fitFcn->SetParLimits(2,0,10);
      //fitFcn->SetParLimits(4,0,1000);
      fitFcn->SetParLimits(5,0,3);
      fitFcn->SetParLimits(6,0,1000);
      TString opt = extraOpt;
      TH1 * h = nullptr;
      if (pass != specialPassNumber) {
         opt  +=   " 0 ";
         if (verbose == 0) opt += "Q";
         if (verbose == 2) opt += "V";
         h = histos[pass];
      } else {
         // last fit or special number do last fit 
         // 
         if (verbose > 0) opt += "V";
         h = histo;
      }
      TStopwatch timerFit; timerFit.Start();
      auto res = h->Fit(fitFcn,opt);
      timerFit.Stop();
      iret = res->Status();
      ntuple->Fill(iret,timerFit.RealTime(),res->NCalls(), res->Chi2(), res->Ndf());
      
      std::cout << "pass " << pass << " fit status  = " << iret << " chi2/ndf " << res->Chi2()/res->Ndf() << " time " << timerFit.RealTime() << " nfcn " << res->NCalls() << std::endl;
      ok &= (iret == 0);
      ret[pass] = (iret!=0);
      double probFit  = fitFcn->GetProb();
      if (iret!=0) Error("DoFit","Fit pass %d failed - iret = %d  chi2prob = %f !",pass,iret,probFit);
      if (probFit < 1.E-9) {
         retChi2[pass] = true;
         if (!ret[pass]) {
            std::cout << "pass " << pass << " - fit is not good (chi2 prob is " << probFit << " ) but status is 0" << std::endl;
            nFailChi2Only++;
         }
      }

      
   }
  

   timer.Stop();

   if (specialPassNumber == npass-1) specialPassNumber = -1;


   auto fitFunction = histo->GetFunction("fitFcn");
   if (fitFunction) fitFunction->SetLineColor(kRed+3);
   gPad->SetFillColor(kYellow-10);

   int nfail = std::accumulate(ret.begin(), ret.end(), 0);

   double cputime = timer.CpuTime();
   printf("%s, npass=%d  : RT=%7.3f s, Cpu=%7.3f s\n",fitter,npass,timer.RealTime(),cputime);
   const char * wtgrad = (usegrad) ? "(G)" : "";
   TPaveLabel *p = new TPaveLabel(0.45,0.7,0.88,0.8,Form("%s %s CPU= %g s \n \t Nfail = %d , %d",fitter,wtgrad,cputime,nfail, nFailChi2Only),"brNDC");
   p->Draw();
   p->SetTextColor(kRed+3);
   p->SetFillColor(kYellow-8);
   std::cout << "*** FitBenchmark summary ****" << std::endl;
   std::cout << "Number of failed fits (total) " << nfail+nFailChi2Only << " flagged "
             << nfail <<   " in total fits:" << npass << std::endl;

   pad->Update();
   TString fileOutName = TString("fitBenchmark.root");
   TString fileStat = (useNewFile) ? "RECREATE" : "UPDATE";
   TFile fileOut(fileOutName,fileStat);
   if (useNewFile) {
      seedValue.Write();  // write also seed */
   }
   ntuple->Write();
   pad->Write(TString("c_") + fitter2);
   histo->Write(TString("h_") + fitter2);
   useNewFile = false;
   return ret;
 }


   int fitBenchmark(int npass=1,  const char * mName = "", const char * aName = "", int seed = 111, int verb=0, double scale = 10.) {

    nbins = scale*nbins0;
    nevts = scale*nevts0;


    // if (initParamMode == InitParamFar) {
    //    initParsMin = std::vector<double>(initParsMin1, initPar
    //    { initParsMin1, initParsMin1 +7} :
    //                                                                 { initParsMin2, initParsMin2 +7};
      
   gRandom->SetSeed(seed);
   seedValue = TParameter<int>("seed",seed);
   
   TCanvas *c1 = new TCanvas("FitBench","Fitting Demo",10,10,900,900);
   c1->SetFillColor(kYellow-9);

   // create convolution signal peak with gaussian
   TF1 * peakFcn = new TF1("peakFcn",lorentzianPeak,0,3,3);
   TF1 * gaus = new TF1("fgaus",[](double *x, double *p){ return ROOT::Math::normal_pdf(x[0],p[0]); },-10,10,1);

   TF1Convolution * sigConv = new TF1Convolution(peakFcn,gaus,0,3, useConvMode == FFTConv);
   sigConv->SetExtraRange(0.1);  
   TF1 * signalFcnNC = new TF1("signalFcnNC",sigConv,0,3,sigConv->GetNpar());
   TF1 * signalFcnVt = new TF1("signalFcnVt",[](double * x, double *p){ return p[0]*TMath::Voigt(x[0]-p[2], p[3],p[1]);} ,0,3,4);

   TF1 * signalFcn = (useNumConv) ? signalFcnNC : signalFcnVt;
   // use '=' in lambda to make lambda valid outside of scope
   auto fitFunction = [=](double *x, double *par) {
      return background(x,par) + signalFcn->EvalPar(x,&par[3]);
   };
   
   TF1 * genFcn = nullptr;
   // in case of no comvolution use formula function
   if (useConvMode == NoConv) {
      auto fbkg = new TF1("fbkg",bkgFormula.c_str(),0,3);
      auto fsig = new TF1("fsig",sigFormula.c_str(),0,3);
      genFcn =  new TF1("genFcn","fbkg+fsig",0,3);
   }
   else
      // create a TF1 from the convolution with the range from 0 to 3 and 7 parameters
      genFcn = new TF1("genFcn",fitFunction,0,3, 3+4);
   
   
   genFcn->SetNpx(10000);

   if (useAD && useConvMode == NoConv) 
      genFcn->GetFormula()->GenerateGradientPar();
   //genFcn->Print("V");
   
   gStyle->SetOptFit();
   gStyle->SetStatY(0.6);

   // to avoid slowdown of integrator plugins
   ROOT::Math::IntegratorOneDimOptions::SetDefaultIntegrator("Gauss");

   verbose = verb;
   bool usegrad = (useAD) ? true : false;   

   std::cout << "Running benchmark"; 
   if (usegrad) std::cout << " with gradients from AD, ";
   if (useConvMode == NumConv) std::cout << " with numerical convolution";
   else if (useConvMode == FFTConv) std::cout << " with FFT  convolution";
   else std::cout << " no convolution";
   std::cout << " " << std::endl;
   
   std::vector<TH1*> histos;
   TFile * f = nullptr;
   // try open a file
   TString hfileName = TString::Format("fitBenchmark_histos_s%d_n%d.root",seed,npass);
   if (!generateHistos) {
      std::cout << "Try reading histograms from file " << hfileName << std::endl;
      f = TFile::Open(hfileName);
      if (f) 
         histos = ReadHistograms(f,npass);
   }
   if (f == nullptr || histos.size() != npass) generateHistos = true;
   if (generateHistos) {
      std::cout << "generating histograms......." << std::endl;
      histos = CreateHistograms(genFcn,npass);
   }

   bool ok = true;

   histos.front()->Draw();
   //genFcn->Draw("SAME");
   gPad->Update();

   if (generateHistos) WriteHistograms(histos,hfileName);
   generateHistos = false;
   
   TString minimName = mName;
   TString algoName = aName;
   if (minimName == "plot") {
      return 0;
   }

   genFcn->SetName("fitFcn");


   std::cout << "do fitting....... " << std::endl;

   
   if (minimName != "") {
      DoFit(histos,genFcn,minimName,aName, gPad, usegrad,""); // try with improve
      //TString fName = "fitBench_" + minimName + ".root";
      //gPad->SaveAs(fName);
      return 0;
   }
      
/*
   c1->Divide(2,4);
   //with Minuit
   c1->cd(1);
//   DoFit(histos,"Minuit",gPad, false,"M"); // try with improve
   DoFit(histos,"Minuit",gPad, false,"M"); // try with improve


   //with Fumili
   c1->cd(2);
   DoFit(histos,"Fumili",gPad);

   //with Minuit2
   c1->cd(3);
   DoFit(histos,"Minuit2",gPad);

   //with Fumili2
   c1->cd(4);
   DoFit(histos,"Fumili2",gPad);

  

   /// using gradients
  //with Minuit2
   c1->cd(5);
   DoFit(histos,"Minuit2",gPad,true);

   //with Fumili2
   c1->cd(6);
   DoFit(histos,"Fumili2",gPad,true);

      //with GSL
   c1->cd(7);
   DoFit(histos,"GSLMultiFit",gPad);

   //with Minuit gradient (failures)
   //c1->cd(8);
//   DoFit(histos,"Minuit",gPad,true);

   

   c1->SaveAs("FitBench.root");
   return 0;//(ok) ? 0 : 1;

*/
   return 0;
}
