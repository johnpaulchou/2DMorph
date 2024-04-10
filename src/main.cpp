#include <iostream>
#include <RooRealVar.h>
#include <RooDataHist.h>
#include <RooGaussian.h>
#include <RooProdPdf.h>
#include <RooPlot.h>
#include <RooWorkspace.h>
#include <RooMomentMorphFuncND.h>
#include <TApplication.h>
#include <RooBinning.h>
#include <RooWrapperPdf.h>

#include <TH2D.h>

using namespace RooFit;
using namespace std;

const double X_MIN=0.;
const double X_MAX=5.;
const double Y_MIN=0.;
const double Y_MAX=4000.;

const double TX_MIN=0.4;
const double TX_MAX=4.0;
const double TY_MIN=200.;
const double TY_MAX=3000.;

const double TX_STEP=0.1;
const double TY_STEP=10;

const string concat(const string& s1, const string& s2)
{
  return s1+s2;
}

RooProdPdf* create2DGauss(string name, RooRealVar& x, RooRealVar& y, double mean1, double mean2, double sigma1, double sigma2)
{
  RooRealVar* mean1v=new RooRealVar(concat(name,"mean1").c_str(), "mean of Gaussian 1", mean1);
  RooRealVar* mean2v=new RooRealVar(concat(name,"mean2").c_str(), "mean of Gaussian 2", mean2);
  RooRealVar* sigma1v=new RooRealVar(concat(name,"sigma1").c_str(), "sigma of Gaussian 1", sigma1,1e-30,1e30);
  RooRealVar* sigma2v=new RooRealVar(concat(name,"sigma2").c_str(), "sigma of Gaussian 2", sigma2,1e-30,1e30);

  // Create the Gaussian PDF in x and y
  RooGaussian* gaussX=new RooGaussian(concat(name,"gaussX").c_str(), "Gaussian in x", x, *mean1v, *sigma1v);
  RooGaussian* gaussY=new RooGaussian(concat(name,"gaussY").c_str(), "Gaussian in y", y, *mean2v, *sigma2v);
  
  // Create a 2D Gaussian PDF by multiplying the individual Gaussians
  RooProdPdf *pdf= new RooProdPdf(concat(name,"gauss2D").c_str(), "2D Gaussian PDF", RooArgList(*gaussX, *gaussY));
  return pdf;
}

int main(int argc, char* argv[])
{
  TApplication app("myApp", nullptr, nullptr);

  // observables
  RooRealVar x("x", "x", X_MIN, X_MAX);
  RooRealVar y("y", "y", Y_MIN, Y_MAX);
  
  // theory parameters
  RooRealVar tx("tx","tx",TX_MIN,TX_MAX);
  RooRealVar ty("ty","ty",TY_MIN,TY_MAX);

  // grid binning
  RooBinning bintx((TX_MAX-TX_MIN)/TX_STEP,TX_MIN,TX_MAX);
  RooBinning binty((TY_MAX-TY_MIN)/TY_STEP,TY_MIN,TY_MAX);
  RooMomentMorphFuncND::Grid2 grid(bintx,binty);

  double m_p, m_o;
  m_o=TX_MIN; m_p=TY_MIN; RooProdPdf* G1=create2DGauss("G1",x,y,m_o,m_p,m_o*0.2,m_p*0.2); grid.addPdf(*G1,bintx.binNumber(m_o),binty.binNumber(m_p));
			     
  m_o=TX_MIN; m_p=TY_MAX; RooProdPdf* G2=create2DGauss("G2",x,y,m_o,m_p,m_o*0.2,m_p*0.2); grid.addPdf(*G2,bintx.binNumber(m_o),binty.binNumber(m_p));
  m_o=TX_MAX; m_p=TY_MIN; RooProdPdf* G3=create2DGauss("G3",x,y,m_o,m_p,m_o*0.2,m_p*0.2); grid.addPdf(*G3,bintx.binNumber(m_o),binty.binNumber(m_p));
  m_o=TX_MAX; m_p=TY_MAX; RooProdPdf* G4=create2DGauss("G4",x,y,m_o,m_p,m_o*0.2,m_p*0.2); grid.addPdf(*G4,bintx.binNumber(m_o),binty.binNumber(m_p));
  
  RooMomentMorphFuncND morph("morph","morph",RooArgList(tx,ty),RooArgList(x,y),grid,RooMomentMorphFuncND::NonLinear);
  morph.setPdfMode();
  RooWrapperPdf pdf("morph_pdf","morph_pdf",morph,true);
  
  // generate at a given point
  tx.setVal(2);
  ty.setVal(300);
  TH2D* hist = new TH2D("hist", "hist", 100, X_MIN, X_MAX, 100, Y_MIN, Y_MAX);
  pdf.generateBinned(RooArgSet(x,y),10000,true)->fillHistogram(hist,RooArgList(x,y));
  //  G4->generateBinned(RooArgSet(x,y),10000,true)->fillHistogram(hist,RooArgList(x,y));
  hist->Draw();
  app.Run();
  
  return 0;
}
