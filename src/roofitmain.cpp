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
#include <TCanvas.h>
#include <TVectorD.h>

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

int getBin(RooBinning& binning, double val)
{
  if(val>=binning.highBound()) return binning.numBins();
  else return binning.binNumber(val);
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
  RooBinning bintx(1,TX_MIN,TX_MAX);
  RooBinning binty(1,TY_MIN,TY_MAX);
  RooMomentMorphFuncND::Grid2 grid(bintx,binty);

  double m_p, m_o;
  m_o=TX_MIN; m_p=TY_MIN; RooProdPdf* G1=create2DGauss("G1",x,y,m_o,m_p,m_o*0.06,m_p*0.03); grid.addPdf(*G1,getBin(bintx,m_o),getBin(binty,m_p));
  m_o=TX_MIN; m_p=TY_MAX; RooProdPdf* G2=create2DGauss("G2",x,y,m_o,m_p,m_o*0.06,m_p*0.03); grid.addPdf(*G2,getBin(bintx,m_o),getBin(binty,m_p));
  m_o=TX_MAX; m_p=TY_MIN; RooProdPdf* G3=create2DGauss("G3",x,y,m_o,m_p,m_o*0.06,m_p*0.03); grid.addPdf(*G3,getBin(bintx,m_o),getBin(binty,m_p));
  m_o=TX_MAX; m_p=TY_MAX; RooProdPdf* G4=create2DGauss("G4",x,y,m_o,m_p,m_o*0.06,m_p*0.03); grid.addPdf(*G4,getBin(bintx,m_o),getBin(binty,m_p));

  RooMomentMorphFuncND morph("morph","morph",RooArgList(tx,ty),RooArgList(x,y),grid,RooMomentMorphFuncND::Linear);
  morph.setPdfMode();
  RooWrapperPdf pdf("morph_pdf","morph_pdf",morph,true);
  
  auto framex = x.frame();
  auto framey = y.frame();

    // generate at a given point
  tx.setVal(2.9);
  ty.setVal(489);
  pdf.plotOn(framex, RooFit::LineColor(kBlue));
  pdf.plotOn(framey, RooFit::LineColor(kBlue));

  // change the point's value
  tx.setVal(3.4);
  ty.setVal(725);
  pdf.plotOn(framex, RooFit::LineColor(kRed));
  pdf.plotOn(framey, RooFit::LineColor(kRed));
  TCanvas c;
  c.Divide(2,1);
  c.cd(1);
  framex->Draw();
  c.cd(2);
  framey->Draw();
  framey->Draw("same");

  app.Run();
  
  return 0;
}
