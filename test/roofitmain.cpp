#include <iostream>
#include <sstream>

#include <RooRealVar.h>
#include <RooDataHist.h>
#include <RooGaussian.h>
#include <RooProdPdf.h>
#include <RooPlot.h>
#include <RooWorkspace.h>
#include <TApplication.h>
#include <RooBinning.h>
#include <RooWrapperPdf.h>
#include <RooHistPdf.h>

#include <RooMomentMorphFuncNDFix.h>

#include <TH2D.h>
#include <TCanvas.h>
#include <TVectorD.h>
#include <TFile.h>

using namespace RooFit;
using namespace std;

TFile *rootfile;

const double X_MIN=0.;
const double X_MAX=5.;
const double Y_MIN=0;
const double Y_MAX=4000.;

const double TX_MIN=0.5;
const double TX_MAX=2.1;
const double TY_MIN=450.;
const double TY_MAX=1050.;

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

RooHistPdf* create2DHistPdf(string name, string filename, string histname, RooRealVar& x, RooRealVar& y)
{
  rootfile=new TFile(filename.c_str());
  if (rootfile->IsZombie()) {
    cout << "Error: Cannot open input file " << filename << endl;
    return 0;
  }
  rootfile->cd();
  TH2* hist2D=dynamic_cast<TH2*>(rootfile->Get(histname.c_str()));
  if(!hist2D) {
    cout << "Error: Cannot get histogram " << histname << endl;
    return 0;
  }
  RooDataHist *dataHist=new RooDataHist((name+"dataHist").c_str(), "dataHist", RooArgList(x, y), hist2D);
  RooHistPdf *histPdf=new RooHistPdf((name+"histPdf").c_str(), "histPdf", RooArgSet(x, y), *dataHist);
  return histPdf;

}

int getBin(RooBinning& binning, double val)
{
  if(val>=binning.highBound()) return binning.numBins();
  else return binning.binNumber(val);
}

int main(int argc, char* argv[])
{
  //  TApplication app("myApp", nullptr, nullptr);
  TFile* newfile=new TFile("../data/hists.root","RECREATE");

  // observables
  RooRealVar x("x", "x", X_MIN, X_MAX);
  RooRealVar y("y", "y", Y_MIN, Y_MAX);
  
  // theory parameters
  RooRealVar tx("tx","tx",TX_MIN,TX_MAX);
  RooRealVar ty("ty","ty",TY_MIN,TY_MAX);

  // grid binning
  RooBinning bintx(1,TX_MIN,TX_MAX);
  RooBinning binty(1,TY_MIN,TY_MAX);
  RooMomentMorphFuncNDFix::Grid2 grid(bintx,binty);

  /*  m_o=TX_MIN; m_p=TY_MIN; RooProdPdf* G1=create2DGauss("G1",x,y,m_o,m_p,m_o*0.06,m_p*0.03); grid.addPdf(*G1,getBin(bintx,m_o),getBin(binty,m_p));
  m_o=TX_MIN; m_p=TY_MAX; RooProdPdf* G2=create2DGauss("G2",x,y,m_o,m_p,m_o*0.06,m_p*0.03); grid.addPdf(*G2,getBin(bintx,m_o),getBin(binty,m_p));
  m_o=TX_MAX; m_p=TY_MIN; RooProdPdf* G3=create2DGauss("G3",x,y,m_o,m_p,m_o*0.06,m_p*0.03); grid.addPdf(*G3,getBin(bintx,m_o),getBin(binty,m_p));
  m_o=TX_MAX; m_p=TY_MAX; RooProdPdf* G4=create2DGauss("G4",x,y,m_o,m_p,m_o*0.06,m_p*0.03); grid.addPdf(*G4,getBin(bintx,m_o),getBin(binty,m_p));*/
  RooAbsPdf* G1=create2DHistPdf("P450O0p5","../data/mass2D_Phi_450_omega_0p5.root","hist_sum_1",x,y);   grid.addPdf(*G1,0,0);
  RooAbsPdf* G2=create2DHistPdf("P450O2p1","../data/mass2D_Phi_450_omega_2p1.root","hist_sum_1",x,y);   grid.addPdf(*G2,1,0);
  RooAbsPdf* G3=create2DHistPdf("P1050O0p5","../data/mass2D_Phi_1050_omega_0p5.root","hist_sum_1",x,y); grid.addPdf(*G3,0,1);
  RooAbsPdf* G4=create2DHistPdf("P1050O2p1","../data/mass2D_Phi_1050_omega_2p1.root","hist_sum_1",x,y); grid.addPdf(*G4,1,1);

  RooMomentMorphFuncNDFix morph("morph","morph",RooArgList(tx,ty),RooArgList(x,y),grid,RooMomentMorphFuncNDFix::Linear);
  morph.setPdfMode();
  RooWrapperPdf pdf("morph_pdf","morph_pdf",morph,true);
  
  auto framex = x.frame();
  auto framey = y.frame();

  for(int i=0; i<=8; i++) {
    tx.setVal(0.5+0.2*i);
    for(int j=0; j<=6; j++) {
      ty.setVal(450.+100.*j);

      std::cout << 0.5+0.2*i << " " << 450.+100.*j << std::endl;
      
      TH2* hist=dynamic_cast<TH2*>(pdf.createHistogram("x,y",100,100));
      ostringstream oss;
      oss << "h" << i << "_" << j;
      TH2* hclone=dynamic_cast<TH2*>(hist->Clone(oss.str().c_str()));
      newfile->cd();
      hclone->Write();
    }
  }

  // app.Run();
  
  return 0;
}
