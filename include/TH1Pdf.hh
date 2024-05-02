#ifndef __THPDF_HH__
#define __THPDF_HH__

#include <TH1.h>


class TH1Pdf
{
public:
  TH1Pdf(const TH1* hist);
  TH1Pdf(const TH1Pdf&) = delete;

  virtual ~TH1Pdf();

  Double_t Integrate(Double_t x1, Double_t x2, Double_t y1=0, Double_t y2=0, Double_t z1=0, Double_t z2=0) const;

  Int_t    getNDim(void) const { return ndim_; }
  Double_t getMeanX(void) const { return xmu_; }
  Double_t getMeanY(void) const { return ymu_; }
  Double_t getMeanZ(void) const { return zmu_; }
  Double_t getSigmaX(void) const { return xsigma_; }
  Double_t getSigmaY(void) const { return ysigma_; }
  Double_t getSigmaZ(void) const { return zsigma_; }

 protected:
  void ComputePdf(const TH1* hist);

  Int_t GetBin(Int_t binx, Int_t biny=0, Int_t binz=0) const {
    return (binx + (xaxis_.GetNbins()+2)*(biny + (yaxis_.GetNbins()+2)*binz));
  }

  
 private:
  TArrayD pdf_;
  TAxis xaxis_, yaxis_, zaxis_;
  int ndim_;
  double xmu_, xsigma_;
  double ymu_, ysigma_;
  double zmu_, zsigma_;


};

void Interpolate1D(TH1D* hist, const TH1Pdf &pdf1, const TH1Pdf &pdf2, double mu1, double mu2, double muint);




#endif

