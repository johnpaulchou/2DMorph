#include "TH1Pdf.hh"
#include <iostream>

TH1Pdf::TH1Pdf(const TH1* hist)
{
  ComputePdf(hist);
}

TH1Pdf::~TH1Pdf() {}

void TH1Pdf::ComputePdf(const TH1* hist)
{
  //copy over the axes
  xaxis_=*hist->GetXaxis();
  yaxis_=*hist->GetYaxis();
  zaxis_=*hist->GetZaxis();

  // get the # of dimensions and set the pdf array size
  ndim_=hist->GetDimension();

  Int_t binx1=1; Int_t binx2=xaxis_.GetNbins();
  Int_t biny1=0; Int_t biny2=0;
  Int_t binz1=0; Int_t binz2=0;
  Int_t arraysize=binx2+2;
  if(ndim_>1) {
    biny1=1;
    biny2=yaxis_.GetNbins();
    arraysize*=(biny2+2);
  }
  if(ndim_>2) {
    binz1=1;
    binz2=zaxis_.GetNbins();
    arraysize*=(binz2+2);
  }
  pdf_.Set(arraysize);
  
  //   Loop over the bins to compute the normalization and means
  xmu_=ymu_=zmu_=0.;
  Double_t norm=0.;
  for (Int_t binx = binx1; binx <= binx2; ++binx) {
    Double_t xcenter=xaxis_.GetBinCenter(binx);
    for(Int_t biny = biny1; biny <= biny2; ++biny) {
      Double_t ycenter=yaxis_.GetBinCenter(biny);
      for (Int_t binz = binz1; binz <= binz2; ++binz) {
	Double_t zcenter=zaxis_.GetBinCenter(binz);
	
	Double_t val=hist->GetBinContent(binx, biny, binz);
	xmu_+=val*xcenter;
	ymu_+=val*ycenter;
	zmu_+=val*zcenter;
	norm+=val;
      }
    }
  }
  xmu_=xmu_/norm;
  if(ndim_>1)
    ymu_=ymu_/norm;
  else
    ymu_=0.0;
  if(ndim_>2)
    zmu_=zmu_/norm;
  else
    zmu_=0.0;
  //  std::cout << "xmu=" << xmu_ << "; ymu=" << ymu_ << "; zmu=" << zmu_ << "; norm=" << norm << std::endl;


  xsigma_=ysigma_=zsigma_=0.;
  //   Loop over the bins to standard deviations and normalized bin contents
  for (Int_t binx = binx1; binx <= binx2; ++binx) {
    Double_t xcenter=xaxis_.GetBinCenter(binx);
    Double_t dx = xaxis_.GetBinWidth(binx);
    for(Int_t biny = biny1; biny <= biny2; ++biny) {
      Double_t ycenter=yaxis_.GetBinCenter(biny);
      Double_t dy = yaxis_.GetBinWidth(biny);
      for (Int_t binz = binz1; binz <= binz2; ++binz) {
	Double_t zcenter=zaxis_.GetBinCenter(binz);
	Double_t dz = zaxis_.GetBinWidth(binz);

	Double_t dv=dx*dy*dz;
	Double_t val=hist->GetBinContent(binx, biny, binz);
	xsigma_+=(xmu_-xcenter)*(xmu_-xcenter)*val;
	ysigma_+=(ymu_-ycenter)*(ymu_-ycenter)*val;
	zsigma_+=(zmu_-zcenter)*(zmu_-zcenter)*val;
	Int_t bin=GetBin(binx,biny,binz);
	pdf_[bin]=val/norm/dv;
	//	std::cout << "binx=" << binx << "; biny=" << biny << "; binz=" << binz << "; bin=" << bin << "; gbin=" << hist->GetBin(binx, biny, binz) << "; pdf[bin]=" << pdf_[bin] << std::endl;
      }
    }
  }
  xsigma_=sqrt(xsigma_/norm);
  if(ndim_>1)
    ysigma_=sqrt(ysigma_/norm);
  else
    ysigma_=0.0;
  if(ndim_>2)
    zsigma_=sqrt(zsigma_/norm);
  else
    zsigma_=0.0;
  
  //  std::cout << "xsigma=" << xsigma_ << "; ysigma=" << ysigma_ << "; zsigma=" << zsigma_ << std::endl;
  return;
}

Double_t TH1Pdf::Integrate(Double_t x1, Double_t x2, Double_t y1, Double_t y2, Double_t z1, Double_t z2) const
{
  if(x2<x1 || (ndim_>1 && y2<y1) || (ndim_>2 && z2<z1)) {
    std::cout << "Boundaries don't make sense" << std::endl;
    return -1.0;
  }
  
  Int_t binx1, binx2, biny1, biny2, binz1, binz2;
  binx1=xaxis_.FindBin(x1); binx2=xaxis_.FindBin(x2);
  if(ndim_>1) { biny1=yaxis_.FindBin(y1); biny2=yaxis_.FindBin(y2); }
  if(ndim_>2) { binz1=zaxis_.FindBin(z1); binz2=zaxis_.FindBin(z2); }

  if(binx1<1) binx1=1;
  if(binx2>xaxis_.GetNbins()) binx2=xaxis_.GetNbins();
  if(ndim_>1 && biny1<1) biny1=1;
  if(ndim_>1 && biny2>yaxis_.GetNbins()) biny2=yaxis_.GetNbins();
  if(ndim_>2 && binz1<1) binz1=1;
  if(ndim_>2 && binz2>zaxis_.GetNbins()) binz2=zaxis_.GetNbins();
  
  double sum=0.;
  for(int binx=binx1; binx<=binx2; binx++) {
    Double_t binxlo=xaxis_.GetBinLowEdge(binx);
    Double_t binxhi=xaxis_.GetBinUpEdge(binx);
    Double_t binxwidth=binxhi-binxlo;
    if(x1>binxlo) binxlo=x1;
    if(x2<binxhi) binxhi=x2;
    double weightx=binxhi-binxlo;
    
    if(ndim_==1) sum+=weightx*pdf_[GetBin(binx)];
    else {
      for(int biny=biny1; biny<=biny2; biny++) {
	Double_t binylo=yaxis_.GetBinLowEdge(biny);
	Double_t binyhi=yaxis_.GetBinUpEdge(biny);
	Double_t binywidth=binyhi-binylo;
	if(y1>binylo) binylo=y1;
	if(y2<binyhi) binyhi=y2;
	double weighty=binyhi-binylo;

	if(ndim_==2) sum+=weightx*weighty*pdf_[GetBin(binx,biny)];
	else {
	  for(int binz=binz1; binz<=binz2; binz++) {
	    Double_t binzlo=zaxis_.GetBinLowEdge(binz);
	    Double_t binzhi=zaxis_.GetBinUpEdge(binz);
	    Double_t binzwidth=binzhi-binzlo;
	    if(z1>binzlo) binzlo=z1;
	    if(z2<binzhi) binzhi=z2;
	    double weightz=binzhi-binzlo;
      
	    sum+=weightx*weighty*weightz*pdf_[GetBin(binx,biny,binz)];
	  }
	}
      }
    }
  }
  return sum;
}


void Interpolate1D(TH1D* hist, const TH1Pdf &pdf1, const TH1Pdf &pdf2, double m1, double m2, double mint)
{
  /*  hist->Reset();

  double c1=(m1-mint)^j(M^-1)
  
  double mu1=pdf1.getMeanX();
  double sig1=pdf1.getSigmaX();
  double mu2=pdf2.getMeanX();
  double sig2=pdf2.getSigmaX();
  for(int binx=1; binx<=hist->GetNbinsX(); binx++) {

  }
  */
  

  return;
}
