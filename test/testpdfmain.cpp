#include "TH1Pdf.hh"
#include <TH1D.h>
#include <TH2D.h>
#include <TF2.h>

int main()
{
  TH1D *h=new TH1D("h","h",20,-50,10);
  h->FillRandom("gaus", 10000);
  TH1Pdf hpdf(h);

  std::cout << hpdf.Integrate(-1000,10000) << std::endl;
  std::cout << hpdf.Integrate(-0.5,0.3) << std::endl;

  
  TH2D *h2=new TH2D("h2","h2",20,-10,10,20,-10,10);
  TF2 *xyg = new TF2("xyg","xygaus",-10,10,-10,10);
  xyg->SetParameters(1,0,2,0,2); //amplitude, meanx,sigmax,meany,sigmay
  h2->FillRandom("xyg",25030);
  TH1Pdf hpdf2(h2);

  std::cout << hpdf2.Integrate(-1,1,.49,.5) << std::endl;
  std::cout << hpdf2.Integrate(-10,10,-10,20) << std::endl;
  
  return 0;
}
