TFile *rootfile=0;

TH1* GetHist(TString filename, TString histname)
{
  cout << "Getting " << histname << " from " << filename << endl;
  rootfile=new TFile(filename);
  if (rootfile->IsZombie()) {
    cout << "Error: Cannot open input file " << filename << endl;
    return 0;
  }
  rootfile->cd();
  TH1* hist=dynamic_cast<TH1*>(rootfile->Get(histname));
  if(!hist) {
    cout << "Error: Cannot get histogram " << histname << endl;
    return 0;
  }
  return hist;
}

void DrawHist1(void)
{
  // Compare the four corners against the generated ones
  
  TH2D* h00=dynamic_cast<TH2D*>(GetHist("../data/mass2D_Phi_450_omega_0p5.root","hist_sum_1"));
  TH2D* h10=dynamic_cast<TH2D*>(GetHist("../data/mass2D_Phi_450_omega_2p1.root","hist_sum_1"));
  TH2D* h01=dynamic_cast<TH2D*>(GetHist("../data/mass2D_Phi_1050_omega_0p5.root","hist_sum_1"));
  TH2D* h11=dynamic_cast<TH2D*>(GetHist("../data/mass2D_Phi_1050_omega_2p1.root","hist_sum_1"));

  TH2F* hGen00=dynamic_cast<TH2F*>(GetHist("../data/hists.root","h0_0"));
  TH2F* hGen10=dynamic_cast<TH2F*>(GetHist("../data/hists.root","h8_0"));
  TH2F* hGen01=dynamic_cast<TH2F*>(GetHist("../data/hists.root","h0_6"));
  TH2F* hGen11=dynamic_cast<TH2F*>(GetHist("../data/hists.root","h8_6"));

  for(int i=0; i<4; i++) {
  
    TCanvas *c=new TCanvas();
    c->Divide(2,1);

    TH2D* h1;
    TH2F* h2;
    if(i==0) { h1=h00; h2=hGen00; }
    if(i==1) { h1=h01; h2=hGen01; }
    if(i==2) { h1=h10; h2=hGen10; }
    if(i==3) { h1=h11; h2=hGen11; }
    
    TH1* hx1=h1->ProjectionX();
    TH1* hx2=h2->ProjectionX();
    TH1* hy1=h1->ProjectionY();
    TH1* hy2=h2->ProjectionY();
    
    hx1->SetLineColor(kCyan+2);
    hx1->SetLineStyle(1);
    hx1->SetLineWidth(2);
    hx1->SetFillColor(kCyan+2);
    hx1->SetFillStyle(1001);
    hx2->SetLineColor(kBlack);
    hx2->SetLineStyle(1);
    hx2->SetLineWidth(2);
    
    hy1->SetLineColor(kCyan+2);
    hy1->SetLineStyle(1);
    hy1->SetLineWidth(2);
    hy1->SetFillColor(kCyan+2);
    hy1->SetFillStyle(1001);
    hy2->SetLineColor(kBlack);
    hy2->SetLineStyle(1);
    hy2->SetLineWidth(2);

    hx1->SetTitle("");
    hx1->GetXaxis()->SetTitle("M_{2p} [GeV]");
    c->cd(1);
    hx1->DrawNormalized();
    hx2->DrawNormalized("same");
    
    TLegend *leg=new TLegend(0.47,0.63,1.0,0.75);
    if(i==0) leg->AddEntry(hx1, "M_{#phi}=450 GeV, M_{#omega}=0.5 GeV", "f");
    if(i==1) leg->AddEntry(hx1, "M_{#phi}=1050 GeV, M_{#omega}=0.5 GeV", "f");
    if(i==2) leg->AddEntry(hx1, "M_{#phi}=450 GeV, M_{#omega}=2.1 GeV", "f");
    if(i==3) leg->AddEntry(hx1, "M_{#phi}=1050 GeV, M_{#omega}=2.1 GeV", "f");
    leg->AddEntry(hx2, "Morphed", "l");
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->Draw();

    hy1->SetTitle("");
    hy1->GetXaxis()->SetTitle("M_{2p#gamma} [GeV]");
    c->cd(2);
    hy1->DrawNormalized();
    hy2->DrawNormalized("same");
    leg->Draw();
  }

  
  return;
}


void DrawHist2(void)
{
  TH2F* hGen[9][7];

  for(int i=0; i<9; i++)
    for(int j=0; j<7; j++) {
      ostringstream oss;
      oss << "h" << i << "_" << j;
      hGen[i][j]=dynamic_cast<TH2F*>(GetHist("../data/hists.root",oss.str().c_str()));
    }

  int colors[9] = {kBlue+4, kBlue+3, kBlue+2, kBlue+1, kBlue, kBlue-4, kBlue-7, kBlue-9, kBlue-10};

  for(int i=0; i<=8; i++) {
    TCanvas *c=new TCanvas();
    TH1* h[7];
    TLegend *leg=new TLegend(0.58,0.65,0.93,0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);

    for(int j=0; j<=6; j++) {
      h[j]=hGen[i][j]->ProjectionY();
    
      h[j]->SetLineWidth(3);
      h[j]->SetLineColor(colors[j]);
      h[j]->SetLineStyle(j+1);
      h[j]->SetMaximum(0.6);
      h[j]->GetYaxis()->SetRangeUser(0,0.6);
      
      if(j==0){
	h[j]->SetStats(0);
	h[j]->GetXaxis()->SetRangeUser(300,1200);
	h[j]->SetTitle("");
	h[j]->GetYaxis()->SetTitle("Events/0.1 GeV");
	h[j]->GetXaxis()->SetTitle("M_{2p#gamma} [GeV]");
	h[j]->DrawNormalized();
      }
      else
	h[j]->DrawNormalized("same");
      ostringstream oss;
      
      oss << "M_{#phi}=" << (450+j*100) <<" GeV, M_{#omega}=" << (0.5+i*0.2) << " GeV";
      leg->AddEntry(h[j],oss.str().c_str(),"l");
    }
    leg->Draw();
    ostringstream fn;
    fn << "Mo" << 5+i*2 << ".pdf";
    c->SaveAs(fn.str().c_str());
  } 

  for(int j=0; j<=6; j++) {
    TCanvas *c=new TCanvas();
    TH1* h[9];
    TLegend *leg=new TLegend(0.52,0.65,0.87,0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);

    for(int i=0; i<=8; i+=2) {
      h[i]=hGen[i][j]->ProjectionX();
      h[i]->Rebin(2);

      h[i]->SetLineWidth(3);
      h[i]->SetLineColor(colors[i/2]);
      h[i]->SetLineStyle(i+1);
      h[i]->SetMaximum(0.4);
      h[i]->GetYaxis()->SetRangeUser(0,0.4);
      
      if(i==0){
	h[i]->SetStats(0);
	h[i]->GetXaxis()->SetRangeUser(0,3);
	h[i]->SetTitle("");
	h[i]->GetYaxis()->SetTitle("Events/0.1 GeV");
	h[i]->GetXaxis()->SetTitle("M_{2p} [GeV]");
	h[i]->DrawNormalized();
      }
      else
	h[i]->DrawNormalized("same");
      ostringstream oss;
      
      oss << "M_{#phi}=" << (450+j*100) <<" GeV, M_{#omega}=" << (0.5+i*0.2) << " GeV";
      leg->AddEntry(h[i],oss.str().c_str(),"l");
    }
    leg->Draw();
    ostringstream fn;
    fn << "Mp" << (450+j*100) << ".pdf";
    c->SaveAs(fn.str().c_str());
  }
  
}

void DrawHist(void)
{
  DrawHist1();
  //  DrawHist2();
  return;
}
