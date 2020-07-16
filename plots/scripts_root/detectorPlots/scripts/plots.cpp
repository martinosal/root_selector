void plot(std::string hist,int x2,TFile* fData=_file0)
{
  TLegend *leg;
  TH1D *h1_proj;
  TH2D* h1 = (TH2D*)fData->Get(("trk_"+hist+"_origin_B").c_str());
  string origin[5] = {"PU","B","C","FRAG","GEANT"};
  float x_d0sig=10.;
//  h1->SetTitle(hist.c_str());

  TCanvas c("c1", "canvas", 1300, 900);
  gStyle->SetOptTitle(0);

  c.Divide(2,1);
  c.cd(1);
  gPad->SetLeftMargin(gPad->GetLeftMargin()*1.3);
  gPad->SetRightMargin(gPad->GetRightMargin()*1.6);
  gStyle->SetOptStat(0);
  gStyle->SetPalette(87);
  TH2D* h2 = (TH2D*)h1->Clone();
  float N=h2->GetEntries();
  h2->Scale(1./N);
  if(x2!=0){
    h2->GetXaxis()->SetRangeUser(0., x2);
  }
  if(x2==0){
    if(hist=="z0sig" || hist=="d0sig")
      h2->GetXaxis()->SetRangeUser(-x_d0sig, x_d0sig);
  }
  h2->GetXaxis()->SetTitle((hist).c_str());
  for(int i=1;i<=5;i++){
    h2->GetYaxis()->SetBinLabel(i,origin[i-1].c_str());
  }

  gStyle->SetPaintTextFormat(".4f");
  h2->Draw("colz");

  c.cd(2);
  gPad->SetLogy();
  leg = new TLegend(0.75,0.75,0.9,0.9);
  h1_proj = new TH1D("a","b",5,-1,4);
//  h1_proj->SetTitle(hist.c_str());
  h1_proj->GetYaxis()->SetRangeUser(0.,0.5);
  std::string idx="projection";

  for(unsigned i=1;i<6;i++){

    TColor *col=gROOT->GetColor(i);
    char n=i;
    h1_proj= (TH1D*) h2->ProjectionX((idx+"_"+std::to_string(n)).c_str(),i,i);
//    h1_proj= (TH1D*) h1->ProjectionX((idx+"_"+std::to_string(n)).c_str(),i,i);
//    float N_proj = h1_proj->GetEntries();
//    h1_proj->Scale(1./N_proj);
    if(x2!=0){
      h1_proj->GetXaxis()->SetRangeUser(0., x2);
    }
    if(x2==0){
      if(hist=="z0sig" || hist=="d0sig")
        h2->GetXaxis()->SetRangeUser(-x_d0sig, x_d0sig);
    }

    if(hist=="d0sig"|| hist=="z0sinthsig" || hist=="logpTfrac"|| hist=="logDR"){
      h1_proj->GetYaxis()->SetRangeUser(1e-5,1e-1);
    }
    else
      h1_proj->GetYaxis()->SetRangeUser(1e-5,1.);

    h1_proj->GetXaxis()->SetTitle((hist).c_str());
//    h1_proj->GetYaxis()->SetRangeUser(1e-3, 1.);
    h1_proj->SetLineWidth(2);
    h1_proj->Draw("same hist PLC");
    leg->AddEntry(h1_proj,origin[i-1].c_str(),"l");

  }
//  leg->SetHeader("The Legend Title");
  leg->Draw();
  c.SaveAs(("/home/salomon/Private/atlas/FTPF/Selector/plots/scripts_root/detectorPlots/"+hist+".pdf").c_str());

  delete h1,h2,c,h1_proj,leg;
}

void plot(std::string hist,std::string alg,int x2,TFile* fData=_file0)
{
  TLegend *leg;
  TH1D *h1_proj;
  TH2D* h1 = (TH2D*)fData->Get(("trk_"+hist+"_origin_"+alg+"_B").c_str());
  string origin[5] = {"PU","B","C","FRAG","GEANT"};
  float x_d0sig=10.;
//  h1->SetTitle((hist+" "+alg).c_str());

  int max_origin=3;
  TCanvas c("c1", "canvas", 1300, 900);
  gStyle->SetOptTitle(0);
  c.Divide(2,1);
  c.cd(1);
  gPad->SetLeftMargin(gPad->GetLeftMargin()*1.3);
  gPad->SetRightMargin(gPad->GetRightMargin()*1.6);
  gStyle->SetOptStat(0);
  gStyle->SetPalette(87);
  TH2D* h2 = (TH2D*)h1->Clone();
  float N=h2->GetEntries();
  h2->Scale(1./N);
  if(x2!=0){
    h2->GetXaxis()->SetRangeUser(0., x2);
  }
  if(x2==0){
    if(hist=="z0sig" || hist=="d0sig")
      h2->GetXaxis()->SetRangeUser(-x_d0sig, x_d0sig);
  }
  for(int i=1;i<=5;i++){
    h2->GetYaxis()->SetBinLabel(i,origin[i-1].c_str());
  }
  h2->GetXaxis()->SetTitle((hist).c_str());
  gStyle->SetPaintTextFormat(".4f");
  h2->Draw("colz");

  c.cd(2);
  gPad->SetLogy();
  leg = new TLegend(0.75,0.75,0.9,0.9);
  h1_proj = new TH1D("a","b",5,-1,4);
//  h1_proj->SetTitle((hist +" "+ alg).c_str());
  h1_proj->GetYaxis()->SetRangeUser(0.,0.5);
  std::string idx="projection";
  //const char *col[4] = { "KBlue", "KRed", "KOrange", "KYellow" };
  for(unsigned i=1;i<max_origin+3;i++){

    TColor *col=gROOT->GetColor(i);
    char n=i;
    h1_proj= (TH1D*) h2->ProjectionX((idx+"_"+std::to_string(n)).c_str(),i,i);
//    h1_proj= (TH1D*) h1->ProjectionX((idx+"_"+std::to_string(n)).c_str(),i,i);
//    float N_proj = h1_proj->GetEntries();
//    h1_proj->Scale(1./N_proj);
    if(x2!=0){
      h1_proj->GetXaxis()->SetRangeUser(0., x2);
    }
    if(x2==0){
      if(hist=="z0sig" || hist=="d0sig")
        h2->GetXaxis()->SetRangeUser(-x_d0sig, x_d0sig);
    }

    if(hist=="d0sig"|| hist=="z0sinthsig" || hist=="logpTfrac"|| hist=="logDR"){
      h1_proj->GetYaxis()->SetRangeUser(1e-5,1e-1);
    }
    else
      h1_proj->GetYaxis()->SetRangeUser(1e-5,1.);

    h1_proj->GetXaxis()->SetTitle((hist).c_str());
//    h1_proj->GetYaxis()->SetRangeUser(1e-3, 1.);
    h1_proj->SetLineWidth(2);
    h1_proj->Draw("same hist PLC");
    leg->AddEntry(h1_proj,origin[i-1].c_str(),"l");

  }
//  leg->SetHeader("The Legend Title");
  leg->Draw();
  c.SaveAs(("/home/salomon/Private/atlas/FTPF/Selector/plots/scripts_root/detectorPlots/"+hist+"_"+alg+".pdf").c_str());

  delete h1,h2,c,h1_proj,leg;
}
