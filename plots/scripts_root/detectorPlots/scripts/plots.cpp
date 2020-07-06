void plot(std::string hist,int x2,TFile* fData=_file0)
{
  TLegend *leg;
  TH1D *h1_proj;
  TH2D* h1 = (TH2D*)fData->Get(("trk_"+hist+"_origin_inB").c_str());

  h1->SetTitle(hist.c_str());

  TCanvas c("c1", "canvas", 1300, 900);
//  gStyle->SetOptTitle(0);
  /*
  TText *t = new TText(0.45,0.97,hist.c_str());
//  t->SetTextAlign(22);
  t->SetTextFont(43);
  t->SetTextSize(20);
  t->Draw();
*/
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
  h2->GetXaxis()->SetTitle(("Number of "+hist).c_str());
  h2->GetYaxis()->SetTitle("track origin");
  gStyle->SetPaintTextFormat(".4f");
  h2->Draw("colz");

  c.cd(2);
  gPad->SetLogy();
  leg = new TLegend(0.7,0.7,0.9,0.9);
  h1_proj = new TH1D("a","b",5,-1,4);
//  h1_proj->SetTitle(hist.c_str());
  h1_proj->GetYaxis()->SetRangeUser(0.,0.5);
  std::string idx="projection";
  const char *label[5]={"PU","B","C","FRAG","GEANT"};
  //const char *col[4] = { "KBlue", "KRed", "KOrange", "KYellow" };

  for(unsigned i=1;i<6;i++){

    TColor *col=gROOT->GetColor(i);
    char n=i;
    h1_proj= (TH1D*) h1->ProjectionX((idx+"_"+std::to_string(n)).c_str(),i,i);
    float N_proj = h1_proj->GetEntries();
    h1_proj->Scale(1./N_proj);
    if(x2!=0){
      h1_proj->GetXaxis()->SetRangeUser(0., x2);
    }
    h1_proj->GetXaxis()->SetTitle(("Number of "+hist).c_str());
    h1_proj->GetYaxis()->SetRangeUser(1e-3, 1.);
    h1_proj->SetLineWidth(2);
    h1_proj->Draw("same hist PLC");
    leg->AddEntry(h1_proj,label[i-1],"l");

  }
//  leg->SetHeader("The Legend Title");
  leg->Draw();
  c.SaveAs(("../plots/scripts_root/detectorPlots/"+hist+".pdf").c_str());

  delete h1,h2,c,h1_proj,leg;
}

void plot(std::string hist,std::string alg,int x2,TFile* fData=_file0)
{
  TLegend *leg;
  TH1D *h1_proj;
  TH2D* h1 = (TH2D*)fData->Get(("trk_"+hist+"_origin_"+alg+"_inB").c_str());

  h1->SetTitle((hist+" "+alg).c_str());

  int max_origin=3;
  TCanvas c("c1", "canvas", 1300, 900);
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
  h2->GetXaxis()->SetTitle(("Number of "+hist).c_str());
  h2->GetYaxis()->SetTitle("track origin");
//  gStyle->SetPaintTextFormat(".4f");
  h2->Draw("colz");

  c.cd(2);
  gPad->SetLogy();
  leg = new TLegend(0.7,0.7,0.9,0.9);
  h1_proj = new TH1D("a","b",5,-1,4);
//  h1_proj->SetTitle((hist +" "+ alg).c_str());
  h1_proj->GetYaxis()->SetRangeUser(0.,0.5);
  std::string idx="projection";
  const char *label[5]={"PU","B","C","FRAG","GEANT"};
  //const char *col[4] = { "KBlue", "KRed", "KOrange", "KYellow" };
  for(unsigned i=1;i<max_origin+3;i++){

    TColor *col=gROOT->GetColor(i);
    char n=i;
    h1_proj= (TH1D*) h1->ProjectionX((idx+"_"+std::to_string(n)).c_str(),i,i);
    float N_proj = h1_proj->GetEntries();
    h1_proj->Scale(1./N_proj);
    if(x2!=0){
      h1_proj->GetXaxis()->SetRangeUser(0., x2);
    }
    h1_proj->GetXaxis()->SetTitle(("Number of "+hist).c_str());
    h1_proj->GetYaxis()->SetRangeUser(1e-3, 1.);
    h1_proj->SetLineWidth(2);
    h1_proj->Draw("same hist PLC");
    leg->AddEntry(h1_proj,label[i-1],"l");

  }
//  leg->SetHeader("The Legend Title");
  leg->Draw();
  c.SaveAs(("../plots/scripts_root/detectorPlots/"+hist+"_"+alg+".pdf").c_str());

  delete h1,h2,c,h1_proj,leg;
}
