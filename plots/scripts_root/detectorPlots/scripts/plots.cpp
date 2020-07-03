
void IBL_plot(std::string hist,int x2,TFile* fData=_file0)
{
  TLegend *leg;
  TH1D *h1_proj;
  TH2D* h1 = (TH2D*)fData->Get((hist+"_inB").c_str());

//  const char *title="IBL Hits";
  std::string title="IBL Hits";

  int max_origin=3;
  TCanvas c("c1", "canvas", 1300, 900);
  c.Divide(2,1);
  c.cd(1);
  gStyle->SetOptStat(0);
  gStyle->SetPalette(87);
  h1->SetTitle((title).c_str());
  float N=h1->GetEntries();
  h1->Scale(1./N);
  h1->GetXaxis()->SetRangeUser(0,x2);
  h1->GetXaxis()->SetTitle("Number of IBL hits");
  h1->GetYaxis()->SetTitle("track origin");
  gStyle->SetPaintTextFormat(".4f");
  h1->Draw("colz text");

  c.cd(2);
  gPad->SetLogy();
  leg = new TLegend(0.7,0.7,0.9,0.9);
  h1_proj = new TH1D("a","b",5,-1,4);
  h1_proj->GetYaxis()->SetRangeUser(0.,0.5);
  std::string idx="projection";
  const char *label[5]={"PU","B","C","FRAG","GEANT"};
  //const char *col[4] = { "KBlue", "KRed", "KOrange", "KYellow" };
  for(unsigned i=1;i<max_origin+3;i++){

    TColor *col=gROOT->GetColor(i);
    char n=i;
    h1_proj= (TH1D*) h1->ProjectionX((idx+"_"+std::to_string(n)).c_str(),i,i);
    h1_proj->GetYaxis()->SetRangeUser(1e-4, 1.);
    h1_proj->SetLineWidth(2);
    h1_proj->Draw("same hist PLC");
    leg->AddEntry(h1_proj,label[i-1],"l");

  }
//  leg->SetHeader("The Legend Title");
  leg->Draw();
  c.SaveAs("../plots/scripts_root/detectorPlots/IBL.pdf");

  delete h1,c,h1_proj,leg;
}

void IBL_plot(std::string hist,std::string alg,int x2,TFile* fData=_file0)
{
  TLegend *leg;
  TH1D *h1_proj;
  TH2D* h1 = (TH2D*)fData->Get((hist+"_"+alg+"_inB").c_str());

//  const char *title="IBL Hits";
  std::string title="IBL Hits";

  int max_origin=3;
  TCanvas c("c1", "canvas", 1300, 900);
  c.Divide(2,1);
  c.cd(1);
  gStyle->SetOptStat(0);
  gStyle->SetPalette(87);
  h1->SetTitle((title+" "+alg).c_str());
  float N=h1->GetEntries();
  h1->Scale(1./N);
  h1->GetXaxis()->SetRangeUser(0,x2);
  h1->GetXaxis()->SetTitle("Number of IBL hits");
  h1->GetYaxis()->SetTitle("track origin");
  gStyle->SetPaintTextFormat(".4f");
  h1->Draw("colz text");

  c.cd(2);
  gPad->SetLogy();
  leg = new TLegend(0.7,0.7,0.9,0.9);
  h1_proj = new TH1D("a","b",5,-1,4);
  h1_proj->GetYaxis()->SetRangeUser(0.,0.5);
  std::string idx="projection";
  const char *label[5]={"PU","B","C","FRAG","GEANT"};
  //const char *col[4] = { "KBlue", "KRed", "KOrange", "KYellow" };
  for(unsigned i=1;i<max_origin+3;i++){

    TColor *col=gROOT->GetColor(i);
    char n=i;
    h1_proj= (TH1D*) h1->ProjectionX((idx+"_"+std::to_string(n)).c_str(),i,i);
    h1_proj->GetYaxis()->SetRangeUser(1e-4, 1.);
    h1_proj->SetLineWidth(2);
    h1_proj->Draw("same hist PLC");
    leg->AddEntry(h1_proj,label[i-1],"l");

  }
//  leg->SetHeader("The Legend Title");
  leg->Draw();
  c.SaveAs(("../plots/scripts_root/detectorPlots/IBL"+alg+".pdf").c_str());

  delete h1,c,h1_proj,leg;
}

void nextToIBL_plot(std::string hist,int x2,TFile* fData=_file0)
{
  TLegend *leg;
  TH1D *h1_proj;
  TH2D* h1 = (TH2D*)fData->Get((hist+"_inB").c_str());

//  const char *title="IBL Hits";
  std::string title="nextToIBL";

  int max_origin=3;
  TCanvas c("c1", "canvas", 1300, 900);
  c.Divide(2,1);
  c.cd(1);
  gStyle->SetOptStat(0);
  gStyle->SetPalette(87);
  h1->SetTitle((title).c_str());
  float N=h1->GetEntries();
  h1->Scale(1./N);
  h1->GetXaxis()->SetRangeUser(0,x2);
  h1->GetXaxis()->SetTitle("Number of nextToIBL hits");
  h1->GetYaxis()->SetTitle("track origin");
  gStyle->SetPaintTextFormat(".4f");
  h1->Draw("colz text");

  c.cd(2);
  gPad->SetLogy();
  leg = new TLegend(0.7,0.7,0.9,0.9);
  h1_proj = new TH1D("a","b",5,-1,4);
  h1_proj->GetYaxis()->SetRangeUser(0.,0.5);
  std::string idx="projection";
  const char *label[5]={"PU","B","C","FRAG","GEANT"};
  //const char *col[4] = { "KBlue", "KRed", "KOrange", "KYellow" };
  for(unsigned i=1;i<max_origin+3;i++){

    TColor *col=gROOT->GetColor(i);
    char n=i;
    h1_proj= (TH1D*) h1->ProjectionX((idx+"_"+std::to_string(n)).c_str(),i,i);
    h1_proj->GetYaxis()->SetRangeUser(1e-4, 1.);
    h1_proj->SetLineWidth(2);
    h1_proj->Draw("same hist PLC");
    leg->AddEntry(h1_proj,label[i-1],"l");

  }
//  leg->SetHeader("The Legend Title");
  leg->Draw();
  c.SaveAs("../plots/scripts_root/detectorPlots/nextToIBL.pdf");

  delete h1,c,h1_proj,leg;
}

void nextToIBL_plot(std::string hist,std::string alg,int x2,TFile* fData=_file0)
{
  TLegend *leg;
  TH1D *h1_proj;
  TH2D* h1 = (TH2D*)fData->Get((hist+"_"+alg+"_inB").c_str());

//  const char *title="IBL Hits";
  std::string title="nextToIBL";

  int max_origin=3;
  TCanvas c("c1", "canvas", 1300, 900);
  c.Divide(2,1);
  c.cd(1);
  gStyle->SetOptStat(0);
  gStyle->SetPalette(87);
  h1->SetTitle((title+" "+alg).c_str());
  float N=h1->GetEntries();
  h1->Scale(1./N);
  h1->GetXaxis()->SetRangeUser(0,x2);
  h1->GetXaxis()->SetTitle("Number of nextToIBL hits");
  h1->GetYaxis()->SetTitle("track origin");
  gStyle->SetPaintTextFormat(".4f");
  h1->Draw("colz text");

  c.cd(2);
  gPad->SetLogy();
  leg = new TLegend(0.7,0.7,0.9,0.9);
  h1_proj = new TH1D("a","b",5,-1,4);
  h1_proj->GetYaxis()->SetRangeUser(0.,0.5);
  std::string idx="projection";
  const char *label[5]={"PU","B","C","FRAG","GEANT"};
  //const char *col[4] = { "KBlue", "KRed", "KOrange", "KYellow" };
  for(unsigned i=1;i<max_origin+3;i++){

    TColor *col=gROOT->GetColor(i);
    char n=i;
    h1_proj= (TH1D*) h1->ProjectionX((idx+"_"+std::to_string(n)).c_str(),i,i);
    h1_proj->GetYaxis()->SetRangeUser(1e-4, 1.);
    h1_proj->SetLineWidth(2);
    h1_proj->Draw("same hist PLC");
    leg->AddEntry(h1_proj,label[i-1],"l");

  }
//  leg->SetHeader("The Legend Title");
  leg->Draw();
  c.SaveAs(("../plots/scripts_root/detectorPlots/nextToIBL"+alg+".pdf").c_str());

  delete h1,c,h1_proj,leg;
}

void sharedIBL_plot(std::string hist,int x2,TFile* fData=_file0)
{
  TLegend *leg;
  TH1D *h1_proj;
  TH2D* h1 = (TH2D*)fData->Get((hist+"_inB").c_str());

//  const char *title="IBL Hits";
  std::string title="sharedIBL";

  int max_origin=3;
  TCanvas c("c1", "canvas", 1300, 900);
  c.Divide(2,1);
  c.cd(1);
  gStyle->SetOptStat(0);
  gStyle->SetPalette(87);
  h1->SetTitle((title).c_str());
  float N=h1->GetEntries();
  h1->Scale(1./N);
  h1->GetXaxis()->SetRangeUser(0,x2);
  h1->GetXaxis()->SetTitle("Number of sharedIBL hits");
  h1->GetYaxis()->SetTitle("track origin");
  gStyle->SetPaintTextFormat(".4f");
  h1->Draw("colz text");

  c.cd(2);
  gPad->SetLogy();
  leg = new TLegend(0.7,0.7,0.9,0.9);
  h1_proj = new TH1D("a","b",5,-1,4);
  h1_proj->GetYaxis()->SetRangeUser(0.,0.5);
  std::string idx="projection";
  const char *label[5]={"PU","B","C","FRAG","GEANT"};
  //const char *col[4] = { "KBlue", "KRed", "KOrange", "KYellow" };
  for(unsigned i=1;i<max_origin+3;i++){

    TColor *col=gROOT->GetColor(i);
    char n=i;
    h1_proj= (TH1D*) h1->ProjectionX((idx+"_"+std::to_string(n)).c_str(),i,i);
    h1_proj->GetYaxis()->SetRangeUser(1e-4, 1.);
    h1_proj->SetLineWidth(2);
    h1_proj->Draw("same hist PLC");
    leg->AddEntry(h1_proj,label[i-1],"l");

  }
//  leg->SetHeader("The Legend Title");
  leg->Draw();
  c.SaveAs("../plots/scripts_root/detectorPlots/sharedIBL.pdf");

  delete h1,c,h1_proj,leg;
}

void sharedIBL_plot(std::string hist,std::string alg,int x2,TFile* fData=_file0)
{
  TLegend *leg;
  TH1D *h1_proj;
  TH2D* h1 = (TH2D*)fData->Get((hist+"_"+alg+"_inB").c_str());

//  const char *title="IBL Hits";
  std::string title="sharedIBL";

  int max_origin=3;
  TCanvas c("c1", "canvas", 1300, 900);
  c.Divide(2,1);
  c.cd(1);
  gStyle->SetOptStat(0);
  gStyle->SetPalette(87);
  h1->SetTitle((title+" "+alg).c_str());
  float N=h1->GetEntries();
  h1->Scale(1./N);
  h1->GetXaxis()->SetRangeUser(0,x2);
  h1->GetXaxis()->SetTitle("Number of sharedIBL hits");
  h1->GetYaxis()->SetTitle("track origin");
  gStyle->SetPaintTextFormat(".4f");
  h1->Draw("colz text");

  c.cd(2);
  gPad->SetLogy();
  leg = new TLegend(0.7,0.7,0.9,0.9);
  h1_proj = new TH1D("a","b",5,-1,4);
  h1_proj->GetYaxis()->SetRangeUser(0.,0.5);
  std::string idx="projection";
  const char *label[5]={"PU","B","C","FRAG","GEANT"};
  //const char *col[4] = { "KBlue", "KRed", "KOrange", "KYellow" };
  for(unsigned i=1;i<max_origin+3;i++){

    TColor *col=gROOT->GetColor(i);
    char n=i;
    h1_proj= (TH1D*) h1->ProjectionX((idx+"_"+std::to_string(n)).c_str(),i,i);
    h1_proj->GetYaxis()->SetRangeUser(1e-4, 1.);
    h1_proj->SetLineWidth(2);
    h1_proj->Draw("same hist PLC");
    leg->AddEntry(h1_proj,label[i-1],"l");

  }
//  leg->SetHeader("The Legend Title");
  leg->Draw();
  c.SaveAs(("../plots/scripts_root/detectorPlots/sharedIBL"+alg+".pdf").c_str());

  delete h1,c,h1_proj,leg;
}

void splitIBLhits(std::string hist,int x2,TFile* fData=_file0)
{
  TLegend *leg;
  TH1D *h1_proj;
  TH2D* h1 = (TH2D*)fData->Get((hist+"_inB").c_str());

//  const char *title="IBL Hits";
  std::string title="splitIBLhits";

  int max_origin=3;
  TCanvas c("c1", "canvas", 1300, 900);
  c.Divide(2,1);
  c.cd(1);
  gStyle->SetOptStat(0);
  gStyle->SetPalette(87);
  h1->SetTitle((title).c_str());
  float N=h1->GetEntries();
  h1->Scale(1./N);
  h1->GetXaxis()->SetRangeUser(0,x2);
  h1->GetXaxis()->SetTitle("Number of splitIBLhits hits");
  h1->GetYaxis()->SetTitle("track origin");
  gStyle->SetPaintTextFormat(".4f");
  h1->Draw("colz text");

  c.cd(2);
  gPad->SetLogy();
  leg = new TLegend(0.7,0.7,0.9,0.9);
  h1_proj = new TH1D("a","b",5,-1,4);
  h1_proj->GetYaxis()->SetRangeUser(0.,0.5);
  std::string idx="projection";
  const char *label[5]={"PU","B","C","FRAG","GEANT"};
  //const char *col[4] = { "KBlue", "KRed", "KOrange", "KYellow" };
  for(unsigned i=1;i<max_origin+3;i++){

    TColor *col=gROOT->GetColor(i);
    char n=i;
    h1_proj= (TH1D*) h1->ProjectionX((idx+"_"+std::to_string(n)).c_str(),i,i);
    h1_proj->GetYaxis()->SetRangeUser(1e-4, 1.);
    h1_proj->SetLineWidth(2);
    h1_proj->Draw("same hist PLC");
    leg->AddEntry(h1_proj,label[i-1],"l");

  }
//  leg->SetHeader("The Legend Title");
  leg->Draw();
  c.SaveAs("../plots/scripts_root/detectorPlots/splitIBLhits.pdf");

  delete h1,c,h1_proj,leg;
}

void splitIBLhits(std::string hist,std::string alg,int x2,TFile* fData=_file0)
{
  TLegend *leg;
  TH1D *h1_proj;
  TH2D* h1 = (TH2D*)fData->Get((hist+"_"+alg+"_inB").c_str());

//  const char *title="IBL Hits";
  std::string title="splitIBLhits";

  int max_origin=3;
  TCanvas c("c1", "canvas", 1300, 900);
  c.Divide(2,1);
  c.cd(1);
  gStyle->SetOptStat(0);
  gStyle->SetPalette(87);
  h1->SetTitle((title+" "+alg).c_str());
  float N=h1->GetEntries();
  h1->Scale(1./N);
  h1->GetXaxis()->SetRangeUser(0,x2);
  h1->GetXaxis()->SetTitle("Number of splitIBLhits hits");
  h1->GetYaxis()->SetTitle("track origin");
  gStyle->SetPaintTextFormat(".4f");
  h1->Draw("colz text");

  c.cd(2);
  gPad->SetLogy();
  leg = new TLegend(0.7,0.7,0.9,0.9);
  h1_proj = new TH1D("a","b",5,-1,4);
  h1_proj->GetYaxis()->SetRangeUser(0.,0.5);
  std::string idx="projection";
  const char *label[5]={"PU","B","C","FRAG","GEANT"};
  //const char *col[4] = { "KBlue", "KRed", "KOrange", "KYellow" };
  for(unsigned i=1;i<max_origin+3;i++){

    TColor *col=gROOT->GetColor(i);
    char n=i;
    h1_proj= (TH1D*) h1->ProjectionX((idx+"_"+std::to_string(n)).c_str(),i,i);
    h1_proj->GetYaxis()->SetRangeUser(1e-4, 1.);
    h1_proj->SetLineWidth(2);
    h1_proj->Draw("same hist PLC");
    leg->AddEntry(h1_proj,label[i-1],"l");

  }
//  leg->SetHeader("The Legend Title");
  leg->Draw();
  c.SaveAs(("../plots/scripts_root/detectorPlots/splitIBLhits"+alg+".pdf").c_str());

  delete h1,c,h1_proj,leg;
}

void nPixhits(std::string hist,int x2,TFile* fData=_file0)
{
  TLegend *leg;
  TH1D *h1_proj;
  TH2D* h1 = (TH2D*)fData->Get((hist+"_inB").c_str());

//  const char *title="IBL Hits";
  std::string title="nPixhits";

  int max_origin=3;
  TCanvas c("c1", "canvas", 1300, 900);
  c.Divide(2,1);
  c.cd(1);
  gStyle->SetOptStat(0);
  gStyle->SetPalette(87);
  h1->SetTitle((title).c_str());
  float N=h1->GetEntries();
  h1->Scale(1./N);
  h1->GetXaxis()->SetRangeUser(0,x2);
  h1->GetXaxis()->SetTitle("Number of nPixhits hits");
  h1->GetYaxis()->SetTitle("track origin");
  gStyle->SetPaintTextFormat(".4f");
  h1->Draw("colz text");

  c.cd(2);
  gPad->SetLogy();
  leg = new TLegend(0.7,0.7,0.9,0.9);
  h1_proj = new TH1D("a","b",5,-1,4);
  h1_proj->GetYaxis()->SetRangeUser(0.,0.5);
  std::string idx="projection";
  const char *label[5]={"PU","B","C","FRAG","GEANT"};
  //const char *col[4] = { "KBlue", "KRed", "KOrange", "KYellow" };
  for(unsigned i=1;i<max_origin+3;i++){

    TColor *col=gROOT->GetColor(i);
    char n=i;
    h1_proj= (TH1D*) h1->ProjectionX((idx+"_"+std::to_string(n)).c_str(),i,i);
    h1_proj->GetYaxis()->SetRangeUser(1e-4, 1.);
    h1_proj->SetLineWidth(2);
    h1_proj->Draw("same hist PLC");
    leg->AddEntry(h1_proj,label[i-1],"l");

  }
//  leg->SetHeader("The Legend Title");
  leg->Draw();
  c.SaveAs("../plots/scripts_root/detectorPlots/nPixhits.pdf");

  delete h1,c,h1_proj,leg;
}

void nPixhits(std::string hist,std::string alg,int x2,TFile* fData=_file0)
{
  TLegend *leg;
  TH1D *h1_proj;
  TH2D* h1 = (TH2D*)fData->Get((hist+"_"+alg+"_inB").c_str());

//  const char *title="IBL Hits";
  std::string title="nPixhits";

  int max_origin=3;
  TCanvas c("c1", "canvas", 1300, 900);
  c.Divide(2,1);
  c.cd(1);
  gStyle->SetOptStat(0);
  gStyle->SetPalette(87);
  h1->SetTitle((title+" "+alg).c_str());
  float N=h1->GetEntries();
  h1->Scale(1./N);
  h1->GetXaxis()->SetRangeUser(0,x2);
  h1->GetXaxis()->SetTitle("Number of nPixhits hits");
  h1->GetYaxis()->SetTitle("track origin");
  gStyle->SetPaintTextFormat(".4f");
  h1->Draw("colz text");

  c.cd(2);
  gPad->SetLogy();
  leg = new TLegend(0.7,0.7,0.9,0.9);
  h1_proj = new TH1D("a","b",5,-1,4);
  h1_proj->GetYaxis()->SetRangeUser(0.,0.5);
  std::string idx="projection";
  const char *label[5]={"PU","B","C","FRAG","GEANT"};
  //const char *col[4] = { "KBlue", "KRed", "KOrange", "KYellow" };
  for(unsigned i=1;i<max_origin+3;i++){

    TColor *col=gROOT->GetColor(i);
    char n=i;
    h1_proj= (TH1D*) h1->ProjectionX((idx+"_"+std::to_string(n)).c_str(),i,i);
    h1_proj->GetYaxis()->SetRangeUser(1e-4, 1.);
    h1_proj->SetLineWidth(2);
    h1_proj->Draw("same hist PLC");
    leg->AddEntry(h1_proj,label[i-1],"l");

  }
//  leg->SetHeader("The Legend Title");
  leg->Draw();
  c.SaveAs(("../plots/scripts_root/detectorPlots/nPixhits"+alg+".pdf").c_str());

  delete h1,c,h1_proj,leg;
}

void sharedPixhits(std::string hist,int x2,TFile* fData=_file0)
{
  TLegend *leg;
  TH1D *h1_proj;
  TH2D* h1 = (TH2D*)fData->Get((hist+"_inB").c_str());

//  const char *title="IBL Hits";
  std::string title="sharedPixhits";

  int max_origin=3;
  TCanvas c("c1", "canvas", 1300, 900);
  c.Divide(2,1);
  c.cd(1);
  gStyle->SetOptStat(0);
  gStyle->SetPalette(87);
  h1->SetTitle((title).c_str());
  float N=h1->GetEntries();
  h1->Scale(1./N);
  h1->GetXaxis()->SetRangeUser(0,x2);
  h1->GetXaxis()->SetTitle("Number of sharedPixhits hits");
  h1->GetYaxis()->SetTitle("track origin");
  gStyle->SetPaintTextFormat(".4f");
  h1->Draw("colz text");

  c.cd(2);
  gPad->SetLogy();
  leg = new TLegend(0.7,0.7,0.9,0.9);
  h1_proj = new TH1D("a","b",5,-1,4);
  h1_proj->GetYaxis()->SetRangeUser(0.,0.5);
  std::string idx="projection";
  const char *label[5]={"PU","B","C","FRAG","GEANT"};
  //const char *col[4] = { "KBlue", "KRed", "KOrange", "KYellow" };
  for(unsigned i=1;i<max_origin+3;i++){

    TColor *col=gROOT->GetColor(i);
    char n=i;
    h1_proj= (TH1D*) h1->ProjectionX((idx+"_"+std::to_string(n)).c_str(),i,i);
    h1_proj->GetYaxis()->SetRangeUser(1e-4, 1.);
    h1_proj->SetLineWidth(2);
    h1_proj->Draw("same hist PLC");
    leg->AddEntry(h1_proj,label[i-1],"l");

  }
//  leg->SetHeader("The Legend Title");
  leg->Draw();
  c.SaveAs("../plots/scripts_root/detectorPlots/sharedPixhits.pdf");

  delete h1,c,h1_proj,leg;
}

void sharedPixhits(std::string hist,std::string alg,int x2,TFile* fData=_file0)
{
  TLegend *leg;
  TH1D *h1_proj;
  TH2D* h1 = (TH2D*)fData->Get((hist+"_"+alg+"_inB").c_str());

//  const char *title="IBL Hits";
  std::string title="sharedPixhits";

  int max_origin=3;
  TCanvas c("c1", "canvas", 1300, 900);
  c.Divide(2,1);
  c.cd(1);
  gStyle->SetOptStat(0);
  gStyle->SetPalette(87);
  h1->SetTitle((title+" "+alg).c_str());
  float N=h1->GetEntries();
  h1->Scale(1./N);
  h1->GetXaxis()->SetRangeUser(0,x2);
  h1->GetXaxis()->SetTitle("Number of sharedPixhits hits");
  h1->GetYaxis()->SetTitle("track origin");
  gStyle->SetPaintTextFormat(".4f");
  h1->Draw("colz text");

  c.cd(2);
  gPad->SetLogy();
  leg = new TLegend(0.7,0.7,0.9,0.9);
  h1_proj = new TH1D("a","b",5,-1,4);
  h1_proj->GetYaxis()->SetRangeUser(0.,0.5);
  std::string idx="projection";
  const char *label[5]={"PU","B","C","FRAG","GEANT"};
  //const char *col[4] = { "KBlue", "KRed", "KOrange", "KYellow" };
  for(unsigned i=1;i<max_origin+3;i++){

    TColor *col=gROOT->GetColor(i);
    char n=i;
    h1_proj= (TH1D*) h1->ProjectionX((idx+"_"+std::to_string(n)).c_str(),i,i);
    h1_proj->GetYaxis()->SetRangeUser(1e-4, 1.);
    h1_proj->SetLineWidth(2);
    h1_proj->Draw("same hist PLC");
    leg->AddEntry(h1_proj,label[i-1],"l");

  }
//  leg->SetHeader("The Legend Title");
  leg->Draw();
  c.SaveAs(("../plots/scripts_root/detectorPlots/sharedPixhits"+alg+".pdf").c_str());

  delete h1,c,h1_proj,leg;
}

void splitPixhits(std::string hist,int x2,TFile* fData=_file0)
{
  TLegend *leg;
  TH1D *h1_proj;
  TH2D* h1 = (TH2D*)fData->Get((hist+"_inB").c_str());

//  const char *title="IBL Hits";
  std::string title="splitPixhits";

  int max_origin=3;
  TCanvas c("c1", "canvas", 1300, 900);
  c.Divide(2,1);
  c.cd(1);
  gStyle->SetOptStat(0);
  gStyle->SetPalette(87);
  h1->SetTitle((title).c_str());
  float N=h1->GetEntries();
  h1->Scale(1./N);
  h1->GetXaxis()->SetRangeUser(0,x2);
  h1->GetXaxis()->SetTitle("Number of splitPixhits hits");
  h1->GetYaxis()->SetTitle("track origin");
  gStyle->SetPaintTextFormat(".4f");
  h1->Draw("colz text");

  c.cd(2);
  gPad->SetLogy();
  leg = new TLegend(0.7,0.7,0.9,0.9);
  h1_proj = new TH1D("a","b",5,-1,4);
  h1_proj->GetYaxis()->SetRangeUser(0.,0.5);
  std::string idx="projection";
  const char *label[5]={"PU","B","C","FRAG","GEANT"};
  //const char *col[4] = { "KBlue", "KRed", "KOrange", "KYellow" };
  for(unsigned i=1;i<max_origin+3;i++){

    TColor *col=gROOT->GetColor(i);
    char n=i;
    h1_proj= (TH1D*) h1->ProjectionX((idx+"_"+std::to_string(n)).c_str(),i,i);
    h1_proj->GetYaxis()->SetRangeUser(1e-4, 1.);
    h1_proj->SetLineWidth(2);
    h1_proj->Draw("same hist PLC");
    leg->AddEntry(h1_proj,label[i-1],"l");

  }
//  leg->SetHeader("The Legend Title");
  leg->Draw();
  c.SaveAs("../plots/scripts_root/detectorPlots/splitPixhits.pdf");

  delete h1,c,h1_proj,leg;
}

void splitPixhits(std::string hist,std::string alg,int x2,TFile* fData=_file0)
{
  TLegend *leg;
  TH1D *h1_proj;
  TH2D* h1 = (TH2D*)fData->Get((hist+"_"+alg+"_inB").c_str());

//  const char *title="IBL Hits";
  std::string title="splitPixhits";

  int max_origin=3;
  TCanvas c("c1", "canvas", 1300, 900);
  c.Divide(2,1);
  c.cd(1);
  gStyle->SetOptStat(0);
  gStyle->SetPalette(87);
  h1->SetTitle((title+" "+alg).c_str());
  float N=h1->GetEntries();
  h1->Scale(1./N);
  h1->GetXaxis()->SetRangeUser(0,x2);
  h1->GetXaxis()->SetTitle("Number of splitPixhits hits");
  h1->GetYaxis()->SetTitle("track origin");
  gStyle->SetPaintTextFormat(".4f");
  h1->Draw("colz text");

  c.cd(2);
  gPad->SetLogy();
  leg = new TLegend(0.7,0.7,0.9,0.9);
  h1_proj = new TH1D("a","b",5,-1,4);
  h1_proj->GetYaxis()->SetRangeUser(0.,0.5);
  std::string idx="projection";
  const char *label[5]={"PU","B","C","FRAG","GEANT"};
  //const char *col[4] = { "KBlue", "KRed", "KOrange", "KYellow" };
  for(unsigned i=1;i<max_origin+3;i++){

    TColor *col=gROOT->GetColor(i);
    char n=i;
    h1_proj= (TH1D*) h1->ProjectionX((idx+"_"+std::to_string(n)).c_str(),i,i);
    h1_proj->GetYaxis()->SetRangeUser(1e-4, 1.);
    h1_proj->SetLineWidth(2);
    h1_proj->Draw("same hist PLC");
    leg->AddEntry(h1_proj,label[i-1],"l");

  }
//  leg->SetHeader("The Legend Title");
  leg->Draw();
  c.SaveAs(("../plots/scripts_root/detectorPlots/splitPixhits"+alg+".pdf").c_str());

  delete h1,c,h1_proj,leg;
}

void nSCThits(std::string hist,int x2,TFile* fData=_file0)
{
  TLegend *leg;
  TH1D *h1_proj;
  TH2D* h1 = (TH2D*)fData->Get((hist+"_inB").c_str());

//  const char *title="IBL Hits";
  std::string title="nSCThits";

  int max_origin=3;
  TCanvas c("c1", "canvas", 1300, 900);
  c.Divide(2,1);
  c.cd(1);
  gStyle->SetOptStat(0);
  gStyle->SetPalette(87);
  h1->SetTitle((title).c_str());
  float N=h1->GetEntries();
  h1->Scale(1./N);
  h1->GetXaxis()->SetRangeUser(0,x2);
  h1->GetXaxis()->SetTitle("Number of nSCThits hits");
  h1->GetYaxis()->SetTitle("track origin");
  gStyle->SetPaintTextFormat(".4f");
  h1->Draw("colz text");

  c.cd(2);
  gPad->SetLogy();
  leg = new TLegend(0.7,0.7,0.9,0.9);
  h1_proj = new TH1D("a","b",5,-1,4);
  h1_proj->GetYaxis()->SetRangeUser(0.,0.5);
  std::string idx="projection";
  const char *label[5]={"PU","B","C","FRAG","GEANT"};
  //const char *col[4] = { "KBlue", "KRed", "KOrange", "KYellow" };
  for(unsigned i=1;i<max_origin+3;i++){

    TColor *col=gROOT->GetColor(i);
    char n=i;
    h1_proj= (TH1D*) h1->ProjectionX((idx+"_"+std::to_string(n)).c_str(),i,i);
    h1_proj->GetYaxis()->SetRangeUser(1e-4, 1.);
    h1_proj->SetLineWidth(2);
    h1_proj->Draw("same hist PLC");
    leg->AddEntry(h1_proj,label[i-1],"l");

  }
//  leg->SetHeader("The Legend Title");
  leg->Draw();
  c.SaveAs("../plots/scripts_root/detectorPlots/nSCThits.pdf");

  delete h1,c,h1_proj,leg;
}

void nSCThits(std::string hist,std::string alg,int x2,TFile* fData=_file0)
{
  TLegend *leg;
  TH1D *h1_proj;
  TH2D* h1 = (TH2D*)fData->Get((hist+"_"+alg+"_inB").c_str());

//  const char *title="IBL Hits";
  std::string title="nSCThits";

  int max_origin=3;
  TCanvas c("c1", "canvas", 1300, 900);
  c.Divide(2,1);
  c.cd(1);
  gStyle->SetOptStat(0);
  gStyle->SetPalette(87);
  h1->SetTitle((title+" "+alg).c_str());
  float N=h1->GetEntries();
  h1->Scale(1./N);
  h1->GetXaxis()->SetRangeUser(0,x2);
  h1->GetXaxis()->SetTitle("Number of nSCThits hits");
  h1->GetYaxis()->SetTitle("track origin");
  gStyle->SetPaintTextFormat(".4f");
  h1->Draw("colz text");

  c.cd(2);
  gPad->SetLogy();
  leg = new TLegend(0.7,0.7,0.9,0.9);
  h1_proj = new TH1D("a","b",5,-1,4);
  h1_proj->GetYaxis()->SetRangeUser(0.,0.5);
  std::string idx="projection";
  const char *label[5]={"PU","B","C","FRAG","GEANT"};
  //const char *col[4] = { "KBlue", "KRed", "KOrange", "KYellow" };
  for(unsigned i=1;i<max_origin+3;i++){

    TColor *col=gROOT->GetColor(i);
    char n=i;
    h1_proj= (TH1D*) h1->ProjectionX((idx+"_"+std::to_string(n)).c_str(),i,i);
    h1_proj->GetYaxis()->SetRangeUser(1e-4, 1.);
    h1_proj->SetLineWidth(2);
    h1_proj->Draw("same hist PLC");
    leg->AddEntry(h1_proj,label[i-1],"l");

  }
//  leg->SetHeader("The Legend Title");
  leg->Draw();
  c.SaveAs(("../plots/scripts_root/detectorPlots/nSCThits"+alg+".pdf").c_str());

  delete h1,c,h1_proj,leg;
}

void sharedSCThits(std::string hist,int x2,TFile* fData=_file0)
{
  TLegend *leg;
  TH1D *h1_proj;
  TH2D* h1 = (TH2D*)fData->Get((hist+"_inB").c_str());

//  const char *title="IBL Hits";
  std::string title="sharedSCThits";

  int max_origin=3;
  TCanvas c("c1", "canvas", 1300, 900);
  c.Divide(2,1);
  c.cd(1);
  gStyle->SetOptStat(0);
  gStyle->SetPalette(87);
  h1->SetTitle((title).c_str());
  float N=h1->GetEntries();
  h1->Scale(1./N);
  h1->GetXaxis()->SetRangeUser(0,x2);
  h1->GetXaxis()->SetTitle("Number of sharedSCThits hits");
  h1->GetYaxis()->SetTitle("track origin");
  gStyle->SetPaintTextFormat(".4f");
  h1->Draw("colz text");

  c.cd(2);
  gPad->SetLogy();
  leg = new TLegend(0.7,0.7,0.9,0.9);
  h1_proj = new TH1D("a","b",5,-1,4);
  h1_proj->GetYaxis()->SetRangeUser(0.,0.5);
  std::string idx="projection";
  const char *label[5]={"PU","B","C","FRAG","GEANT"};
  //const char *col[4] = { "KBlue", "KRed", "KOrange", "KYellow" };
  for(unsigned i=1;i<max_origin+3;i++){

    TColor *col=gROOT->GetColor(i);
    char n=i;
    h1_proj= (TH1D*) h1->ProjectionX((idx+"_"+std::to_string(n)).c_str(),i,i);
    h1_proj->GetYaxis()->SetRangeUser(1e-4, 1.);
    h1_proj->SetLineWidth(2);
    h1_proj->Draw("same hist PLC");
    leg->AddEntry(h1_proj,label[i-1],"l");

  }
//  leg->SetHeader("The Legend Title");
  leg->Draw();
  c.SaveAs("../plots/scripts_root/detectorPlots/sharedSCThits.pdf");

  delete h1,c,h1_proj,leg;
}

void sharedSCThits(std::string hist,std::string alg,int x2,TFile* fData=_file0)
{
  TLegend *leg;
  TH1D *h1_proj;
  TH2D* h1 = (TH2D*)fData->Get((hist+"_"+alg+"_inB").c_str());

//  const char *title="IBL Hits";
  std::string title="sharedSCThits";

  int max_origin=3;
  TCanvas c("c1", "canvas", 1300, 900);
  c.Divide(2,1);
  c.cd(1);
  gStyle->SetOptStat(0);
  gStyle->SetPalette(87);
  h1->SetTitle((title+" "+alg).c_str());
  float N=h1->GetEntries();
  h1->Scale(1./N);
  h1->GetXaxis()->SetRangeUser(0,x2);
  h1->GetXaxis()->SetTitle("Number of sharedSCThits hits");
  h1->GetYaxis()->SetTitle("track origin");
  gStyle->SetPaintTextFormat(".4f");
  h1->Draw("colz text");

  c.cd(2);
  gPad->SetLogy();
  leg = new TLegend(0.7,0.7,0.9,0.9);
  h1_proj = new TH1D("a","b",5,-1,4);
  h1_proj->GetYaxis()->SetRangeUser(0.,0.5);
  std::string idx="projection";
  const char *label[5]={"PU","B","C","FRAG","GEANT"};
  //const char *col[4] = { "KBlue", "KRed", "KOrange", "KYellow" };
  for(unsigned i=1;i<max_origin+3;i++){

    TColor *col=gROOT->GetColor(i);
    char n=i;
    h1_proj= (TH1D*) h1->ProjectionX((idx+"_"+std::to_string(n)).c_str(),i,i);
    h1_proj->GetYaxis()->SetRangeUser(1e-4, 1.);
    h1_proj->SetLineWidth(2);
    h1_proj->Draw("same hist PLC");
    leg->AddEntry(h1_proj,label[i-1],"l");

  }
//  leg->SetHeader("The Legend Title");
  leg->Draw();
  c.SaveAs(("../plots/scripts_root/detectorPlots/sharedSCThits"+alg+".pdf").c_str());

  delete h1,c,h1_proj,leg;
}
