#include <TCanvas.h>
#include <TFile.h>
#include <TH2F.h>
#include <string>
#include <vector>
#include "AtlasUtils.C"
#include "AtlasLabels.C"

TCanvas* show2Dplot(TH2F* h,std::string ss="", bool absolute=true)
{
  ss = std::string(h->GetName())+ss;
  TCanvas* cLinTot = new TCanvas(ss.c_str(),ss.c_str(),800,700);
  gPad->SetLogz();
  //cLinTot->SetLeftMargin(cLinTot->GetLeftMargin()*0.9);
  cLinTot->SetTopMargin(cLinTot->GetTopMargin()*2.);
  cLinTot->SetLeftMargin(cLinTot->GetLeftMargin()*1.);
  cLinTot->SetRightMargin(cLinTot->GetRightMargin()*3.8);
  cLinTot->SetBottomMargin(cLinTot->GetBottomMargin()*1.3);
  std::string labX[6] = {"1B", "1D, 0B", "0B 0D", ">1B", ">1D, 0B", ">1B, >1D"};
  for (int i=0; i<6; ++i)
    {
      h->GetXaxis()->SetBinLabel(i+1,labX[i].c_str());
      h->GetYaxis()->SetBinLabel(i+1,labX[i].c_str());
    }
  h->LabelsDeflate("X");
  h->LabelsDeflate("Y");
  h->LabelsOption("v");
  h->GetXaxis()->SetLabelSize(0.04);
  h->GetYaxis()->SetLabelSize(0.04);
  if (absolute)  h->GetZaxis()->SetTitle("Number of jets");
  else  h->GetZaxis()->SetTitle("Fraction of jets [%]");
  h->GetXaxis()->SetTitle("Cone labeling");
  h->GetYaxis()->SetTitle("Ghost labeling");
  h->GetZaxis()->SetTitleOffset(1.25);
  h->GetYaxis()->SetTitleOffset(1.5);
  h->GetXaxis()->SetTitleOffset(1.9);
  h->Draw("ZCOLTEXT45");
  /*
  float x=0.15;
  float y=0.05;
  TLatex l;
  l.SetNDC();
  l.SetTextFont(42);
  l.SetTextSize(0.04);
  l.DrawLatex(x,y,"Cone labeling");
  //l.SetAngle(90.);
  //x = 0.5;
  //y=0.5;
  //l.PaintLatex(x,y,90.,0.04,"Ghost labeling");
  */
  return cLinTot;
}

void show_labelMatrix()
{

  bool absolute=true;
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
    //
  std::string loc="/afs/le.infn.it/user/s/spagnolo/atlas/Athena/FTAGmyfork/root_selector/DAOD_selector/finalHistos/";
  TFile *_file0 = TFile::Open((loc+"debug_bTag_AntiKtVR30Rmax4Rmin02TrackJets_BTagging201903_ConeIncl.root").c_str());
  TFile *_file1 = TFile::Open((loc+"debug_bTag_AntiKtVR30Rmax4Rmin02TrackGhostTagJets_GhostIncl_12GeV.root").c_str());
  TFile *_file2 = TFile::Open((loc+"debug_bTag_AntiKt4EMPFlowJets_BTagging201903_ConeIncl.root").c_str());
  TH2F* hmatVR20 = (TH2F*)_file0->Get("jetFlavorLabelMatrix");
  TH2F* hmatVR12 = (TH2F*)_file1->Get("jetFlavorLabelMatrix");
  TH2F* hmatEMPf = (TH2F*)_file2->Get("jetFlavorLabelMatrix");
  TCanvas* c = new TCanvas("c","c",800,600);
  hmatVR20->Draw();
  hmatVR12->Draw();
  hmatEMPf->Draw();
  
  double x = 0.15;
  double y = 0.915;
  TLatex l;
  l.SetNDC();
  l.SetTextFont(42);
  l.SetTextSize(0.03);
  std::string collection = "AntiKt4EMPFlowJets";

  std::vector<TCanvas*> vCanvas;

  //  if (absolute)
  //   {
  //      gStyle->SetPaintTextFormat("6.0f");
      
  vCanvas.push_back(show2Dplot((TH2F*)hmatVR20,"_VR20",absolute));
      ATLASLabel(x,y,"Internal");
      collection = "AntiKtVR30Rmax4Rmin02TrackJets";
      l.DrawLatex(x+0.35,y+0.035,("mc16d, t#bar{t}, jet p_{T} > 20 GeV"));
      l.DrawLatex(x+0.35,y,(collection).c_str());
      vCanvas.back()->SaveAs("labelmatrix_AntiKtVR30Rmax4Rmin02TrackJets_20GeV.pdf");
      
      vCanvas.push_back(show2Dplot((TH2F*)hmatVR12,"_VR12",absolute));
      ATLASLabel(x,y,"Internal");
      collection = "AntiKtVR30Rmax4Rmin02TrackJets";
      l.DrawLatex(x+0.35,y+0.035,("mc16d, t#bar{t}, jet p_{T} > 12 GeV"));
      l.DrawLatex(x+0.35,y,(collection).c_str());
      vCanvas.back()->SaveAs("labelmatrix_AntiKtVR30Rmax4Rmin02TrackJets_12GeV.pdf");
      
      vCanvas.push_back(show2Dplot((TH2F*)hmatEMPf,"_EMPf",absolute)); 
      ATLASLabel(x,y,"Internal");
      collection = "AntiKt4EMPFlowJets";
      l.DrawLatex(x+0.35,y+0.035,("mc16d, t#bar{t}, jet p_{T} > 20 GeV"));
      l.DrawLatex(x+0.35,y,(collection).c_str());
      vCanvas.back()->SaveAs("labelmatrix_AntiKt4EMPFlowJets_20GeV.pdf");

      //return;
      //    }
      //  else
      //    {
      gStyle->SetPaintTextFormat("6.4f");
      hmatVR12N=hmatVR12->Clone();
      ((TH2F*)hmatVR12N)->Scale(100./((TH2F*)hmatVR12N)->GetEntries());
      hmatEMPfN=hmatEMPf->Clone();
      ((TH2F*)hmatEMPfN)->Scale(100./((TH2F*)hmatEMPfN)->GetEntries());
      hmatVR20N=hmatVR20->Clone();
      ((TH2F*)hmatVR20N)->Scale(100./((TH2F*)hmatVR20N)->GetEntries());
      c->cd();
      hmatVR12N->Draw();
      hmatVR20N->Draw();
      hmatEMPfN->Draw();

      absolute = false;
      
      vCanvas.push_back(show2Dplot((TH2F*)hmatVR20N,"_VR20N",absolute));
      ATLASLabel(x,y,"Internal");
      collection = "AntiKtVR30Rmax4Rmin02TrackJets";
      l.DrawLatex(x+0.35,y+0.035,("mc16d, t#bar{t}, jet p_{T} > 20 GeV"));
      l.DrawLatex(x+0.35,y,(collection).c_str());
      vCanvas.back()->SaveAs("labelmatrix_AntiKtVR30Rmax4Rmin02TrackJets_20GeV_rel.pdf");
      
      vCanvas.push_back(show2Dplot((TH2F*)hmatVR12N,"_VR12N",absolute)); 
      ATLASLabel(x,y,"Internal");
      collection = "AntiKtVR30Rmax4Rmin02TrackJets";
      l.DrawLatex(x+0.35,y+0.035,("mc16d, t#bar{t}, jet p_{T} > 12 GeV"));
      l.DrawLatex(x+0.35,y,(collection).c_str());
      vCanvas.back()->SaveAs("labelmatrix_AntiKtVR30Rmax4Rmin02TrackJets_12GeV_rel.pdf");
      
      vCanvas.push_back(show2Dplot((TH2F*)hmatEMPfN,"_EMPfN",absolute)); 
      ATLASLabel(x,y,"Internal");
      collection = "AntiKt4EMPFlowJets";
      l.DrawLatex(x+0.35,y+0.035,("mc16d, t#bar{t}, jet p_{T} > 20 GeV"));
      l.DrawLatex(x+0.35,y,(collection).c_str());
      vCanvas.back()->SaveAs("labelmatrix_AntiKt4EMPFlowJets_20GeV_rel.pdf");

      //    }
}

