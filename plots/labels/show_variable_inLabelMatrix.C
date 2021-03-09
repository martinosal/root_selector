#include <TCanvas.h>
#include <TFile.h>
#include <TH2F.h>
#include <string>
#include <vector>
#include <TLegend.h>
#include "AtlasUtils.C"
#include "AtlasLabels.C"


void show_variable_inLabelMatrix(std::string var, std::string ConeLab)
{

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  TCanvas* cLinTot = new TCanvas("a","a",800,600);
  gPad->SetLogy();
    //
  std::string loc="/afs/le.infn.it/user/s/spagnolo/atlas/Athena/FTAGmyfork/root_selector/DAOD_selector/finalHistos/forLabelMatrixPlots/";
  TFile *_file0 = TFile::Open((loc+"debug_Labels_bTag_AntiKtVR30Rmax4Rmin02TrackJets_BTagging201903_Cone_10GeV.root").c_str());
  //  TFile *_file0 = TFile::Open((loc+"debug_Labels_bTag_AntiKtVR30Rmax4Rmin02TrackJets_BTagging201903_Cone_12GeV.root").c_str());
  //  TFile *_file0 = TFile::Open((loc+"debug_Labels_bTag_AntiKtVR30Rmax4Rmin02TrackJets_BTagging201903_Cone_20GeV.root").c_str());
  //  TFile *_file0 = TFile::Open((loc+"debug_Labels_bTag_AntiKt4EMPFlowJets_BTagging201903_Cone.root").c_str());


  std::string HistoName[6] = {"1B", "1D0B", "0B0D", "2B", "2D0B", "2B2D"};
  std::string hNameBin; //= HistoName[iX]+"_"+HistoName[iY];
  std::string hSel;
  unsigned int index=99;
  for (unsigned int i=0; i<6; ++i)
    {
      if (ConeLab==HistoName[i])
	{
	  index=i;
	  break;
	}
    }
  if (index==99) {
    std::cout<<" ConeLab "<<ConeLab<<" not found ... stop here "<<std::endl;
    return;
  }
  TH1D* h[6];
  TH1D* hN[6];
  _file0->ls();
  for (unsigned int i=0; i<6; ++i)
    {
      hNameBin=HistoName[index]+"_"+HistoName[i]+"_";
      hSel="Labels_"+hNameBin+var;
      std::cout<<"Looking for histogram named <"<<hSel<<">"<<std::endl;
      h[i] = (TH1D*)_file0->Get(hSel.c_str());
      h[i]->Draw();
      //      if (h[i]->GetEntries()<100) continue;      
    }
  std::cout<<"Files found ! "<<std::endl;
  
  double ymax = 0.;
  unsigned int imax = 0;
  for (unsigned int i=0; i<6; ++i)
    {
      if (h[i]->GetEntries()<100) continue;
      std::cout<<"Histogram "<<h[i]->GetName()<<" has "<<h[i]->GetEntries()<<" and it will be used"<<std::endl;

      hN[i]=(TH1D*)h[i]->Clone();
      hN[i]->Scale(1./h[i]->GetEntries());


      hN[i]->SetLineColor(i+1);
      hN[i]->SetMarkerColor(i+1);
      if (ymax<hN[i]->GetMaximum())
	{
	  ymax = hN[i]->GetMaximum();
	  imax=i;
	}
    }
  std::cout<<"Maximum for histo "<<h[imax]->GetName()<<std::endl;


  TLegend* legend = new TLegend(0.4,0.7,0.98,0.9);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
    

  hN[imax]->SetMarkerSize(0.8);
  hN[imax]->Draw();
  //legend->AddEntry(hN[imax],(HistoName[index]+" cone Lab, "+HistoName[imax]+" ghost Lab").c_str(),"P");
  //return;
  for (unsigned int i=0; i<6; ++i)
    {
      if (h[i]->GetEntries()<100) continue;
      //if (i==imax) continue; 
      if (i!=index) hN[i]->SetMarkerSize(0.8);
      else hN[i]->SetMarkerSize(1.2);
      hN[i]->Draw("SAME");
      legend->AddEntry(hN[i],(HistoName[index]+" cone Lab, "+HistoName[i]+" ghost Lab").c_str(),"PE");
    }
  legend->Draw();
  
  cLinTot->SaveAs(("VR30_10GeV_Labels_"+HistoName[index]+"_"+var+".pdf").c_str());
  //cLinTot->SaveAs(("VR30_12GeV_Labels_"+HistoName[index]+"_"+var+".pdf").c_str());
  //cLinTot->SaveAs(("VR30_20GeV_Labels_"+HistoName[index]+"_"+var+".pdf").c_str());
  //cLinTot->SaveAs(("EMPf_20GeV_Labels_"+HistoName[index]+"_"+var+".pdf").c_str());
  
}

