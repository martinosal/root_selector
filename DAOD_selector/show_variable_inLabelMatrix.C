#include <TCanvas.h>
#include <TFile.h>
#include <TH2F.h>
#include <string>
#include <vector>
#include <TLegend.h>
#include "atlasstyle-00-04-02/AtlasUtils.C"
#include "atlasstyle-00-04-02/AtlasLabels.C"
#include "atlasstyle-00-04-02/AtlasStyle.C"

void show_variable_inLabelMatrix(std::string coll, std::string var, std::string ConeLab)
{
  SetAtlasStyle();
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  TCanvas* cLinTot = new TCanvas("a","a",800,600);
  if(var=="jetPt" || var=="jetEta" || var=="nTrkAlgo_IP3D" || var=="nTrkAlgo_RNNIP" || var=="nTrkAlgo_JF" || var=="nTrkAlgo_SV1" || var=="JFdR" || var=="JFeFc" || var=="JFmVtx" || var=="JFnVtx" || var=="JFsig3d" || var=="SV1eFc" || var=="SV1mVtx" || var=="SV1sig3d" || var=="bHPt" || var=="bHPtFraction" || var=="bHjetDR") gPad->SetLogy();
  std::string loc="/afs/le.infn.it/user/s/spagnolo/atlas/Athena/FTAGmyfork/root_selector/DAOD_selector/finalHistos/forLabelMatrixPlots/";
  std::string fileName;
  std::string jetCollection;
  std::string outputName;

  if(coll.find("Cone10")!=string::npos){
    fileName = "debug_Labels_bTag_AntiKtVR30Rmax4Rmin02TrackJets_BTagging201903_Cone_10GeV.root";
    jetCollection = "AntiKtVR30Rmax4Rmin02TrackJets, p_{T}>10 GeV";
    outputName = "VR30_10GeV_Labels_";
  }
  else if(coll.find("Cone12")!=string::npos){
    fileName = "debug_Labels_bTag_AntiKtVR30Rmax4Rmin02TrackJets_BTagging201903_Cone_12GeV.root";
    jetCollection = "AntiKtVR30Rmax4Rmin02TrackJets, p_{T}>12 GeV";
    outputName = "VR30_12GeV_Labels_";
  }
  else if(coll.find("Cone20")!=string::npos){
    fileName = "debug_Labels_bTag_AntiKtVR30Rmax4Rmin02TrackJets_BTagging201903_Cone_20GeV.root";
    jetCollection = "AntiKtVR30Rmax4Rmin02TrackJets, p_{T}>20 GeV";
    outputName = "VR30_20GeV_Labels_";
  }
  else if(coll.find("EMPFlow")!=string::npos){
    fileName = "debug_Labels_bTag_AntiKt4EMPFlowJets_BTagging201903_Cone.root";
    jetCollection = "AntiKt4EMPFlowJets";
    outputName = "EMPf_20GeV_Labels_";
  }
  else {cout << "Collection not found. Possible collection: Cone10, Cone12, Cone20, EMPFlow"; return;}

  TFile *_file0 = TFile::Open((loc+fileName).c_str());

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
  //_file0->ls();
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
      hN[i]->SetTitle((HistoName[i]).c_str());
      if (ymax<hN[i]->GetMaximum())
	{
	  ymax = hN[i]->GetMaximum();
	  imax=i;
	}
    }
  std::cout<<"Maximum for histo "<<h[imax]->GetName()<<" maximum: "<<ymax<<std::endl;

  TLegend* legend;
  if(var=="jetPt" || var=="JFdR" || var=="JFeFc" || var=="JFmVtx" || var=="JFnVtx" || var=="JFsig3d" || var=="SV1mVtx" || var=="SV1sig3d" || var=="bHPt" || var=="bHPtFraction" || var=="bHjetDR") legend = new TLegend(0.40,0.58,0.95,0.78);
  else if(var=="jetEta")  legend = new TLegend(0.40,0.67,0.95,0.82);
  else if(var=="nTrkAlgo_IP2D") legend = new TLegend(0.15,0.64,0.70,0.79);
  else if(var=="nTrkAlgo_IP3D" || var=="nTrkAlgo_RNNIP" || var=="nTrkAlgo_JF") legend = new TLegend(0.37,0.2,0.92,0.37);
  //else if(var=="nTrkAlgo_SV1" || ConeLab!="1D0B") legend = new TLegend(0.12,0.2,0.67,0.37);
  else if(var=="nTrkAlgo_SV1") legend = new TLegend(0.40,0.64,0.95,0.79);
  else if(var=="SV1eFc") legend = new TLegend(0.40,0.62,0.95,0.82);
  else if(var=="SV1nVtx") legend = new TLegend(0.42,0.57,0.97,0.77);

  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->SetTextFont(42);
  legend->SetTextSize(0.035);
    
  //hN[imax]->SetMarkerSize(0.8);
  hN[imax]->Draw();
  //legend->AddEntry(hN[imax],(HistoName[index]+" cone Lab, "+HistoName[imax]+" ghost Lab").c_str(),"P");
  //return;
  std::vector<int> color = {2,3,6,7,9,94};
  for (unsigned int i=0; i<6; ++i)
    {
      if (h[i]->GetEntries()<100) continue;
      //if (i==imax) continue; 
      std::string title = hN[i]->GetTitle();
      if(title!=ConeLab) {
	hN[i]->SetMarkerSize(0.8);
	hN[i]->SetLineColor(color[i]);
	hN[i]->SetMarkerColor(color[i]);
      }
      else {
	hN[i]->SetMarkerSize(1.2);
	hN[i]->SetLineColor(1);
	hN[i]->SetMarkerColor(1);
      }
      if(var=="jetEta") hN[i]->GetXaxis()->SetTitle("Jet #eta");
      if(var=="jetPt" || var=="SV1mVtx" || var=="SV1sig3d" || var=="bHPt") hN[i]->SetMaximum(ymax*10);
      else if (var=="jetEta") hN[i]->SetMaximum(ymax*15);
      else if(var=="nTrkAlgo_IP2D") {
	if(ConeLab=="1B" || ConeLab=="1D0B") hN[i]->SetMaximum(ymax+0.19);
	else hN[i]->SetMaximum(ymax+0.13);
      }
      else if(var=="nTrkAlgo_IP2D") hN[i]->SetMaximum(ymax+0.13);
      else if(var=="SV1nVtx") hN[i]->SetMaximum(ymax+0.22);
      else if(var=="nTrkAlgo_SV1" || var=="SV1eFc") hN[i]->SetMaximum(ymax*20);
      hN[i]->GetYaxis()->SetDecimals();
      hN[i]->Draw("SAME");
      legend->AddEntry(hN[i],(HistoName[index]+" cone Lab, "+HistoName[i]+" ghost Lab").c_str(),"P");
    }
  legend->Draw();

  double x = 0.2;
  double y = 0.88;
  ATLASLabel(x,y,"Internal");
  TLatex l;
  l.SetNDC();
  l.SetTextFont(42);
  l.SetTextSize(0.025);
  l.DrawLatex(x,y-0.03,("mc16d, t#bar{t}, "+jetCollection).c_str());

  string outputDir = "/afs/le.infn.it/project/itk/forWEB/FTag/DOADsel_finalHistos/hForLabels/plots/";
  cLinTot->SaveAs((outputDir+coll+"_"+HistoName[index]+"_"+var+".pdf").c_str());
}
