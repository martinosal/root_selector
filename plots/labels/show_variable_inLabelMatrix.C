#include <TCanvas.h>
#include <TFile.h>
#include <TH2F.h>
#include <string>
#include <vector>
#include <TLegend.h>
#include "atlasstyle-00-04-02/AtlasUtils.C"
#include "atlasstyle-00-04-02/AtlasLabels.C"
//#include "atlasstyle-00-04-02/AtlasStyle.C"

void show_variable_inLabelMatrix(std::string coll, std::string var, std::string ConeLab, bool savePlot=true)
{
  std::string loc="/afs/le.infn.it/user/s/spagnolo/atlas/Athena/FTAGmyfork/root_selector/DAOD_selector/finalHistos/hForLabels/";
  //std::string loc="/afs/le.infn.it/project/itk/forWEB/FTag/DOADsel_finalHistos/hForLabels/";
  std::string fileName;
  std::string jetCollection;
  std::string outputName;

  if(coll.find("VR10Cone")!=string::npos){
    if(var=="JFdR" || var=="JFeFc" || var=="JFmVtx" || var=="JFnVtx" || var=="JFsig3d" || var=="nTrkAlgo_JF") {cout<<"Variable "<< var << " not needed; return."<<endl; return;}
    fileName = "debug_Labels_bTag_AntiKtVR30Rmax4Rmin02TrackJets_BTagging201903_Cone_10GeV.root";
    jetCollection = "AntiKtVR30Rmax4Rmin02TrackJets, p_{T}>10 GeV";
    outputName = "VR30Cone_10GeV_Labels_";
  }
  else if(coll.find("EMPFlow")!=string::npos){
    fileName = "debug_Labels_bTag_AntiKt4EMPFlowJets_BTagging201903_Cone.root";
    jetCollection = "AntiKt4EMPFlowJets";
    outputName = "EMPf_20GeV_Labels_";
  }
  else if(coll.find("VR10Ghost")!=string::npos){
    if(var=="JFdR" || var=="JFeFc" || var=="JFmVtx" || var=="JFnVtx" || var=="JFsig3d" || var=="nTrkAlgo_JF") {cout<<"Variable "<<var<<" not needed; return."<<endl; return;}
    fileName = "debug_Labels_bTag_AntiKtVR30Rmax4Rmin02TrackGhostTagJets_Ghost_10GeV.root";
    jetCollection = "AntiKtVR30Rmax4Rmin02TrackGhost, p_{T}>10 GeV";
    outputName = "VR30Ghost_10GeV_Labels_";
}
  else {cout << "Collection not found. Possible collection: VR10Cone, EMPFlow, VR10Ghost"; return;}

  if(ConeLab=="0B0D" || ConeLab=="1D0B" || ConeLab=="2D0B"){
    if(var=="bHPt" || var=="bHPtFraction" || var=="bHjetDR") {
      cout << "0 b-hadron found, var "<<var<<" not needed; return. " << endl; 
      return;
    }
}

  SetAtlasStyle();
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  TCanvas* cLinTot = new TCanvas("a","a",800,600);
  if(var=="jetPt" || var=="jetEta" || var=="nTrkAlgo_IP3D" || var=="nTrkAlgo_JF" || var=="nTrkAlgo_SV1" || var=="JFdR" || var=="JFeFc" || var=="JFmVtx" || var=="JFnVtx" || var=="JFsig3d" || var=="SV1eFc" || var=="SV1mVtx" || var=="SV1sig3d" || var=="bHPt" || var=="bHPtFraction" || var=="bHjetDR" || var=="SV1dR") gPad->SetLogy();


  TFile *_file0;
  _file0 = TFile::Open((loc+fileName).c_str());
  if (_file0 != NULL) std::cout<<"File open"<<std::endl;

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
    std::cout<<" Bad matrix coloumn (cone label) requested "<<ConeLab<<" ... stop here "<<std::endl;
    return;
  }
  TH1D* h[6];
  TH1D* hN[6];
  //_file0->ls();
  for (unsigned int i=0; i<6; ++i)
    {
      h[i]=NULL;
      hNameBin=HistoName[index]+"_"+HistoName[i]+"_";
      hSel="Labels_"+hNameBin+var;
      std::cout<<"Looking for histogram named <"<<hSel<<">"<<std::endl;
      h[i] = (TH1D*)_file0->Get(hSel.c_str());
      if (h[i]== NULL)
	{
	  std::cout<<" Histogram <"
	    <<hSel<<"> not found !!!!"
	    <<std::endl;
	      return;					       
	}
      h[i]->Draw();
      //      if (h[i]->GetEntries()<100) continue;
      std::cout<<"Found with n. of entries = "<<h[i]->GetEntries()<<" under/overflow="<< h[i]->GetBinContent(0)<<" / "<< h[i]->GetBinContent(h[i]->GetNbinsX()+1)<<" entries in range ="<<h[i]->GetEntries()-h[i]->GetBinContent(h[i]->GetNbinsX()+1)-h[i]->GetBinContent(0)<<std::endl;
    }
  std::cout<<"Files found ! "<<std::endl;
  
  double ymax = 0.;
  unsigned int imax = 0;

  for (unsigned int i=0; i<6; ++i)
    {
      if (h[i]->GetEntries()<100) continue;
      std::cout<<"Before scaling Histogram "<<h[i]->GetName()<<" has "<<h[i]->GetEntries()<<" entries and it will be used - integral "<<h[i]->Integral()<<std::endl;
      std::cout<<"Before scaling  Histogram "<<h[i]->GetName()<<" has "<<h[i]->GetEntries()<<" entries and it will be used - integral (sum of bin content) "<<h[i]->Integral()<<" integral (mat style) "<<h[i]->Integral("width")<<std::endl;

      if (h[i]->GetSumw2N()  == 0) h[i]->Sumw2(kTRUE);
      hN[i]=(TH1D*)h[i]->Clone();
      if (hN[i]->GetSumw2N() == 0) hN[i]->Sumw2(kTRUE);
      hN[i]->Scale(1./h[i]->Integral(1,h[i]->GetNbinsX()));
      //hN[i]->Scale(1.);
      std::cout<<"After scaling  Histogram "<<hN[i]->GetName()<<" has "<<hN[i]->GetEntries()<<" entries and it will be used - integral (sum of bin content) "<<hN[i]->Integral()<<" integral (mat style) "<<hN[i]->Integral("width")<<std::endl;
      hN[i]->SetTitle((HistoName[i]).c_str());
      if (ymax<hN[i]->GetMaximum())
	{
	  ymax = hN[i]->GetMaximum();
	  imax=i;
	}
    }
  std::cout<<"Maximum for histo "<<hN[imax]->GetName()<<" maximum: "<<ymax<<std::endl;

  TLegend* legend;
  if(var=="jetPt" || var=="JFdR" || var=="JFeFc" || var=="JFmVtx" || var=="JFnVtx" || var=="JFsig3d" || var=="SV1mVtx" || var=="SV1sig3d" || var=="bHPt" || var=="bHjetDR" || var=="nTrkAlgo_JF" || var=="SV1dR") legend = new TLegend(0.40,0.58,0.95,0.78);
  else if(var=="bHPtFraction") legend = new TLegend(0.50,0.62,0.95,0.78);
  else if(var=="jetEta")  legend = new TLegend(0.40,0.67,0.95,0.82);
  else if(var=="nTrkAlgo_IP3D") legend = new TLegend(0.42,0.66,0.97,0.81);
  else if(var=="nTrkAlgo_SV1") legend = new TLegend(0.40,0.64,0.95,0.79);
  else if(var=="SV1eFc") legend = new TLegend(0.40,0.62,0.95,0.82);
  else if(var=="SV1nVtx") legend = new TLegend(0.42,0.57,0.97,0.77);
  else legend = new TLegend(0.42,0.66,0.97,0.81);
  
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->SetTextFont(42);
  legend->SetTextSize(0.035);
    
  hN[imax]->Draw();
  //return;
  //std::vector<int> color = {2,3,6,7,9,94};
  std::vector<int> color = {63,2,8,9,6,94};
  for (unsigned int i=0; i<6; ++i)
    {
      if (h[i]->GetEntries()<100) continue;
      std::string title = hN[i]->GetTitle();
      if(title!=ConeLab) {
	hN[i]->SetMarkerSize(1.);
	hN[i]->SetMarkerStyle(23);
	hN[i]->SetLineColor(color[i]);
	hN[i]->SetMarkerColor(color[i]);
      }
      else {
	hN[i]->SetMarkerSize(1.5);
	hN[i]->SetLineColor(1);
	hN[i]->SetMarkerColor(1);
      }
      if(var=="nTrkAlgo_IP3D") hN[i]->GetXaxis()->SetTitle("# of tracks used by IPxD and RNNIP");
      if(var=="JFdR") hN[i]->GetXaxis()->SetTitle("JF vertex #DeltaR");
      if(var=="bHjetDR") hN[i]->GetXaxis()->SetTitle("#DeltaR(b-hadron,jet)");
      if(var=="jetPt" || var=="SV1mVtx" || var=="SV1sig3d" || var=="bHPt") hN[i]->SetMaximum(ymax*10);
      else if (var=="jetEta") hN[i]->SetMaximum(ymax*15);
      else if(var=="SV1nVtx") hN[i]->SetMaximum(ymax+0.22);
      else if(var=="nTrkAlgo_SV1" || var=="SV1eFc" || var=="nTrkAlgo_IP3D" || var=="SV1dR") hN[i]->SetMaximum(ymax*20);
      hN[i]->GetXaxis()->SetDecimals();
      hN[i]->GetYaxis()->SetDecimals();
      hN[i]->Draw("SAME");
      if(title==ConeLab) hN[i]->Draw("HISTSAME");
      legend->AddEntry(hN[i],(HistoName[index]+" cone Lab, "+HistoName[i]+" ghost Lab").c_str(),"PL");
    }
  legend->Draw();

  double x = 0.2;
  double y = 0.88;
  ATLASLabel(x,y,"Internal");
  TLatex l;
  l.SetNDC();
  l.SetTextFont(42);
  l.SetTextSize(0.03);
  l.DrawLatex(x,y-0.03,("mc16d, t#bar{t}, "+jetCollection).c_str());

  //string outputDir = "/afs/le.infn.it/project/itk/forWEB/FTag/DOADsel_finalHistos/hForLabels/plots/";
  string outputDir = "./tmp/";
  //cLinTot->SaveAs((outputDir+coll+"/"+coll+"_"+HistoName[index]+"_"+var+".pdf").c_str());
  if (savePlot) cLinTot->SaveAs((outputDir+"/"+coll+"_"+HistoName[index]+"_"+var+".pdf").c_str());
}
