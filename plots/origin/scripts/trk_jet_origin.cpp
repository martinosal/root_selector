
// std::string locForOutputFiles()
// {
//   std::string myPath = "./plots/";
//   return myPath;
// }

void jet_orig_total(std::string hist,TFile* fData=_file0, int DpT=50, double pt_max=1000){

  TH2F* h = (TH2F*)fData->Get(hist.c_str());
  const Int_t n_Xbins = h->GetNbinsX();
  const Int_t n_Ybins = h->GetNbinsY();
  const double y_max = h->GetYaxis()->GetXmax();

//  int DpT=100;
  double bin_to_pT=(double) n_Ybins/y_max;
  int rebin=(int) DpT/bin_to_pT;
  int Nbins=(int) pt_max/DpT;

  double M[Nbins][n_Xbins];
  double M_der[Nbins][7];


  TCanvas c("c", "canvas", 1300, 900);
  c.SetGrid();

  gPad->SetLogy();

  TH1D** px = new TH1D*[Nbins];
  for(int k=0;k<Nbins;k++){
    TString sx=Form("hj%d",k);
    px[k] = h->ProjectionX(sx,bin_to_pT*DpT*k,bin_to_pT*DpT*(k+1));
    px[k]->SetStats(0);
//    float N=px[k]->GetEntries();
//    px[k]->Scale(1./N);
    std::string s = std::to_string(k);
    std::string s1 = std::to_string((int) DpT*k);
    std::string s2 = std::to_string((int) DpT*(k+1));
    px[k]->SetTitle(("origin trk in ["+s1+","+s2+"] Gev jet pT").c_str());
//    px[k]->Draw("hist");
//    c.SaveAs(("/home/salomon/Private/atlas/FTPF/Selector/plots/origin/"+hist+"Xproj_"+s+".pdf").c_str());
    for(int i=1;i<=n_Xbins;i++){
      M[k][i]=px[k]->GetBinContent(i);
//      std::cout<<M[k][i]<<"\t";
    }
//    std::cout<<"\n";
  }

  std::vector<int> truthlabel={0,1,2,3,4,5,6,11,12,13,14,15,16,101,102,103,104,105,106,111,112,113,114,115,116};//25 categories
  std::vector<int> derived_truthlabel={0,1,2,3,4,5};

  std::vector<int> map_0={0};
  std::vector<int> map_1={2,4,6};
  std::vector<int> map_2={12,14,16,102,104,106,112,114,116};
  std::vector<int> map_3={1,11,101,111};
  std::vector<int> map_4={3,5};
  std::vector<int> map_5={13,103,113,15,105,115};

  std::vector<std::vector<int>> list;
  list.push_back(map_0);
    list.push_back(map_1);
      list.push_back(map_2);
        list.push_back(map_3);
          list.push_back(map_4);
            list.push_back(map_5);

  float y=0;

  for(int k=0;k<Nbins;k++){
    std::cout<<"\ntrk in ["<<DpT*k<<","<<DpT*(k+1)<<"] Gev jet pT";
    for(int i=0;i<list.size();i++){
      y=0;
      std::cout<<"\norigin: "<<i<<"\t";
      for(std::vector<int>::iterator it = list.at(i).begin(); it != list.at(i).end(); ++it){
//        std::cout<<*it<<" "<<M[k][*it+2]<<"\t";
        y+=M[k][*it+2];
      }
      std::cout<<y<<"\t";
      M_der[k][i+1]=y;
    }
  }

  double s=0.;
  for(int i=0;i<=6;i++){
//    std::cout<<"\n";
    s=0.;
    for(int k=0;k<Nbins;k++){
      s+=M_der[k][i];
    }
    std::cout<<s<<"\n";
  }

    std::cout<<"\n";
    delete h,M,c,px;
}

void jet_orig(std::string hist,TFile* fData=_file0, int DpT=50, double pt_max=1000){

  TH2F* h = (TH2F*)fData->Get(hist.c_str());
  const Int_t n_Xbins = h->GetNbinsX();
  const Int_t n_Ybins = h->GetNbinsY();
//  const double pt_max = h->GetYaxis()->GetXmax();

//  int DpT=100;
  int rebin=(int) 2*DpT;
  int Nbins=(int) pt_max/(DpT);

  double M[Nbins][5];

  std::cout<<n_Ybins<<", "<<DpT<<", "<<Nbins<<"\n";

  string origin[5] = {"PU","B","C","FRAG","GEANT"};
  int colors[5] = {2,4,8,12,44};
  TCanvas c("c", "canvas", 1300, 900);
  c.SetGrid();

  gPad->SetLogy();

  TH1D** px = new TH1D*[Nbins];
  for(int k=0;k<Nbins;k++){
    TString sx=Form("hj%d",k);
    px[k] = h->ProjectionX(sx,(int) DpT*k/2,(int) DpT*(k+1)/2);
    px[k]->SetStats(0);
    float N=px[k]->GetEntries();
    px[k]->Scale(1./N);
    std::string s = std::to_string(k);
    std::string s1 = std::to_string((int) DpT*k/2);
    std::string s2 = std::to_string((int) DpT*(k+1)/2);
    px[k]->SetTitle(("origin trk in ["+s1+","+s2+"] Gev jet pT").c_str());
//    px[k]->Draw("hist");
//    c.SaveAs(("/home/salomon/Private/atlas/FTPF/Selector/plots/origin/"+hist+"Xproj_"+s+".pdf").c_str());
    for(int i=1;i<=5;i++){
      M[k][i]=px[k]->GetBinContent(i);
      std::cout<<M[k][i]<<"\t";
    }
    std::cout<<"\n";
  }


  auto leg = new TLegend(0.2,0.2,0.35,0.4);//x1,y1,x2,y2
  TH1D** py_2 = new TH1D*[5];
  for(int i=1;i<=5;i++){
    py_2[i] = new TH1D(Form("syj%d",i),"",Nbins,0,pt_max);
    for(int k=0;k<Nbins;k++){
      py_2[i]->AddBinContent(k+1,M[k][i]);
      py_2[i]->SetStats(0);
      py_2[i]->GetYaxis()->SetRangeUser(1e-3, 1.);
      py_2[i]->GetXaxis()->SetTitle("jet pT [GeV]");
      py_2[i]->SetLineWidth(2);
      py_2[i]->SetLineColor(colors[i-1]);
      py_2[i]->SetOption("col");
      py_2[i]->SetMarkerStyle(5);
      py_2[i]->SetMarkerSize(2);
    }
    leg->AddEntry(py_2[i],origin[i-1].c_str(),"l");
      py_2[i]->Draw("same hist");
      py_2[i]->Draw("hist p same");
//    c.SaveAs(("/home/salomon/Private/atlas/FTPF/Selector/plots/origin/"+hist+"_"+origin[i-1]+"_binN.pdf").c_str());
  }
    leg->Draw();
    c.SaveAs((locForOutputFiles()+hist+"_binN.pdf").c_str());
    delete h,M,c,px,py_2;
}

void trk_orig(std::string hist,TFile* fData=_file0, int DpT=50){//fractional composition

  TH2F* h = (TH2F*)fData->Get(hist.c_str());
  const Int_t n_Xbins = h->GetNbinsX();
  const Int_t n_Ybins = h->GetNbinsY();
  const double pt_max = h->GetYaxis()->GetXmax();

//  int DpT=100;
  int rebin=(int) 2*DpT;
  int Nbins=(int) pt_max/(DpT);

  double M[Nbins][5];

  std::cout<<n_Ybins<<", "<<DpT<<", "<<Nbins<<"\n";

  string origin[5] = {"PU","B","C","FRAG","GEANT"};
  int colors[5] = {2,4,8,12,44};
  TCanvas c("c", "canvas", 1300, 900);
  c.SetGrid();

  gPad->SetLogy();

  TH1D** px = new TH1D*[Nbins];
  for(int k=0;k<Nbins;k++){
    TString sx=Form("ht%d",k);
    px[k] = h->ProjectionX(sx,(int) DpT*k/2,(int) DpT*(k+1)/2);
    px[k]->SetStats(0);
    float N=px[k]->GetEntries();
    px[k]->Scale(1./N);
    std::string s = std::to_string(k);
    std::string s1 = std::to_string((int) DpT*k/2);
    std::string s2 = std::to_string((int) DpT*(k+1)/2);
    px[k]->SetTitle(("origin trk pT in ["+s1+","+s2+"] Gev").c_str());
//    px[k]->Draw("hist");
//    c.SaveAs(("/home/salomon/Private/atlas/FTPF/Selector/plots/origin/"+hist+"Xproj_"+s+".pdf").c_str());
    for(int i=1;i<=5;i++){
      M[k][i]=px[k]->GetBinContent(i);
      std::cout<<M[k][i]<<"\t";
    }
    std::cout<<"\n";
  }


  auto leg = new TLegend(0.2,0.2,0.35,0.4);//x1,y1,x2,y2
  TH1D** py_2 = new TH1D*[5];
  for(int i=1;i<=5;i++){
    py_2[i] = new TH1D(Form("syt%d",i),"",Nbins,0,pt_max);
    for(int k=0;k<Nbins;k++){
      py_2[i]->AddBinContent(k+1,M[k][i]);
      py_2[i]->SetStats(0);
      py_2[i]->GetYaxis()->SetRangeUser(1e-3, 1.);
      py_2[i]->GetXaxis()->SetTitle("trk pT [GeV]");
      py_2[i]->SetLineWidth(2);
      py_2[i]->SetLineColor(colors[i-1]);
      py_2[i]->SetOption("col");
      py_2[i]->SetMarkerStyle(5);
      py_2[i]->SetMarkerSize(2);
    }
    leg->AddEntry(py_2[i],origin[i-1].c_str(),"l");
    py_2[i]->Draw("same hist");
    py_2[i]->Draw("hist p same");
//    c.SaveAs(("/home/salomon/Private/atlas/FTPF/Selector/plots/origin/"+hist+"_"+origin[i-1]+"_binN.pdf").c_str());
  }
    leg->Draw();
    c.SaveAs((locForOutputFiles()+hist+"_binN.pdf").c_str());
    delete h,M,c,px,py_2;
}

void avtrk_orig(std::string hist,std::string hist_jet_pt,TFile* fData=_file0, int DpT=50, double pt_max=1000){

  TH2F* h = (TH2F*)fData->Get(hist.c_str());
  TH1F* h_jet = (TH1F*)fData->Get(hist_jet_pt.c_str());
  const Int_t n_Xbins = h->GetNbinsX();
  const Int_t n_Ybins = h->GetNbinsY();
//  const double pt_max = h->GetYaxis()->GetXmax();

  const Int_t n_bins_jet = h_jet->GetNbinsX();
  const double pt_max_jet = h_jet->GetXaxis()->GetXmax();

  if(n_bins_jet!=n_Ybins || pt_max!=pt_max_jet)  std::cout<<"ERROR\n";
//  int DpT=100;
  int rebin=(int) DpT/2;
  int Nbins=(int) pt_max/DpT;

  TCanvas c("c", "canvas", 1300, 900);
  c.SetGrid();
  gPad->SetLogy();

  TH1* jet_pt=h_jet->RebinX(rebin,"jet_pt");
//  jet_pt->Draw("hist");

//  c.SaveAs(("/home/salomon/Private/atlas/FTPF/Selector/plots/origin/"+hist_jet_pt+".pdf").c_str());

  double M[Nbins][5];

  std::cout<<n_Ybins<<", "<<pt_max<<", "<<Nbins<<"\n";
  std::cout<<n_bins_jet<<", "<<pt_max_jet<<", "<<Nbins<<"\n";

  string origin[5] = {"PU","B","C","FRAG","GEANT"};
  int colors[5] = {2,4,8,12,44};


  TH1D** px = new TH1D*[Nbins];

  for(int k=0;k<Nbins;k++){
    TString sx=Form("ha%d",k);
    px[k] = h->ProjectionX(sx,(int) DpT*k/2,(int) DpT*(k+1)/2);
    px[k]->SetStats(0);
//    float N=px[k]->GetEntries();
//    px[k]->Scale(1./N);
    std::string s = std::to_string(k);
    std::string s1 = std::to_string((int) DpT*k/2);
    std::string s2 = std::to_string((int) DpT*(k+1)/2);
    px[k]->SetTitle(("origin trk pT in ["+s1+","+s2+"] Gev").c_str());
//    px[k]->Draw("hist");
//    c.SaveAs(("/home/salomon/Private/atlas/FTPF/Selector/plots/origin/"+hist+"Xproj_"+s+".pdf").c_str());
    for(int i=1;i<=5;i++){
      M[k][i]=px[k]->GetBinContent(i);
      std::cout<<M[k][i]<<"\t";
    }
    std::cout<<"\n";
  }
  std::cout<<"\n";

  auto leg = new TLegend(0.75,0.75,0.9,0.9);//x1,y1,x2,y2
//  auto leg = new TLegend(0.75,0.15,0.85,0.35);//x1,y1,x2,y2
  TH1D** py_2 = new TH1D*[5];
  for(int i=1;i<=5;i++){
    py_2[i] = new TH1D(Form("sya%d",i),"",Nbins,0,pt_max);
    for(int k=0;k<Nbins;k++){
      py_2[i]->AddBinContent(k+1,M[k][i]);
    }
//    c.SaveAs(("/home/salomon/Private/atlas/FTPF/Selector/plots/origin/"+hist+"_"+origin[i-1]+"_binN.pdf").c_str());
  }

//  TH1D *h2 = new TH1D("h","efficiency_vs_jetpT",Nbins,0,pt_max);
//  h2=(TH1D*)jet_pt->Clone();
  TH1D** h2 = new TH1D*[5];
  double jetpt_bin=0;
  double avtrk=0;

  for(int i=1;i<=5;i++){
    h2[i]= new TH1D(Form("h2a%d",i),"",Nbins,0,pt_max);
    for(int k=0;k<Nbins;k++){
      jetpt_bin=(double) jet_pt->GetBinContent(k+1);
      avtrk= (double) M[k][i]/jetpt_bin;
      h2[i]->AddBinContent(k+1,avtrk);
//      std::cout<<M[k][i]<<"\t"<<jetpt_bin<<"\t"<<avtrk<<"\n";
    }
  }

  for(int i=1;i<=5;i++){
    h2[i]->SetStats(0);
    h2[i]->GetYaxis()->SetRangeUser(1e-1, 500);
    h2[i]->GetXaxis()->SetTitle("jet pT [GeV]");
    h2[i]->GetYaxis()->SetTitle("<trk>");
    h2[i]->SetLineWidth(2);
    h2[i]->SetLineColor(colors[i-1]);
    h2[i]->SetOption("col");
    h2[i]->SetMarkerStyle(5);
    h2[i]->SetMarkerSize(2);
    leg->AddEntry(h2[i],origin[i-1].c_str(),"l");
    h2[i]->Draw("same hist");
    h2[i]->Draw("hist p same");
  }

  leg->Draw();
  c.SaveAs((locForOutputFiles()+hist+"_"+hist_jet_pt+"_avtracks.pdf").c_str());
  delete h,h2,M,c,px,py_2;
}
