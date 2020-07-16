void sig(std::string hist,TFile* fData=_file0){

//  TH1F* h = (TH1F*)fData->Get(hist.c_str());
  TH1F* h = (TH1F*)fData->Get(("trk_"+hist+"_B").c_str());
  TCanvas c("c", "canvas", 1300, 900);

  gPad->SetLogy();
  float N=h->GetEntries();
  h->Scale(1./N);
  h->SetStats(0);
  h->SetTitle((hist).c_str());
  h->Draw("hist");

  c.SaveAs(("/home/salomon/Private/atlas/FTPF/Selector/plots/scripts_root/detectorPlots/"+hist+".pdf").c_str());
  delete h;
}

void orig(std::string hist,TFile* fData=_file0){

  TH1F* h = (TH1F*)fData->Get(("trk_"+hist+"_origin_B").c_str());
  const Int_t n_bins = h->GetNbinsX();
  string origin[5] = {"PU","B","C","FRAG","GEANT"};
  TCanvas c("c", "canvas", 1300, 900);

  gPad->SetLogy();
  float N=h->GetEntries();
  h->Scale(1./N);
//  auto h2 = (TH1F*)h->Clone();

  for(int i=1;i<=n_bins;i++){
    h->GetXaxis()->SetBinLabel(i,origin[i-1].c_str());
  }
  h->GetYaxis()->SetRangeUser(1e-2, 1.);
  gStyle->SetPaintTextFormat(".4f");
  h->SetStats(0);
  h->SetTitle((hist).c_str());
  h->Draw("hist TEXT0 min0");

  c.SaveAs(("/home/salomon/Private/atlas/FTPF/Selector/plots/scripts_root/detectorPlots/"+hist+"_origin.pdf").c_str());
  delete h;
}
