
void assoc_tracks(std::string hist,std::string hist_pt,TFile* fData=_file0, string flag="jet"){
  TH2F* h = (TH2F*)fData->Get(hist.c_str());//truth label origin
  //hist=n_tracks_jetpt_B,n_tracks_jetpt_C,n_tracks_bHpt_B,n_tracks_cHpt_C,
  TH1F* h_pt = (TH1F*)fData->Get(hist_pt.c_str());

  h->RebinX(25);
  h_pt->RebinX(25);
  const Int_t n_binsX_h=h->GetNbinsX();
  const Int_t n_binsY_h=h->GetNbinsY();
  const Int_t n_binsX_hpt=h_pt->GetNbinsX();
  if(n_binsX_h!=n_binsX_hpt)  std::cout<<"Warning\n";


  TH1D* h0 = new TH1D("h0","h0",n_binsX_h,0.,1000.);
  float s=0;
  for(unsigned i=1;i<n_binsX_h;i++){
    s=0;
    for(unsigned j=1;j<n_binsY_h;j++){
      s=s+j*double(h->GetBinContent(i,j));
    }
    h0->AddBinContent(i,s);

    std::cout<<i<<"\t"<<s<<"\t"<<h_pt->GetBinContent(i+1)<<std::endl;
  }



//  TH1D* h_proj=h->ProjectionX("h_projected");

  TH1D* h2;
  h2=(TH1D*)h0->Clone();
  h2->Divide(h_pt);



  TCanvas c("c", "canvas", 1300, 900);
  c.SetGrid();
  gPad->SetLogy();
//  h_proj->SetTitle("");
  if(flag=="jet") h2->GetXaxis()->SetTitle("jet pT [GeV]");
  if(flag=="bH") h2->GetXaxis()->SetTitle("bH pT [GeV]");
  if(flag=="cH") h2->GetXaxis()->SetTitle("cH pT [GeV]");
  h2->SetStats(0);
  h2->GetYaxis()->SetTitle("average number of assoc. tracks");
  h2->GetXaxis()->SetRangeUser(0., 300.);
  h2->GetYaxis()->SetRangeUser(1e-1, 100.);
//  h_proj->GetYaxis()->SetRangeUser(0., 0.3);

  h2->Draw();
  c.SaveAs(("/home/salomon/Private/atlas/FTPF/Selector/plots/origin/tmp_QT/"+hist+"_assoctrks.pdf").c_str());

}
