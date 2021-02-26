std::string locForOutputFiles()
{
  std::string myPath = "./plots/";
  return myPath;
}

void eff_DR(string flag="jet",string flav="b"){
//TTree* tree = (TTree*) file.Get("treename");
//TCanvas c1("c1", "canvas", 1300, 900);
  std::string num="tmp_num";
  std::string den="tmp_den";

  if(flag=="jet" && flav=="b"){
    num="matched_origin_trk_DR_jetpt_B";//500, 0., 1000.,200,0.,2.
    den="child_DR_jetpt_B";//500, 0., 1000.,200,0.,2.
  }
  if(flag=="jet" && flav=="c"){
    num="matched_origin_trk_DR_jetpt_C";//500, 0., 1000.,200,0.,2.
    den="child_DR_jetpt_C";//500, 0., 1000.,200,0.,2.
  }
  if(flag=="bH" && flav=="b"){
    num="matched_origin_trk_DR_bHpt_B";//500, 0., 1000.,200,0.,2.
    den="child_DR_bHpt_B";//500, 0., 1000.,200,0.,2.
  }
  if(flag=="cH" && flav=="c"){
    num="matched_origin_trk_DR_cHpt_C";//500, 0., 1000.,200,0.,2.
    den="child_DR_cHpt_C";//500, 0., 1000.,200,0.,2.
  }
  TH2F* h_num = (TH2F*)_file0->Get(num.c_str());
  TH2F* h_den = (TH2F*)_file0->Get(den.c_str());

//  TH1D* h_num_1 = h_num->ProjectionX();
//  TH1D* h_den_1 = h_den->ProjectionX();
  h_num->RebinX(25);
  h_num->RebinY(5);
  h_den->RebinX(25);
  h_den->RebinY(5);

  TH2F* h=(TH2F*)h_num->Clone();
//  h_num->Sumw2();
//  h_den->Sumw2();
  h->Divide(h_den);

  const Int_t n_binsX_h=h_num->GetNbinsX();
  const Int_t n_binsY_h=h_num->GetNbinsY();
  float cont=0.,num_cont=0.,den_cont=0.;
  for(unsigned i=1;i<n_binsX_h;i++){
    for(unsigned j=1;j<n_binsY_h;j++){
      cont=h->GetBinContent(i,j);
      num_cont=h_num->GetBinContent(i,j);
      den_cont=h_den->GetBinContent(i,j);
      if(den_cont<30 && den_cont>0){
        std::cout<<i<<","<<j<<"\t"<<cont<<"\t"<<den_cont<<std::endl;
        h->SetBinContent(i,j,0);
      }
    }
  }


  TCanvas*  c = new TCanvas("c1", "canvas", 1300, 900);
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);

  //h->SetTitle("efficiency_jet_pT");
  if(flag=="jet") h->GetXaxis()->SetTitle("jet pT [GeV]");
  if(flag=="bH")  h->GetXaxis()->SetTitle("bH pT [GeV]");
  if(flag=="cH")  h->GetXaxis()->SetTitle("cH pT [GeV]");
  if(flav=="b")  h->GetYaxis()->SetTitle("DR_bH,jet");
  if(flav=="c")  h->GetYaxis()->SetTitle("DR_cH,jet");

  h->GetXaxis()->SetRangeUser(0., 300.);
  h->GetYaxis()->SetRangeUser(0., 0.3);
  h->GetZaxis()->SetRangeUser(0., 1.);
  h->Draw("COLZ");

  std::string path =locForOutputFiles();//"/home/salomon/Private/atlas/FTPF/Selector/plots/origin/";
  std::cout<<" path = "<<path<<std::endl;

  c->SaveAs((path+num+den+".pdf").c_str());

}

void eff_pt(string flag="jet",string flav="b"){
//TTree* tree = (TTree*) file.Get("treename");
//TCanvas c1("c1", "canvas", 1300, 900);
  std::string num="tmp_num";
  std::string den="tmp_den";

  if(flag=="jet" && flav=="b"){
    num="matched_origin_trk_jetpT_B";//500, 0., 1000.,200,0.,2.
    den="child_jetpT_B";//500, 0., 1000.,200,0.,2.
  }
  if(flag=="jet" && flav=="c"){
    num="matched_origin_trk_jetpT_C";//500, 0., 1000.,200,0.,2.
    den="child_jetpT_C";//500, 0., 1000.,200,0.,2.
  }
  if(flag=="bH" && flav=="b"){
    num="matched_origin_trk_bHpT_B";//500, 0., 1000.,200,0.,2.
    den="child_bHpT_B";//500, 0., 1000.,200,0.,2.
  }
  if(flag=="cH" && flav=="c"){
    num="matched_origin_trk_cHpT_C";//500, 0., 1000.,200,0.,2.
    den="child_cHpT_C";//500, 0., 1000.,200,0.,2.
  }

  TH1F* h_num = (TH1F*)_file0->Get(num.c_str());
  TH1F* h_den = (TH1F*)_file0->Get(den.c_str());

//  TH1D* h_num_1 = h_num->ProjectionX();
//  TH1D* h_den_1 = h_den->ProjectionX();
  h_num->Rebin(25);
  h_den->Rebin(25);
  TH1D* h=(TH1D*)h_num->Clone();
  h_num->Sumw2();
  h_den->Sumw2();
  h->Divide(h_num, h_den, 1., 1., "B");

  TCanvas* c = new TCanvas("c1", "canvas", 1300, 900);
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);

//  h->SetTitle("efficiency_jet_pT");
  if(flag=="jet") h->GetXaxis()->SetTitle("jet pT [GeV]");
  if(flag=="bH")  h->GetXaxis()->SetTitle("bH pT [GeV]");
  if(flag=="cH")  h->GetXaxis()->SetTitle("cH pT [GeV]");
  h->GetYaxis()->SetTitle("efficiency");
  h->GetXaxis()->SetRangeUser(0., 300.);
  h->GetYaxis()->SetRangeUser(0., 1.);
  h->Draw();

  std::string path =locForOutputFiles();//"/home/salomon/Private/atlas/FTPF/Selector/plots/origin/";
  std::cout<<" path = "<<path<<std::endl;
  c->SaveAs((path+num+den+".pdf").c_str());

}
