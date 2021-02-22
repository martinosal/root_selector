/*
root [1] .L ../plots/origin/scripts/truth_label_d0_origin.cpp
root [2] derivedtruthlabel("trk_origin_truth_label_B",_file0)
*/
void jet_composition(TFile* fData=_file0){
  TCanvas c("c", "canvas", 1300, 900);
  gStyle->SetOptStat(0);
  trk_pT_jet_DR_B->GetYaxis()->SetRangeUser(0, 0.7);
  trk_pT_jet_DR_B->GetXaxis()->SetRangeUser(0, 300.);
  trk_pT_jet_DR_B->Draw("COLZ");
  c.SaveAs("/home/salomon/Private/atlas/FTPF/Selector/plots/origin/trk_pT_jet_DR_B.pdf");
}

void jet_truthlabel_composition(std::string hist, TFile* fData=_file0){
  TH2F* h = (TH2F*)fData->Get(hist.c_str());//truth label origin
  const Int_t nbins_x = h->GetNbinsX();
  const Int_t nbins_y= h->GetNbinsY();
  const TAxis *axis_x=h->GetXaxis();
  const Int_t x_min=axis_x->GetXmin(); // should be -1 relative to i = 1;

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



  int y=0;

  for(int i=0;i<list.size();i++){
    y=0;
    for(std::vector<int>::iterator it = list.at(i).begin(); it != list.at(i).end(); ++it){
      for(int k=0;k<nbins_y;k++){
        y+=h->GetBinContent(k,*it+2);
      }
    }
    std::cout<<"derived truth label "<<i<<":\t"<<y<<"\n";
  }
}
