{
int r=5;
TCanvas c1("c", "canvas", 1200, 900);
c1.Divide(2,2);
c1.cd(1);
TH1F *h1 = new TH1F("h1","1",100,0,10.);
h1=(TH1F*)matched_pTfraction_inB->Clone();
h1->Rebin(r);
TH1F *h1b = new TH1F("h1b","1b",100,0,10.);
h1b=(TH1F*)h1->Clone();
TH1F *h2 = new TH1F("h2","2",100,0,10.);
h2=(TH1F*)single_matched_pTfraction_inB->Clone();
h2->Rebin(r);
h1->Add(h2,-1);
h1->Divide(h1b);
h1->SetTitle("(matched-singlematched)/matched_vs_DpT/pT");
h1->GetXaxis()->SetTitle("DpT/pT");
h1->GetYaxis()->SetTitle("fraction");
h1->GetXaxis()->SetRangeUser(0., 10);
h1->GetYaxis()->SetRangeUser(0., 1.);
h1->Draw();

//TCanvas c2("c", "canvas", 1200, 900);
c1.cd(2);
TH1F *h3 = new TH1F("h3","h3",100,0,10.);
h3=(TH1F*)single_matched_pTfraction_inB->Clone();
h3->Rebin(r);
TH1F *h4 = new TH1F("h4","h4",100,0,10.);
h4=(TH1F*)matched_pTfraction_inB->Clone();
h4->Rebin(r);
h3->Divide(h4);
h3->SetTitle("singlematched/matched_vs_DpT/pT");
h3->GetXaxis()->SetTitle("DpT/pT");
h3->GetYaxis()->SetTitle("fraction");
h3->GetXaxis()->SetRangeUser(0., 10);
h3->GetYaxis()->SetRangeUser(0., 1.);
h3->Draw();    


//TCanvas c3("c", "canvas", 1200, 900);
c1.cd(3);
TH1F *h5 = new TH1F("h5","5",100,0,1.);
h5=(TH1F*)matched_DR_trk_pTfraction->Clone();
h5->Rebin(r);
TH1F *h5b = new TH1F("h5b","5b",100,0,1.);
h5b=(TH1F*)h5->Clone();
TH1F *h6 = new TH1F("h6","6",100,0,1.);
h6=(TH1F*)single_matched_DR_trk_inB->Clone();
h6->Rebin(r);
h5->Add(h6,-1);
h5->Divide(h5b);
h5->SetTitle("(matched-singlematched)/matched_vs_DR_child,track");
h5->GetXaxis()->SetTitle("DR_child,track");
h5->GetYaxis()->SetTitle("fraction");
h5->GetXaxis()->SetRangeUser(0., 1.);
h5->GetYaxis()->SetRangeUser(0., 1.);
h5->Draw();

//TCanvas c4("c", "canvas", 1200, 900);
c1.cd(4);
TH1F *h7 = new TH1F("h7","h7",100,0,1.);
h7=(TH1F*)single_matched_DRfraction_inB->Clone();
h7->Rebin(r);
TH1F *h8 = new TH1F("h8","h8",100,0,1.);
h8=(TH1F*)matched_DRfraction_inB->Clone();
h8->Rebin(r);
h7->Divide(h8);
h7->SetTitle("singlematched/matched_vs_DR_child,track");
h7->GetXaxis()->SetTitle("DR_child,track");
h7->GetYaxis()->SetTitle("fraction");
h7->GetXaxis()->SetRangeUser(0., 1.);
h7->GetYaxis()->SetRangeUser(0., 1.);
h7->Draw();    

}
