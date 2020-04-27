{
TCanvas c1("c", "canvas", 1200, 900);
TH1F *h1 = new TH1F("h","efficiency_vs_DR",100,0,2.);
h1=(TH1F*)matched_DR_400_inB->Clone();
h1->Divide(child_DR_400_inB);
h1->SetTitle("efficiency_vs_DR");
h1->GetXaxis()->SetTitle("DeltaR");
h1->GetYaxis()->SetTitle("efficiency");
h1->GetXaxis()->SetRangeUser(0., 0.6);
h1->GetYaxis()->SetRangeUser(0., 1.1);
h1->Draw();

TCanvas c2("c", "canvas", 1200, 900);
TH1F *h2 = new TH1F("h","efficiency_vs_pT",100,0,200.);
h2=(TH1F*)matched_pT_400_inB->Clone();
h2->Divide(child_pT_400_inB);
h2->SetTitle("efficiency_vs_pT");
h2->GetXaxis()->SetTitle("pT [GeV]");
h2->GetYaxis()->SetTitle("efficiency");
h2->GetXaxis()->SetRangeUser(0., 100.);
h2->GetYaxis()->SetRangeUser(0., 1.1);
h2->Draw();

TCanvas c3("c", "canvas", 1200, 900);
TH1F *h3 = new TH1F("h","efficiency_vs_eta",100,0,2.);
h3=(TH1F*)matched_eta_400_inB->Clone();
h3->Divide(child_eta_400_inB);
h3->SetTitle("efficiency_vs_eta");
h3->GetXaxis()->SetTitle("eta");
h3->GetYaxis()->SetTitle("efficiency");
h3->GetXaxis()->SetRangeUser(-2.6, 2.6);
h3->GetYaxis()->SetRangeUser(0., 1.1);
h3->Draw();

TCanvas c4("c", "canvas", 1200, 900);
TH1F *h4 = new TH1F("h","efficiency_vs_phi",100,0,200.);
h4=(TH1F*)matched_phi_400_inB->Clone();
h4->Divide(child_phi_400_inB);
h4->SetTitle("efficiency_vs_phi");
h4->GetXaxis()->SetTitle("phi");
h4->GetYaxis()->SetTitle("efficiency");
h4->GetXaxis()->SetRangeUser(-4., 4.);
h4->GetYaxis()->SetRangeUser(0., 1.1);
h4->Draw();
}
