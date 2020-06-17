{
TCanvas c1("c", "canvas", 1200, 900);
c1.Divide(2,2);
c1.cd(1);
TH1F *h1 = new TH1F("h","efficiency_vs_DR",100,0,2.);
h1=(TH1F*)matched_child_DR_inB->Clone();
matched_child_DR_inB->Sumw2();
child_DR_inB->Sumw2();
h1->Divide(matched_child_DR_inB, child_DR_inB, 1., 1., "B");
h1->SetTitle("efficiency_vs_DR");
h1->GetXaxis()->SetTitle("Delta R_child,jet");
h1->GetYaxis()->SetTitle("efficiency");
h1->GetXaxis()->SetRangeUser(0., 0.6);
h1->GetYaxis()->SetRangeUser(0., 1.1);
h1->Draw();

c1.cd(2);
TH1F *h2 = new TH1F("h","efficiency_vs_pT",100,0,200.);
h2=(TH1F*)matched_child_pT_inB->Clone();
matched_child_pT_inB->Sumw2();
child_pT_inB->Sumw2();
h2->Divide(matched_child_pT_inB,child_pT_inB,1., 1., "B");
h2->SetTitle("efficiency_vs_pT");
h2->GetXaxis()->SetTitle("pT_child [GeV]");
h2->GetYaxis()->SetTitle("efficiency");
h2->GetXaxis()->SetRangeUser(0., 100.);
h2->GetYaxis()->SetRangeUser(0., 1.1);
h2->Draw();

c1.cd(3);
TH1F *h3 = new TH1F("h","efficiency_vs_eta",100,0,2.);
h3=(TH1F*)matched_child_eta_inB->Clone();
matched_child_eta_inB->Sumw2();
child_eta_inB->Sumw2();
h3->Divide(matched_child_eta_inB,child_eta_inB,1., 1., "B");
h3->SetTitle("efficiency_vs_eta");
h3->GetXaxis()->SetTitle("eta_child");
h3->GetYaxis()->SetTitle("efficiency");
h3->GetXaxis()->SetRangeUser(-2.6, 2.6);
h3->GetYaxis()->SetRangeUser(0., 1.1);
h3->Draw();

c1.cd(4);
TH1F *h4 = new TH1F("h","efficiency_vs_phi",100,0,200.);
h4=(TH1F*)matched_child_phi_inB->Clone();
matched_child_phi_inB->Sumw2();
child_phi_inB->Sumw2();
h4->Divide(matched_child_phi_inB,child_phi_inB,1., 1., "B");
h4->SetTitle("efficiency_vs_phi");
h4->GetXaxis()->SetTitle("phi_child");
h4->GetYaxis()->SetTitle("efficiency");
h4->GetXaxis()->SetRangeUser(-4., 4.);
h4->GetYaxis()->SetRangeUser(0., 1.1);
h4->Draw();
}
