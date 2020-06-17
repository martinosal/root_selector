{
Lxy_max=40.;
TCanvas c1("c", "canvas", 1200, 900);
gStyle->SetOptStat(0);
TH1F *h1 = new TH1F("h","efficiency_vs_DR",100,0,2.);
h1=(TH1F*)matched_child_pi_Lxy_inB->Clone();
matched_child_pi_Lxy_inB->Sumw2();
child_pi_Lxy_inB->Sumw2();
h1->Divide(matched_child_pi_Lxy_inB, child_pi_Lxy_inB, 1., 1., "B");
h1->SetTitle("AktVR30Rmax4Rmin02Tr_BTagging201903: pi - efficiency_vs_Lxy");
h1->GetXaxis()->SetTitle("Lxy");
h1->GetYaxis()->SetTitle("efficiency");
h1->GetXaxis()->SetRangeUser(0., Lxy_max);
h1->GetYaxis()->SetRangeUser(0., 1.);
h1->Draw();

TCanvas c2("c", "canvas", 1200, 900);
gStyle->SetOptStat(0);
TH1F *h2 = new TH1F("h","efficiency_vs_Lxy",100,0,200.);
h2=(TH1F*)matched_child_K_Lxy_inB->Clone();
matched_child_K_Lxy_inB->Sumw2();
child_K_Lxy_inB->Sumw2();
h2->Divide(matched_child_K_Lxy_inB, child_K_Lxy_inB, 1., 1., "B");
h2->SetTitle("AktVR30Rmax4Rmin02Tr_BTagging201903: K - efficiency_vs_Lxy");
h2->GetXaxis()->SetTitle("Lxy");
h2->GetYaxis()->SetTitle("efficiency");
h2->GetXaxis()->SetRangeUser(0., Lxy_max);
h2->GetYaxis()->SetRangeUser(0., 1.);
h2->Draw();

TCanvas c3("c", "canvas", 1200, 900);
gStyle->SetOptStat(0);
TH1F *h3 = new TH1F("h","efficiency_vs_Lxy",100,0,200.);
h3=(TH1F*)matched_child_mu_Lxy_inB->Clone();
matched_child_mu_Lxy_inB->Sumw2();
child_mu_Lxy_inB->Sumw2();
h3->Divide(matched_child_mu_Lxy_inB, child_mu_Lxy_inB, 1., 1., "B");
h3->SetTitle("AktVR30Rmax4Rmin02Tr_BTagging201903: mu - efficiency_vs_Lxy");
h3->GetXaxis()->SetTitle("Lxy");
h3->GetYaxis()->SetTitle("efficiency");
h3->GetXaxis()->SetRangeUser(0., Lxy_max);
h3->GetYaxis()->SetRangeUser(0., 1.);
h3->Draw();

TCanvas c4("c", "canvas", 1200, 900);
gStyle->SetOptStat(0);
TH1F *h4 = new TH1F("h","efficiency_vs_Lxy",100,0,200.);
h4=(TH1F*)matched_child_p_Lxy_inB->Clone();
matched_child_p_Lxy_inB->Sumw2();
child_p_Lxy_inB->Sumw2();
h4->Divide(matched_child_p_Lxy_inB, child_p_Lxy_inB, 1., 1., "B");
h4->SetTitle("AktVR30Rmax4Rmin02Tr_BTagging201903: p - efficiency_vs_Lxy");
h4->GetXaxis()->SetTitle("Lxy");
h4->GetYaxis()->SetTitle("efficiency");
h4->GetXaxis()->SetRangeUser(0., Lxy_max);
h4->GetYaxis()->SetRangeUser(0., 1.);
h4->Draw();

TCanvas c5("c", "canvas", 1200, 900);
gStyle->SetOptStat(0);
TH1F *h5 = new TH1F("h","efficiency_vs_Lxy",100,0,200.);
h5=(TH1F*)matched_child_e_Lxy_inB->Clone();
matched_child_e_Lxy_inB->Sumw2();
child_e_Lxy_inB->Sumw2();
h5->Divide(matched_child_e_Lxy_inB, child_e_Lxy_inB, 1., 1., "B");
h5->SetTitle("AktVR30Rmax4Rmin02Tr_BTagging201903: e - efficiency_vs_Lxy");
h5->GetXaxis()->SetTitle("Lxy");
h5->GetYaxis()->SetTitle("efficiency");
h5->GetXaxis()->SetRangeUser(0., Lxy_max);
h5->GetYaxis()->SetRangeUser(0., 1.);
h5->Draw();

c1.SaveAs("../plots/pdgID/pi_Lxy_eff.pdf");
c2.SaveAs("../plots/pdgID/K_Lxy_eff.pdf");
c3.SaveAs("../plots/pdgID/mu_Lxy_eff.pdf");
c4.SaveAs("../plots/pdgID/p_Lxy_eff.pdf");
c4.SaveAs("../plots/pdgID/e_Lxy_eff.pdf");
}
