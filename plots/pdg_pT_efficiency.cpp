{
pt_max=40.;
TCanvas c1("c", "canvas", 1200, 900);
gStyle->SetOptStat(0);
TH1F *h1 = new TH1F("h","efficiency_vs_DR",100,0,2.);
h1=(TH1F*)matched_child_pi_pT_inB->Clone();
matched_child_pi_pT_inB->Sumw2();
child_pi_pT_inB->Sumw2();
h1->Divide(matched_child_pi_pT_inB, child_pi_pT_inB, 1., 1., "B");
h1->SetTitle("bTag_AntiKt4EMPFlowJets_201903: pi - efficiency_vs_pT");
h1->GetXaxis()->SetTitle("pT");
h1->GetYaxis()->SetTitle("efficiency");
h1->GetXaxis()->SetRangeUser(0., pt_max);
h1->GetYaxis()->SetRangeUser(0., 1.);
h1->Draw();

TCanvas c2("c", "canvas", 1200, 900);
gStyle->SetOptStat(0);
TH1F *h2 = new TH1F("h","efficiency_vs_pT",100,0,200.);
h2=(TH1F*)matched_child_K_pT_inB->Clone();
matched_child_K_pT_inB->Sumw2();
child_K_pT_inB->Sumw2();
h2->Divide(matched_child_K_pT_inB, child_K_pT_inB, 1., 1., "B");
h2->SetTitle("bTag_AntiKt4EMPFlowJets_201903: K - efficiency_vs_pT");
h2->GetXaxis()->SetTitle("pT");
h2->GetYaxis()->SetTitle("efficiency");
h2->GetXaxis()->SetRangeUser(0., pt_max);
h2->GetYaxis()->SetRangeUser(0., 1.);
h2->Draw();

TCanvas c3("c", "canvas", 1200, 900);
gStyle->SetOptStat(0);
TH1F *h3 = new TH1F("h","efficiency_vs_pT",100,0,200.);
h3=(TH1F*)matched_child_mu_pT_inB->Clone();
matched_child_mu_pT_inB->Sumw2();
child_mu_pT_inB->Sumw2();
h3->Divide(matched_child_mu_pT_inB, child_mu_pT_inB, 1., 1., "B");
h3->SetTitle("bTag_AntiKt4EMPFlowJets_201903: mu - efficiency_vs_pT");
h3->GetXaxis()->SetTitle("pT");
h3->GetYaxis()->SetTitle("efficiency");
h3->GetXaxis()->SetRangeUser(0., pt_max);
h3->GetYaxis()->SetRangeUser(0., 1.);
h3->Draw();

TCanvas c4("c", "canvas", 1200, 900);
gStyle->SetOptStat(0);
TH1F *h4 = new TH1F("h","efficiency_vs_pT",100,0,200.);
h4=(TH1F*)matched_child_p_pT_inB->Clone();
matched_child_p_pT_inB->Sumw2();
child_p_pT_inB->Sumw2();
h4->Divide(matched_child_p_pT_inB, child_p_pT_inB, 1., 1., "B");
h4->SetTitle("bTag_AntiKt4EMPFlowJets_201903: p - efficiency_vs_pT");
h4->GetXaxis()->SetTitle("pT");
h4->GetYaxis()->SetTitle("efficiency");
h4->GetXaxis()->SetRangeUser(0., pt_max);
h4->GetYaxis()->SetRangeUser(0., 1.);
h4->Draw();

TCanvas c5("c", "canvas", 1200, 900);
gStyle->SetOptStat(0);
TH1F *h5 = new TH1F("h","efficiency_vs_pT",100,0,200.);
h5=(TH1F*)matched_child_e_pT_inB->Clone();
matched_child_e_pT_inB->Sumw2();
child_e_pT_inB->Sumw2();
h5->Divide(matched_child_e_pT_inB, child_e_pT_inB, 1., 1., "B");
h5->SetTitle("bTag_AntiKt4EMPFlowJets_201903: e - efficiency_vs_pT");
h5->GetXaxis()->SetTitle("pT");
h5->GetYaxis()->SetTitle("efficiency");
h5->GetXaxis()->SetRangeUser(0., pt_max);
h5->GetYaxis()->SetRangeUser(0., 1.);
h5->Draw();

c1.SaveAs("../plots/pdgID/pi_eff.pdf");
c2.SaveAs("../plots/pdgID/K_eff.pdf");
c3.SaveAs("../plots/pdgID/mu_eff.pdf");
c4.SaveAs("../plots/pdgID/p_eff.pdf");
c5.SaveAs("../plots/pdgID/e_eff.pdf");
}
