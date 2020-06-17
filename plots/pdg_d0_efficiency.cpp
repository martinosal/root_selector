{
d0_min=-2.;d0_max=4.;
TCanvas c1("c", "canvas", 1200, 900);
gStyle->SetOptStat(0);
TH1F *h1 = new TH1F("h","efficiency_vs_DR",100,0,2.);
h1=(TH1F*)matched_child_pi_d0_truth_inB->Clone();
matched_child_pi_d0_truth_inB->Sumw2();
child_pi_d0_truth_inB->Sumw2();
h1->Divide(matched_child_pi_d0_truth_inB, child_pi_d0_truth_inB, 1., 1., "B");
h1->SetTitle("bTag_AntiKt4EMPFlowJets_201903: pi - efficiency_vs_d0");
h1->GetXaxis()->SetTitle("d0");
h1->GetYaxis()->SetTitle("efficiency");
h1->GetXaxis()->SetRangeUser(d0_min, d0_max);
h1->GetYaxis()->SetRangeUser(d0_min, 1.);
h1->Draw();

TCanvas c2("c", "canvas", 1200, 900);
gStyle->SetOptStat(0);
TH1F *h2 = new TH1F("h","efficiency_vs_d0",100,0,200.);
h2=(TH1F*)matched_child_K_d0_truth_inB->Clone();
matched_child_K_d0_truth_inB->Sumw2();
child_K_d0_truth_inB->Sumw2();
h2->Divide(matched_child_K_d0_truth_inB, child_K_d0_truth_inB, 1., 1., "B");
h2->SetTitle("bTag_AntiKt4EMPFlowJets_201903: K - efficiency_vs_d0");
h2->GetXaxis()->SetTitle("d0");
h2->GetYaxis()->SetTitle("efficiency");
h2->GetXaxis()->SetRangeUser(d0_min, d0_max);
h2->GetYaxis()->SetRangeUser(d0_min, 1.);
h2->Draw();

TCanvas c3("c", "canvas", 1200, 900);
gStyle->SetOptStat(0);
TH1F *h3 = new TH1F("h","efficiency_vs_d0",100,0,200.);
h3=(TH1F*)matched_child_mu_d0_truth_inB->Clone();
matched_child_mu_d0_truth_inB->Sumw2();
child_mu_d0_truth_inB->Sumw2();
h3->Divide(matched_child_mu_d0_truth_inB, child_mu_d0_truth_inB, 1., 1., "B");
h3->SetTitle("bTag_AntiKt4EMPFlowJets_201903: mu - efficiency_vs_d0");
h3->GetXaxis()->SetTitle("d0");
h3->GetYaxis()->SetTitle("efficiency");
h3->GetXaxis()->SetRangeUser(d0_min, d0_max);
h3->GetYaxis()->SetRangeUser(0., 1.);
h3->Draw();

TCanvas c4("c", "canvas", 1200, 900);
gStyle->SetOptStat(0);
TH1F *h4 = new TH1F("h","efficiency_vs_d0",100,0,200.);
h4=(TH1F*)matched_child_p_d0_truth_inB->Clone();
matched_child_p_d0_truth_inB->Sumw2();
child_p_d0_truth_inB->Sumw2();
h4->Divide(matched_child_p_d0_truth_inB, child_p_d0_truth_inB, 1., 1., "B");
h4->SetTitle("bTag_AntiKt4EMPFlowJets_201903: p - efficiency_vs_d0");
h4->GetXaxis()->SetTitle("d0");
h4->GetYaxis()->SetTitle("efficiency");
h4->GetXaxis()->SetRangeUser(d0_min, d0_max);
h4->GetYaxis()->SetRangeUser(0., 1.);
h4->Draw();

TCanvas c5("c", "canvas", 1200, 900);
gStyle->SetOptStat(0);
TH1F *h5 = new TH1F("h","efficiency_vs_d0",100,0,200.);
h5=(TH1F*)matched_child_e_d0_truth_inB->Clone();
matched_child_e_d0_truth_inB->Sumw2();
child_e_d0_truth_inB->Sumw2();
h5->Divide(matched_child_e_d0_truth_inB, child_e_d0_truth_inB, 1., 1., "B");
h5->SetTitle("bTag_AntiKt4EMPFlowJets_201903: e - efficiency_vs_d0");
h5->GetXaxis()->SetTitle("d0");
h5->GetYaxis()->SetTitle("efficiency");
h5->GetXaxis()->SetRangeUser(d0_min, d0_max);
h5->GetYaxis()->SetRangeUser(0., 1.);
h5->Draw();

c1.SaveAs("../plots/pdgID/pi_d0_eff.pdf");
c2.SaveAs("../plots/pdgID/K_d0_eff.pdf");
c3.SaveAs("../plots/pdgID/mu_d0_eff.pdf");
c4.SaveAs("../plots/pdgID/p_d0_eff.pdf");
c5.SaveAs("../plots/pdgID/e_d0_eff.pdf");
}
