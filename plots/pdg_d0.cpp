{
  double a=-2.,b=4;
TCanvas b1("c", "canvas", 1200, 900);
//b1.Setlogy();
gPad->SetLogy();
gStyle->SetOptStat(1);
child_pi_d0_truth_inB->SetTitle("bTag_AntiKt4EMPFlowJets_201903: pi");
child_pi_d0_truth_inB->GetXaxis()->SetTitle("d0");
child_pi_d0_truth_inB->GetXaxis()->SetRangeUser(a, b);
child_pi_d0_truth_inB->Draw();

TCanvas b2("c", "canvas", 1200, 900);
//b2.Setlogy();
gPad->SetLogy();
gStyle->SetOptStat(1);
child_mu_d0_truth_inB->SetTitle("bTag_AntiKt4EMPFlowJets_201903: mu");
child_mu_d0_truth_inB->GetXaxis()->SetTitle("d0");
child_mu_d0_truth_inB->GetXaxis()->SetRangeUser(a, b);
child_mu_d0_truth_inB->Draw();

TCanvas b3("c", "canvas", 1200, 900);
//b3.Setlogy();
gPad->SetLogy();
gStyle->SetOptStat(1);
child_K_d0_truth_inB->SetTitle("bTag_AntiKt4EMPFlowJets_201903: K");
child_K_d0_truth_inB->GetXaxis()->SetTitle("d0");
child_K_d0_truth_inB->GetXaxis()->SetRangeUser(a, b);
child_K_d0_truth_inB->Draw();

TCanvas b4("c", "canvas", 1200, 900);
//b4.Setlogy();
gPad->SetLogy();
gStyle->SetOptStat(1);
child_p_d0_truth_inB->SetTitle("bTag_AntiKt4EMPFlowJets_201903: proton");
child_p_d0_truth_inB->GetXaxis()->SetTitle("d0");
child_p_d0_truth_inB->GetXaxis()->SetRangeUser(a, b);
child_p_d0_truth_inB->Draw();

TCanvas b5("c", "canvas", 1200, 900);
//b4.Setlogy();
gPad->SetLogy();
gStyle->SetOptStat(1);
child_e_d0_truth_inB->SetTitle("bTag_AntiKt4EMPFlowJets_201903: electron");
child_e_d0_truth_inB->GetXaxis()->SetTitle("d0");
child_e_d0_truth_inB->GetXaxis()->SetRangeUser(a, b);
child_e_d0_truth_inB->Draw();

TCanvas a1("c", "canvas", 1200, 900);
//b1.Setlogy();
gPad->SetLogy();
gStyle->SetOptStat(1);
TH1F *h1 = new TH1F("h","h",300,0,300.);
h1=(TH1F*)child_pi_d0_truth_inB->Clone();
h1->Add(matched_child_pi_d0_truth_inB,-1);
h1->SetTitle("bTag_AntiKt4EMPFlowJets_201903: not-matched_pi");
h1->GetXaxis()->SetTitle("d0");
h1->GetXaxis()->SetRangeUser(a, b);
h1->Draw();

TCanvas a2("c", "canvas", 1200, 900);
//b2.Setlogy();
gPad->SetLogy();
gStyle->SetOptStat(1);
TH1F *h2 = new TH1F("h","h",300,0,300.);
h2=(TH1F*)child_mu_d0_truth_inB->Clone();
h2->Add(matched_child_mu_d0_truth_inB,-1);
h2->SetTitle("bTag_AntiKt4EMPFlowJets_201903: not-matched_mu");
h2->GetXaxis()->SetTitle("d0");
h2->GetXaxis()->SetRangeUser(a, b);
h2->Draw();

TCanvas a3("c", "canvas", 1200, 900);
//b3.Setlogy();
gPad->SetLogy();
gStyle->SetOptStat(1);
TH1F *h3 = new TH1F("h","h",300,0,300.);
h3=(TH1F*)child_K_d0_truth_inB->Clone();
h3->Add(matched_child_K_d0_truth_inB,-1);
h3->SetTitle("bTag_AntiKt4EMPFlowJets_201903: not-matched_K");
h3->GetXaxis()->SetTitle("d0");
h3->GetXaxis()->SetRangeUser(a, b);
h3->Draw();

TCanvas a4("c", "canvas", 1200, 900);
//b4.Setlogy();
gPad->SetLogy();
gStyle->SetOptStat(1);
TH1F *h4 = new TH1F("h","h",300,0,300.);
h4=(TH1F*)child_p_d0_truth_inB->Clone();
h4->Add(matched_child_p_d0_truth_inB,-1);
h4->SetTitle("bTag_AntiKt4EMPFlowJets_201903: not-matched_p");
h4->GetXaxis()->SetTitle("d0");
h4->GetXaxis()->SetRangeUser(a, b);
h4->Draw();

TCanvas a5("c", "canvas", 1200, 900);
//b4.Setlogy();
gPad->SetLogy();
gStyle->SetOptStat(1);
TH1F *h5 = new TH1F("h","h",300,0,300.);
h5=(TH1F*)child_e_d0_truth_inB->Clone();
h5->Add(matched_child_e_d0_truth_inB,-1);
h5->SetTitle("bTag_AntiKt4EMPFlowJets_201903: not-matched_e");
h5->GetXaxis()->SetTitle("d0");
h5->GetXaxis()->SetRangeUser(a, b);
h5->Draw();


b1.SaveAs("../plots/pdgID/pi_d0.pdf");
b2.SaveAs("../plots/pdgID/mu_d0.pdf");
b3.SaveAs("../plots/pdgID/K_d0.pdf");
b4.SaveAs("../plots/pdgID/p_d0.pdf");
b5.SaveAs("../plots/pdgID/e_d0.pdf");

a1.SaveAs("../plots/pdgID/notmatched_pi_d0.pdf");
a2.SaveAs("../plots/pdgID/notmatched_mu_d0.pdf");
a3.SaveAs("../plots/pdgID/notmatched_K_d0.pdf");
a4.SaveAs("../plots/pdgID/notmatched_p_d0.pdf");
a5.SaveAs("../plots/pdgID/notmatched_e_d0.pdf");


}
