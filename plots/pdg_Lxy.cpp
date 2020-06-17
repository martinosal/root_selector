{
  double a=0.,b=40;
TCanvas b1("c", "canvas", 1200, 900);
//b1.Setlogy();
gPad->SetLogy();
gStyle->SetOptStat(1);
child_pi_Lxy_inB->SetTitle("pi");
child_pi_Lxy_inB->GetXaxis()->SetTitle("Lxy");
child_pi_Lxy_inB->GetXaxis()->SetRangeUser(a, b);
child_pi_Lxy_inB->Draw();

TCanvas b2("c", "canvas", 1200, 900);
//b2.Setlogy();
gPad->SetLogy();
gStyle->SetOptStat(1);
child_mu_Lxy_inB->SetTitle("mu");
child_mu_Lxy_inB->GetXaxis()->SetTitle("Lxy");
child_mu_Lxy_inB->GetXaxis()->SetRangeUser(a, b);
child_mu_Lxy_inB->Draw();

TCanvas b3("c", "canvas", 1200, 900);
//b3.Setlogy();
gPad->SetLogy();
gStyle->SetOptStat(1);
child_K_Lxy_inB->SetTitle("K");
child_K_Lxy_inB->GetXaxis()->SetTitle("Lxy");
child_K_Lxy_inB->GetXaxis()->SetRangeUser(a, b);
child_K_Lxy_inB->Draw();

TCanvas b4("c", "canvas", 1200, 900);
//b4.Setlogy();
gPad->SetLogy();
gStyle->SetOptStat(1);
child_p_Lxy_inB->SetTitle("proton");
child_p_Lxy_inB->GetXaxis()->SetTitle("Lxy");
child_p_Lxy_inB->GetXaxis()->SetRangeUser(a, b);
child_p_Lxy_inB->Draw();

TCanvas b5("c", "canvas", 1200, 900);
//b4.Setlogy();
gPad->SetLogy();
gStyle->SetOptStat(1);
child_e_Lxy_inB->SetTitle("electron");
child_e_Lxy_inB->GetXaxis()->SetTitle("Lxy");
child_e_Lxy_inB->GetXaxis()->SetRangeUser(a, b);
child_e_Lxy_inB->Draw();

TCanvas a1("c", "canvas", 1200, 900);
//b1.Setlogy();
gPad->SetLogy();
gStyle->SetOptStat(1);
TH1F *h1 = new TH1F("h","h",300,0,300.);
h1=(TH1F*)child_pi_Lxy_inB->Clone();
h1->Add(matched_child_pi_Lxy_inB,-1);
h1->SetTitle("not-matched_pi");
h1->GetXaxis()->SetTitle("Lxy");
h1->GetXaxis()->SetRangeUser(a, b);
h1->Draw();

TCanvas a2("c", "canvas", 1200, 900);
//b2.Setlogy();
gPad->SetLogy();
gStyle->SetOptStat(1);
TH1F *h2 = new TH1F("h","h",300,0,300.);
h2=(TH1F*)child_mu_Lxy_inB->Clone();
h2->Add(matched_child_mu_Lxy_inB,-1);
h2->SetTitle("not-matched_mu");
h2->GetXaxis()->SetTitle("Lxy");
h2->GetXaxis()->SetRangeUser(a, b);
h2->Draw();

TCanvas a3("c", "canvas", 1200, 900);
//b3.Setlogy();
gPad->SetLogy();
gStyle->SetOptStat(1);
TH1F *h3 = new TH1F("h","h",300,0,300.);
h3=(TH1F*)child_K_Lxy_inB->Clone();
h3->Add(matched_child_K_Lxy_inB,-1);
h3->SetTitle("not-matched_K");
h3->GetXaxis()->SetTitle("Lxy");
h3->GetXaxis()->SetRangeUser(a, b);
h3->Draw();

TCanvas a4("c", "canvas", 1200, 900);
//b4.Setlogy();
gPad->SetLogy();
gStyle->SetOptStat(1);
TH1F *h4 = new TH1F("h","h",300,0,300.);
h4=(TH1F*)child_p_Lxy_inB->Clone();
h4->Add(matched_child_p_Lxy_inB,-1);
h4->SetTitle("not-matched_p");
h4->GetXaxis()->SetTitle("Lxy");
h4->GetXaxis()->SetRangeUser(a, b);
h4->Draw();

TCanvas a5("c", "canvas", 1200, 900);
//b4.Setlogy();
gPad->SetLogy();
gStyle->SetOptStat(1);
TH1F *h5 = new TH1F("h","h",300,0,300.);
h5=(TH1F*)child_e_Lxy_inB->Clone();
h5->Add(matched_child_e_Lxy_inB,-1);
h5->SetTitle("not-matched_e");
h5->GetXaxis()->SetTitle("Lxy");
h5->GetXaxis()->SetRangeUser(a, b);
h5->Draw();


b1.SaveAs("../plots/pdgID/pi_Lxy.pdf");
b2.SaveAs("../plots/pdgID/mu_Lxy.pdf");
b3.SaveAs("../plots/pdgID/K_Lxy.pdf");
b4.SaveAs("../plots/pdgID/p_Lxy.pdf");
b4.SaveAs("../plots/pdgID/e_Lxy.pdf");

a1.SaveAs("../plots/pdgID/notmatched_pi_Lxy.pdf");
a2.SaveAs("../plots/pdgID/notmatched_mu_Lxy.pdf");
a3.SaveAs("../plots/pdgID/notmatched_K_Lxy.pdf");
a4.SaveAs("../plots/pdgID/notmatched_p_Lxy.pdf");
a4.SaveAs("../plots/pdgID/notmatched_e_Lxy.pdf");


}
