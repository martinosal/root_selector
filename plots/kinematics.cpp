{
TCanvas c1("c", "canvas", 1200, 900);
gStyle->SetOptStat(0);
c1.Divide(2,2);
c1.cd(1);
TH1F *h1 = new TH1F("h","DR",100,0,2.);
h1=(TH1F*)child_DR_inB->Clone();
h1->Add(matched_child_DR_inB,-1);
h1->GetXaxis()->SetTitle("Delta R_child,jet");
h1->GetXaxis()->SetRangeUser(0., 1.5);
//h1->GetYaxis()->SetRangeUser(0., 1.1);
h1->Draw();

c1.cd(2);
TH1F *h2 = new TH1F("h","pT",100,0,200.);
h2=(TH1F*)child_pT_inB->Clone();
h2->Add(matched_child_pT_inB,-1);
h2->GetXaxis()->SetTitle("pT_child [GeV]");
h2->GetXaxis()->SetRangeUser(0., 100.);
//h2->GetYaxis()->SetRangeUser(0., 1.1);
h2->Draw();

c1.cd(3);
TH1F *h3 = new TH1F("h","eta",100,0,2.);
h3=(TH1F*)child_eta_inB->Clone();
h3->Add(matched_child_eta_inB,-1);
h3->GetXaxis()->SetTitle("eta_child");
h3->GetXaxis()->SetRangeUser(-2.6, 2.6);
//h3->GetYaxis()->SetRangeUser(0., 1.1);
h3->Draw();

c1.cd(4);
TH1F *h4 = new TH1F("h","phi",100,0,200.);
h4=(TH1F*)child_phi_inB->Clone();
h4->Add(matched_child_phi_inB,-1);
h4->GetXaxis()->SetTitle("phi_child");
h4->GetXaxis()->SetRangeUser(-4., 4.);
//h4->GetYaxis()->SetRangeUser(0., 1.1);
h4->Draw();
}
