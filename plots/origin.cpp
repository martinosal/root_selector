{
a1=-1.;a2=4.;

TCanvas c1("c", "canvas", 1200, 900);
c1.SetLogy();
gStyle->SetOptStat(0);
TH1F *h1 = new TH1F("h","pT fraction",100,a1,a2);
h1=(TH1F*)matched_trk_origin_inB->Clone();
h1->SetTitle("Akt4EMPf_retagT: matched origin");
h1->GetXaxis()->SetTitle("d0");
//h1->GetYaxis()->SetTitle("efficiency");
h1->GetXaxis()->SetRangeUser(a1,a2);
//h1->GetYaxis()->SetRangeUser(0., 1);
h1->Draw();
leg = new TLegend(0.7,0.7,0.9,0.9);
leg->AddEntry(h1,"Pile Up: -1","l");
leg->AddEntry(h1,"From B: 0","l");
leg->AddEntry(h1,"From C: 1","l");
leg->AddEntry(h1,"FRAG: 2","l");
leg->AddEntry(h1,"GEANT: 1","l");
leg->Draw();

TCanvas c3("c", "canvas", 1200, 900);
gStyle->SetOptStat(0);
trk_origin_inB->GetXaxis()->SetRangeUser(a1,a2);
trk_origin_inB->SetTitle("Akt4EMPf_retagT: origin");
trk_origin_inB->Draw();
leg = new TLegend(0.7,0.7,0.9,0.9);
leg->AddEntry(h1,"Pile Up: -1","l");
leg->AddEntry(h1,"From B: 0","l");
leg->AddEntry(h1,"From C: 1","l");
leg->AddEntry(h1,"FRAG: 2","l");
leg->AddEntry(h1,"GEANT: 1","l");
leg->Draw();

TCanvas c2("c", "canvas", 1200, 900);
c2.SetLogy();
TH1F *h2 = new TH1F("h2","efficiency_vs_origin",100,a1,a2);
gStyle->SetOptStat(0);
h2=(TH1F*)matched_trk_origin_inB->Clone();
h2->Divide(trk_origin_inB);
h2->SetTitle("Akt4EMPf_retagT: efficiency_vs_origin");
h2->GetXaxis()->SetTitle("origin");
h2->GetYaxis()->SetTitle("efficiency");
h2->GetXaxis()->SetRangeUser(a1,a2);
h2->GetYaxis()->SetRangeUser(1e-4, 1.1);
h2->Draw("");
leg = new TLegend(0.7,0.7,0.9,0.9);
leg->AddEntry(h1,"Pile Up: -1","l");
leg->AddEntry(h1,"From B: 0","l");
leg->AddEntry(h1,"From C: 1","l");
leg->AddEntry(h1,"FRAG: 2","l");
leg->AddEntry(h1,"GEANT: 1","l");
leg->Draw();
}
