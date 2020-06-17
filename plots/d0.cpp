{
a1=-2;a2=6.;
TCanvas c1("c", "canvas", 1200, 900);
gStyle->SetOptStat(0);
TH1F *h1 = new TH1F("h1","efficiency_vs_d0",100,a1,a2);
h1=(TH1F*)matched_child_d0_inB->Clone();
h1->Divide(child_d0_inB);
h1->SetTitle("efficiency_vs_d0");
h1->GetXaxis()->SetTitle("d0");
h1->GetYaxis()->SetTitle("efficiency");
h1->GetXaxis()->SetRangeUser(a1,a2);
h1->GetYaxis()->SetRangeUser(0., 1.1);
//leg = new TLegend(0.6,0.8,0.9,0.9);
//leg->AddEntry("origin=-1","PUFAKE");
h1->Draw("");
/*
TCanvas c2("c2", "canvas", 1200, 900);
TH1F *h2 = new TH1F("h1","single_matched/matched_vs_d0",100,a1,a2);
gStyle->SetOptStat(0);
h1=(TH1F*)single_matched_child_d0_inB->Clone();
h1->Divide(matched_child_d0_inB);
h1->SetTitle("single_matched/matched_vs_d0");
h1->GetXaxis()->SetTitle("d0");
h1->GetYaxis()->SetTitle("single_matched/matched");
h1->GetXaxis()->SetRangeUser(a1,a2);
h1->GetYaxis()->SetRangeUser(0., 1.1);
//leg = new TLegend(0.6,0.8,0.9,0.9);
//leg->AddEntry("origin=-1","PUFAKE");
h1->Draw("");
*/
}
