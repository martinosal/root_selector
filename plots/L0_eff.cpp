{
TCanvas c1("c", "canvas", 1200, 900);
TH1F *h1 = new TH1F("h1","efficiency_vs_origin_2",100,0,300);
h1=(TH1F*)matched_child_Lxy_inB->Clone();
h1->Divide(child_Lxy_inB);
h1->SetTitle("efficiency_vs_Lxy");
h1->GetXaxis()->SetTitle("Lxy");
h1->GetYaxis()->SetTitle("efficiency");
h1->GetXaxis()->SetRangeUser(0, 300);
h1->GetYaxis()->SetRangeUser(0., 1.1);
//leg = new TLegend(0.6,0.8,0.9,0.9);
//leg->AddEntry("origin=-1","PUFAKE");
h1->Draw("");

TCanvas c2("c", "canvas", 1200, 900);
TH1F *h2 = new TH1F("h2","efficiency_vs_origin_2",100,0,300);
h2=(TH1F*)matched_child_Lxyz_inB->Clone();
h2->Divide(child_Lxyz_inB);
h2->SetTitle("efficiency_vs_Lxyz");
h2->GetXaxis()->SetTitle("Lxyz");
h2->GetYaxis()->SetTitle("efficiency");
h2->GetXaxis()->SetRangeUser(0, 300);
h2->GetYaxis()->SetRangeUser(0., 1.1);
//leg = new TLegend(0.6,0.8,0.9,0.9);
//leg->AddEntry("origin=-1","PUFAKE");
h2->Draw("");
}
