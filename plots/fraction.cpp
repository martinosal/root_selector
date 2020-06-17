{
  gStyle->SetOptStat(0);
TCanvas c1("c", "canvas", 1200, 600);
//c1.SetLogy();
c1.Divide(2,1);
c1.cd(1);
TH1D *h1 = matched_pT_child_pTfraction_inB->ProjectionY();
//h1=(TH1F*)matched_pT_child_pTfraction_inB->Clone();
h1->SetTitle("pT fraction");
h1->GetXaxis()->SetTitle("DpT/pT");
//h1->GetYaxis()->SetTitle("efficiency");
h1->GetXaxis()->SetRangeUser(0.01, 10);
//h1->GetYaxis()->SetRangeUser(0., 1);
h1->Draw();

c1.cd(2);
TH1F *h2 = new TH1F("h","DR",100,0,1.);
h2=(TH1F*)matched_DR_trk_inB->Clone();
h2->SetTitle("DR");
h2->GetXaxis()->SetTitle("DR_child,track");
//h2->GetYaxis()->SetTitle("efficiency");
h2->GetXaxis()->SetRangeUser(0., 1);
//h1->GetYaxis()->SetRangeUser(0., 1);
h2->Draw();
}
