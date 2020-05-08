{
TCanvas c0("c", "canvas", 1200, 900);
gStyle->SetPalette(55);
TH2F *h0 = new TH2F("h","heat map",100,0.,1.,100,0.,10.);
h0=(TH2F*)matched_DR_trk_pTfraction->Clone();
//h->Rebin(3);
h0->Divide(single_matched_DR_trk_pTfraction);
//h->Rebin(3);
h0->SetTitle("heat map efficiency");
h0->GetXaxis()->SetTitle("DR_trk");
h0->GetYaxis()->SetTitle("DpT/pT");
h0->GetXaxis()->SetRangeUser(0., 1.);
h0->GetYaxis()->SetRangeUser(0., 10.);
h0->Draw("COLZ");

TCanvas c1("c", "canvas", 1200, 900);
gStyle->SetPalette(55);
TH2F *h1 = new TH2F("h","heat map",100,-1.,1.,100,-1.,1.);
h1=(TH2F*)matched_child_Dphi_Deta_400->Clone();
//h->Rebin(3);
h1->Divide(single_matched_child_Dphi_Deta_400);
//h->Rebin(3);
h1->SetTitle("heat map efficiency");
h1->GetXaxis()->SetTitle("Dphi");
h1->GetYaxis()->SetTitle("Deta");
h1->GetXaxis()->SetRangeUser(-1.,1.);
h1->GetYaxis()->SetRangeUser(-1.,1.);
h1->Draw("COLZ");
}
