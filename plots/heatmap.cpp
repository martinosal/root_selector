{
TCanvas c0("c", "canvas", 1200, 900);
gStyle->SetPalette(55);
TH2F *h0 = new TH2F("h","heat map",100,0.,150.,100,0.,0.6);
h0=(TH2F*)matched_child_pT_DR_400_inB->Clone();
//h->Rebin(3);
h0->Divide(child_pT_DR_400_inB);
//h->Rebin(3);
h0->SetTitle("heat map efficiency");
h0->GetXaxis()->SetTitle("pT_child [GeV]");
h0->GetYaxis()->SetTitle("Delta R_child,jet");
h0->GetXaxis()->SetRangeUser(0., 150.);
h0->GetYaxis()->SetRangeUser(0., 0.6);
h0->Draw("COLZ");
}
