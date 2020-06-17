{
a1=0.;a2=300.;
TCanvas c1("c", "canvas", 1200, 900);
gStyle->SetOptStat(0);
trk_pT_jet_DR_inB->SetTitle("AktVR30Rmax4Rmin02Tr_BTagging201903_inB");
trk_pT_jet_DR_inB->GetXaxis()->SetTitle("pT_jet [GeV]");
trk_pT_jet_DR_inB->GetYaxis()->SetTitle("DR");
trk_pT_jet_DR_inB->GetXaxis()->SetRangeUser(a1,a2);
trk_pT_jet_DR_inB->GetYaxis()->SetRangeUser(0., 1.);
//leg = new TLegend(0.6,0.8,0.9,0.9);
//leg->AddEntry("origin=-1","PUFAKE");
trk_pT_jet_DR_inB->Draw("BOX");
}
