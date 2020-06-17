{
int r=2;
TCanvas c8("c", "canvas", 1200, 900);
c8.SetLogy();
gStyle->SetOptStat(0);
Double_t norm_8 = DL1_light->GetEntries();
DL1_light->Scale(1/norm_8);
DL1_light->Rebin(r);
DL1_light->SetLineColor(kRed);
DL1_light->SetLineWidth(3);
DL1_light->GetYaxis()->SetRangeUser(1e-4, 1.);
DL1_light->SetTitle("DL1");
DL1_light->Draw("hist");
Double_t norm_8b = DL1_inC->GetEntries();
DL1_inC->Scale(1/norm_8b);
DL1_inC->Rebin(r);
DL1_inC->SetLineColor(kGreen+2);
DL1_inC->SetLineWidth(3);
DL1_inC->Draw("same hist");
Double_t norm_8c = DL1_inB->GetEntries();
DL1_inB->Scale(1/norm_8c);
DL1_inB->Rebin(r);
DL1_inB->SetLineColor(kBlue);
DL1_inB->SetLineWidth(3);
DL1_inB->Draw("same hist");
leg = new TLegend(0.6,0.8,0.9,0.9);
//leg->SetHeader("The Legend Title");
leg->AddEntry(DL1_light,"light jets","l");
leg->AddEntry(DL1_inC,"c jets","l");
leg->AddEntry(DL1_inB,"b jets","l");
leg->Draw();

}
