{

//IP3D LogLikelihoodRatio inclusive
int r=2;
TCanvas c2("c", "canvas", 1200, 900);
c2.SetLogy();
gStyle->SetOptStat(0);
Double_t norm_2 = ip3d_llr_l->GetEntries();
ip3d_llr_l->Scale(1/norm_2);
ip3d_llr_l->Rebin(r);
ip3d_llr_l->SetLineColor(kRed);
ip3d_llr_l->SetLineWidth(3);
ip3d_llr_l->GetYaxis()->SetRangeUser(1e-4, 1.);
ip3d_llr_l->SetTitle("IP3D LogLikelihoodRatio");
ip3d_llr_l->Draw("hist");
Double_t norm_2b = ip3d_llr_inC->GetEntries();
ip3d_llr_inC->Scale(1/norm_2b);
ip3d_llr_inC->Rebin(r);
ip3d_llr_inC->SetLineColor(kBlue);
ip3d_llr_inC->SetLineWidth(3);
ip3d_llr_inC->Draw("same hist");
Double_t norm_2c = ip3d_llr_inB->GetEntries();
ip3d_llr_inB->Scale(1/norm_2c);
ip3d_llr_inB->Rebin(r);
ip3d_llr_inB->SetLineColor(kGreen+2);
ip3d_llr_inB->SetLineWidth(3);
ip3d_llr_inB->Draw("same hist");
leg = new TLegend(0.6,0.8,0.9,0.9);
//leg->SetHeader("The Legend Title");
leg->AddEntry(ip3d_llr_l,"light jets","l");
leg->AddEntry(ip3d_llr_inC,"c jets","l");
leg->AddEntry(ip3d_llr_inB,"b jets","l");
leg->Draw();
}
