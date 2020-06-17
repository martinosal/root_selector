{
//IP3D LogLikelihoodRatio exclusive

TCanvas c4("c", "canvas", 1200, 900);
c4.SetLogy();
gStyle->SetOptStat(0);
Double_t norm_4 = ip2d_llr_l->GetEntries();
ip3d_llr_l->Scale(1/norm_4);
ip3d_llr_l->Rebin(r);
ip3d_llr_l->SetLineColor(kRed);
ip3d_llr_l->SetLineWidth(3);
ip3d_llr_l->GetYaxis()->SetRangeUser(1e-4, 1.);
ip3d_llr_l->SetTitle("IP3D LogLikelihoodRatio");
ip3d_llr_l->Draw("hist");
Double_t norm_4b = ip3d_llr_exC->GetEntries();
ip3d_llr_exC->Scale(1/norm_4b);
ip3d_llr_exC->Rebin(r);
ip3d_llr_exC->SetLineColor(kBlue);
ip3d_llr_exC->SetLineWidth(3);
ip3d_llr_exC->Draw("same hist");
Double_t norm_4c = ip3d_llr_exB->GetEntries();
ip3d_llr_exB->Scale(1/norm_4c);
ip3d_llr_exB->Rebin(r);
ip3d_llr_exB->SetLineColor(kGreen+2);
ip3d_llr_exB->SetLineWidth(3);
ip3d_llr_exB->Draw("same hist");
leg = new TLegend(0.6,0.8,0.9,0.9);
//leg->SetHeader("The Legend Title");
leg->AddEntry(ip3d_llr_l,"light jets","l");
leg->AddEntry(ip3d_llr_exC,"c jets","l");
leg->AddEntry(ip3d_llr_exB,"b jets","l");
leg->Draw();

}
