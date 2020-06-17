//IP2D LogLikelihoodRatio inclusive
{
  int r=2;
TCanvas c1("c", "canvas", 1200, 900);
c1.SetLogy();
gStyle->SetOptStat(0);
Double_t norm_1 = ip2d_llr_l->GetEntries();
ip2d_llr_l->Scale(1/norm_1);
ip2d_llr_l->Rebin(r);
ip2d_llr_l->SetLineColor(kRed);
ip2d_llr_l->SetLineWidth(3);
ip2d_llr_l->GetYaxis()->SetRangeUser(1e-4, 1.);
ip2d_llr_l->SetTitle("IP2D LogLikelihoodRatio");
ip2d_llr_l->Draw("hist");
Double_t norm_1b = ip2d_llr_inC->GetEntries();
ip2d_llr_inC->Scale(1/norm_1b);
ip2d_llr_inC->Rebin(r);
ip2d_llr_inC->SetLineColor(kBlue);
ip2d_llr_inC->SetLineWidth(3);
ip2d_llr_inC->Draw("same hist");
Double_t norm_1c = ip2d_llr_inB->GetEntries();
ip2d_llr_inB->Scale(1/norm_1c);
ip2d_llr_inB->Rebin(r);
ip2d_llr_inB->SetLineColor(kGreen+2);
ip2d_llr_inB->SetLineWidth(3);
ip2d_llr_inB->Draw("same hist");
leg = new TLegend(0.6,0.8,0.9,0.9);
//leg->SetHeader("The Legend Title");
leg->AddEntry(ip2d_llr_l,"light jets","l");
leg->AddEntry(ip2d_llr_inC,"c jets","l");
leg->AddEntry(ip2d_llr_inB,"b jets","l");
leg->Draw();

}
