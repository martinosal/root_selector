{

//IP2D LogLikelihoodRatio C_ex_B_in

TCanvas c5("c", "canvas", 1200, 900);
c5.SetLogy();
gStyle->SetOptStat(0);
Double_t norm_5 = ip2d_llr_l->GetEntries();
ip2d_llr_l->Scale(1/norm_5);
ip2d_llr_l->Rebin(r);
ip2d_llr_l->SetLineColor(kRed);
ip2d_llr_l->SetLineWidth(3);
ip2d_llr_l->GetYaxis()->SetRangeUser(1e-3, 1.);
ip2d_llr_l->SetTitle("IP2D LogLikelihoodRatio");
ip2d_llr_l->Draw("hist");
Double_t norm_5b = ip2d_llr_exC->GetEntries();
ip2d_llr_exC->Scale(1/norm_5b);
ip2d_llr_exC->Rebin(r);
ip2d_llr_exC->SetLineColor(kBlue);
ip2d_llr_exC->SetLineWidth(3);
ip2d_llr_exC->Draw("same hist");
Double_t norm_5c = ip2d_llr_inB->GetEntries();
ip2d_llr_inB->Scale(1/norm_5c);
ip2d_llr_inB->Rebin(r);
ip2d_llr_inB->SetLineColor(kGreen+2);
ip2d_llr_inB->SetLineWidth(3);
ip2d_llr_inB->Draw("same hist");
leg = new TLegend(0.6,0.8,0.9,0.9);
//leg->SetHeader("The Legend Title");
leg->AddEntry(ip2d_llr_l,"light jets","l");
leg->AddEntry(ip2d_llr_exC,"c jets","l");
leg->AddEntry(ip2d_llr_inB,"b jets","l");
leg->Draw();

}
