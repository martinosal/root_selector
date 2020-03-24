void plots(){
    TCanvas c("c", "canvas", 1200, 900);
    gStyle->SetOptStat(0);
    Double_t norm_2 = ip2d_llr_inC->GetEntries();
    ip2d_llr_inC->Scale(1/norm_2);
    ip2d_llr_inC->SetLineColor(kBlue);
    ip2d_llr_inC->SetLineWidth(3);
    ip2d_llr_inC->Draw("hist");
    Double_t norm_1 = ip2d_llr_inB->GetEntries();
    ip2d_llr_inB->Scale(1/norm_1);
    ip2d_llr_inB->SetLineColor(kGreen+2);
    ip2d_llr_inB->SetLineWidth(3);
    ip2d_llr_inB->Draw("same hist title");
    gPad->BuildLegend();
    
}


