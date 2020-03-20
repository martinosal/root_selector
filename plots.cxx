void plots(){
    TCanvas c("c", "canvas", 1200, 900);
    gStyle->SetOptStat(0);
    ip2d_llr_B->SetLineColor(kBlue);
    ip2d_llr->Draw();
    ip2d_llr_C->SetLineColor(kRed);
    ip2d_llr_C->Draw("same");
    gPad->BuildLegend();
}
