{
float P_1=0.239,P_2=-1.220,P_3=-1.64*1e-5;    
TCanvas c1("c", "canvas", 900, 500);
//c1.Divide(2,1);
//c1.cd(1);
TH1F *h1 = new TH1F("h","heat map",100,0,2.);
h1=(TH1F*)nomatched_child_pT_jet_DR_400_inB->Clone();
h1->SetTitle("");
h1->GetXaxis()->SetTitle("pT_jet");
h1->GetYaxis()->SetTitle("DR_child,jet");
h1->GetXaxis()->SetRangeUser(0., 200);
h1->GetYaxis()->SetRangeUser(0., 1.5);
h1->Draw("BOX");
TF1 f1("f","0.239+exp(-1.220 + -1.64*1e-2 * x)",0.,200.);
f1.Draw("same");
}
