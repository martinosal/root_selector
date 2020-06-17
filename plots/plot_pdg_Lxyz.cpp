{
  Lxyz_max=10000.;
  TH1F *h1 = new TH1F("h","h",300,0,300.);
  h1=(TH1F*)child_pi_Lxyz_inB->Clone();
  TH1F *h2 = new TH1F("h","h",300,0,300.);
  h2=(TH1F*)child_mu_Lxyz_inB->Clone();
  TH1F *h3 = new TH1F("h","h",300,0,300.);
  h3=(TH1F*)child_K_Lxyz_inB->Clone();
  TH1F *h4 = new TH1F("h","h",300,0,300.);
  h4=(TH1F*)child_p_Lxyz_inB->Clone();
  TH1F *h5 = new TH1F("h","h",300,0,300.);
  h5=(TH1F*)child_e_Lxyz_inB->Clone();
  TCanvas cc("c", "canvas", 1200, 900);
  //b4.Setlogy();
  gPad->SetLogy();
  gStyle->SetOptStat(0);
  Double_t n1 = h1->GetEntries();
  h1->Scale(1/n1);
  h1->GetXaxis()->SetRangeUser(0., Lxyz_max);
  h1->SetLineColor(kRed);
  h1->SetLineWidth(2);
  h1->SetTitle("bTag_AntiKt4EMPFlowJets_201903: children_Lxyz");
  h1->GetXaxis()->SetTitle("Lxyz");
  Double_t n2 = h2->GetEntries();
  h2->Scale(1/n2);
  h2->SetLineColor(kGreen+1);
  h2->SetLineWidth(2);
  Double_t n3 = h3->GetEntries();
  h3->Scale(1/n3);
  h3->SetLineColor(kBlue);
  h3->SetLineWidth(2);
  Double_t n4 = h4->GetEntries();
  h4->Scale(1/n4);
  h4->SetLineColor(kOrange);
  h4->SetLineWidth(2);
  Double_t n5 = h5->GetEntries();
  h5->Scale(1/n5);
  h5->SetLineColor(kBlack);
  h5->SetLineWidth(2);

  h1->Draw("hist");
  h2->Draw("same hist");
  h3->Draw("same hist");
  h4->Draw("same hist");
  h5->Draw("same hist");

  leg = new TLegend(0.8,0.8,0.9,0.9);
  //leg->SetHeader("The Legend Title");
  leg->AddEntry(h1,"pi","l");
  leg->AddEntry(h2,"mu","l");
  leg->AddEntry(h3,"K","l");
  leg->AddEntry(h4,"p","l");
  leg->AddEntry(h5,"e","l");
  leg->Draw();

  cc.SaveAs("../plots/pdgID/Lxyz.pdf");
}
