void makePlots()
{

  gROOT->ProcessLine(".L ../plots/efficiency/eff.cpp");
  gROOT->ProcessLine(".L ../plots/origin/scripts/truth_label_origin.cpp");
  gROOT->ProcessLine(".L ../plots/origin/scripts/assoc_tracks.cpp");
		      
		      
  gROOT->ProcessLine("eff_pt(\"jet\",\"b\")");
  gROOT->ProcessLine("eff_pt(\"jet\",\"c\")");
  gROOT->ProcessLine("eff_pt(\"bH\" ,\"b\")");
  gROOT->ProcessLine("eff_pt(\"cH\" ,\"c\")");
		      
  gROOT->ProcessLine("eff_DR(\"jet\",\"b\")");
  gROOT->ProcessLine("eff_DR(\"jet\",\"c\")");
  gROOT->ProcessLine("eff_DR(\"bH\" ,\"b\")");
  gROOT->ProcessLine("eff_DR(\"cH\" ,\"c\")");

  gROOT->ProcessLine("jet_orig(\"Bjet_cut_origin_truth_label_pT\",_file0,50,250,\"jet\") ");
  gROOT->ProcessLine("jet_orig(\"Bjet_cut_origin_truth_label_bHpT\",_file0,50,250,\"bH\")");
  gROOT->ProcessLine("jet_orig(\"Cjet_cut_origin_truth_label_pT\",_file0,50,250,\"jet\") ");
  gROOT->ProcessLine("jet_orig(\"Cjet_cut_origin_truth_label_cHpT\",_file0,50,250,\"cH\")");
  gROOT->ProcessLine("jet_orig(\"ljet_cut_origin_truth_label_pT\",_file0,50,250,\"jet\") ");

  gROOT->ProcessLine("jet_orig(\"trk_jet_pT_origin_truth_label_IP3D_B\",_file0,50,250,\"jet\") ");
  gROOT->ProcessLine("jet_orig(\"trk_bH_pT_origin_truth_label_IP3D_B\",_file0,50,250,\"bH\")   ");
  gROOT->ProcessLine("jet_orig(\"trk_jet_pT_origin_truth_label_SV1_B\",_file0,50,250,\"jet\")  ");
  gROOT->ProcessLine("jet_orig(\"trk_bH_pT_origin_truth_label_SV1_B\",_file0,50,250,\"bH\")    ");
 											       
  gROOT->ProcessLine("avtrk_orig(\"Bjet_cut_origin_truth_label_pT\",\"pT_B\",_file0,50,250,\"jet\")    ");
  gROOT->ProcessLine("avtrk_orig(\"Bjet_cut_origin_truth_label_bHpT\",\"bHpT_B\",_file0,50,250,\"bH\") ");
  gROOT->ProcessLine("avtrk_orig(\"Cjet_cut_origin_truth_label_pT\",\"pT_C\",_file0,50,250,\"jet\")    ");
  gROOT->ProcessLine("avtrk_orig(\"Cjet_cut_origin_truth_label_cHpT\",\"cHpT_C\",_file0,50,250,\"cH\") ");
		      
  gROOT->ProcessLine("assoc_tracks(\"n_tracks_bHpt_B\",\"bHpT_B\",_file0,\"bH\") ");
  gROOT->ProcessLine("assoc_tracks(\"n_tracks_cHpt_C\",\"cHpT_C\",_file0,\"cH\") ");
  gROOT->ProcessLine("assoc_tracks(\"n_tracks_jetpt_B\",\"pT_B\",_file0,\"jet\") ");
  gROOT->ProcessLine("assoc_tracks(\"n_tracks_jetpt_C\",\"pT_C\",_file0,\"jet\") ");

}
