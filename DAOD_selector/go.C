void go()
{

  // to compile and load the library 
  gROOT->ProcessLine(".L DAOD_selector.C++");
  // to load the library already compiled 
  //gROOT->ProcessLine(".L DAOD_selector.C++");

  //gROOT->ProcessLine(".x launch_selector.C");
  //gROOT->ProcessLine(".x launch_selector_allSamples.C");
  //gROOT->ProcessLine(".x launch_selector_oneSample.C");

  //// to produce plots histograms with all variants with go.sh (using lets_go.sh)
  ////
  //gROOT->ProcessLine(".x launch_selector_1Sample.C(\"Lecce\",\"bTag_AntiKt4EMPFlowJets_BTagging201903\",\"Cone\"     ,false,3000)");


  ////// to produce plots for jet labeling studies, run one after the other the following
  //////
  //gROOT->ProcessLine(".x launch_selector_1Sample_labelingStudyOnly.C(\"Lecce\",\"bTag_AntiKt4EMPFlowJets_BTagging201903\",            \"Cone\"              ,false)");
  //gROOT->ProcessLine(".x launch_selector_1Sample_labelingStudyOnly.C(\"Lecce\",\"bTag_AntiKtVR30Rmax4Rmin02TrackJets_BTagging201903\",\"Cone\"              ,false)");
  //gROOT->ProcessLine(".x launch_selector_1Sample_labelingStudyOnly.C(\"Lecce\",\"bTag_AntiKtVR30Rmax4Rmin02TrackGhostTagJets\",       \"Ghost\"             ,false)");
  //gROOT->ProcessLine(".x launch_selector_1Sample_labelingStudyOnly.C(\"Lecce\",\"bTag_AntiKtVR30Rmax4Rmin02TrackJets_BTagging201903\",\"Cone\"              ,true)");
  //gROOT->ProcessLine(".x launch_selector_1Sample_labelingStudyOnly.C(\"Lecce\",\"bTag_AntiKtVR30Rmax4Rmin02TrackGhostTagJets\",       \"Ghost\"             ,true)");
  // change in launch_selector_1Sample_labelingStudyOnly.C the pt_min cut from 10 GeV to 12 GeV or viceversa and 
  //gROOT->ProcessLine(".x launch_selector_1Sample_labelingStudyOnly.C(\"Lecce\",\"bTag_AntiKtVR30Rmax4Rmin02TrackJets_BTagging201903\",\"Cone\"              ,true)");
  //gROOT->ProcessLine(".x launch_selector_1Sample_labelingStudyOnly.C(\"Lecce\",\"bTag_AntiKtVR30Rmax4Rmin02TrackGhostTagJets\",       \"Ghost\"             ,true)");


}
