void go()
{

  gROOT->ProcessLine(".L DAOD_selector.C++");

  //gROOT->ProcessLine(".x launch_selector.C");
  //gROOT->ProcessLine(".x launch_selector_allSamples.C");
  //gROOT->ProcessLine(".x launch_selector_oneSample.C");

  gROOT->ProcessLine(".x launch_selector_1Sample.C(\"Lecce\",\"bTag_AntiKt4EMPFlowJets_BTagging201903\",\"Cone\"     ,false,3000)");

}
