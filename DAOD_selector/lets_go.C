void letsgo()
{

  gROOT->ProcessLine(".L DAOD_selector_C.so");

//gROOT->ProcessLine(".x launch_selector_1Sample.C(\"Lecce\",\"Jcoll\",\"Label\" ,false,3000)");
  gROOT->ProcessLine(".x launch_selector_1Sample.C(\"Lecce\",\"Jcoll\",\"Label\" ,false)");
  gROOT->ProcessLine(".q");
}
