void go()
{

  gROOT->ProcessLine(".L DAOD_selector.C++");

  gROOT->ProcessLine(".x launch_selector.C");

}
