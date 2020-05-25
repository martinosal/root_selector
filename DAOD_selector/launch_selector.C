#include "TChain.h"
#include "laptop_files.C"
#include "lxplus_files.C"

void launch_selector()
{
  bool laptop=true;
  bool lxplus=false;
  bool debug=true;
  bool selections=true;
  bool discriminants=false;
  bool shrinking_cone=false;
  bool selection_alg=true;
  bool cut=true;
  bool retagT=false;

  TChain *f = new TChain("bTag_AntiKt4EMPFlowJets_BTagging201903");

  if(laptop)  {laptop_files(f);}
  if(lxplus)  {lxplus_files(f);}

  DAOD_selector a;

  a.setFlags(lxplus,debug,selections,discriminants,shrinking_cone,selection_alg,cut,retagT);

  f->Process(&a);

//  f->Process("DAOD_selector.C");
//  f->Process("DAOD_selector.C","",100,1); //for developing
}
