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
  bool cut=false;
  bool retagT=false;

  TChain *f = new TChain("bTag_AntiKt4EMTopoJets");

  if(laptop)  {laptop_files(f);}
  if(lxplus)  {lxplus_files(f);}

  selector_1 a;

  a.setFlags(lxplus,debug,selections,discriminants,shrinking_cone,selection_alg,cut,retagT);

  f->Process(&a);

//  f->Process("selector_1.C");
//  f->Process("selector_1.C","",100,1); //for developing
}
