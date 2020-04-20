#include "TChain.h"
#include "list_files.C"

void launch_selector()
{
  bool list=true;

  TChain *f = new TChain("bTag_AntiKt4EMTopoJets");

  if(list)  {list_files(f);}

  f->Process("selector_1.C");
//  f->Process("selector_1.C","",10,1); //for developing
}
