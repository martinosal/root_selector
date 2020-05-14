#include "TChain.h"
#include "list_files.C"

void launch_selector()
{
  bool list=true;

  TChain *f = new TChain("bTag_AntiKt4EMTopoJets");

  if(list)  {list_files(f);}

//  selector_1 a;
//  f->Process(&a);
  f->Process("selector_1.C");
//  f->Process("selector_1.C","",100,1); //for developing
}
