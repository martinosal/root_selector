#include "TChain.h"
#include "laptop_files.C"
#include "lxplus_files.C"

void launch_selector()
{
  bool laptop=true;
  bool lxplus=false;
  bool debug=true;
  bool selections=true;
  bool discriminants=true;
  bool shrinking_cone=false;
  bool selection_alg=true;
  bool origin_selection=true;
  bool geometric_selection=true;
  bool cut=false;
  bool retagT=false;

//  const char *jetcollection="bTag_AntiKt4EMPFlowJets_BTagging201903";
  const char *jetcollection="bTag_AntiKt4EMPFlowJets";
//  const char *jetcollection="bTag_AntiKtVR30Rmax4Rmin02TrackJets_BTagging201903";

  std::cout<<"\nJet Collection: " << jetcollection << "\n";

  TChain *f = new TChain(jetcollection);

  if(laptop)  {laptop_files(f);}
  if(lxplus)  {lxplus_files(f);}

  DAOD_selector a;

  a.setFlags(lxplus,debug,selections,discriminants,shrinking_cone,selection_alg,origin_selection,geometric_selection,cut,retagT);

  f->Process(&a);

//  f->Process("DAOD_selector.C");
//  f->Process("DAOD_selector.C","",100,1); //for developing

}
