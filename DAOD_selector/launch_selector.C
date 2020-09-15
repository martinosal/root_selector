#include "TChain.h"
#include "laptop_files.C"
#include "lxplus_files.C"
#include "lecce_files.C"

void launch_selector()
{
  bool laptop=true;
  bool lxplus=false;
  bool lecce=false;


  bool debug=true;
  bool selections=true;
  bool discriminants=true;
  bool shrinking_cone=false;
  bool selection_alg=true;
  bool origin_selection=true;
  bool geometric_selection=true;
  bool cut=false;
  bool retagT=false;
  string decay_mode="false";//can be "leptonic" or "hadronic", set "false" or any other value for decay_mode=false

  float jet_pT_infcut=20*1e3,jet_pT_supcut=1000*1e3,jet_eta_cut=2.5,jet_JVT_cut=0.5;
  //float jet_pT_infcut=20*1e3,jet_pT_supcut=300*1e3,jet_eta_cut=2.5,jet_JVT_cut=0.5;
  float DR_bcH_cut=0.3,pT_bcH_cut=5*1e3;
  float trk_pT_cut=0.4*1e3,trk_eta_cut=2.4,trk_d0_cut=1e3;
  //float trk_pT_cut=1e3,trk_eta_cut=2.4,trk_d0_cut=1.;

//  const char *jetcollection="bTag_AntiKt4EMPFlowJets_BTagging201903";
  const char *jetcollection="bTag_AntiKt4EMPFlowJets";
//  const char *jetcollection="bTag_AntiKtVR30Rmax4Rmin02TrackJets_BTagging201903";

  std::cout<<"\nJet Collection: " << jetcollection << "\n";

  TChain *f = new TChain(jetcollection);

  if(lecce)   {lecce_files(f);}
  if(laptop)  {laptop_files(f);}
  if(lxplus)  {lxplus_files(f);}

  DAOD_selector a;

  a.setFlags(lxplus,debug,selections,discriminants,shrinking_cone,selection_alg,origin_selection,geometric_selection,cut,retagT,decay_mode);
  a.setCuts(jet_pT_infcut,jet_pT_supcut,jet_eta_cut,jet_JVT_cut,DR_bcH_cut,pT_bcH_cut,trk_pT_cut,trk_eta_cut,trk_d0_cut);

  f->Process(&a);
//  f->Process(&a,"",1,10); //for developing

//  f->Process("DAOD_selector.C");

}
