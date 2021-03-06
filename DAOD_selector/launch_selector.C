#include "TChain.h"
#include "laptop_files.C"
#include "lxplus_files.C"
#include "lecce_files.C"

void launch_selector()
{
  bool laptop=false;
  bool lxplus=true;
  bool lecce=false;

  bool debug=true;
  bool derived_origin=true;
  bool selections=true;
  bool discriminants=true;
  bool shrinking_cone=false;
  bool selection_alg=true;
  bool origin_selection=true;
  bool geometric_selection=true;
  bool cut=true;
  bool retag=false;
  double m_p1=0.,m_p2=0.,m_p3=0.;

  std::string decay_mode="false";//can be "leptonic" or "hadronic", set "false" or any other value for decay_mode=false
  std::string jetlabeling="DIAGjetlab";//can be "DIAGjetlab", "CONEjetlab" or "GHOSTjetlab"
  std::string collection="EMPFlow";


  float jet_pT_infcut,jet_pT_supcut=1e8,jet_eta_cut=2.5,jet_JVT_cut=0.2;//was 0.5
  
  float DR_bcH_cut=0.3,pT_bcH_cut=5*1e3;
  float trk_pT_cut=0.4*1e3,trk_eta_cut=2.4,trk_d0_cut=50.*1e3,trk_z0sinth_cut=50.*1e3;
  //float trk_pT_cut=1e3,trk_eta_cut=2.4,trk_d0_cut=1.;


    const char *jetcollection="";

    if(collection=="VR30ghost"){
      jet_pT_infcut=12*1e3;
      jetcollection="bTag_AntiKtVR30Rmax4Rmin02TrackGhostTagJets";
    }
    if(collection=="VR30BTag"){
      jet_pT_infcut=12*1e3;
      jetcollection="bTag_AntiKtVR30Rmax4Rmin02TrackJets_BTagging201903";
    }
    if(collection=="EMPFlow"){
      jet_pT_infcut=20*1e3;
      jetcollection="bTag_AntiKt4EMPFlowJets_BTagging201903";
    }

    //const char *jetcollection="bTag_AntiKtVR30Rmax4Rmin02TrackJets_BTagging201903";



  std::cout<<"\nJet Collection: " << jetcollection << "\n";
  std::cout<<"Jet labeling: "<<jetlabeling<<std::endl;

  TChain *f = new TChain(jetcollection);

  if(lecce)   {lecce_files(f);}
  if(laptop)  {laptop_files(f);}
  if(lxplus)  {lxplus_files(f,std::string(jetcollection));}

  DAOD_selector a;

  a.setFlags(lxplus,debug,derived_origin,selections,discriminants,shrinking_cone,selection_alg,origin_selection,geometric_selection,cut,retag,m_p1,m_p2,m_p3,decay_mode, jetlabeling);
  a.setCuts(jet_pT_infcut,jet_pT_supcut,jet_eta_cut,jet_JVT_cut,DR_bcH_cut,pT_bcH_cut,trk_pT_cut,trk_eta_cut,trk_d0_cut,trk_z0sinth_cut);

  if(!decay_mode.compare("leptonic") || !decay_mode.compare("hadronic")){
    a.setOutputFNameString(std::string(jetcollection)+"_"+jetlabeling+"_"+decay_mode);
  }
  else{
    a.setOutputFNameString(std::string(jetcollection)+"_"+jetlabeling);
  }


      f->Process(&a);
  //    f->Process(&a,"",10000,1); //for developing

}
