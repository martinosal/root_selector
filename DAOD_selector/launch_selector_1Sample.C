#include "TChain.h"
#include "laptop_files.C"
#include "lxplus_files.C"
#include "lecce_files.C"

void launch_selector_1Sample(std::string location="Lecce", 
			     std::string jcoll="bTag_AntiKt4EMPFlowJets_BTagging201903", 
			     std::string labScheme="Cone", 
			     bool lowPt=false, 
			     int nEv=10000000)
{

  //////////=====================================

  //int nEventsToProcess = 10000000;
  //int nEventsToProcess = 3000;
  int nEventsToProcess = nEv;
  //
  bool laptop=false;
  bool lxplus=false;
  bool lecce =true;
  lecce  = (location=="Lecce"  || location=="lecce"  || location=="LECCE"  );
  lxplus = (location=="Lxplus" || location=="lxplus" || location=="LXPLUS" );
  laptop = (location=="Laptop" || location=="laptop" || location=="LAPTOP" );
  //
  bool lowPtForVR = lowPt;
  //
  //jetcollection="bTag_AntiKt4EMPFlowJets_BTagging201903";
  //jetcollection="bTag_AntiKtVR30Rmax4Rmin02TrackJets";
  //jetcollection="bTag_AntiKtVR30Rmax4Rmin02TrackJets_BTagging201903";
  //jetcollection="bTag_AntiKtVR30Rmax4Rmin02TrackGhostTagJets";
  const char *jetcollection=jcoll.c_str();
  //
  //Labeling ="Cone";
  //Labeling ="Ghost";
  //Labeling ="GhostCone";
  //Labeling ="ConeIncl";
  //Labeling ="GhostIncl";
  //Labeling ="GhostConeIncl";
  std::string Labeling=labScheme;

  //////////=====================================

  
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
  string decay_mode="false";//can be "leptonic" or "hadronic", set "false" or any other value for decay_mode=false

  float jet_pT_infcut=20*1e3,jet_pT_supcut=1000*1e4,jet_eta_cut=2.5,jet_JVT_cut=0.2;//was 0.5
  //float jet_pT_infcut=20*1e3,jet_pT_supcut=300*1e3,jet_eta_cut=2.5,jet_JVT_cut=0.5;
  float DR_bcH_cut=0.3,pT_bcH_cut=5*1e3;
  float trk_pT_cut=0.4*1e3,trk_eta_cut=2.4,trk_d0_cut=50.*1e3,trk_z0sinth_cut=50.*1e3;
  //float trk_pT_cut=1e3,trk_eta_cut=2.4,trk_d0_cut=1.;


  if (lowPtForVR){
  //////// changing cuts
  // jet pt > 12 GeV 
    jet_pT_infcut=12*1e3;
  //////// changing cuts
  }

  

  DAOD_selector a;

  a.setFlags(lxplus,debug,derived_origin,selections,discriminants,shrinking_cone,selection_alg,origin_selection,geometric_selection,cut,retag,m_p1,m_p2,m_p3,decay_mode);
  a.setCuts(jet_pT_infcut,jet_pT_supcut,jet_eta_cut,jet_JVT_cut,DR_bcH_cut,pT_bcH_cut,trk_pT_cut,trk_eta_cut,trk_d0_cut,trk_z0sinth_cut);


  TChain *f = new TChain(jetcollection);

  std::cout<<"\nJet Collection: " << jetcollection << "\n";
  if(laptop)  {laptop_files(f);}
  if(lecce)   {lecce_files( f,std::string(jetcollection));}
  if(lxplus)  {lxplus_files(f,std::string(jetcollection));}

  if (lowPtForVR)
    {
      if (Labeling=="Cone"      || Labeling=="ConeIncl"     || 
	  Labeling=="Ghost"     || Labeling=="GhostIncl"    ||
	  Labeling=="GhostCone" || Labeling=="GhostConeIncl")
	{
	  a.setOutputFNameString(std::string(jetcollection)+"_"+Labeling+"_12GeV");
	  a.setJetLabeling(Labeling);
	}
      /*
      if (Labeling=="Cone")
	{
	  a.setOutputFNameString(std::string(jetcollection)+"_ConeLab_12GeV");
	  a.setJetLabeling("Cone");
	}
      else if (Labeling=="Ghost")
	{
	  a.setOutputFNameString(std::string(jetcollection)+"_GhostLab_12GeV");
	  a.setJetLabeling("Ghost");
	}
      else if (Labeling=="GhostCone")
	{
	  a.setOutputFNameString(std::string(jetcollection)+"_GhostConeLab_12GeV");
	  a.setJetLabeling("GhostCone");
	}
      */
    }
  else
    {
      if (Labeling=="Cone"      || Labeling=="ConeIncl"     || 
	  Labeling=="Ghost"     || Labeling=="GhostIncl"    ||
	  Labeling=="GhostCone" || Labeling=="GhostConeIncl")
	{
	  a.setOutputFNameString(std::string(jetcollection)+"_"+Labeling);
	  a.setJetLabeling(Labeling);
	}
    }
    
  //f->Process(&a);
  f->Process(&a,"",nEventsToProcess); //to runprocess only the first 3000 events

  
}
