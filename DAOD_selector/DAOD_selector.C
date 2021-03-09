#define CHECK_BIT(var,pos) ((var) & (1<<(pos)))
#define DAOD_selector_cxx
// The class definition in selector.h has been generated automatically
// by the ROOT utility TTree::MakeSelector(). This class is derived
// from the ROOT class TSelector. For more information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.


// The following methods are defined in this file:
//    Begin():        called every time a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.
//
// To use this file, try the following session on your Tree T:
//
// root> T->Process("selector.C")
// root> T->Process("selector.C","some options")
// root> T->Process("selector.C+")
//


#include "DAOD_selector.h"
#include <TStyle.h>
#include "TTree.h"

DAOD_selector::DAOD_selector()
{
   fOutputString="";
   fJetLabeling ="Cone";
   fDoFlavorLabelMatrix=false;
}


void DAOD_selector::setFlags(bool lxplus_flag, bool debug_flag, bool derived_origin_flag, bool selections_flag, bool discriminants_flag, bool shrinking_cone_flag, bool selection_alg_flag, bool origin_selection_flag, bool geometric_selection_flag, bool cut_flag, bool retag_flag, double p1, double p2, double p3, string decay_mode_flag)
{
   std::cout<<"\n In DAOD_selector::setFlags"<<std::endl;

   lxplus=lxplus_flag;
   debug=debug_flag;
   derived_origin=derived_origin_flag;
   selections=selections_flag;
   discriminants=discriminants_flag;
   shrinking_cone=shrinking_cone_flag;
   selection_alg=selection_alg_flag;
   origin_selection=origin_selection_flag;
   geometric_selection=geometric_selection_flag;
   cut=cut_flag;
   retag=retag_flag;
   if(retag==false){
     m_p1=0.239;
     m_p2=-1.220;
     m_p3=-1.64*1e-2;
   }
   if(retag==true){
     m_p1=p1;
     m_p2=p2;
     m_p3=p3;
   }
   decay_mode=decay_mode_flag;
   return;
}

void DAOD_selector::setCuts(float m_jet_pT_infcut, float m_jet_pT_supcut, float m_jet_eta_cut, float m_jet_JVT_cut, float m_DR_bcH_cut, float m_pT_bcH_cut, float m_trk_pT_cut, float m_trk_eta_cut, float m_trk_d0_cut, float m_trk_z0sinth_cut)
{
   std::cout<<"\n In DAOD_selector::setCuts"<<std::endl;
   jet_pT_infcut=m_jet_pT_infcut;
   jet_pT_supcut=m_jet_pT_supcut;
   jet_eta_cut=m_jet_eta_cut;
   jet_JVT_cut=m_jet_JVT_cut;
   jet_Isol_cut=1.;
   m_DR_bcH_truth_cut=m_DR_bcH_cut;
   m_pT_bcH_truth_cut=m_pT_bcH_cut;
   trk_pT_cut=m_trk_pT_cut;
   trk_eta_cut=m_trk_eta_cut;
   trk_d0_cut=m_trk_d0_cut;
   trk_z0sinth_cut=m_trk_z0sinth_cut;


   //// print cuts [
   //// ]
   return;
}

void DAOD_selector::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).


   std::cout<<"\n In DOAD_selector::Begin ... "<<std::endl;


   list.push_back(map_0);
     list.push_back(map_1);
       list.push_back(map_2);
         list.push_back(map_3);
           list.push_back(map_4);
             list.push_back(map_5);



   initFlagsAndCuts();
   initCounters();

   //tmp_pTfraction=0.,tmp_DR=0.,tmp_min_pTfraction=1.,tmp_min_DR=1.;

   m1=0,m2=0,m1_ex=0,m2_ex=0,mm1_ex=0,mm2_ex=0,m1_ov=0,m2_ov=0,mm=0;



   JF_ntrk=0,SV1_ntrk=0,SV0_ntrk=0,IP2D_ntrk=0,IP3D_ntrk=0;
   b_cnt=0,b_trkcut_cnt=0,trkcut_cnt=0;
   SV1input_trks=0,nSV1jets=0,nSV1trk=0,nSV1outputjets=0;
   nIPxDtrk=0,nIPxDoutputjets=0;
   nJFtrk=0,nJFoutputjets=0;

   //   std::cout<<"\n";
   openOutputFile(fOutputString);

   bookHistosForSelections();

   bookHistosForFlavorLabelStudies();

   bookHistosForDiscriminants();

   bookHistosForShrinkingCone();

   bookHistosForSelectionAlgos();


   TString option = GetOption();
   std::cout<<"Out of DOAD_selector::Begin ... "<<std::endl;
}

void DAOD_selector::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).
   std::cout<<"\n In  DOAD_selector::SlaveBegin ... "<<std::endl;
   TString option = GetOption();
   std::cout<<"Out of DOAD_selector::SlaveBegin ... "<<std::endl;
}

Bool_t DAOD_selector::Process(Long64_t entry)
{
   // The Process() function is called for each entry in the tree (or possibly
   // keyed object in the case of PROOF) to be processed. The entry argument
   // specifies which entry in the currently loaded tree is to be processed.
   // When processing keyed objects with PROOF, the object is already loaded
   // and is available via the fObject pointer.
   //
   // This function should contain the \"body\" of the analysis. It can contain
   // simple or elaborate selection criteria, run algorithms on the data
   // of the event and typically fill histograms.
   //
   // The processing can be stopped by calling Abort().
   //
   // Use fStatus to set the return value of TTree::Process().
   //
   // The return value is currently not used.

   fReader.SetLocalEntry(entry);

//   all TTreeReaderArray<int> and TTreeReaderArray<float> have the same dim, namely dim=*njets, the nuber of jets in each event.

   int nBjets=0,nCjets=0,nljets=0,nBjets_2=0,nCjets_2=0,nljets_2=0,nBcheck=0,nCcheck=0,nlcheck=0;
   int isSV1tagged=0;
   std::vector<int> isJet,isJet_OR,isBcheck,isCcheck,islcheck;

   //if (m_Ntot%1000==0) std::cout<<".... in Process ... events processed so far "<<m_Ntot<<"/"<<fReader.GetCurrentEntry()<<"  out of "<<fReader.GetEntries(true)<<" Run / event nb = "<<*runnb<<"/"<<*eventnb<<std::endl;
   if (m_Ntot%1000==0) std::cout<<".... in Process ... events processed so far "<<m_Ntot<<"/"<<fReader.GetCurrentEntry()<<" Run / event nb = "<<*runnb<<"/"<<*eventnb<<std::endl;

   m_Ntot++;

   double D_phi=0.,D_eta=0.,DR=0.,DR_bhadrjet=0,DR_chadrjet=0,px=0.,py=0.,pz=0.,Dx_1=0.,Dy_1=0.,Dz_1=0.,Dx_2=0.,Dy_2=0.,Dz_2=0,Lxy=0,Lxyz=0,Dxy_1=0.,x0=0.,y0=0.,Dx_3=0.,Dy_3=0.,Dxy_3=0.,rand_n=0.,R0=0.,child_z0=0.,trk_d0_signed=0.,child_d0_signed=0,trk_z0_signed=0.,child_z0_signed=0.,trk_z0sinth_signed=0.,child_z0sinth_signed=0;

   double D_phi_trk=0.,D_eta_trk=0.,DR_trk=0.,DpT_trk=0.;
   double tmp_pTfraction=0.,tmp_DR=0.,tmp_min_pTfraction=1.,tmp_min_DR=1.;
   int match=0,max_size=0;

   isJet.clear();
   isJet_OR.clear();

   //return kTRUE;
   // here loop over jets
   std::vector<int> isB_1,isB_2,isC_1,isC_2;
   Int_t nMyJets = *njets;
   //std::cout<<" N. of jets = "<<nMyJets<<std::endl; 
   for(int i=0;i<nMyJets;i++) {

     if(jet_nBHadr[i]>0){
       nBjets++;
       isB_1.push_back(i);
     }

     //     continue;
     if(jet_nCHadr[i]>0){
       nCjets++;
       isC_1.push_back(i);
     }

     if(jet_nBHadr[i]==0 && jet_nCHadr[i]==0){
       nljets++;
     }
     //   return kTRUE;



     if (jet_pt[i]<jet_pT_infcut) continue; //failing the pt cut (lower limit)
     m_njets_2_passPtMin++;
     if (jet_pt[i]>jet_pT_supcut) continue; //failing the pt cut (upper limit)
     m_njets_2_passPtMax++;
     if (abs(jet_eta[i])>jet_eta_cut) continue; //failing the pt cut
     m_njets_2_passEtaRange++;
     if (jet_isBadMedium[i] !=0 ) continue;  //failing is bad cut
     m_njets_2_passBadMedium++;
     if(jet_aliveAfterOR[i] <1  ) continue; // failing overlap removal
     m_njets_2_passOR++;
     if(jet_aliveAfterORmu[i]<1 ) continue; // failing overlap removal
     m_njets_2_passORmu++;

     // JVT cut
     if(jet_pt[i]<60*1e3 && abs(jet_eta[i])<2.4){
       if(jet_JVT[i]<jet_JVT_cut) continue;      // failing JVT
     }
     m_njets_2_passJVT++;
     isJet.push_back(i);         //// Here's the vector of jets passing the jet-selection




     if(jet_nBHadr[i]>0){
       nBjets_2++;
       isB_2.push_back(i);      //// Here's the vector of jets passing the jet-selection with >=1b hadron
     }

     if(jet_nCHadr[i]>0){
       nCjets_2++;
       isC_2.push_back(i);       //// Here's the vector of jets passing the jet-selection with >=1c hadron
     }
     if(jet_nBHadr[i]==0 && jet_nCHadr[i]==0){   //// Here's the vector of jets passing the jet-selection without any b and c hadron
       nljets_2++;
     }
   }


//number of B-C-l jets before jet cut
   m_njets+=*njets;
   m_nBjets+=nBjets;
   m_nCjets+=nCjets;
   m_nljets+=nljets;

//number of B-C-l jets after jet cut
   m_njets_2+=isJet.size();
   m_nBjets_2+=nBjets_2;
   m_nCjets_2+=nCjets_2;
   m_nljets_2+=nljets_2;

// number of jets with >=1 b hadron and >=1 c hadron
   m_nJetBCoverlap+=overlap(isB_1,isC_1);
   m_nJetBCoverlap_postJetSel+=overlap(isB_2,isC_2);


   OverlapRemoval(isJet, isJet_OR);//checks isolated jets with DR>1
   m_njets_2_passIsol=m_njets_2_passIsol+isJet_OR.size();


   if (fJetLabeling=="ConeLocalImpl")
     {
       // jet labelling (b-, c-, l-) based on truth (with pt, Deta cuts), Clusive samples
       getTrueJetFlavourLabel(isJet_OR, isBcheck, isCcheck, islcheck);
     }
   else if (fJetLabeling=="Cone")
     {
       // jet labelling (b-, c-, l-) fron DAOD flag HadronConeExclusiveExtendedTruthLabelID ===>>> get pure b/c-jets (no double hadrons) and light=no b/c hadrons (includes tau jets)
       getHadronConeExtFlavourLabel(isJet_OR, isBcheck, isCcheck, islcheck);
     }
   else if (fJetLabeling=="ConeIncl")
     {
       // jet labelling (b-, c-, l-) fron DAOD flag HadronConeExclusiveTruthLabelID ===>>> b/c-jets contain b/c hadrons (double hadrons allowed) and light=no b/c hadrons (includes tau jets)
       getHadronConeFlavourLabel(isJet_OR, isBcheck, isCcheck, islcheck);
     }
   else if (fJetLabeling=="Ghost")
     {
       // jet labelling (b-, c-, l-) fron DAOD flags  ===>>> ghostB/ChadronCount  ==> can select pure b/c-jets (no double hadrons) and light == 0 b/c hadrons
       getGhostExtJetFlavourLabel(isJet_OR, isBcheck, isCcheck, islcheck, false);
     }
   else if (fJetLabeling=="GhostIncl")
     {
       // jet labelling (b-, c-, l-) fron DAOD flags  ===>>> ghostB/ChadronCount  ==> can select pure b/c-jets (no double hadrons) and light == 0 b/c hadrons
       getGhostJetFlavourLabel(isJet_OR, isBcheck, isCcheck, islcheck, false);
     }
   else if (fJetLabeling=="GhostCone")
     {
       // jet labelling (b-, c-, l-) fron DAOD flags - matching labels --- exclusive labeling schemes 
       getGhostExtJetFlavourLabel(isJet_OR, isBcheck, isCcheck, islcheck, true);//diag_trms=false: not diagonal terms
     }
   else if (fJetLabeling=="GhostConeIncl")
     {
       // jet labelling (b-, c-, l-) fron DAOD flags - matching labels --- exclusive labeling schemes 
       getGhostJetFlavourLabel(isJet_OR, isBcheck, isCcheck, islcheck, true);//diag_trms=false: not diagonal terms
     }

   getFlavorLabelMatrix(isJet_OR);
   if (fDoFlavorLabelMatrix) getJetFeaturesInFlavorLabelMatrix();
   
   m_nBcheck+=isBcheck.size();
   m_nCcheck+=isCcheck.size();
   m_nlcheck+=islcheck.size();

   ov_check+=overlap(isBcheck,isCcheck);
   nBcheck=isBcheck.size();
   nCcheck=isCcheck.size();
   nlcheck=islcheck.size();

/*
   if(nB==0)  {m_noB++;}

   if (n1+2*n2+3*n3!=nB) {
        std::cout << "Warning: n1+n2+n3!=nB\t" << n1<<"\t"<<n2<<"\t"<<n3<<"\t"<<nB<<"\n";
    }
*/


   if(selections){
     unsigned size_jet=0;
     double RNNIP=0.,DL1=0.;

     if(nBcheck>0){

       for(std::vector<int>::iterator it = isBcheck.begin(); it != isBcheck.end(); ++it){
         b_cnt++;
         nSV1trk=0,nIPxDtrk=0,nJFtrk=0;

         hist_nBjets->Fill(1);
         size_jet=jet_trk_pt[*it].size();

         hist_pt_B->Fill(jet_pt[*it]*0.001);
         if(jet_bH_pt[*it].size()==1) hist_bHpt_B->Fill(jet_bH_pt[*it].at(0)*0.001);
         hist_eta_B->Fill(jet_eta[*it]);
         hist_phi_B->Fill(jet_phi[*it]);
         hist_E_B->Fill(jet_E[*it]);

         if(jet_ip2d_pb[*it]!=-99){
           hist_ip2d_llr_B->Fill(jet_ip2d_llr[*it]); //llr is computed as log(pb/pu)
           hist_ip2d_llr_jetpt_B->Fill(jet_ip2d_llr[*it],jet_pt[*it]*0.001);
           if(jet_bH_pt[*it].size()==1)
             hist_ip2d_llr_jetpt_singleB->Fill(jet_ip2d_llr[*it],jet_pt[*it]*0.001);
         }
         if(jet_ip3d_pb[*it]!=-99){
           hist_ip3d_llr_B->Fill(jet_ip3d_llr[*it]); //llr is computed as log(pb/pu)
           hist_ip3d_llr_jetpt_B->Fill(jet_ip3d_llr[*it],jet_pt[*it]*0.001);
           if(jet_bH_pt[*it].size()==1)
             hist_ip3d_llr_jetpt_singleB->Fill(jet_ip3d_llr[*it],jet_pt[*it]*0.001);
         }
         if(jet_rnnip_pb[*it]!=-99){
           RNNIP=log(jet_rnnip_pb[*it]/(m_fcRNNIP*jet_rnnip_pc[*it]+(1.-m_fcRNNIP)*jet_rnnip_pu[*it]));
           hist_rnnip_llr_B->Fill(RNNIP); //llr is computed as log(pb/pu)
           hist_rnnip_llr_jetpt_B->Fill(RNNIP,jet_pt[*it]*0.001);
           if(jet_bH_pt[*it].size()==1)
             hist_rnnip_llr_jetpt_singleB->Fill(RNNIP,jet_pt[*it]*0.001);
         }
         if(sv1_llr[*it]!=-99){
           hist_sv1_llr_B->Fill(sv1_llr[*it]); //llr is computed as log(pb/pu)
           hist_sv1_llr_jetpt_B->Fill(sv1_llr[*it],jet_pt[*it]*0.001);
           if(jet_bH_pt[*it].size()==1)
             hist_sv1_llr_jetpt_singleB->Fill(sv1_llr[*it],jet_pt[*it]*0.001);

           hist_jet_sv1_Nvtx_B->Fill(jet_sv1_Nvtx[*it]);
           hist_jet_sv1_ntrkv_B->Fill(jet_sv1_ntrkv[*it]);
           hist_jet_sv1_n2t_B->Fill(jet_sv1_n2t[*it]);
           hist_jet_sv1_m_B->Fill(jet_sv1_m[*it]);
           hist_jet_sv1_efc_B->Fill(jet_sv1_efc[*it]);
           hist_jet_sv1_sig3d_B->Fill(jet_sv1_sig3d[*it]);
           hist_jet_sv1_deltaR_B->Fill(jet_sv1_deltaR[*it]);
           hist_jet_sv1_Lxy_B->Fill(jet_sv1_Lxy[*it]);
           hist_jet_sv1_L3d_B->Fill(jet_sv1_L3d[*it]);
         }
         if(jet_jf_llr[*it]!=-99){
           hist_jf_llr_B->Fill(jet_jf_llr[*it]); //llr is computed as log(pb/pu)
           hist_jf_llr_jetpt_B->Fill(jet_jf_llr[*it],jet_pt[*it]*0.001);
           if(jet_bH_pt[*it].size()==1)
             hist_jf_llr_jetpt_singleB->Fill(jet_jf_llr[*it],jet_bH_pt[*it][0]*0.001);
         }
         if(jet_dl1_pb[*it]!=-99){
           DL1=log(jet_dl1_pb[*it]/(m_fc*jet_dl1_pc[*it]+(1-m_fc)*jet_dl1_pu[*it]));
           hist_dl1_B->Fill(DL1);
           hist_dl1_llr_jetpt_B->Fill(DL1,jet_pt[*it]*0.001);
           if(jet_bH_pt[*it].size()==1)
             hist_dl1_llr_jetpt_singleB->Fill(DL1,jet_bH_pt[*it][0]*0.001);
         }

         hist_n_tracks_jetpt_B->Fill(0.001*jet_pt[*it],size_jet);
         if(jet_bH_pt[*it].size()==1){
            hist_n_tracks_bHpt_B->Fill(0.001*jet_bH_pt[*it].at(0),size_jet);
         }
         int nBCcount=0;
         for(unsigned i=0;i<size_jet;i++){
           int origin=jet_trk_truth_derivedlabel[*it].at(i);
           if(origin==5){
             nBCcount++;
           }
         }
         hist_n_BCtracks_jetpt_B->Fill(0.001*jet_pt[*it],nBCcount);
         if(jet_bH_pt[*it].size()==1){
           hist_n_BCtracks_bHpt_B->Fill(0.001*jet_bH_pt[*it].at(0),nBCcount);
         }



         trkcut_cnt=0;
         for(unsigned i=0;i<size_jet;i++){
           trkcut_cnt++;

//           int origin=jet_trk_orig[*it].at(i);
           int origin_truth_label=jet_trk_truth_label[*it].at(i);
           int origin_derived_label=jet_trk_truth_derivedlabel[*it].at(i);
           int origin=origin_derived_label;

           hist_Bjet_origin->Fill(origin_derived_label);
           if(abs(jet_trk_eta[*it].at(i))<trk_eta_cut && jet_trk_pt[*it].at(i)>trk_pT_cut){// && abs(jet_trk_d0[*it].at(i))<trk_d0_cut && abs(jet_trk_z0[*it].at(i)*sin(jet_trk_theta[*it].at(i)))<trk_z0sinth_cut){//
             n_trk_pT_cut++;

             hist_Bjet_cut_origin->Fill(origin);
             hist_Bjet_cut_origin_pT->Fill(origin,1e-3*jet_trk_pt[*it].at(i));
             hist_Bjet_cut_origin_jetpT->Fill(origin,1e-3*jet_pt[*it]);
             hist_Bjet_cut_origin_truth_label->Fill(jet_trk_truth_label[*it].at(i));
             hist_Bjet_cut_origin_truth_label_pT->Fill(jet_trk_truth_label[*it].at(i),1e-3*jet_pt[*it]);
             if(jet_bH_pt[*it].size()==1){
                hist_Bjet_cut_origin_truth_label_bHpT->Fill(jet_trk_truth_label[*it].at(i),1e-3*jet_bH_pt[*it].at(0));
             }

             D_eta=jet_trk_eta[*it][i]-jet_eta[*it];
             if(abs(jet_trk_phi[*it].at(i)-jet_phi[*it])>M_PI){
               D_phi=2*M_PI-abs(jet_trk_phi[*it].at(i)-jet_phi[*it]);
             }
             if(abs(jet_trk_phi[*it].at(i)-jet_phi[*it])<M_PI){
               D_phi=jet_trk_phi[*it].at(i)-jet_phi[*it];
             }
             DR=sqrt(D_eta*D_eta+D_phi*D_phi);



             hist_trk_pT_B->Fill(1e-3*jet_trk_pt[*it].at(i));
             hist_trk_eta_B->Fill(jet_trk_eta[*it].at(i));
             hist_trk_phi_B->Fill(jet_trk_phi[*it].at(i));
             hist_trk_Deta_B->Fill(D_eta);
             hist_trk_Dphi_B->Fill(D_phi);
             hist_trk_Dphi_Deta_B->Fill(D_phi,D_eta);
             hist_trk_DR_B->Fill(DR);
             hist_trk_pT_DR_B->Fill(1e-3*jet_trk_pt[*it].at(i),DR);
             hist_trk_pT_jet_DR_B->Fill(1e-3*jet_pt[*it],DR);
             hist_trk_pdgId_B->Fill(jet_trk_pdg_id[*it].at(i));
             hist_trk_origin_B->Fill(origin);
             hist_trk_origin_truth_label_B->Fill(origin_truth_label);

             A=sin(jet_phi[*it]-jet_trk_phi[*it].at(i))*jet_trk_d0[*it].at(i);
             sgn_d0=A/abs(A);
             A=(jet_eta[*it]-jet_trk_eta[*it].at(i))*jet_trk_z0[*it].at(i);
             sgn_z0=A/abs(A);
             trk_d0_signed=sgn_d0*abs(jet_trk_d0[*it].at(i));
             trk_z0_signed=sgn_z0*abs(jet_trk_z0[*it].at(i));
             trk_z0sinth_signed=trk_z0_signed*sin(jet_trk_theta[*it].at(i));

             hist_trk_d0_B->Fill(trk_d0_signed);
             hist_trk_d0_jet_pT_BC_B->Fill(1e-3*jet_pt[*it],trk_d0_signed);
             hist_trk_DR_jet_pt_BC->Fill(1e-3*jet_pt[*it],DR);
             hist_trk_sigd0_B->Fill(jet_trk_d0Err[*it].at(i));
             hist_trk_z0_B->Fill(trk_z0_signed);
             hist_trk_sigz0_B->Fill(jet_trk_z0Err[*it].at(i));
             hist_trk_z0sinth_B->Fill(trk_z0sinth_signed);
             hist_trk_d0sig_B->Fill(trk_d0_signed/jet_trk_d0Err[*it].at(i));
             hist_trk_z0sinthsig_B->Fill(trk_z0sinth_signed/jet_trk_z0Err[*it].at(i));
             hist_trk_d0sig_origin_B->Fill(trk_d0_signed/jet_trk_d0Err[*it].at(i),origin);
             hist_trk_d0sig_truthlabel_B->Fill(trk_d0_signed/jet_trk_d0Err[*it].at(i),origin_truth_label);
             hist_trk_z0sinthsig_origin_B->Fill(trk_z0sinth_signed/jet_trk_z0Err[*it].at(i),origin);
             hist_trk_logpTfrac_origin_B->Fill(log(jet_trk_pt[*it].at(i)/jet_pt[*it]),origin);
             hist_trk_logDR_origin_B->Fill(log(DR),origin);
             hist_trk_IBLhits_origin_B->Fill(jet_trk_nInnHits[*it].at(i),origin);
             hist_trk_NextToIBLhits_origin_B->Fill(jet_trk_nNextToInnHits[*it].at(i),origin);
             hist_trk_sharedIBLhits_origin_B->Fill(jet_trk_nsharedBLHits[*it].at(i),origin);
             hist_trk_splitIBLhits_origin_B->Fill(jet_trk_nsplitBLHits[*it].at(i),origin);//jet_trk_nsplitPixHits
             hist_trk_nPixhits_origin_B->Fill(jet_trk_nPixHits[*it].at(i),origin);
             hist_trk_sharedPixhits_origin_B->Fill(jet_trk_nsharedPixHits[*it].at(i),origin);
             hist_trk_splitPixhits_origin_B->Fill(jet_trk_nsplitPixHits[*it].at(i),origin);
             hist_trk_nSCThits_origin_B->Fill(jet_trk_nSCTHits[*it].at(i),origin);
             hist_trk_sharedSCThits_origin_B->Fill(jet_trk_nsharedSCTHits[*it].at(i),origin);

//               std::cout<<jet_trk_algo[*it].at(i)<<"\n";
             std::vector<int> bin_alg(5);
             int x=jet_trk_algo[*it].at(i);
//SV1 START
             if(jet_trk_nSCTHits[*it].at(i)>=4 && jet_trk_nPixHits[*it].at(i)>=1 && jet_trk_nPixHits[*it].at(i)+jet_trk_nSCTHits[*it].at(i)>=7 && jet_trk_nBLHits[*it].at(i)>=0 && jet_trk_pt[*it].at(i)>=700 && jet_trk_z0[*it].at(i)<=25 && jet_trk_z0Err[*it].at(i)<=5 && abs(jet_trk_d0[*it].at(i))<=5 && jet_trk_d0Err[*it].at(i)<=1 && DR<=0.4 && jet_trk_chi2[*it].at(i)>=3){
                isSV1tagged=1;
                SV1input_trks+=1;
                hist_trk_SV1input_origin_B->Fill(origin);

             }
//END SV1


             int idx=0,p=0;
             while(idx<5){
               p=x%2;
               bin_alg.at(idx)=p;
               if(idx==4) {
                 JF_ntrk=JF_ntrk+p;
                 nJFtrk+=p;
                 if(p==1){
                   hist_trk_pT_JF_B->Fill(1e-3*jet_trk_pt[*it].at(i));
                   hist_jet_pT_origin_JF_B->Fill(origin,1e-3*jet_pt[*it]);
                   hist_jet_pT_origin_truth_label_JF_B->Fill(origin_truth_label,1e-3*jet_pt[*it]);

                   hist_trk_eta_JF_B->Fill(jet_trk_eta[*it].at(i));
                   hist_trk_pT_jet_DR_JF_B->Fill(1e-3*jet_pt[*it],DR);
                   hist_trk_d0_JF_B->Fill(trk_d0_signed);
                   hist_trk_z0sinth_JF_B->Fill(trk_z0sinth_signed);

                   hist_trk_origin_JF_B->Fill(origin);
                   hist_trk_origin_truth_label_JF_B->Fill(origin_truth_label);
                   hist_trk_d0sig_JF_B->Fill(trk_d0_signed/jet_trk_d0Err[*it].at(i));
                   hist_trk_z0sinthsig_JF_B->Fill(trk_z0sinth_signed/jet_trk_z0Err[*it].at(i));
                   hist_trk_d0sig_origin_JF_B->Fill(trk_d0_signed/jet_trk_d0Err[*it].at(i),origin);
                   hist_trk_z0sinthsig_origin_JF_B->Fill(trk_z0sinth_signed/jet_trk_z0Err[*it].at(i),origin);
                   hist_trk_logpTfrac_origin_JF_B->Fill(log(jet_trk_pt[*it].at(i)/jet_pt[*it]),origin);
                   hist_trk_logDR_origin_JF_B->Fill(log(DR),origin);
                   hist_trk_IBLhits_origin_JF_B->Fill(jet_trk_nInnHits[*it].at(i),origin);
                   hist_trk_NextToIBLhits_origin_JF_B->Fill(jet_trk_nNextToInnHits[*it].at(i),origin);
                   hist_trk_sharedIBLhits_origin_JF_B->Fill(jet_trk_nsharedBLHits[*it].at(i),origin);
                   hist_trk_splitIBLhits_origin_JF_B->Fill(jet_trk_nsplitBLHits[*it].at(i),origin);//jet_trk_nsplitPixHits
                   hist_trk_nPixhits_origin_JF_B->Fill(jet_trk_nPixHits[*it].at(i),origin);
                   hist_trk_sharedPixhits_origin_JF_B->Fill(jet_trk_nsharedPixHits[*it].at(i),origin);
                   hist_trk_splitPixhits_origin_JF_B->Fill(jet_trk_nsplitPixHits[*it].at(i),origin);
                   hist_trk_nSCThits_origin_JF_B->Fill(jet_trk_nSCTHits[*it].at(i),origin);
                   hist_trk_sharedSCThits_origin_JF_B->Fill(jet_trk_nsharedSCTHits[*it].at(i),origin);
                 }
               }
               if(idx==3) {
                 SV1_ntrk=SV1_ntrk+p;
                 nSV1trk+=p;
                 if(p==1){
                   hist_trk_pT_SV1_B->Fill(1e-3*jet_trk_pt[*it].at(i));
                   hist_jet_pT_origin_SV1_B->Fill(origin,1e-3*jet_pt[*it]);
                   hist_jet_pT_origin_truth_label_SV1_B->Fill(origin_truth_label,1e-3*jet_pt[*it]);
                   if(jet_bH_pt[*it].size()==1) hist_bH_pT_origin_truth_label_SV1_B->Fill(origin_truth_label,1e-3*jet_bH_pt[*it].at(0));

                   hist_trk_eta_SV1_B->Fill(jet_trk_eta[*it].at(i));
                   hist_trk_pT_jet_DR_SV1_B->Fill(1e-3*jet_pt[*it],DR);
                   hist_trk_d0_SV1_B->Fill(trk_d0_signed);
                   hist_trk_z0sinth_SV1_B->Fill(trk_z0sinth_signed);

                   hist_trk_origin_SV1_B->Fill(origin);
                   hist_trk_origin_truth_label_SV1_B->Fill(origin_truth_label);
                   hist_trk_d0sig_SV1_B->Fill(trk_d0_signed/jet_trk_d0Err[*it].at(i));
                   hist_trk_z0sinthsig_SV1_B->Fill(trk_z0sinth_signed/jet_trk_z0Err[*it].at(i));
                   hist_trk_d0sig_origin_SV1_B->Fill(trk_d0_signed/jet_trk_d0Err[*it].at(i),origin);
                   hist_trk_z0sinthsig_origin_SV1_B->Fill(trk_z0sinth_signed/jet_trk_z0Err[*it].at(i),origin);
                   hist_trk_logpTfrac_origin_SV1_B->Fill(log(jet_trk_pt[*it].at(i)/jet_pt[*it]),origin);
                   hist_trk_logDR_origin_SV1_B->Fill(log(DR),origin);
                   hist_trk_IBLhits_origin_SV1_B->Fill(jet_trk_nInnHits[*it].at(i),origin);
                   hist_trk_NextToIBLhits_origin_SV1_B->Fill(jet_trk_nNextToInnHits[*it].at(i),origin);
                   hist_trk_sharedIBLhits_origin_SV1_B->Fill(jet_trk_nsharedBLHits[*it].at(i),origin);
                   hist_trk_splitIBLhits_origin_SV1_B->Fill(jet_trk_nsplitBLHits[*it].at(i),origin);//jet_trk_nsplitPixHits
                   hist_trk_nPixhits_origin_SV1_B->Fill(jet_trk_nPixHits[*it].at(i),origin);
                   hist_trk_sharedPixhits_origin_SV1_B->Fill(jet_trk_nsharedPixHits[*it].at(i),origin);
                   hist_trk_splitPixhits_origin_SV1_B->Fill(jet_trk_nsplitPixHits[*it].at(i),origin);
                   hist_trk_nSCThits_origin_SV1_B->Fill(jet_trk_nSCTHits[*it].at(i),origin);
                   hist_trk_sharedSCThits_origin_SV1_B->Fill(jet_trk_nsharedSCTHits[*it].at(i),origin);
                   hist_trk_chi2_SV1_B->Fill(jet_trk_chi2[*it].at(i));


                 }
               }
               if(idx==2) {
                 SV0_ntrk=SV0_ntrk+p;
                 if(p==1){
                   hist_trk_pT_SV0_B->Fill(1e-3*jet_trk_pt[*it].at(i));
                   hist_trk_eta_SV0_B->Fill(jet_trk_eta[*it].at(i));
                   hist_trk_pT_jet_DR_SV0_B->Fill(1e-3*jet_pt[*it],DR);
                   hist_trk_d0_SV0_B->Fill(trk_d0_signed);
                   hist_trk_z0sinth_SV0_B->Fill(trk_z0sinth_signed);

                   hist_trk_origin_SV0_B->Fill(origin);
                   hist_trk_d0sig_SV0_B->Fill(trk_d0_signed/jet_trk_d0Err[*it].at(i));
                   hist_trk_z0sinthsig_SV0_B->Fill(trk_z0sinth_signed/jet_trk_z0Err[*it].at(i));
                   hist_trk_d0sig_origin_SV0_B->Fill(trk_d0_signed/jet_trk_d0Err[*it].at(i),origin);
                   hist_trk_z0sinthsig_origin_SV0_B->Fill(trk_z0sinth_signed/jet_trk_z0Err[*it].at(i),origin);
                   hist_trk_logpTfrac_origin_SV0_B->Fill(log(jet_trk_pt[*it].at(i)/jet_pt[*it]),origin);
                   hist_trk_logDR_origin_SV0_B->Fill(log(DR),origin);
                   hist_trk_IBLhits_origin_SV0_B->Fill(jet_trk_nInnHits[*it].at(i),origin);
                   hist_trk_NextToIBLhits_origin_SV0_B->Fill(jet_trk_nNextToInnHits[*it].at(i),origin);
                   hist_trk_sharedIBLhits_origin_SV0_B->Fill(jet_trk_nsharedBLHits[*it].at(i),origin);
                   hist_trk_splitIBLhits_origin_SV0_B->Fill(jet_trk_nsplitBLHits[*it].at(i),origin);//jet_trk_nsplitPixHits
                   hist_trk_nPixhits_origin_SV0_B->Fill(jet_trk_nPixHits[*it].at(i),origin);
                   hist_trk_sharedPixhits_origin_SV0_B->Fill(jet_trk_nsharedPixHits[*it].at(i),origin);
                   hist_trk_splitPixhits_origin_SV0_B->Fill(jet_trk_nsplitPixHits[*it].at(i),origin);
                   hist_trk_nSCThits_origin_SV0_B->Fill(jet_trk_nSCTHits[*it].at(i),origin);
                   hist_trk_sharedSCThits_origin_SV0_B->Fill(jet_trk_nsharedSCTHits[*it].at(i),origin);
                 }
               }
               if(idx==1) {
                 IP3D_ntrk=IP3D_ntrk+p;
                 nIPxDtrk+=p;
                 if(p==1){
                   hist_trk_pT_IP3D_B->Fill(1e-3*jet_trk_pt[*it].at(i));
                   hist_jet_pT_origin_IP3D_B->Fill(origin,1e-3*jet_pt[*it]);
                   hist_jet_pT_origin_truth_label_IP3D_B->Fill(origin_truth_label,1e-3*jet_pt[*it]);
                   if(jet_bH_pt[*it].size()==1) hist_bH_pT_origin_truth_label_IP3D_B->Fill(origin_truth_label,1e-3*jet_bH_pt[*it].at(0));

                   hist_trk_eta_IP3D_B->Fill(jet_trk_eta[*it].at(i));
                   hist_trk_pT_jet_DR_IP3D_B->Fill(1e-3*jet_pt[*it],DR);
                   hist_trk_d0_IP3D_B->Fill(trk_d0_signed);
                   hist_trk_z0sinth_IP3D_B->Fill(trk_z0sinth_signed);

                   hist_trk_origin_IP3D_B->Fill(origin);
                   hist_trk_origin_truth_label_IP3D_B->Fill(origin_truth_label);
                   hist_trk_d0sig_IP3D_B->Fill(trk_d0_signed/jet_trk_d0Err[*it].at(i));
                   hist_trk_z0sinthsig_IP3D_B->Fill(trk_z0sinth_signed/jet_trk_z0Err[*it].at(i));
                   hist_trk_d0sig_origin_IP3D_B->Fill(trk_d0_signed/jet_trk_d0Err[*it].at(i),origin);
                   hist_trk_z0sinthsig_origin_IP3D_B->Fill(trk_z0sinth_signed/jet_trk_z0Err[*it].at(i),origin);
                   hist_trk_logpTfrac_origin_IP3D_B->Fill(log(jet_trk_pt[*it].at(i)/jet_pt[*it]),origin);
                   hist_trk_logDR_origin_IP3D_B->Fill(log(DR),origin);
                   hist_trk_IBLhits_origin_IP3D_B->Fill(jet_trk_nInnHits[*it].at(i),origin);
                   hist_trk_NextToIBLhits_origin_IP3D_B->Fill(jet_trk_nNextToInnHits[*it].at(i),origin);
                   hist_trk_sharedIBLhits_origin_IP3D_B->Fill(jet_trk_nsharedBLHits[*it].at(i),origin);
                   hist_trk_splitIBLhits_origin_IP3D_B->Fill(jet_trk_nsplitBLHits[*it].at(i),origin);//jet_trk_nsplitPixHits
                   hist_trk_nPixhits_origin_IP3D_B->Fill(jet_trk_nPixHits[*it].at(i),origin);
                   hist_trk_sharedPixhits_origin_IP3D_B->Fill(jet_trk_nsharedPixHits[*it].at(i),origin);
                   hist_trk_splitPixhits_origin_IP3D_B->Fill(jet_trk_nsplitPixHits[*it].at(i),origin);
                   hist_trk_nSCThits_origin_IP3D_B->Fill(jet_trk_nSCTHits[*it].at(i),origin);
                   hist_trk_sharedSCThits_origin_IP3D_B->Fill(jet_trk_nsharedSCTHits[*it].at(i),origin);
                 }
               }
               if(idx==0) {
                 IP2D_ntrk=IP2D_ntrk+p;
                 if(p==1){
                   hist_trk_pT_IP2D_B->Fill(1e-3*jet_trk_pt[*it].at(i));
                   hist_trk_eta_IP2D_B->Fill(jet_trk_eta[*it].at(i));
                   hist_trk_pT_jet_DR_IP2D_B->Fill(1e-3*jet_pt[*it],DR);
                   hist_trk_d0_IP2D_B->Fill(trk_d0_signed);
                   hist_trk_z0sinth_IP2D_B->Fill(trk_z0sinth_signed);

                   hist_trk_origin_IP2D_B->Fill(origin);
                   hist_trk_d0sig_IP2D_B->Fill(trk_d0_signed/jet_trk_d0Err[*it].at(i));
                   hist_trk_z0sinthsig_IP2D_B->Fill(trk_z0sinth_signed/jet_trk_z0Err[*it].at(i));
                   hist_trk_d0sig_origin_IP2D_B->Fill(trk_d0_signed/jet_trk_d0Err[*it].at(i),origin);
                   hist_trk_z0sinthsig_origin_IP2D_B->Fill(trk_z0sinth_signed/jet_trk_z0Err[*it].at(i),origin);
                   hist_trk_logpTfrac_origin_IP2D_B->Fill(log(jet_trk_pt[*it].at(i)/jet_pt[*it]),origin);
                   hist_trk_logDR_origin_IP2D_B->Fill(log(DR),origin);
                   hist_trk_IBLhits_origin_IP2D_B->Fill(jet_trk_nInnHits[*it].at(i),origin);
                   hist_trk_NextToIBLhits_origin_IP2D_B->Fill(jet_trk_nNextToInnHits[*it].at(i),origin);
                   hist_trk_sharedIBLhits_origin_IP2D_B->Fill(jet_trk_nsharedBLHits[*it].at(i),origin);
                   hist_trk_splitIBLhits_origin_IP2D_B->Fill(jet_trk_nsplitBLHits[*it].at(i),origin);//jet_trk_nsplitPixHits
                   hist_trk_nPixhits_origin_IP2D_B->Fill(jet_trk_nPixHits[*it].at(i),origin);
                   hist_trk_sharedPixhits_origin_IP2D_B->Fill(jet_trk_nsharedPixHits[*it].at(i),origin);
                   hist_trk_splitPixhits_origin_IP2D_B->Fill(jet_trk_nsplitPixHits[*it].at(i),origin);
                   hist_trk_nSCThits_origin_IP2D_B->Fill(jet_trk_nSCTHits[*it].at(i),origin);
                   hist_trk_sharedSCThits_origin_IP2D_B->Fill(jet_trk_nsharedSCTHits[*it].at(i),origin);
                 }
               }
               x=x/2;
               idx++;
             }
/*
             std::cout<<"\n";
             for (unsigned i=0;i<5;i++)
                std::cout << bin_alg.at(i) << ' ';
             */
/*
             std::cout<<"\n"<<jet_trk_algo[*it].at(i)<<"\t(";
             for(int z=bin_alg.size()-1;z>=0;z--){
               if(z>0)
                 std::cout<<bin_alg.at(z)<<",";
               if(z==0)
                 std::cout<<bin_alg.at(z);
             }
             std::cout<<")";
*/
/*
             if(origin==-1){
               hist_trk_d0_PUB->Fill(trk_d0_signed);
               n_trk_PU_pT_cut++;
             }
             if(origin==0){
               hist_trk_d0_BB->Fill(trk_d0_signed);
               n_trk_B++;
             }
             if(origin==1){
               hist_trk_d0_CB->Fill(trk_d0_signed);
               n_trk_C++;
             }
             if(origin==2){
               hist_trk_d0_FRAGB->Fill(trk_d0_signed);
               n_trk_FRAG_pT_cut++;
             }
             if(origin==3){
               hist_trk_d0_GEANTB->Fill(trk_d0_signed);
               n_trk_GEANT_pT_cut++;
             }
*/

             if(origin_selection){
               if(origin==5){
                 hist_matched_origin_pT_B->Fill(1e-3*jet_trk_pt[*it].at(i));
                 hist_matched_origin_eta_B->Fill(jet_trk_eta[*it].at(i));
                 hist_matched_origin_phi_B->Fill(jet_trk_phi[*it].at(i));
                 hist_matched_origin_Deta_B->Fill(D_eta);
                 hist_matched_origin_Dphi_B->Fill(D_phi);
                 hist_matched_origin_Dphi_Deta_B->Fill(D_phi,D_eta);
                 hist_matched_origin_DR_B->Fill(DR);
                 hist_matched_origin_pT_DR_B->Fill(1e-3*jet_trk_pt[*it].at(i),DR);
                 hist_matched_origin_pT_jet_DR_B->Fill(1e-3*jet_pt[*it],DR);
                 hist_matched_origin_pdgId_B->Fill(jet_trk_pdg_id[*it].at(i));
                 hist_matched_origin_origin_B->Fill(origin);
                 hist_matched_origin_d0_B->Fill(trk_d0_signed);
//                   hist_matched_origin_Lxy_B->Fill();
//                   hist_matched_origin_Lxyz_B->Fill();
                hist_matched_origin_jetpT_B->Fill(1e-3*jet_pt[*it]);
                if(jet_bH_pt[*it].size()==1){

                  D_eta=jet_bH_eta[*it].at(0)-jet_eta[*it];
                  if(abs(jet_bH_phi[*it].at(0)-jet_phi[*it])>M_PI){
                    D_phi=2*M_PI-abs(jet_bH_phi[*it].at(0)-jet_phi[*it]);
                  }
                  if(abs(jet_bH_phi[*it].at(0)-jet_phi[*it])<M_PI){
                    D_phi=jet_bH_phi[*it].at(0)-jet_phi[*it];
                  }
                  DR_bhadrjet=sqrt(D_eta*D_eta+D_phi*D_phi);
                  hist_trk_DR_jetpt_B->Fill(1e-3*jet_pt[*it],DR_bhadrjet);
                  hist_trk_DR_bHpt_B->Fill(1e-3*jet_bH_pt[*it].at(0),DR_bhadrjet);
                  hist_matched_origin_bHpT_B->Fill(1e-3*jet_bH_pt[*it].at(0));
                }
               }

             }


           }
         }
         if(trkcut_cnt>0) b_trkcut_cnt+=1;
         if(isSV1tagged>0) nSV1jets+=1;

         if(nSV1trk>0){
           nSV1outputjets++;
           hist_jet_pt_SV1_B->Fill(1e-3*jet_pt[*it]);
         }
         if(nIPxDtrk>0){
           nIPxDoutputjets++;
           hist_jet_pt_IP3D_B->Fill(1e-3*jet_pt[*it]);
         }
         if(nJFtrk>0){
           nJFoutputjets++;
           hist_jet_pt_JF_B->Fill(1e-3*jet_pt[*it]);
         }
       }
       //hist_n_trackstrkcut_cnt

     }

     if(nCcheck>0){
       for(std::vector<int>::iterator it = isCcheck.begin(); it != isCcheck.end(); ++it){
         size_jet=jet_trk_pt[*it].size();
         hist_nCjets->Fill(1);

         hist_pt_C->Fill(jet_pt[*it]*0.001);
         if(jet_cH_pt[*it].size()==1) hist_cHpt_C->Fill(jet_cH_pt[*it].at(0)*0.001);
         hist_eta_C->Fill(jet_eta[*it]);
         hist_phi_C->Fill(jet_phi[*it]);
         hist_E_C->Fill(jet_E[*it]);

         if(jet_ip2d_pb[*it]!=-99){
           hist_ip2d_llr_C->Fill(jet_ip2d_llr[*it]); //llr is computed as log(pb/pu)
           hist_ip2d_llr_jetpt_C->Fill(jet_ip2d_llr[*it],jet_pt[*it]*0.001);
             if(jet_cH_pt[*it].size()==1)
               hist_ip2d_llr_jetpt_singleC->Fill(jet_ip2d_llr[*it],jet_cH_pt[*it][0]*0.001);
         }

         if(jet_ip3d_pb[*it]!=-99){
           hist_ip3d_llr_C->Fill(jet_ip3d_llr[*it]); //llr is computed as log(pb/pu)
           hist_ip3d_llr_jetpt_C->Fill(jet_ip3d_llr[*it],jet_pt[*it]*0.001);
           if(jet_cH_pt[*it].size()==1)
             hist_ip3d_llr_jetpt_singleC->Fill(jet_ip3d_llr[*it],jet_cH_pt[*it][0]*0.001);
         }
         if(jet_rnnip_pb[*it]!=-99){
           RNNIP=log(jet_rnnip_pb[*it]/(m_fcRNNIP*jet_rnnip_pc[*it]+(1.-m_fcRNNIP)*jet_rnnip_pu[*it]));
           hist_rnnip_llr_C->Fill(RNNIP); //llr is computed as log(pb/pu)
           hist_rnnip_llr_jetpt_C->Fill(RNNIP,jet_pt[*it]*0.001);
           if(jet_cH_pt[*it].size()==1)
             hist_rnnip_llr_jetpt_singleC->Fill(RNNIP,jet_cH_pt[*it][0]*0.001);
         }
         if(sv1_llr[*it]!=-99){
           hist_sv1_llr_C->Fill(sv1_llr[*it]); //llr is computed as log(pb/pu)
           hist_sv1_llr_jetpt_C->Fill(sv1_llr[*it],jet_pt[*it]*0.001);
           if(jet_cH_pt[*it].size()==1)
             hist_sv1_llr_jetpt_singleC->Fill(sv1_llr[*it],jet_cH_pt[*it][0]*0.001);
         }
         if(jet_jf_llr[*it]!=-99){
           hist_jf_llr_C->Fill(jet_jf_llr[*it]); //llr is computed as log(pb/pu)
           hist_jf_llr_jetpt_C->Fill(jet_jf_llr[*it],jet_pt[*it]*0.001);
           if(jet_cH_pt[*it].size()==1)
             hist_jf_llr_jetpt_singleC->Fill(jet_jf_llr[*it],jet_cH_pt[*it][0]*0.001);
         }
         if(jet_dl1_pb[*it]!=-99){
           DL1=log(jet_dl1_pb[*it]/(m_fc*jet_dl1_pc[*it]+(1-m_fc)*jet_dl1_pu[*it]));
           hist_dl1_C->Fill(DL1);
           hist_dl1_llr_jetpt_C->Fill(DL1,jet_pt[*it]*0.001);
           if(jet_cH_pt[*it].size()==1)
             hist_dl1_llr_jetpt_singleC->Fill(DL1,jet_cH_pt[*it][0]*0.001);
         }


         hist_n_tracks_jetpt_C->Fill(0.001*jet_pt[*it],size_jet);
         if(jet_cH_pt[*it].size()==1){
            hist_n_tracks_cHpt_C->Fill(0.001*jet_cH_pt[*it].at(0),size_jet);
         }
         int nBCcount=0;
         for(unsigned i=0;i<size_jet;i++){
           int origin=jet_trk_truth_derivedlabel[*it].at(i);
           if(origin==5){
             nBCcount++;
           }
         }
         hist_n_BCtracks_jetpt_C->Fill(0.001*jet_pt[*it],nBCcount);
         if(jet_cH_pt[*it].size()==1){
           hist_n_BCtracks_cHpt_C->Fill(0.001*jet_cH_pt[*it].at(0),nBCcount);
         }



         for(unsigned i=0;i<size_jet;i++){
//           int origin=jet_trk_orig[*it].at(i);
          int origin_truth_label=jet_trk_truth_label[*it].at(i);
          int origin_derived_label=jet_trk_truth_derivedlabel[*it].at(i);
          int origin=origin_derived_label;

           hist_Cjet_origin->Fill(origin);
           if(abs(jet_trk_eta[*it].at(i))<trk_eta_cut && jet_trk_pt[*it].at(i)>trk_pT_cut){// && abs(jet_trk_d0[*it].at(i))<trk_d0_cut && abs(jet_trk_z0[*it].at(i)*sin(jet_trk_theta[*it].at(i)))<trk_z0sinth_cut){//
             hist_Cjet_cut_origin->Fill(origin);
             hist_Cjet_cut_origin_pT->Fill(origin,1e-3*jet_trk_pt[*it].at(i));
             hist_Cjet_cut_origin_jetpT->Fill(origin,1e-3*jet_pt[*it]);
             hist_Cjet_cut_origin_truth_label->Fill(jet_trk_truth_label[*it].at(i));
             hist_Cjet_cut_origin_truth_label_pT->Fill(jet_trk_truth_label[*it].at(i),1e-3*jet_pt[*it]);
             if(jet_cH_pt[*it].size()==1){
                hist_Cjet_cut_origin_truth_label_cHpT->Fill(jet_trk_truth_label[*it].at(i),1e-3*jet_cH_pt[*it].at(0));
             }
           }

           if(origin==5){
             hist_matched_origin_jetpT_C->Fill(1e-3*jet_pt[*it]);
             if(jet_cH_pt[*it].size()==1){

               D_eta=jet_cH_eta[*it].at(0)-jet_eta[*it];
               if(abs(jet_cH_phi[*it].at(0)-jet_phi[*it])>M_PI){
                 D_phi=2*M_PI-abs(jet_cH_phi[*it].at(0)-jet_phi[*it]);
               }
               if(abs(jet_cH_phi[*it].at(0)-jet_phi[*it])<M_PI){
                 D_phi=jet_cH_phi[*it].at(0)-jet_phi[*it];
               }
               DR_chadrjet=sqrt(D_eta*D_eta+D_phi*D_phi);
               hist_trk_DR_jetpt_C->Fill(1e-3*jet_pt[*it],DR_chadrjet);
               hist_trk_DR_cHpt_C->Fill(1e-3*jet_cH_pt[*it].at(0),DR_chadrjet);
               hist_matched_origin_cHpT_C->Fill(1e-3*jet_cH_pt[*it].at(0));
             }
            }

         }
       }
     }

     if(nlcheck>0){
       for(std::vector<int>::iterator it = islcheck.begin(); it != islcheck.end(); ++it){
         if(jet_nBHadr[*it]>0)
           std::cout<<"b hadron in light jet\n";
         size_jet=jet_trk_pt[*it].size();
         for(unsigned i=0;i<size_jet;i++){
           if(jet_trk_truth_label[*it].at(i)==13 || jet_trk_truth_label[*it].at(i)==103 || jet_trk_truth_label[*it].at(i)==113 || jet_trk_truth_label[*it].at(i)==15 || jet_trk_truth_label[*it].at(i)==105 || jet_trk_truth_label[*it].at(i)==115){
             if(jet_trk_truth_label[*it].at(i)==13) cnt_13+=1;
             if(jet_trk_truth_label[*it].at(i)==103) cnt_103+=1;
             if(jet_trk_truth_label[*it].at(i)==113) cnt_113+=1;
             if(jet_trk_truth_label[*it].at(i)==15) cnt_15+=1;
             if(jet_trk_truth_label[*it].at(i)==105) cnt_105+=1;
             if(jet_trk_truth_label[*it].at(i)==115) cnt_115+=1;
           }
         }
         hist_nljets->Fill(1);

         hist_pt_l->Fill(jet_pt[*it]*0.001);
         hist_eta_l->Fill(jet_eta[*it]);
         hist_phi_l->Fill(jet_phi[*it]);
         hist_E_l->Fill(jet_E[*it]);

         if(jet_ip2d_pb[*it]!=-99){
           hist_ip2d_llr_l->Fill(jet_ip2d_llr[*it]); //llr is computed as log(pb/pu)
           if(jet_ip2d_llr[*it]!=jet_ip2[*it])
             std::cout<<"W\n";
         }
         if(jet_ip3d_pb[*it]!=-99){
           hist_ip3d_llr_l->Fill(jet_ip3d_llr[*it]); //llr is computed as log(pb/pu)
         }
         if(jet_rnnip_pb[*it]!=-99){
           RNNIP=log(jet_rnnip_pb[*it]/(m_fcRNNIP*jet_rnnip_pc[*it]+(1.-m_fcRNNIP)*jet_rnnip_pu[*it]));
           hist_rnnip_llr_l->Fill(RNNIP); //llr is computed as log(pb/pu)
         }
         if(sv1_llr[*it]!=-99){
           hist_sv1_llr_l->Fill(sv1_llr[*it]); //llr is computed as log(pb/pu)
         }
         if(jet_jf_llr[*it]!=-99){
           hist_jf_llr_l->Fill(jet_jf_llr[*it]); //llr is computed as log(pb/pu)
         }
         if(jet_dl1_pb[*it]!=-99){
           DL1=log(jet_dl1_pb[*it]/(m_fc*jet_dl1_pc[*it]+(1-m_fc)*jet_dl1_pu[*it]));
           hist_dl1_l->Fill(DL1);
         }


         hist_n_tracks_jetpt_l->Fill(0.001*jet_pt[*it],size_jet);


         for(unsigned i=0;i<size_jet;i++){
//           int origin=jet_trk_orig[*it].at(i);
          int origin_truth_label=jet_trk_truth_label[*it].at(i);
          int origin_derived_label=jet_trk_truth_derivedlabel[*it].at(i);
          int origin=origin_derived_label;
           hist_ljet_origin->Fill(origin);
           if(abs(jet_trk_eta[*it].at(i))<trk_eta_cut && jet_trk_pt[*it].at(i)>trk_pT_cut){// && abs(jet_trk_d0[*it].at(i))<trk_d0_cut && abs(jet_trk_z0[*it].at(i)*sin(jet_trk_theta[*it].at(i)))<trk_z0sinth_cut){//
             hist_ljet_cut_origin->Fill(origin);
             hist_ljet_cut_origin_pT->Fill(origin,1e-3*jet_trk_pt[*it].at(i));
             hist_ljet_cut_origin_jetpT->Fill(origin,1e-3*jet_pt[*it]);
             hist_ljet_cut_origin_truth_label->Fill(jet_trk_truth_label[*it].at(i));
             hist_ljet_cut_origin_truth_label_pT->Fill(jet_trk_truth_label[*it].at(i),1e-3*jet_pt[*it]);
           }
         }


       }
     }

   }



   if(shrinking_cone){


 //select by: n1==1 && n2==0 && n3==0
     if(nBcheck>0){

       for(std::vector<int>::iterator it = isBcheck.begin(); it != isBcheck.end(); ++it){
         if(jet_pt[*it]*0.001<m_pt_max_shrCone){

           std::vector<float> eta=jet_trk_eta[*it];
           std::vector<float> phi=jet_trk_phi[*it];
           if(eta.size()!=phi.size()) std::cout<< "ERROR"<< "\n";
           int size=eta.size();//size is the number of tracks, on which we make a cut
   //        if(size>=m_track_cut){
           std::vector<float> R(size);
           for(int j=0;j<size;j++){
             R.at(j)=sqrt(eta.at(j)*eta.at(j)+phi.at(j)*phi.at(j));
           }
           double R_m;
   //          R_m=std::accumulate(R.begin(), R.end(), 0.0)/R.size();
           R_m=sqrt(jet_eta[*it]*jet_eta[*it]+jet_phi[*it]*jet_phi[*it]);

           float D_phi=0;
           float DR_Max=0,tmp_M=0,sq_sum=0,std_dev=0;
           for(int j=0;j<size;j++){
   //            tmp_M=abs(R.at(j)-R_m);//<---------------------------------------WRONG!!!!!!!!!!!!!
             if(abs(jet_trk_phi[*it].at(j)-jet_phi[*it])>M_PI){
               D_phi=2*M_PI-abs(jet_trk_phi[*it].at(j)-jet_phi[*it]);
             }
             if(abs(jet_trk_phi[*it].at(j)-jet_phi[*it])<M_PI){
               D_phi=jet_trk_phi[*it].at(j)-jet_phi[*it];
             }
               //if(jet_trk_phi[*it].at(j)>0 && jet_phi[*it]<0) std::cout<<D_phi<<"\n";
             tmp_M=sqrt((eta.at(j)-jet_eta[*it])*(eta.at(j)-jet_eta[*it])+D_phi*D_phi);
   //            sq_sum+=tmp_M*tmp_M;
             if(tmp_M>DR_Max){
               DR_Max=tmp_M;
             }
           }

   //          std_dev=sqrt(sq_sum/(size-1));
   //          hist_std_dev_DR_1->Fill(std_dev);
           hist_n_tracks->Fill(size);

   //        g->SetPoint(g->GetN(), jet_pt[*it]*0.001, DR_Max);
           jet_DR_pT->Fill(jet_pt[*it]*0.001, DR_Max);
           hist_tracks_DR->Fill(DR_Max, size);
           hist_DR_1->Fill(DR_Max);

           int quot=(int) jet_pt[*it]*0.001/m_Delta_pt_shrCone;
           bin_v.at(quot).push_back(DR_Max);

   //        }
         }
       }
     }
   }

   if(derived_origin){

     if(nBcheck>0){

       for(std::vector<int>::iterator it = isBcheck.begin(); it != isBcheck.end(); ++it){
         int y=0;
         unsigned size_jet=jet_trk_pt[*it].size();

         for(unsigned i=0;i<list.size();i++){
           y=0;
           for(std::vector<int>::iterator list_it = list.at(i).begin(); list_it != list.at(i).end(); ++list_it){
             for(unsigned k=0;k<size_jet;k++){
               if((int) jet_trk_truth_label[*it].at(k)==(int) *list_it){
                 y+=1;
                 M[i]=M[i]+1;
                 continue;
               }
             }
           }
         }
       }
     }

   }//END DERIVED ORIGIN





   if(selection_alg){

     unsigned size_jet=0,size_child=0;
     if(nCcheck>0){
       for(std::vector<int>::iterator it = isCcheck.begin(); it != isCcheck.end(); ++it){//isBcheck isB
          size_child=jet_cH_child_px[*it].size();
          std::vector<double> child_Pt,child_Eta,child_Phi;

          for(unsigned j=0;j<size_child;j++){
            TLorentzVector v(
              jet_cH_child_px[*it].at(j),
              jet_cH_child_py[*it].at(j),
              jet_cH_child_pz[*it].at(j),
              jet_cH_child_E[*it].at(j)
            );
            child_Pt.push_back(v.Pt());
            child_Eta.push_back(v.Eta());
            child_Phi.push_back(v.Phi());

            if(abs(child_Eta.at(j))<trk_eta_cut && child_Pt.at(j)>trk_pT_cut){// && abs(jet_cH_child_d0[*it].at(j))<trk_d0_cut && abs((jet_cH_child_z0[*it].at(j)-(*PVz))*sin(jet_cH_child_theta[*it].at(j)))<trk_z0sinth_cut){//CHILD SELECTION CRITERIA
              hist_child_jetpT_C->Fill(1e-3*jet_pt[*it]);
              if(jet_cH_pt[*it].size()==1){

                D_eta=jet_cH_eta[*it].at(0)-jet_eta[*it];
                if(abs(jet_cH_phi[*it].at(0)-jet_phi[*it])>M_PI){
                  D_phi=2*M_PI-abs(jet_cH_phi[*it].at(0)-jet_phi[*it]);
                }
                if(abs(jet_cH_phi[*it].at(0)-jet_phi[*it])<M_PI){
                  D_phi=jet_cH_phi[*it].at(0)-jet_phi[*it];
                }
                DR_chadrjet=sqrt(D_eta*D_eta+D_phi*D_phi);
                hist_child_DR_jetpt_C->Fill(1e-3*jet_pt[*it],DR_chadrjet);
                hist_child_DR_cHpt_C->Fill(1e-3*jet_cH_pt[*it].at(0),DR_chadrjet);
                hist_child_cHpT_C->Fill(1e-3*jet_cH_pt[*it].at(0));
              }
            }
          }
        }


      }

     size_jet=0,size_child=0;
     if(nBcheck>0){
       for(std::vector<int>::iterator it = isBcheck.begin(); it != isBcheck.end(); ++it){//isBcheck isB
         std::vector<int> matching_trk,matching_child;
  //         int b=(jet_bH_pdgId[*it][i]/100)%10;
         size_jet=jet_trk_pt[*it].size();

         if(size_jet!=jet_trk_eta[*it].size() || size_jet!=jet_trk_phi[*it].size() || size_jet!=jet_trk_pdg_id[*it].size()){
           std::cout<<"WARNING\n";
         }
         size_child=jet_bH_child_px[*it].size();
         if(size_child!=jet_bH_child_py[*it].size() || size_child!=jet_bH_child_pz[*it].size() || size_child!=jet_bH_child_E[*it].size() || size_child!=jet_bH_child_pdg_id[*it].size()){
           std::cout<<"WARNING\n";
         }

  //       if(size_jet!=jet_trk_orig[*it].size()) std::cout<<"W\n";

         if(size_jet!=0 && size_child!=0){

           if(size_jet<=size_child) max_size=size_child;
           else max_size=size_jet;

           hist_n_trk->Fill(size_jet);
           hist_n_child->Fill(size_child);

           jet.SetPtEtaPhiE(jet_pt[*it],jet_eta[*it],jet_phi[*it],jet_E[*it]);
           double jet_v[]={jet.Px(),jet.Py(),jet.Pz()};


           std::vector<double> child_Pt,child_Eta,child_Phi;
           std::vector<int> child_idx;

           double child_IP[size_child],child_Lxy[size_child],child_Lxyz[size_child],child_z0sinth_signed_v[size_child],child_d0_signed_v[size_child];

           for(unsigned j=0;j<size_child;j++){
             TLorentzVector v(
               jet_bH_child_px[*it].at(j),
               jet_bH_child_py[*it].at(j),
               jet_bH_child_pz[*it].at(j),
               jet_bH_child_E[*it].at(j)
             );
             child_Pt.push_back(v.Pt());
             child_Eta.push_back(v.Eta());
             child_Phi.push_back(v.Phi());

             if(abs(child_Eta.at(j))<trk_eta_cut && child_Pt.at(j)>trk_pT_cut){// && abs(jet_bH_child_d0[*it].at(j))<trk_d0_cut && abs((jet_bH_child_z0[*it].at(j)-(*PVz))*sin(jet_bH_child_theta[*it].at(j)))<trk_z0sinth_cut){//CHILD SELECTION CRITERIA
               if(jet_bH_child_prod_x[*it].size()!=size_child || jet_bH_child_decay_x[*it].size()!=size_child){
                 std::cout<<"WARNING\n";
               }

               child_idx.push_back(j);

               px=jet_bH_child_px[*it].at(j);
               py=jet_bH_child_py[*it].at(j);
               pz=jet_bH_child_pz[*it].at(j);

               if(jet_bH_child_decay_x[*it].at(j)==-999){//NO DECAY CHILD
                 if(jet_bH_child_charge[*it].at(j)>1.){
                   c=-1;
                 }
                 if(jet_bH_child_charge[*it].at(j)<1.){
                   c=+1;
                 }

                 Lxy=800;
                 Lxyz=800;

                 hist_child_Lxy_B->Fill(Lxy);
                 hist_child_Lxyz_B->Fill(Lxyz);
                 child_Lxy[j]=Lxy;
                 child_Lxyz[j]=Lxyz;

                 if(abs(jet_bH_child_pdg_id[*it].at(j))==211){
                   hist_child_pi_notD->Fill(1e-3*child_Pt[j]);
                 }
                 if(abs(jet_bH_child_pdg_id[*it].at(j))==321){
                   hist_child_K_notD->Fill(1e-3*child_Pt[j]);
                 }
/*
                 TRandom3 *rand = new TRandom3(0);
                 rand_n=rand->Gaus(0.59,0.05943);
                 R0=abs(child_Pt.at(j)/rand_n);//mm
                 x0=jet_bH_child_prod_x[*it].at(j)-c*R0*(py/(sqrt(px*px+py*py)));
                 y0=jet_bH_child_prod_y[*it].at(j)+c*R0*(px/(sqrt(px*px+py*py)));
                 Dx_3=x0-(*PVx);
                 Dy_3=y0-(*PVy);
                 Dxy_3=sqrt((Dx_3*Dx_3)+(Dy_3*Dy_3));

                 gamma=1.-abs(R0)/(Dxy_3);
                 A=(Dx_3*gamma*py-Dy_3*gamma*px)*(jet_v[0]*py-jet_v[1]*px);
                 sgn=A/abs(A);

                 d0=sgn*abs(Dxy_3-abs(R0));
*/
//                 hist_child_nodecay_IP->Fill(d0);


               }//NO DECAY CHILD

               if(jet_bH_child_decay_x[*it].at(j)!=-999 && jet_bH_child_prod_x[*it].at(j)!=-99){//DECAY CHILD
//                 std::cout<<jet_bH_child_decay_x[*it].at(j)<<"\n";
                 Dx_2=jet_bH_child_decay_x[*it].at(j)-jet_bH_child_prod_x[*it].at(j);//x_DV - x_SV
                 Dy_2=jet_bH_child_decay_y[*it].at(j)-jet_bH_child_prod_y[*it].at(j);
                 Dz_2=jet_bH_child_decay_z[*it].at(j)-jet_bH_child_prod_z[*it].at(j);

                 Lxy=sqrt( jet_bH_child_E[*it].at(j)*jet_bH_child_E[*it].at(j)/(px*px+py*py)-(px*px+py*py+pz*pz)/(px*px+py*py) )*sqrt( Dx_2*Dx_2+Dy_2*Dy_2 );
                 Lxyz=sqrt( jet_bH_child_E[*it].at(j)*jet_bH_child_E[*it].at(j)/(px*px+py*py+pz*pz)-1. )*sqrt( Dx_2*Dx_2+Dy_2*Dy_2+Dz_2*Dz_2 );

                 hist_child_Lxy_B->Fill(Lxy);
                 hist_child_Lxyz_B->Fill(Lxyz);
                 child_Lxy[j]=Lxy;
                 child_Lxyz[j]=Lxyz;
//                 Dxy_2=sqrt((Dx_1*Dx_1)+(Dy_1*Dy_1));

//                 vx=1e3*c*px/jet_bH_child_E[*it].at(j);
//                 vy=1e3*c*py/jet_bH_child_E[*it].at(j);
/*
                 R0=0.5*sqrt(px*px+py*py)*(Dx_2*Dx_2+Dy_2*Dy_2)/(Dy_2*px-Dx_2*py);
                 x0=jet_bH_child_prod_x[*it].at(j)-R0*py/sqrt(px*px+py*py);
                 y0=jet_bH_child_prod_y[*it].at(j)+R0*px/sqrt(px*px+py*py);
                 Dx_3=x0-(*PVx);
                 Dy_3=y0-(*PVy);
                 Dxy_3=sqrt((Dx_3*Dx_3)+(Dy_3*Dy_3));


                 gamma=1.-abs(R0)/(Dxy_3);
                 A=(Dx_3*gamma*py-Dy_3*gamma*px)*(jet_v[0]*py-jet_v[1]*px);
                 sgn=A/abs(A);

                 d0=sgn*abs(Dxy_3-abs(R0));
*/
//                 hist_child_decay_IP->Fill(d0);

                 hist_pT_vs_R0_ratio_B->Fill(child_Pt.at(j)/abs(R0));

//                 ncorr+=jet_bH_child_charge[*it].at(j)*d0/abs(d0);
//                 std::cout<<jet_bH_child_charge[*it].at(j)*R0/abs(jet_bH_child_charge[*it].at(j)*R0)<<"\n";
//                 std::cout<<jet_bH_child_charge[*it].at(j)<<"\t"<<d0/abs(d0)<<"\n";
               }//DECAY CHILD

               child_IP[j]=jet_bH_child_d0[*it].at(j);
               hist_child_d0->Fill(child_IP[j]);
               hist_child_d0_pT->Fill(child_IP[j],1e-3*sqrt(px*px+py*py));

               A=sin(jet_phi[*it]-child_Phi.at(j))*jet_bH_child_d0[*it].at(j);
               sgn_d0=A/abs(A);
               child_d0_signed=sgn_d0*abs(jet_bH_child_d0[*it].at(j));
               child_z0=jet_bH_child_z0[*it].at(j)-(*PVz);
               A=(jet_eta[*it]-child_Eta.at(j))*child_z0;
               sgn_z0=A/abs(A);
               child_z0_signed=sgn_z0*abs(child_z0);
//               double sintheta=sqrt((jet_bH_child_px[*it].at(j)*jet_bH_child_px[*it].at(j)+jet_bH_child_py[*it].at(j)*jet_bH_child_py[*it].at(j))/(jet_bH_child_px[*it].at(j)*jet_bH_child_px[*it].at(j)+jet_bH_child_py[*it].at(j)*jet_bH_child_py[*it].at(j)+jet_bH_child_pz[*it].at(j)*jet_bH_child_pz[*it].at(j)));
//               if((sintheta-sin(jet_bH_child_theta[*it].at(j)))/sin(jet_bH_child_theta[*it].at(j))>1e-2)
//                std::cout<<"child sintheta error: "<<(sintheta-sin(jet_bH_child_theta[*it].at(j)))/sin(jet_bH_child_theta[*it].at(j))<<"\n";
               child_z0sinth_signed=child_z0_signed*sin(jet_bH_child_theta[*it].at(j));
//               child_z0sinth_signed=child_z0_signed*sintheta;
               child_d0_signed_v[j]=child_d0_signed;
               child_z0sinth_signed_v[j]=child_z0sinth_signed;

               hist_child_z0->Fill(child_z0_signed);
               hist_child_d0_truth->Fill(child_d0_signed);
               if(jet_bH_child_theta[*it].at(j)==-99)
               std::cout<<"W\n";
               hist_child_z0sinth_B->Fill(child_z0sinth_signed);

               D_eta=child_Eta[j]-jet_eta[*it];
               if(abs(child_Phi[j]-jet_phi[*it])>M_PI){
                 D_phi=2*M_PI-abs(child_Phi[j]-jet_phi[*it]);
               }
               if(abs(child_Phi[j]-jet_phi[*it])<M_PI){
                 D_phi=child_Phi[j]-jet_phi[*it];
               }
               DR=sqrt(D_eta*D_eta+D_phi*D_phi);
               hist_child_pT_B->Fill(1e-3*child_Pt[j]);
               hist_child_eta_B->Fill(child_Eta[j]);
               hist_child_phi_B->Fill(child_Phi[j]);
               hist_child_Deta_B->Fill(D_eta);
               hist_child_Dphi_B->Fill(D_phi);
               hist_child_Dphi_Deta_B->Fill(D_phi,D_eta);
               hist_child_DR_B->Fill(DR);
               hist_child_pT_DR_B->Fill(1e-3*child_Pt[j],DR);
               hist_child_pT_jet_DR_B->Fill(1e-3*jet_pt[*it],DR);
               hist_child_pdgID_B->Fill(jet_bH_child_pdg_id[*it].at(j));

               hist_child_jetpT_B->Fill(1e-3*jet_pt[*it]);
               if(jet_bH_pt[*it].size()==1){

                 D_eta=jet_bH_eta[*it].at(0)-jet_eta[*it];
                 if(abs(jet_bH_phi[*it].at(0)-jet_phi[*it])>M_PI){
                   D_phi=2*M_PI-abs(jet_bH_phi[*it].at(0)-jet_phi[*it]);
                 }
                 if(abs(jet_bH_phi[*it].at(0)-jet_phi[*it])<M_PI){
                   D_phi=jet_bH_phi[*it].at(0)-jet_phi[*it];
                 }
                 DR_bhadrjet=sqrt(D_eta*D_eta+D_phi*D_phi);
                 hist_child_DR_jetpt_B->Fill(1e-3*jet_pt[*it],DR_bhadrjet);
                 hist_child_DR_bHpt_B->Fill(1e-3*jet_bH_pt[*it].at(0),DR_bhadrjet);
                 hist_child_bHpT_B->Fill(1e-3*jet_bH_pt[*it].at(0));
               }


               if(abs(jet_bH_child_pdg_id[*it].at(j))==211){
                 hist_child_pi->Fill(1e-3*child_Pt[j]);
                 hist_child_pi_Lxy_B->Fill(Lxy);
                 hist_child_pi_Lxyz_B->Fill(Lxyz);//CHECK
//                 hist_child_pi_d0_truth_B->Fill(child_d0_signed);
//                 hist_child_pi_z0_truth_B->Fill(jet_bH_child_z0[*it].at(j));
               }
               if(abs(jet_bH_child_pdg_id[*it].at(j))==321){
                 hist_child_K->Fill(1e-3*child_Pt[j]);
                 hist_child_K_Lxy_B->Fill(Lxy);
                 hist_child_K_Lxyz_B->Fill(Lxyz);
//                 hist_child_K_d0_truth_B->Fill(child_d0_signed);
//                 hist_child_K_z0_truth_B->Fill(jet_bH_child_z0[*it].at(j));
               }
               if(abs(jet_bH_child_pdg_id[*it].at(j))==13){
                 hist_child_mu->Fill(1e-3*child_Pt[j]);
                 hist_child_mu_Lxy_B->Fill(Lxy);
                 hist_child_mu_Lxyz_B->Fill(Lxyz);
//                 hist_child_mu_d0_truth_B->Fill(child_d0_signed);
//                 hist_child_mu_z0_truth_B->Fill(jet_bH_child_z0[*it].at(j));
               }
               if(abs(jet_bH_child_pdg_id[*it].at(j))==2212){
                 hist_child_p->Fill(1e-3*child_Pt[j]);
                 hist_child_p_Lxy_B->Fill(Lxy);
                 hist_child_p_Lxyz_B->Fill(Lxyz);
//                 hist_child_p_d0_truth_B->Fill(child_d0_signed);
//                 hist_child_p_z0_truth_B->Fill(jet_bH_child_z0[*it].at(j));
               }
               if(abs(jet_bH_child_pdg_id[*it].at(j))==11){
                 hist_child_e->Fill(1e-3*child_Pt[j]);
                 hist_child_e_Lxy_B->Fill(Lxy);
                 hist_child_e_Lxyz_B->Fill(Lxyz);
//                 hist_child_e_d0_truth_B->Fill(child_d0_signed);
//                 hist_child_e_z0_truth_B->Fill(jet_bH_child_z0[*it].at(j));
               }
             }
           }



           if(geometric_selection){
   	         int den=0;

             if(child_Pt.size()!=size_child || child_Eta.size()!=size_child || child_Phi.size()!=size_child){
               std::cout<<"WARNING\n";
             }

             int bool_matrix[size_jet][size_child];
             for(unsigned l=0;l<size_jet;l++){
               for(unsigned k=0;k<size_child;k++){
                  bool_matrix[l][k]=0;
               }
             }

             for(unsigned j=0;j<size_child;j++){
               if(abs(child_Eta.at(j))<trk_eta_cut && child_Pt.at(j)>trk_pT_cut){// && abs(jet_bH_child_d0[*it].at(j))<trk_d0_cut && abs((jet_bH_child_z0[*it].at(j)-(*PVz))*sin(jet_bH_child_theta[*it].at(j)))<trk_z0sinth_cut){//
                 den++;

                 if(jet_bH_child_charge[*it].at(j)==0){
                   std::cout<<"CHARGE 0\n";
                 }
                 for(unsigned i=0;i<size_jet;i++){
                   if(abs(jet_trk_eta[*it].at(i))<trk_eta_cut && jet_trk_pt[*it].at(i)>trk_pT_cut){// && abs(jet_trk_d0[*it].at(i))<trk_d0_cut && abs(jet_trk_z0[*it].at(i)*sin(jet_trk_theta[*it].at(i)))<trk_z0sinth_cut){//
                     if( (jet_trk_parent_pdgid[*it].at(i)-jet_bH_child_parent_pdg_id[*it].at(j))==0 && abs(jet_trk_parent_pdgid[*it].at(i))!=999 ){
                       if( (jet_trk_pdg_id[*it].at(i)-jet_bH_child_pdg_id[*it].at(j))==0 && abs(jet_trk_pdg_id[*it].at(i))!=999 ){
  //                       std::cout<<"("<<i<<","<<j<<") ";
                         bool_matrix[i][j] = (int) 1;
                       }
                     }
                   }
                 }
               }
             }



/*
    std::cout<<"\nCORRESPONDENCE MATRIX\t"<<size_jet<<"x"<<size_child<<"\tevent:"<<m_Ntot<<"\n";
             std::cout<<"\t";
             for(int j=0;j<size_child;j++)
                std::cout<<jet_bH_child_pdg_id[*it].at(j)<<"|"<<jet_bH_child_parent_pdg_id[*it].at(j)<<" ";
             for(int i=0;i<size_jet;i++){
               std::cout<<"\n"<<jet_trk_pdg_id[*it].at(i)<<"|"<<jet_trk_parent_pdgid[*it].at(i)<<"   \t";
               for(int j=0;j<size_child;j++){
                 std::cout << bool_matrix[i][j]<<"\t";
               }
             }

             std::cout<<"\n";
*/
             for(q=0;q<max_size;q++){
               m_qc=(int) q%size_child;
               m_qj=(int) q%size_jet;
               sc=0;
/*
                 if(m_qc>=size_child || m_qj>=size_jet){
                   std::cout<<m_Ntot<<"\t"<<q<<"\t"<<size_jet<<","<<m_qj<<"\t"<<size_child<<","<<m_qc<<"\n";
                 }
*/

               if(cut){
                 tmp_min_pTfraction=m_pTfraction_cut;//CUT
                 tmp_min_DR=m_DRcut;//CUT
               }
               if(!cut){
                 tmp_min_pTfraction=m_pTfraction_nocut;//NOCUT
                 tmp_min_DR=m_DRnocut;//NOCUT
               }

               a=-2;b=-2;

               for(unsigned j=0;j<size_child;j++){
                 if(bool_matrix[m_qj][j]==1){
    //                 std::cout<<"row\n";
                   sc=sc+1;

                   D_eta=child_Eta.at(j)-jet_trk_eta[*it].at(m_qj);
                   if(abs(child_Phi.at(j)-jet_trk_phi[*it].at(m_qj))>M_PI){
                     D_phi=2*M_PI-abs(child_Phi.at(j)-jet_trk_phi[*it].at(m_qj));
                   }
                   if(abs(child_Phi.at(j)-jet_trk_phi[*it].at(m_qj))<M_PI){
                     D_phi=child_Phi.at(j)-jet_trk_phi[*it].at(m_qj);
                   }
                   tmp_DR=sqrt(D_eta*D_eta+D_phi*D_phi);
                   tmp_pTfraction=1e-3*abs(jet_trk_pt[*it].at(m_qj)-child_Pt.at(j))/(1e-3*child_Pt.at(j));

                   if((tmp_pTfraction<tmp_min_pTfraction) && (tmp_DR<tmp_min_DR)){//NOCUT
                     tmp_min_pTfraction=tmp_pTfraction;
                     tmp_min_DR=tmp_DR;
                     a=m_qj;
                     b=j;
                   }
                 }
               }
               for(unsigned i=0;i<size_jet;i++){
                 if(bool_matrix[i][m_qc]==1){
    //                 std::cout<<"col\n";
                   sc=sc+1;

                   D_eta=child_Eta.at(m_qc)-jet_trk_eta[*it].at(i);
                   if(abs(child_Phi.at(m_qc)-jet_trk_phi[*it].at(i))>M_PI){
                     D_phi=2*M_PI-abs(child_Phi.at(m_qc)-jet_trk_phi[*it].at(i));
                   }
                   if(abs(child_Phi.at(m_qc)-jet_trk_phi[*it].at(i))<M_PI){
                     D_phi=child_Phi.at(m_qc)-jet_trk_phi[*it].at(i);
                   }
                   tmp_DR=sqrt(D_eta*D_eta+D_phi*D_phi);
                   tmp_pTfraction=1e-3*abs(jet_trk_pt[*it].at(i)-child_Pt.at(m_qc))/(1e-3*child_Pt.at(m_qc));

                   if((tmp_pTfraction<tmp_min_pTfraction) && (tmp_DR<tmp_min_DR)){//CUT    tmp_min_DR=0.1
                     tmp_min_pTfraction=tmp_pTfraction;
                     tmp_min_DR=tmp_DR;
                     a=i;
                     b=m_qc;
                   }
                 }
               }

    /*
                 if(a==-2 && b==-2 && sc>0){
                   std::cout<<"CLUSED CHILDS:\t"<<m_Ntot<<"\t"<<tmp_DpT<<"\t"<<tmp_DR<<"\n";
    //               den--;
                 }
    */

               if(a!=-2 && b!=-2){
                 match++;
                 matching_trk.push_back(a);
                 matching_child.push_back(b);
//                 std::cout<<"matchings: "<<a<<","<<b<<"\t";

                 if(jet_trk_orig[*it].at(a)==0 || jet_trk_orig[*it].at(a)==1){
                   m_match_overlap++;
                 }
                 if(jet_trk_orig[*it].at(a)!=0 && jet_trk_orig[*it].at(a)!=1){
                   m_match_notoverlap++;
                 }

                 for(unsigned i=0;i<child_idx.size();i++){
                   if(child_idx.at(i)==b){
                     child_idx.erase(child_idx.begin()+i);
                   }
                 }

                 px=jet_bH_child_px[*it].at(b);
                 py=jet_bH_child_py[*it].at(b);
                 pz=jet_bH_child_pz[*it].at(b);

                 if(jet_bH_child_decay_x[*it].at(b)!=-999 && jet_bH_child_prod_x[*it].at(b)!=-99){//DECAY CHILD
  //                 std::cout<<jet_bH_child_decay_x[*it].at(b)<<"\n";
                   Dx_2=jet_bH_child_decay_x[*it].at(b)-jet_bH_child_prod_x[*it].at(b);//x_DV - x_SV
                   Dy_2=jet_bH_child_decay_y[*it].at(b)-jet_bH_child_prod_y[*it].at(b);
                   Dz_2=jet_bH_child_decay_z[*it].at(b)-jet_bH_child_prod_z[*it].at(b);
                   Lxy=sqrt( jet_bH_child_E[*it].at(b)*jet_bH_child_E[*it].at(b)/(px*px+py*py)-(px*px+py*py+pz*pz)/(px*px+py*py) )*sqrt( Dx_2*Dx_2+Dy_2*Dy_2 );
                   Lxyz=sqrt( jet_bH_child_E[*it].at(b)*jet_bH_child_E[*it].at(b)/(px*px+py*py+pz*pz)-1. )*sqrt( Dx_2*Dx_2+Dy_2*Dy_2+Dz_2*Dz_2 );
                 }
                 if(jet_bH_child_decay_x[*it].at(b)==-999){//NO DECAY CHILD
                   Lxy=800;
                   Lxyz=800;
                 }

                   //CINEMATICA RISPETTO AL JET
                 D_eta=child_Eta.at(b)-jet_eta[*it];
                 if(abs(child_Phi.at(b)-jet_phi[*it])>M_PI){
                   D_phi=2*M_PI-abs(child_Phi.at(b)-jet_phi[*it]);
                 }
                 if(abs(child_Phi.at(b)-jet_phi[*it])<M_PI){
                   D_phi=child_Phi.at(b)-jet_phi[*it];
                 }
                 DR=sqrt(D_eta*D_eta+D_phi*D_phi);
                 hist_matched_pT_B->Fill(1e-3*child_Pt.at(b));
                 hist_matched_eta_B->Fill(child_Eta.at(b));
                 hist_matched_phi_B->Fill(child_Phi.at(b));
                 hist_matched_Deta_B->Fill(D_eta);
                 hist_matched_Dphi_B->Fill(D_phi);
                 hist_matched_Dphi_Deta_B->Fill(D_phi,D_eta);
                 hist_matched_DR_B->Fill(DR);
                 hist_matched_pT_DR_B->Fill(1e-3*child_Pt.at(b),DR);
                 hist_matched_pT_jet_DR_B->Fill(1e-3*jet_pt[*it],DR);
                 hist_matched_pdgId_B->Fill(jet_bH_child_pdg_id[*it].at(b));

                 if(jet_bH_child_decay_x[*it].at(b)==-999){//NO DECAY MATCHED CHILD
                   if(abs(jet_bH_child_pdg_id[*it].at(b))==211){
                     hist_matched_child_pi_notD->Fill(1e-3*child_Pt[b]);
                   }
                   if(abs(jet_bH_child_pdg_id[*it].at(b))==321){
                     hist_matched_child_K_notD->Fill(1e-3*child_Pt[b]);
                   }
                 }

                 if(abs(jet_bH_child_pdg_id[*it].at(b))==211){
                   hist_matched_child_pi->Fill(1e-3*child_Pt.at(b));
                   hist_matched_child_pi_Lxy_B->Fill(Lxy);
                   hist_matched_child_pi_Lxyz_B->Fill(Lxyz);
//                   hist_matched_child_pi_d0_truth_B->Fill(sgn*abs(jet_bH_child_d0[*it].at(b)));
//                   hist_matched_child_pi_z0_truth_B->Fill(jet_bH_child_z0[*it].at(b));
                 }
                 if(abs(jet_bH_child_pdg_id[*it].at(b))==321){
                   hist_matched_child_K->Fill(1e-3*child_Pt.at(b));
                   hist_matched_child_K_Lxy_B->Fill(Lxy);
                   hist_matched_child_K_Lxyz_B->Fill(Lxyz);
//                   hist_matched_child_K_d0_truth_B->Fill(sgn*abs(jet_bH_child_d0[*it].at(b)));
//                   hist_matched_child_K_z0_truth_B->Fill(jet_bH_child_z0[*it].at(b));
                 }
                 if(abs(jet_bH_child_pdg_id[*it].at(b))==13){
                   hist_matched_child_mu->Fill(1e-3*child_Pt.at(b));
                   hist_matched_child_mu_Lxy_B->Fill(Lxy);
                   hist_matched_child_mu_Lxyz_B->Fill(Lxyz);
//                   hist_matched_child_mu_d0_truth_B->Fill(sgn*abs(jet_bH_child_d0[*it].at(b)));
//                   hist_matched_child_mu_z0_truth_B->Fill(jet_bH_child_z0[*it].at(b));
                 }
                 if(abs(jet_bH_child_pdg_id[*it].at(b))==2212){
                   hist_matched_child_p->Fill(1e-3*child_Pt.at(b));
                   hist_matched_child_p_Lxy_B->Fill(Lxy);
                   hist_matched_child_p_Lxyz_B->Fill(Lxyz);
//                   hist_matched_child_p_d0_truth_B->Fill(sgn*abs(jet_bH_child_d0[*it].at(b)));
//                   hist_matched_child_p_z0_truth_B->Fill(jet_bH_child_z0[*it].at(b));
                 }
                 if(abs(jet_bH_child_pdg_id[*it].at(b))==11){
                   hist_matched_child_e->Fill(1e-3*child_Pt.at(b));
                   hist_matched_child_e_Lxy_B->Fill(Lxy);
                   hist_matched_child_e_Lxyz_B->Fill(Lxyz);
//                   hist_matched_child_e_d0_truth_B->Fill(sgn*abs(jet_bH_child_d0[*it].at(b)));
//                   hist_matched_child_e_z0_truth_B->Fill(jet_bH_child_z0[*it].at(b));
                 }

                 hist_matched_origin_B->Fill(jet_trk_orig[*it].at(a));
                 hist_matched_d0_B->Fill(child_d0_signed);
                 hist_matched_Lxy_B->Fill(child_Lxy[b]);
                 hist_matched_Lxyz_B->Fill(child_Lxyz[b]);

                 //CINEMATICA RISPETTO ALLA TRACCIA
                 DpT_trk=1e-3*(child_Pt.at(b)-jet_trk_pt[*it].at(a));
                 D_eta_trk=child_Eta.at(b)-jet_trk_eta[*it].at(a);
                 if(abs(child_Phi.at(b)-jet_trk_phi[*it].at(a))>M_PI){
                   D_phi_trk=2*M_PI-abs(child_Phi.at(b)-jet_trk_phi[*it].at(a));
                 }
                 if(abs(child_Phi.at(b)-jet_trk_phi[*it].at(a))<M_PI){
                   D_phi_trk=child_Phi.at(b)-jet_trk_phi[*it].at(a);
                 }
                 DR_trk=sqrt(D_eta_trk*D_eta_trk+D_phi_trk*D_phi_trk);
                 hist_matched_pT_child_pTfraction_B->Fill(1e-3*child_Pt.at(b),-1.*DpT_trk/(1e-3*child_Pt.at(b)));
                 hist_matched_DR_trk_B->Fill(DR_trk);
                 hist_matched_DR_trk_pTfraction->Fill(DR_trk,-1.*DpT_trk/(1e-3*child_Pt.at(b)));

                 if(sc==1){
                   m_sc+=1;
                   /*
                   hist_single_matched_pT_B->Fill(1e-3*child_Pt.at(b));
                   hist_single_matched_eta_B->Fill(child_Eta.at(b));
                   hist_single_matched_phi_B->Fill(child_Phi.at(b));
                   hist_single_matched_Deta_B->Fill(D_eta);
                   hist_single_matched_Dphi_B->Fill(D_phi);
                   hist_single_matched_Dphi_Deta_B->Fill(D_phi,D_eta);
                   hist_single_matched_DR_B->Fill(DR);//RISPETTO AL JET
                   hist_single_matched_pT_DR_B->Fill(1e-3*child_Pt.at(b),DR);
                   hist_single_matched_pT_jet_DR_B->Fill(1e-3*jet_pt[*it],DR);
                   hist_single_matched_pdgId_B->Fill(jet_bH_child_pdg_id[*it].at(b));
                   hist_single_matched_origin_B->Fill(jet_trk_orig[*it].at(a));

                   hist_single_matched_pT_child_pTfraction_B->Fill(1e-3*child_Pt.at(b),-1.*DpT_trk/(1e-3*child_Pt.at(b)));
                   hist_single_matched_DR_trk_B->Fill(DR_trk);
                   hist_single_matched_DR_trk_pTfraction->Fill(DR_trk,-1.*DpT_trk/(1e-3*child_Pt.at(b)));

                   hist_single_matched_d0_B->Fill(child_IP[b]);
                   */
                 }
                 if(sc==2){
                   m_sc2++;
                 }
                 if(sc>=3){
                   m_sc3++;
                 }
               }

//  std::cout<<"\n"<<"q="<<q<<"\t(a,b)=("<<a+1<<","<<b+1<<")"<<"\t"<<sc<<"\n";

               if(sc>0){
                 for(unsigned j=0;j<size_child;j++){
                   bool_matrix[a][j]=0;
                 }
                 for(unsigned i=0;i<size_jet;i++){
                   bool_matrix[i][b]=0;
                 }
               }

/*
    std::cout<<"ELIMINATION STARTS:\n";
               for(int i=0;i<size_jet;i++){
                 std::cout<<"\n";
                 for(int j=0;j<size_child;j++){
                   std::cout << bool_matrix[i][j]<<"\t";
                 }
               }
             std::cout<<"\n";
*/

             }

             m_nomatch+=child_idx.size();

             int j=0;
             double DR_shrink;
             for(unsigned l=0;l<child_idx.size();l++){
               j=child_idx.at(l);
               //CINEMATICA RISPETTO AL JET
//               if(abs(child_Eta.at(j))<trk_eta_cut && child_Pt.at(j)>trk_pT_cut){
                 D_eta=child_Eta.at(j)-jet_eta[*it];
                 if(abs(child_Phi.at(j)-jet_phi[*it])>M_PI){
                   D_phi=2*M_PI-abs(child_Phi.at(j)-jet_phi[*it]);
                 }
                 if(abs(child_Phi.at(j)-jet_phi[*it])<M_PI){
                   D_phi=child_Phi.at(j)-jet_phi[*it];
                 }
                 DR=sqrt(D_eta*D_eta+D_phi*D_phi);
                 /*
                 hist_nomatched_pT_B->Fill(1e-3*child_Pt.at(j));
                 hist_nomatched_eta_B->Fill(child_Eta.at(j));
                 hist_nomatched_phi_B->Fill(child_Phi.at(j));
                 hist_nomatched_Deta_B->Fill(D_eta);
                 hist_nomatched_Dphi_B->Fill(D_phi);
                 hist_nomatched_Dphi_Deta_B->Fill(D_phi,D_eta);
                 hist_nomatched_DR_B->Fill(DR);
                 hist_nomatched_pT_DR_B->Fill(1e-3*child_Pt.at(j),DR);
                 */
                 hist_nomatched_pT_jet_DR_B->Fill(1e-3*jet_pt[*it],DR);
//                 hist_nomatched_pdgId_B->Fill(jet_bH_child_pdg_id[*it].at(j));

                 DR_shrink=shrinking_cone_DR(1e-3*jet_pt[*it],m_p1,m_p2,m_p3);
                 if(DR<=DR_shrink){
                   hist_nomatchedIN_pT_B->Fill(1e-3*child_Pt.at(j));
                   hist_nomatchedIN_eta_B->Fill(child_Eta.at(j));
                   hist_nomatchedIN_phi_B->Fill(child_Phi.at(j));
                   hist_nomatchedIN_DR_B->Fill(DR);
                   hist_nomatchedIN_pT_jet_DR_B->Fill(1e-3*jet_pt[*it],DR);
                   hist_nomatchedIN_d0_B->Fill(child_d0_signed_v[j]);
                   hist_nomatchedIN_z0sinth_B->Fill(child_z0sinth_signed_v[j]);
                 }
                 if(DR>DR_shrink){
                   hist_nomatchedOUT_pT_B->Fill(1e-3*child_Pt.at(j));
                   hist_nomatchedOUT_eta_B->Fill(child_Eta.at(j));
                   hist_nomatchedOUT_phi_B->Fill(child_Phi.at(j));
                   hist_nomatchedOUT_DR_B->Fill(DR);
                   hist_nomatchedOUT_pT_jet_DR_B->Fill(1e-3*jet_pt[*it],DR);
                   hist_nomatchedOUT_d0_B->Fill(child_d0_signed_v[j]);
                   hist_nomatchedOUT_z0sinth_B->Fill(child_z0sinth_signed_v[j]);
                 }
//               }
             }


    //std::cout<<"event:\t"<<m_Ntot<< "\tn of matched tracks:\t"<< match << "\tn of child tracks(denominator):\t" << den << "\tratio:\t" <<(float) match/den<<"\t";
             m_match+=match;
             m_den+=den;
             hist_n_match->Fill(match);
             hist_efficiency_B->Fill((float) match/den);
             match=0;

           }

         }
         std::vector<int> origin_trk;
         int origin=0;
         for(unsigned i=0;i<size_jet;i++){
           if(abs(jet_trk_eta[*it].at(i))<trk_eta_cut && jet_trk_pt[*it].at(i)>trk_pT_cut){// && abs(jet_trk_d0[*it].at(i))<trk_d0_cut && abs(jet_trk_z0[*it].at(i)*sin(jet_trk_theta[*it].at(i)))<trk_z0sinth_cut){//
             origin=jet_trk_truth_derivedlabel[*it].at(i);
             if(origin==5){
               origin_trk.push_back(i);
             }
           }
         }
         //GEOMETRIC AND ORIGIN SELECTION DONE
//         std::cout<<"\n";
         for(unsigned i=0;i<origin_trk.size();i++){//TRK - CUT
               m1_ov=0;
               for(unsigned a=0;a<matching_trk.size();a++){
                 if(origin_trk.at(i)==matching_trk.at(a)){
                   m1++;
                   m1_ov++;
                 }
               }
               if(m1_ov==0){//ORIGIN - NOT-GEOMETRY
                 m1_ex++;

//                 std::cout<<"ORIGIN - NOT-GEOMETRY: track: "<<jet_trk_pdg_id[*it].at(origin_trk.at(i))<<"|"<<jet_trk_parent_pdgid[*it].at(origin_trk.at(i))<<"   \t"<<origin_trk.at(i)<<"   \tpT: "<<jet_trk_pt[*it].at(origin_trk.at(i))<<"   eta: "<<jet_trk_eta[*it].at(origin_trk.at(i))<<"   \tEvent: "<<m_Ntot<<"\n";
/*
                 for(unsigned k=0;k<jet_bH_child_pdg_id[*it].size();k++){
                   std::cout<<jet_bH_child_pdg_id[*it].at(k)<<"\t"<<jet_bH_child_parent_pdg_id[*it].at(k)<<"\t"<<k<<"\n";
                 }
*/
               }
         }
         for(unsigned a=0;a<matching_trk.size();a++){//CHILD CUT
               m2_ov=0;
               for(unsigned i=0;i<origin_trk.size();i++){
                 if(origin_trk.at(i)==matching_trk.at(a)){
                   m2++;
                   m2_ov++;
                 }
               }
               if(m2_ov==0){//GEOMETRY - NOT-ORIGIN
                 m2_ex++;
                 if(jet_trk_orig[*it].at(matching_trk.at(a))==-1){
                   m_GeomNOr_PU++;
                 }
                 if(jet_trk_orig[*it].at(matching_trk.at(a))==2){
                   m_GeomNOr_F++;
                 }
                 if(jet_trk_orig[*it].at(matching_trk.at(a))==3){
                   m_GeomNOr_G++;
                 }


//                 std::cout<<"GEOMETRY - NOT-ORIGIN: track: "<<jet_trk_pdg_id[*it].at(matching_trk.at(a))<<"|"<<jet_trk_parent_pdgid[*it].at(matching_trk.at(a))<<"   \t"<<matching_trk.at(a)<<"   \tpT: "<<jet_trk_pt[*it].at(matching_trk.at(a))<<"   eta: "<<jet_trk_eta[*it].at(matching_trk.at(a))<<"   \tEvent: "<<m_Ntot<<"\n";
//                 matching_ex.push_back(a);
               }
         }
         if(m1!=m2){
           std::cout<<"W"<<"\n";
         }
         mm+=m1;
         mm1_ex+=m1_ex;
         mm2_ex+=m2_ex;
         m1=0;m2=0;m1_ex=0;m2_ex=0;m1_ov=0;m2_ov=0;



       }// end of for(isBcheck) (b-jets in a given event)

     }// end of if (nBcheck>0)

   }// end of if (selection_alg)


   return kTRUE;
}

void DAOD_selector::SlaveTerminate()
{
   std::cout<<"\n In DAOD_selector::SlaveTerminate"<<std::endl;
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.

   std::cout<<"Out of DAOD_selector::SlaveTerminate"<<std::endl;
}

void DAOD_selector::Terminate()
{
   std::cout<<"\n===\n==="<<std::endl;
   std::cout<<" In DAOD_selector::Terminate"<<std::endl;
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.

   if(derived_origin){
     hist_jet_trk_truth_label->Write();
   }

   if(selections){
/*
     hist_pt_1->Write();
     hist_eta_1->Write();
     hist_phi_1->Write();
     hist_E_1->Write();

     hist_pt_2->Write();
     hist_eta_2->Write();
     hist_phi_2->Write();
     hist_E_2->Write();

     hist_pt_2b->Write();
     hist_eta_2b->Write();
     hist_phi_2b->Write();
     hist_E_2b->Write();

     hist_pt_3a->Write();
     hist_eta_3a->Write();
     hist_phi_3a->Write();
     hist_E_3a->Write();

     hist_pt_3b->Write();
     hist_eta_3b->Write();
     hist_phi_3b->Write();
     hist_E_3b->Write();

     hist_pt_4->Write();
     hist_eta_4->Write();
     hist_phi_4->Write();
     hist_E_4->Write();
*/
     hist_pt_B->Write();
     hist_bHpt_B->Write();
     hist_eta_B->Write();
     hist_phi_B->Write();
     hist_E_B->Write();
     hist_Bjet_origin->Write();
     hist_Bjet_cut_origin->Write();
     hist_Bjet_cut_origin_pT->Write();
     hist_Bjet_cut_origin_jetpT->Write();
     hist_Bjet_cut_origin_truth_label->Write();
     hist_Bjet_cut_origin_truth_label_pT->Write();
     hist_Bjet_cut_origin_truth_label_bHpT->Write();

     hist_pt_C->Write();
     hist_cHpt_C->Write();
     hist_eta_C->Write();
     hist_phi_C->Write();
     hist_E_C->Write();
     hist_Cjet_origin->Write();
     hist_Cjet_cut_origin->Write();
     hist_Cjet_cut_origin_pT->Write();
     hist_Cjet_cut_origin_jetpT->Write();
     hist_Cjet_cut_origin_truth_label->Write();
     hist_Cjet_cut_origin_truth_label_pT->Write();
     hist_Cjet_cut_origin_truth_label_cHpT->Write();

     hist_pt_l->Write();
     hist_eta_l->Write();
     hist_phi_l->Write();
     hist_E_l->Write();
     hist_ljet_origin->Write();
     hist_ljet_cut_origin->Write();
     hist_ljet_cut_origin_pT->Write();
     hist_ljet_cut_origin_jetpT->Write();
     hist_ljet_cut_origin_truth_label->Write();
     hist_ljet_cut_origin_truth_label_pT->Write();

     hist_nBjets->Write();
     hist_n_tracks_jetpt_B->Write();
     hist_n_tracks_bHpt_B->Write();
     hist_n_BCtracks_jetpt_B->Write();
     hist_n_BCtracks_bHpt_B->Write();

     hist_nCjets->Write();
     hist_n_tracks_jetpt_C->Write();
     hist_n_tracks_cHpt_C->Write();
     hist_n_BCtracks_jetpt_C->Write();
     hist_n_BCtracks_cHpt_C->Write();

     hist_nljets->Write();
     hist_n_tracks_jetpt_l->Write();

     hist2_jetFlavorMatrix->Write();

   }

   
   if (fDoFlavorLabelMatrix)
     {

       saveHistosInMaps();
       /*
       //getPtFromFlavorLabelMatrix();
       for(unsigned int i1=0; i1<6; i1++){
         for(unsigned int i2=0; i2<6; i2++){
           HistopT[i1][i2]->Write();
         }
       }
       */       
     }


   if(discriminants){

//     hist_jet_IP2_B->Write();
/*
     hist_ip2d_pb->Write();
     hist_ip2d_pc->Write();
     hist_ip2d_pu->Write();
     hist_ip2d_llr->Write();
     hist_ip3d_pb->Write();
     hist_ip3d_pc->Write();
     hist_ip3d_pu->Write();
     hist_ip3d_llr->Write();
     hist_dl1_pb->Write();
     hist_dl1_pc->Write();
     hist_dl1_pu->Write();
*/

     hist_ip2d_llr_B->Write();
     hist_ip2d_llr_jetpt_B->Write();
     hist_ip2d_llr_jetpt_singleB->Write();
     hist_ip2d_llr_C->Write();
     hist_ip2d_llr_jetpt_C->Write();
     hist_ip2d_llr_jetpt_singleC->Write();
     hist_ip2d_llr_l->Write();

     hist_ip3d_llr_B->Write();
     hist_ip3d_llr_jetpt_singleB->Write();
     hist_ip3d_llr_jetpt_B->Write();
     hist_ip3d_llr_C->Write();
     hist_ip3d_llr_jetpt_singleC->Write();
     hist_ip3d_llr_jetpt_C->Write();
     hist_ip3d_llr_l->Write();

     hist_rnnip_llr_B->Write();
     hist_rnnip_llr_jetpt_singleB->Write();
     hist_rnnip_llr_jetpt_B->Write();
     hist_rnnip_llr_C->Write();
     hist_rnnip_llr_jetpt_singleC->Write();
     hist_rnnip_llr_jetpt_C->Write();
     hist_rnnip_llr_l->Write();

     hist_sv1_llr_B->Write();
     hist_sv1_llr_jetpt_B->Write();
     hist_sv1_llr_jetpt_singleB->Write();
     hist_sv1_llr_C->Write();
     hist_sv1_llr_jetpt_C->Write();
     hist_sv1_llr_jetpt_singleC->Write();
     hist_sv1_llr_l->Write();
     hist_jet_sv1_Nvtx_B->Write();
     hist_jet_sv1_ntrkv_B->Write();
     hist_jet_sv1_n2t_B->Write();
     hist_jet_sv1_m_B->Write();
     hist_jet_sv1_efc_B->Write();
     hist_jet_sv1_sig3d_B->Write();
     hist_jet_sv1_deltaR_B->Write();
     hist_jet_sv1_Lxy_B->Write();
     hist_jet_sv1_L3d_B->Write();

     hist_jf_llr_B->Write();
     hist_jf_llr_jetpt_B->Write();
     hist_jf_llr_jetpt_singleB->Write();
     hist_jf_llr_C->Write();
     hist_jf_llr_jetpt_C->Write();
     hist_jf_llr_jetpt_singleC->Write();
     hist_jf_llr_l->Write();
     hist_dl1_B->Write();
     hist_dl1_llr_jetpt_B->Write();
     hist_dl1_llr_jetpt_singleB->Write();
     hist_dl1_C->Write();
     hist_dl1_llr_jetpt_C->Write();
     hist_dl1_llr_jetpt_singleC->Write();
     hist_dl1_l->Write();

   }


    if(shrinking_cone){

      hist_n_tracks->Write();
      hist_tracks_DR->SetMarkerStyle(kFullCircle);
 //    hist_tracks_DR->SetMarkerSize(10);
      hist_tracks_DR->Write();
      hist_DR_1->Write();
 //    hist_std_dev_DR_1->Write();
      jet_DR_pT->SetMarkerStyle(kFullCircle);
      jet_DR_pT->Write();
    }


    if(selection_alg){

      hist_trk_pT_B->Write();
      hist_trk_eta_B->Write();
      hist_trk_phi_B->Write();
      hist_trk_Deta_B->Write();
      hist_trk_Dphi_B->Write();
      hist_trk_Dphi_Deta_B->Write();
      hist_trk_DR_B->Write();
      hist_trk_pT_DR_B->Write();
      hist_trk_pT_jet_DR_B->Write();
      hist_trk_pdgId_B->Write();
      hist_trk_origin_B->Write();
      hist_trk_origin_truth_label_B->Write();

      hist_trk_d0_B->Write();
      hist_trk_d0_jet_pT_BC_B->Write();
      hist_trk_DR_jet_pt_BC->Write();
      hist_trk_sigd0_B->Write();
      hist_trk_z0_B->Write();
      hist_trk_sigz0_B->Write();
      hist_trk_z0sinth_B->Write();
      hist_trk_d0sig_B->Write();
      hist_trk_z0sinthsig_B->Write();
      hist_trk_d0sig_origin_B->Write();
      hist_trk_d0sig_truthlabel_B->Write();
      hist_trk_z0sinthsig_origin_B->Write();
      hist_trk_logpTfrac_origin_B->Write();
      hist_trk_logDR_origin_B->Write();
      hist_trk_IBLhits_origin_B->Write();
      hist_trk_NextToIBLhits_origin_B->Write();
      hist_trk_sharedIBLhits_origin_B->Write();
      hist_trk_splitIBLhits_origin_B->Write();
      hist_trk_nPixhits_origin_B->Write();
      hist_trk_sharedPixhits_origin_B->Write();
      hist_trk_splitPixhits_origin_B->Write();
      hist_trk_nSCThits_origin_B->Write();
      hist_trk_sharedSCThits_origin_B->Write();

      hist_trk_pT_JF_B->Write();
      hist_jet_pT_origin_JF_B->Write();
      hist_jet_pT_origin_truth_label_JF_B->Write();
      hist_jet_pt_JF_B->Write();
      hist_trk_eta_JF_B->Write();
      hist_trk_pT_jet_DR_JF_B->Write();
      hist_trk_origin_JF_B->Write();
      hist_trk_origin_truth_label_JF_B->Write();
      hist_trk_d0_JF_B->Write();
      hist_trk_z0sinth_JF_B->Write();
      hist_trk_d0sig_JF_B->Write();
      hist_trk_z0sinthsig_JF_B->Write();
      hist_trk_d0sig_origin_JF_B->Write();
      hist_trk_z0sinthsig_origin_JF_B->Write();
      hist_trk_logpTfrac_origin_JF_B->Write();
      hist_trk_logDR_origin_JF_B->Write();
      hist_trk_IBLhits_origin_JF_B->Write();
      hist_trk_NextToIBLhits_origin_JF_B->Write();
      hist_trk_sharedIBLhits_origin_JF_B->Write();
      hist_trk_splitIBLhits_origin_JF_B->Write();
      hist_trk_nPixhits_origin_JF_B->Write();
      hist_trk_sharedPixhits_origin_JF_B->Write();
      hist_trk_splitPixhits_origin_JF_B->Write();
      hist_trk_nSCThits_origin_JF_B->Write();
      hist_trk_sharedSCThits_origin_JF_B->Write();

      hist_trk_pT_SV1_B->Write();
      hist_jet_pT_origin_SV1_B->Write();
      hist_jet_pT_origin_truth_label_SV1_B->Write();
      hist_bH_pT_origin_truth_label_SV1_B->Write();
      hist_jet_pt_SV1_B->Write();
      hist_trk_eta_SV1_B->Write();
      hist_trk_pT_jet_DR_SV1_B->Write();
      hist_trk_origin_SV1_B->Write();
      hist_trk_origin_truth_label_SV1_B->Write();
      hist_trk_d0_SV1_B->Write();
      hist_trk_z0sinth_SV1_B->Write();
      hist_trk_d0sig_SV1_B->Write();
      hist_trk_z0sinthsig_SV1_B->Write();
      hist_trk_d0sig_origin_SV1_B->Write();
      hist_trk_z0sinthsig_origin_SV1_B->Write();
      hist_trk_logpTfrac_origin_SV1_B->Write();
      hist_trk_logDR_origin_SV1_B->Write();
      hist_trk_IBLhits_origin_SV1_B->Write();
      hist_trk_NextToIBLhits_origin_SV1_B->Write();
      hist_trk_sharedIBLhits_origin_SV1_B->Write();
      hist_trk_splitIBLhits_origin_SV1_B->Write();
      hist_trk_nPixhits_origin_SV1_B->Write();
      hist_trk_sharedPixhits_origin_SV1_B->Write();
      hist_trk_splitPixhits_origin_SV1_B->Write();
      hist_trk_nSCThits_origin_SV1_B->Write();
      hist_trk_sharedSCThits_origin_SV1_B->Write();
      hist_trk_chi2_SV1_B->Write();
      hist_trk_SV1input_origin_B->Write();

      hist_trk_pT_SV0_B->Write();
      hist_trk_eta_SV0_B->Write();
      hist_trk_pT_jet_DR_SV0_B->Write();
      hist_trk_origin_SV0_B->Write();
      hist_trk_d0_SV0_B->Write();
      hist_trk_z0sinth_SV0_B->Write();
      hist_trk_d0sig_SV0_B->Write();
      hist_trk_z0sinthsig_SV0_B->Write();
      hist_trk_d0sig_origin_SV0_B->Write();
      hist_trk_z0sinthsig_origin_SV0_B->Write();
      hist_trk_logpTfrac_origin_SV0_B->Write();
      hist_trk_logDR_origin_SV0_B->Write();
      hist_trk_IBLhits_origin_SV0_B->Write();
      hist_trk_NextToIBLhits_origin_SV0_B->Write();
      hist_trk_sharedIBLhits_origin_SV0_B->Write();
      hist_trk_splitIBLhits_origin_SV0_B->Write();
      hist_trk_nPixhits_origin_SV0_B->Write();
      hist_trk_sharedPixhits_origin_SV0_B->Write();
      hist_trk_splitPixhits_origin_SV0_B->Write();
      hist_trk_nSCThits_origin_SV0_B->Write();
      hist_trk_sharedSCThits_origin_SV0_B->Write();

      hist_trk_pT_IP3D_B->Write();
      hist_jet_pT_origin_IP3D_B->Write();
      hist_jet_pT_origin_truth_label_IP3D_B->Write();
      hist_bH_pT_origin_truth_label_IP3D_B->Write();
      hist_trk_eta_IP3D_B->Write();
      hist_jet_pt_IP3D_B->Write();
      hist_trk_pT_jet_DR_IP3D_B->Write();
      hist_trk_origin_IP3D_B->Write();
      hist_trk_origin_truth_label_IP3D_B->Write();
      hist_trk_d0_IP3D_B->Write();
      hist_trk_z0sinth_IP3D_B->Write();
      hist_trk_d0sig_IP3D_B->Write();
      hist_trk_z0sinthsig_IP3D_B->Write();
      hist_trk_d0sig_origin_IP3D_B->Write();
      hist_trk_z0sinthsig_origin_IP3D_B->Write();
      hist_trk_logpTfrac_origin_IP3D_B->Write();
      hist_trk_logDR_origin_IP3D_B->Write();
      hist_trk_IBLhits_origin_IP3D_B->Write();
      hist_trk_NextToIBLhits_origin_IP3D_B->Write();
      hist_trk_sharedIBLhits_origin_IP3D_B->Write();
      hist_trk_splitIBLhits_origin_IP3D_B->Write();
      hist_trk_nPixhits_origin_IP3D_B->Write();
      hist_trk_sharedPixhits_origin_IP3D_B->Write();
      hist_trk_splitPixhits_origin_IP3D_B->Write();
      hist_trk_nSCThits_origin_IP3D_B->Write();
      hist_trk_sharedSCThits_origin_IP3D_B->Write();

      hist_trk_pT_IP2D_B->Write();
      hist_trk_eta_IP2D_B->Write();
      hist_trk_pT_jet_DR_IP2D_B->Write();
      hist_trk_origin_IP2D_B->Write();
      hist_trk_d0_IP2D_B->Write();
      hist_trk_z0sinth_IP2D_B->Write();
      hist_trk_d0sig_IP2D_B->Write();
      hist_trk_z0sinthsig_IP2D_B->Write();
      hist_trk_d0sig_origin_IP2D_B->Write();
      hist_trk_z0sinthsig_origin_IP2D_B->Write();
      hist_trk_logpTfrac_origin_IP2D_B->Write();
      hist_trk_logDR_origin_IP2D_B->Write();
      hist_trk_IBLhits_origin_IP2D_B->Write();
      hist_trk_NextToIBLhits_origin_IP2D_B->Write();
      hist_trk_sharedIBLhits_origin_IP2D_B->Write();
      hist_trk_splitIBLhits_origin_IP2D_B->Write();
      hist_trk_nPixhits_origin_IP2D_B->Write();
      hist_trk_sharedPixhits_origin_IP2D_B->Write();
      hist_trk_splitPixhits_origin_IP2D_B->Write();
      hist_trk_nSCThits_origin_IP2D_B->Write();
      hist_trk_sharedSCThits_origin_IP2D_B->Write();

      hist_trk_d0_PUB->Write();
      hist_trk_d0_BB->Write();
      hist_trk_d0_CB->Write();
      hist_trk_d0_FRAGB->Write();
      hist_trk_d0_GEANTB->Write();

      hist_child_pT_B->Write();
      hist_child_eta_B->Write();
      hist_child_phi_B->Write();
      hist_child_Deta_B->Write();
      hist_child_Dphi_B->Write();
      hist_child_Dphi_Deta_B->Write();
      hist_child_DR_B->Write();
      hist_child_pT_DR_B->Write();
      hist_child_pT_jet_DR_B->Write();
      hist_child_pdgID_B->Write();

      hist_child_jetpT_B->Write();
      hist_child_DR_jetpt_B->Write();
      hist_child_DR_bHpt_B->Write();
      hist_child_bHpT_B->Write();

      hist_child_jetpT_C->Write();
      hist_child_DR_jetpt_C->Write();
      hist_child_DR_cHpt_C->Write();
      hist_child_cHpT_C->Write();

      hist_child_pi_notD->Write();
      hist_child_K_notD->Write();
      hist_child_pi->Write();
      hist_child_pi_Lxy_B->Write();
      hist_child_pi_Lxyz_B->Write();
      hist_child_pi_d0_truth_B->Write();
      hist_child_pi_z0_truth_B->Write();
      hist_child_K->Write();
      hist_child_K_Lxy_B->Write();
      hist_child_K_Lxyz_B->Write();
      hist_child_K_d0_truth_B->Write();
      hist_child_K_z0_truth_B->Write();
      hist_child_mu->Write();
      hist_child_mu_Lxy_B->Write();
      hist_child_mu_Lxyz_B->Write();
      hist_child_mu_d0_truth_B->Write();
      hist_child_mu_z0_truth_B->Write();
      hist_child_p->Write();
      hist_child_p_Lxy_B->Write();
      hist_child_p_Lxyz_B->Write();
      hist_child_p_d0_truth_B->Write();
      hist_child_p_z0_truth_B->Write();
      hist_child_e->Write();
      hist_child_e_Lxy_B->Write();
      hist_child_e_Lxyz_B->Write();
      hist_child_e_d0_truth_B->Write();
      hist_child_e_z0_truth_B->Write();

      hist_child_Lxy_B->Write();
      hist_child_Lxyz_B->Write();
//      hist_child_decay_IP->Write();
//      hist_child_nodecay_IP->Write();
//      hist_child_linear_IP->Write();
      hist_child_d0_truth->Write();
      hist_child_d0->Write();
      hist_child_z0->Write();
      hist_child_d0_pT->Write();
      hist_child_z0sinth_B->Write();
      hist_pT_vs_R0_ratio_B->Write();

      if(origin_selection){

        hist_matched_origin_pT_B->Write();
        hist_matched_origin_eta_B->Write();
        hist_matched_origin_phi_B->Write();
        hist_matched_origin_Deta_B->Write();
        hist_matched_origin_Dphi_B->Write();
        hist_matched_origin_Dphi_Deta_B->Write();
        hist_matched_origin_DR_B->Write();
        hist_matched_origin_pT_DR_B->Write();
        hist_matched_origin_pT_jet_DR_B->Write();
        hist_matched_origin_d0_B->Write();
        hist_matched_origin_pdgId_B->Write();
        hist_matched_origin_origin_B->Write();
//        hist_matched_origin_Lxy_B->Write();
//        hist_matched_origin_Lxyz_B->Write();

        hist_matched_origin_jetpT_B->Write();
        hist_trk_DR_jetpt_B->Write();
        hist_trk_DR_bHpt_B->Write();
        hist_matched_origin_bHpT_B->Write();

        hist_matched_origin_jetpT_C->Write();
        hist_trk_DR_jetpt_C->Write();
        hist_trk_DR_cHpt_C->Write();
        hist_matched_origin_cHpT_C->Write();

      }

      if(geometric_selection){

        hist_efficiency_B->Write();
        hist_n_trk->Write();
        hist_n_child->Write();
        hist_n_match->Write();

        hist_matched_pT_B->Write();
        hist_matched_eta_B->Write();
        hist_matched_phi_B->Write();
        hist_matched_Deta_B->Write();
        hist_matched_Dphi_B->Write();
        hist_matched_Dphi_Deta_B->Write();
        hist_matched_DR_B->Write();
        hist_matched_pT_DR_B->Write();
        hist_matched_pT_jet_DR_B->Write();
        hist_matched_pdgId_B->Write();

        hist_matched_child_pi_notD->Write();
        hist_matched_child_K_notD->Write();
        hist_matched_child_pi->Write();
        hist_matched_child_pi_Lxy_B->Write();
        hist_matched_child_pi_Lxyz_B->Write();
        hist_matched_child_pi_d0_truth_B->Write();
        hist_matched_child_pi_z0_truth_B->Write();
        hist_matched_child_K->Write();
        hist_matched_child_K_Lxy_B->Write();
        hist_matched_child_K_Lxyz_B->Write();
        hist_matched_child_K_d0_truth_B->Write();
        hist_matched_child_K_z0_truth_B->Write();
        hist_matched_child_mu->Write();
        hist_matched_child_mu_Lxy_B->Write();
        hist_matched_child_mu_Lxyz_B->Write();
        hist_matched_child_mu_d0_truth_B->Write();
        hist_matched_child_mu_z0_truth_B->Write();
        hist_matched_child_p->Write();
        hist_matched_child_p_Lxy_B->Write();
        hist_matched_child_p_Lxyz_B->Write();
        hist_matched_child_p_d0_truth_B->Write();
        hist_matched_child_p_z0_truth_B->Write();
        hist_matched_child_e->Write();
        hist_matched_child_e_Lxy_B->Write();
        hist_matched_child_e_Lxyz_B->Write();
        hist_matched_child_e_d0_truth_B->Write();
        hist_matched_child_e_z0_truth_B->Write();

        hist_matched_origin_B->Write();
        hist_matched_pT_child_pTfraction_B->Write();
        hist_matched_DR_trk_B->Write();
        hist_matched_DR_trk_pTfraction->Write();
        hist_matched_Lxy_B->Write();
        hist_matched_Lxyz_B->Write();
        hist_matched_d0_B->Write();
  /*
        hist_nomatched_pT_B->Write();
        hist_nomatched_eta_B->Write();
        hist_nomatched_phi_B->Write();
        hist_nomatched_Deta_B->Write();
        hist_nomatched_Dphi_B->Write();
        hist_nomatched_Dphi_Deta_B->Write();
        hist_nomatched_DR_B->Write();
        hist_nomatched_pT_DR_B->Write();
        */
        hist_nomatched_pT_jet_DR_B->Write();
        hist_nomatchedIN_pT_B->Write();
        hist_nomatchedIN_eta_B->Write();
        hist_nomatchedIN_phi_B->Write();
        hist_nomatchedIN_DR_B->Write();
        hist_nomatchedIN_pT_jet_DR_B->Write();
        hist_nomatchedIN_d0_B->Write();
        hist_nomatchedIN_z0sinth_B->Write();
        hist_nomatchedOUT_pT_B->Write();
        hist_nomatchedOUT_eta_B->Write();
        hist_nomatchedOUT_phi_B->Write();
        hist_nomatchedOUT_DR_B->Write();
        hist_nomatchedOUT_pT_jet_DR_B->Write();
        hist_nomatchedOUT_d0_B->Write();
        hist_nomatchedOUT_z0sinth_B->Write();
//        hist_nomatched_pdgId_B->Write();

  /*
        hist_single_matched_pT_B->Write();
        hist_single_matched_eta_B->Write();
        hist_single_matched_phi_B->Write();
        hist_single_matched_Deta_B->Write();
        hist_single_matched_Dphi_B->Write();
        hist_single_matched_Dphi_Deta_B->Write();
        hist_single_matched_DR_B->Write();
        hist_single_matched_pT_DR_B->Write();
        hist_single_matched_pT_jet_DR_B->Write();
        hist_single_matched_pdgId_B->Write();
        hist_single_matched_origin_B->Write();
        hist_single_matched_pT_child_pTfraction_B->Write();
        hist_single_matched_DR_trk_B->Write();
        hist_single_matched_DR_trk_pTfraction->Write();
        hist_single_matched_d0_B->Write();
  */
      }
    }


    file->Close();
/*
    if(selections){
      std::cout<< "\nSELECTIONS\n\n";
      std::cout<< "fraction of events with one single b:\t" << (double) m_b/m_Ntot << "\n";
      std::cout<< "fraction of events with two single b:\t" << (double) m_bb/m_Ntot << "\n";
      std::cout<< "fraction of events without b:\t\t" << (double) m_noB/m_Ntot << "\n";
    }
*/
    std::cout<<"=================================================================="<<std::endl;
    std::cout<<"====== Output file Name : "<<getOutputFNameString()<<std::endl;
    if(discriminants && selections){
      std::cout<< "\nDISCRIMINANTS\n\n";
    }

    if(shrinking_cone){
      std::cout<<"\nSHRINKING CONE\n";
      float R_bin[bin_1]{0},dev_R[bin_1]{0},std_dev_R[bin_1]{0};
      float x[bin_1]{0},y[bin_1]{0},ex[bin_1]{0},ey[bin_1]{0};

      int n=0;
      std::cout <<"\n";
      std::vector<float> tmp(tracksize,0.),max(tracksize,0.);

      int sup=0;
      if(tracksize<=m_track_cut)  sup=tracksize;
      if(tracksize>m_track_cut)  sup=m_track_cut;

      for(int i=0;i<bin_1;i++){

        std::fill(tmp.begin(), tmp.end(), 0);
        std::fill(max.begin(), max.end(), 0);

        if(bin_v.at(i).size()>=m_track_cut){


          for(int l=0;l<sup;l++){
            for(unsigned k=0;k<bin_v.at(i).size();k++){
              if(l==0){
                tmp[l]=bin_v.at(i).at(k);
                if(tmp[l]>max[l]) max[l]=tmp[l];
              }
              if(l>0){
                tmp[l]=bin_v.at(i).at(k);
                if(tmp[l]>max[l] && tmp[l]<max[l-1]) max[l]=tmp[l];
              }
            }
          }

 //        for(int a=0;a<max.size();a++){
 //          std::cout << std::fixed << std::setprecision(5)<<max[a]<<"\t";
 //        }
 //         std::cout<<"\n";

 //        R_bin[i]=std::accumulate(bin_v.at(i).begin(), bin_v.at(i).end(), 0.0)/bin_v.at(i).size();
          R_bin[i]=std::accumulate(max.begin(), max.end(), 0.0)/max.size();

 //        for(int j=0;j<bin_v.at(i).size();j++){
 //          dev_R[i]+=(bin_v.at(i).at(j)-R_bin[i])*(bin_v.at(i).at(j)-R_bin[i]);
 //        }

          for(unsigned j=0;j<max.size();j++){
            dev_R[i]+=(max.at(j)-R_bin[i])*(max.at(j)-R_bin[i]);
          }
 //      std_dev_R[i]=sqrt(dev_R[i]/(bin_v.at(i).size()-1));
          std_dev_R[i]=sqrt(dev_R[i]/(max.size()-1));


 //      std::cout << std::fixed << std::setprecision(5) << R_bin[i] << "\t+-\t" << std_dev_R[i] << "\t\t" << "with\t" << bin_v.at(i).size() << "\tpoints with pT in\t" << "["  << std::fixed << std::setprecision(1) << i*Delta << ","  << (i+1)*Delta << "]\tGeV" << "\n";
          std::cout << std::fixed << std::setprecision(5) << R_bin[i] << "\t+- " << std_dev_R[i] << "\t" << "with\t" << sup << "\tpoints with pT in\t" << "["  << std::fixed << std::setprecision(1) << i*m_Delta_pt_shrCone << ","  << (i+1)*m_Delta_pt_shrCone << "]\tGeV" << "\n";

          x[i]=i*m_Delta_pt_shrCone+m_Delta_pt_shrCone*0.5;
          y[i]=R_bin[i];
          ey[i]=std_dev_R[i];
          n++;

        }
      }
      TMultiGraph *mg = new TMultiGraph();
      TGraphErrors* g_E = new TGraphErrors(n-1, x, y, ex, ey);
 //    g_E->RemovePoint(0);
      g_E->SetMarkerColor(1);
      g_E->SetMarkerSize(1.);
      g_E->SetMarkerStyle(20);
      mg->Add(g_E, "PL");

      double x2[30]{0},y2[30]{0},p1=0.239,p2=-1.220,p3=-1.64*1e-2;
      for(int i=0;i<30;i++){
        x2[i]=200*i/30;
        y2[i]=p1+exp(p2+p3*x2[i]);
      }
      TGraph* g_E2 = new TGraph(30, x2, y2);
      g_E2->SetMarkerColor(2);
 //    g_E2->SetMarkerSize(1.);
      g_E2->SetMarkerStyle(20);
      mg->Add(g_E2,"A");

      mg->Draw("A");
    }

    //if(selection_alg){
    std::cout<<"\nTRUTH LEVEL COMPOSITION\n# of jets: "<< m_njets << " of which:\n" << "b jets: " << m_nBjets << ", c jets: " << m_nCjets << ", light jets: "<< m_nljets << "\nB-C overlap: " << m_nJetBCoverlap<<" ("<< (float) 100*m_nJetBCoverlap/m_nBjets<<" %)\n";
    
    std::cout<<"\nB-C JET CUTS\n"<<1e-3*jet_pT_infcut<<" < Jet pT [GeV] < "<<1e-3*jet_pT_supcut<<"\nJet JVT > "<<jet_JVT_cut<<" (for jet pT in [20, 60] GeV, |eta| < 2.4)\nJet |eta| < "<<jet_eta_cut<<"\nJet Isolation DR > "<<jet_Isol_cut<<", OverlapRemoval (e+-/mu+-), NOT isBad\n";
    
    std::cout<<" CutFlow: passing PtMin           "<<m_njets_2_passPtMin<<std::endl;
    std::cout<<" CutFlow: passing PtMax           "<<m_njets_2_passPtMax<<std::endl;
    std::cout<<" CutFlow: passing Eta Range       "<<m_njets_2_passEtaRange<<std::endl;
    std::cout<<" CutFlow: passing BadMedium       "<<m_njets_2_passBadMedium<<std::endl;
    std::cout<<" CutFlow: passing OR              "<<m_njets_2_passOR<<std::endl;
    std::cout<<" CutFlow: passing ORmu            "<<m_njets_2_passORmu<<std::endl;
    std::cout<<" CutFlow: passing JVT             "<<m_njets_2_passJVT<<std::endl;
    std::cout<<" CutFlow: passing isolation       "<<m_njets_2_passIsol<<std::endl;
    
    
    std::cout<<"\n# of jets: "<< m_njets_2 << " of which:\n" << "b jets: " << m_nBjets_2 << ", c jets: " << m_nCjets_2 << ", light jets: "<< m_nljets_2 << "\nB-C overlap: " <<m_nJetBCoverlap_postJetSel<<" ("<< (float) 100*m_nJetBCoverlap_postJetSel/m_nBjets_2<<" %)\n";
    
    
    std::cout<<"\nB-C HADRONS CUTS\n"<<"b-c Hadr pT > "<<1e-3*m_pT_bcH_truth_cut<<" GeV\n"<<"b-c Hadr DR < "<<m_DR_bcH_truth_cut<<"\n";
    if(!decay_mode.compare("leptonic") || !decay_mode.compare("hadronic")){
      std::cout<<"\n"<<decay_mode<<" mode:";
    }
    
    std::cout<<"\nApplying jet labeling of type <"<<getJetLabeling()<<"> "<<std::endl;
    std::cout<<"b jets: " << m_nBcheck << ", c jets: "<< m_nCcheck << ", light jets "<<m_nlcheck<<"\t\t Total Labeled = "<<m_nBcheck+m_nCcheck+m_nlcheck<<"\t Total Lab/Total "<<Double_t(m_nBcheck+m_nCcheck+m_nlcheck)/Double_t(m_njets_2_passIsol)
	     <<"\nB-C overlap: "<<ov_check<<" ("<< (float) 100*ov_check/m_nBcheck<<" %)\n";
    
    if(selection_alg){
      std::cout<<"\nTRACK/CHILDREN CUTS\npT > "<<1e-3*trk_pT_cut<<" GeV\n"<<"|eta| < "<<trk_eta_cut<<"\n";//<<"|d0| < "<<trk_d0_cut<<" mm"<<"\n"<<"|z0*sin(theta)| < "<<trk_z0sinth_cut<<" mm"<<"\n";
      std::cout<<"\n(JF_ntrk, SV1_ntrk, SV0_ntrk, IP2D_ntrk, IP3D_ntrk)\n"<<JF_ntrk<<", \t"<<SV1_ntrk<<", \t"<<SV0_ntrk<<", \t"<<IP2D_ntrk<<", \t"<<IP3D_ntrk<<"\n";
      std::cout<<"Number of jets: "<<b_cnt<<"\tCut b-jets: "<<b_trkcut_cnt<<"\n";
      std::cout<<"Number of total Input tracks from SV1: "<<SV1input_trks<<".\tNumber of SV1 input jets: "<<nSV1jets<<"\n";
      std::cout<<"output SV1 jets: "<<nSV1outputjets<<" ("<<(float) 100*nSV1outputjets/b_cnt<<" %)\n";
      std::cout<<"output IPxD jets: "<<nIPxDoutputjets<<" ("<<(float) 100*nIPxDoutputjets/b_cnt<<" %)\n";
      std::cout<<"output JF jets: "<<nJFoutputjets<<" ("<<(float) 100*nJFoutputjets/b_cnt<<" %)\n";
      if(origin_selection){
        std::cout<<"\nORIGIN SELECTION\n";
        std::cout<< std::fixed << std::setprecision(3) << "\n";
        std::cout<< "tracks from BC:   \t" << n_trk_B+n_trk_C << "\t(" << (float) 100*(n_trk_B+n_trk_C)/n_trk_pT_cut << " %)\n";
        std::cout<< "tracks from PU:   \t" << n_trk_PU_pT_cut << "\t(" << (float) 100*n_trk_PU_pT_cut/n_trk_pT_cut << " %)\n";
        std::cout<< "tracks from FRAG:\t" << n_trk_FRAG_pT_cut << "\t(" << (float) 100*n_trk_FRAG_pT_cut/n_trk_pT_cut << " %)\n";
        std::cout<< "tracks from GEANT:\t" << n_trk_GEANT_pT_cut << "\t(" << (float) 100*n_trk_GEANT_pT_cut/n_trk_pT_cut << " %)\n";
        std::cout<< "Total number of tracks: "<<n_trk_B+n_trk_C+n_trk_PU_pT_cut+n_trk_FRAG_pT_cut+n_trk_GEANT_pT_cut<<"\n";
        std::cout<< "\nmatches:\t\t"<< n_trk_B+n_trk_C <<"\t("<<mm+mm1_ex<<")"<<"\n";
//        std::cout<< "no matches:\t\t"<< m_den-(n_trk_B+n_trk_C) <<"\n";
        std::cout<< "\naverage efficiency:\t" << (float) (n_trk_B+n_trk_C)/m_den << "\n";
      }
      if(derived_origin==true){
        std::cout<<"\nDERIVED ORIGIN SELECTION\n";
        std::cout<< std::fixed << std::setprecision(3) << "\n";
        for(int i=0;i<6;i++)
          std::cout<< "tracks from "<<der_origin[i].c_str()<<":\t"<<(int) M[i]<<"\n";
        std::cout<< "Total number of tracks: "<<(int) (M[0]+M[1]+M[2]+M[3]+M[4]+M[5])<<"\n";
        std::cout<< "\nmatches:\t\t"<< (int) (M[2]+M[5]) <<"\n";
        std::cout<< "\naverage efficiency:\t" << (float) (M[2]+M[5])/m_den << "\n";
      }
      std::cout<<"\nLight Jets tracks from B/C content:\nGood leptons from B/C: "<<cnt_13<<"\t"<<cnt_103<<"\t"<<cnt_113<<"\nGood hadrons from B/C: "<<cnt_15<<"\t"<<cnt_105<<"\t"<<cnt_115<<"\n";

      if(geometric_selection){
        if(cut && selection_alg)  std::cout<<"\nGEOMETRICAL SELECTION ALGORITHM - CUT\n";
        if(!cut && selection_alg) std::cout<<"\nGEOMETRICAL SELECTION ALGORITHM - NO CUT\n";
        std::cout<< std::fixed << std::setprecision(3) << "\n";
        std::cout<< "number of children:\t" << m_den << "\n";
        std::cout<< "matches:\t\t"<< m_match <<"\t("<<mm+mm2_ex<<")"<<"\n";
        std::cout<< "no matches:\t\t"<< m_nomatch <<"\n";
        std::cout<< "single matches:\t\t" << m_sc << "\t(" <<(float) 100*m_sc/m_match << " %)\n";
        std::cout<< "double matches:\t\t" << m_sc2 << "\t(" <<(float) 100*m_sc2/m_match << " %)\n";
        std::cout<< "triple+ matches:\t" << m_sc3 << "\t(" <<(float) 100*m_sc3/m_match << " %)\n";
        std::cout<< "\naverage efficiency:\t" << (float) m_match/m_den << "\n";
      }

      if(origin_selection || geometric_selection){
        std::cout<<"\nORIGIN/GEOMETRIC OVERLAP\n";
        std::cout<< "\norigin-geometric overlap:\t" << m_match_overlap <<"\t("<<mm<<")"<<"\n";
        std::cout<< "overlap/origin_matches:\t\t" << (float) m_match_overlap/(M[2]+M[5]) <<"\n";
        std::cout<< "overlap/geometric_matches:\t" << (float) m_match_overlap/m_match <<"\n";
        std::cout<<"origin - not-geometry: "<<mm1_ex<<"\n";
        std::cout<<"geometry - not-origin: "<<mm2_ex<<"\n";
        std::cout<<"geometry - not-origin from PU: "<<m_GeomNOr_PU<<"\t("<<(float) 100*m_GeomNOr_PU/mm2_ex<<" %)"<<"\n";
        std::cout<<"geometry - not-origin from Frag: "<<m_GeomNOr_F<<"\t("<<(float) 100*m_GeomNOr_F/mm2_ex<<" %)"<<"\n";
        std::cout<<"geometry - not-origin from Geant: "<<m_GeomNOr_G<<"\t("<<(float) 100*m_GeomNOr_G/mm2_ex<<" %)"<<"\n";
      }
    }


   std::cout<<"=================================================================="<<std::endl;
   std::cout<<"Out of DAOD_selector::Terminate for "<<getOutputFNameString()<<std::endl;
   std::cout<<"===\n===\n"<<std::endl;
   
   return;
}
void DAOD_selector::openOutputFile(std::string fileNameStringID)
{
   std::cout<<"\n In DAOD_selector::openOutputFile"<<std::endl;
   if(lxplus && !debug){
     if(retag){
       std::cout<<"RETAG TRUE\n";
       if(cut){
         std::cout<<"CUT\n";
         file = new TFile(("../output_files/lxplus_output_doretag_cut"+fileNameStringID+".root").c_str(),"RECREATE");
       }
       if(!cut){
         std::cout<<"NO CUT\n";
         file = new TFile(("../output_files/lxplus_output_doretag_nocut"+fileNameStringID+".root").c_str(),"RECREATE");
       }
     }
     if(!retag){
       std::cout<<"RETAG FALSE\n";
       if(cut){
         std::cout<<"CUT\n";
         file = new TFile(("../output_files/lxplus_output_doRetagF_cut"+fileNameStringID+".root").c_str(),"RECREATE");
       }
       if(!cut){
         std::cout<<"NO CUT\n";
         file = new TFile(("../output_files/lxplus_output_doRetagF_nocut"+fileNameStringID+".root").c_str(),"RECREATE");
       }
     }
   }
   if(debug){
     std::cout<<"DEBUG MODE"<<std::endl;
     if(!decay_mode.compare("leptonic") || !decay_mode.compare("hadronic"))
       file = new TFile(("debug_"+decay_mode+fileNameStringID+".root").c_str(),"RECREATE");
     else
       file = new TFile(("debug_"+fileNameStringID+".root").c_str(),"RECREATE");
   }
   if(!lxplus && !debug){
     if(retag){
       std::cout<<"RETAG TRUE\n";
       if(cut){
         std::cout<<"CUT\n";
         file = new TFile(("../output_files/output_doretag_cut"+fileNameStringID+".root").c_str(),"RECREATE");
       }
       if(!cut){
         std::cout<<"NO CUT\n";
         file = new TFile(("../output_files/output_doretag_nocut"+fileNameStringID+".root").c_str(),"RECREATE");
       }
     }
     if(!retag){
       std::cout<<"RETAG FALSE\n";
       if(cut){
         std::cout<<"CUT\n";
         file = new TFile(("../output_files/output_doRetagF_cut"+fileNameStringID+".root").c_str(),"RECREATE");
       }
       if(!cut){
         std::cout<<"NO CUT\n";
         file = new TFile(("../output_files/output_doRetagF_nocut"+fileNameStringID+".root").c_str(),"RECREATE");
       }
     }
   }

   std::cout<<"Output fileName is "<<file->GetName()<<std::endl;
   std::cout<<"Out of DAOD_selector::openOutputFile"<<std::endl;
}

void DAOD_selector::bookHistosForSelections()
{

   if(derived_origin){
     hist_jet_trk_truth_label = new TH1F("derived_jet_trk_truth_label ","derived_jet_trk_truth_label",8,0,7);
   }
   if(selections){
/*
     hist_pt_1 = new TH1F("pT_1", "n1==1", 100, 0., 1000.);
     hist_eta_1 = new TH1F("eta_1", "n1==1", 100, -5., 5.);
     hist_phi_1 = new TH1F("phi_1", "n1==1", 100, -4.,4.);
     hist_E_1 = new TH1F("E_1", "n1==1", 100, 0., 10e5);
     hist_pt_2 = new TH1F("pT_n1==2", "n1==2", 100, 0., 1000.);
     hist_eta_2 = new TH1F("eta_n1==2", "n1==2", 100, -5., 5.);
     hist_phi_2 = new TH1F("phi_n1==2", "n1==2", 100, -4.,4.);
     hist_E_2 = new TH1F("E_n1==2", "n1==2", 100, 0., 10e5);
     hist_pt_2b = new TH1F("pT_2b", "n2==1", 100, 0., 1000.);
     hist_eta_2b = new TH1F("eta_2b", "n2==1", 100, -5., 5.);
     hist_phi_2b = new TH1F("phi_2b", "n2==1", 100, -4.,4.);
     hist_E_2b = new TH1F("E_2b", "n2==1", 100, 0., 10e5);
     hist_pt_3a = new TH1F("pT_3a", "n1==3", 100, 0., 1000.);
     hist_eta_3a = new TH1F("eta_3a", "n1==3", 100, -5., 5.);
     hist_phi_3a = new TH1F("phi_3a", "n1==3", 100, -4.,4.);
     hist_E_3a = new TH1F("E_3a", "n1==3", 100, 0., 10e5);
     hist_pt_3b = new TH1F("pT_3b", "n3==1", 100, 0., 1000.);
     hist_eta_3b = new TH1F("eta_3b", "n3==1", 100, -5., 5.);
     hist_phi_3b = new TH1F("phi_3b", "n3==1", 100, -4.,4.);
     hist_E_3b = new TH1F("E_3b", "n3==1", 100, 0., 10e5);
     hist_pt_4 = new TH1F("pT_4", "n1==4", 100, 0., 1000.);
     hist_eta_4 = new TH1F("eta_4", "n1==4", 100, -5., 5.);
     hist_phi_4 = new TH1F("phi_4", "n1==4", 100,  -4.,4.);
     hist_E_4 = new TH1F("E_4", "n1==4", 100, 0., 10e5);
*/
     hist_pt_B = new TH1F("pT_B", "pT_B", 500, 0., 1000.);
     hist_bHpt_B = new TH1F("bHpT_B", "bHpT_B", 500, 0., 1000.);
     hist_eta_B = new TH1F("eta_B", "eta_B", 100, -5., 5.);
     hist_phi_B = new TH1F("phi_B", "phi_B", 100, -4.,4.);
     hist_E_B = new TH1F("E_B", "E_B", 100, 0., 10e5);
     hist_Bjet_origin = new TH1F("Bjet_origin","Bjet_origin",6,0,6);
     hist_Bjet_cut_origin = new TH1F("Bjet_cut_origin","Bjet_cut_origin",6,0,6);
     hist_Bjet_cut_origin_pT = new TH2F("Bjet_cut_origin_pT","Bjet_cut_origin_pT",6,0,6,300,0.,150.);
     hist_Bjet_cut_origin_jetpT = new TH2F("Bjet_cut_origin_jetpT","Bjet_cut_origin_jetpT",6,0,6,500,0.,1000.);
     hist_Bjet_cut_origin_truth_label = new TH1F("Bjet_cut_origin_truth_label","Bjet_cut_origin_truth_label",121,-1,120);
     hist_Bjet_cut_origin_truth_label_pT = new TH2F("Bjet_cut_origin_truth_label_pT","Bjet_cut_origin_truth_label_pT",121,-1,120,500,0.,1000.);
     hist_Bjet_cut_origin_truth_label_bHpT = new TH2F("Bjet_cut_origin_truth_label_bHpT","Bjet_cut_origin_truth_label_bHpT",121,-1,120,500,0.,1000.);


     hist_pt_C = new TH1F("pT_C", "C", 500, 0., 1000.);
     hist_cHpt_C = new TH1F("cHpT_C", "cHpT_C", 500, 0., 1000.);
     hist_eta_C = new TH1F("eta_C", "C", 100, -5., 5.);
     hist_phi_C = new TH1F("phi_C", "C", 100, -4.,4.);
     hist_E_C = new TH1F("E_C", "C", 100, 0., 10e5);
     hist_Cjet_origin = new TH1F("Cjet_origin","Cjet_origin",6,0,6);
     hist_Cjet_cut_origin = new TH1F("Cjet_cut_origin","Cjet_cut_origin",6,0,6);
     hist_Cjet_cut_origin_pT = new TH2F("Cjet_cut_origin_pT","Cjet_cut_origin_pT",6,0,6,300,0.,150.);
     hist_Cjet_cut_origin_jetpT = new TH2F("Cjet_cut_origin_jetpT","Cjet_cut_origin_jetpT",6,0,6,500,0.,1000.);
     hist_Cjet_cut_origin_truth_label = new TH1F("Cjet_cut_origin_truth_label","Cjet_cut_origin_truth_label",121,-1,120);
     hist_Cjet_cut_origin_truth_label_pT = new TH2F("Cjet_cut_origin_truth_label_pT","Cjet_cut_origin_truth_label_pT",121,-1,120,500,0.,1000.);
     hist_Cjet_cut_origin_truth_label_cHpT = new TH2F("Cjet_cut_origin_truth_label_cHpT","Cjet_cut_origin_truth_label_cHpT",121,-1,120,500,0.,1000.);

     hist_pt_l = new TH1F("pT_l", "l", 500, 0., 1000.);
     hist_eta_l = new TH1F("eta_l", "l", 100, -5., 5.);
     hist_phi_l = new TH1F("phi_l", "l", 100, -4.,4.);
     hist_E_l = new TH1F("E_l", "l", 100, 0., 10e5);
     hist_ljet_origin = new TH1F("ljet_origin","ljet_origin",6,0,6);
     hist_ljet_cut_origin = new TH1F("ljet_cut_origin","ljet_cut_origin",6,0,6);
     hist_ljet_cut_origin_pT = new TH2F("ljet_cut_origin_pT","ljet_cut_origin_pT",6,0,6,300,0.,150.);
     hist_ljet_cut_origin_jetpT = new TH2F("ljet_cut_origin_jetpT","ljet_cut_origin_jetpT",6,0,6,500,0.,1000.);
     hist_ljet_cut_origin_truth_label = new TH1F("ljet_cut_origin_truth_label","ljet_cut_origin_truth_label",121,-1,120);
     hist_ljet_cut_origin_truth_label_pT = new TH2F("ljet_cut_origin_truth_label_pT","ljet_cut_origin_truth_label_pT",121,-1,120,500,0.,1000.);

     hist_nBjets = new TH1F("nBjets","nBjets",1,1,2);
     hist_n_tracks_jetpt_B = new TH2F("n_tracks_jetpt_B", "n_tracks_jetpt_B", 500,0.,1000.,100, 0, 100);
     hist_n_tracks_bHpt_B = new TH2F("n_tracks_bHpt_B", "n_tracks_bHpt_B", 500,0.,1000.,100, 0, 100);
     hist_n_BCtracks_jetpt_B = new TH2F("n_BCtracks_jetpt_B", "n_BCtracks_jetpt_B", 500,0.,1000.,100, 0, 100);
     hist_n_BCtracks_bHpt_B = new TH2F("n_BCtracks_bHpt_B", "n_BCtracks_bHpt_B", 500,0.,1000.,100, 0, 100);

     hist_nCjets = new TH1F("nCjets","nCjets",1,1,2);
     hist_n_tracks_jetpt_C = new TH2F("n_tracks_jetpt_C", "n_tracks_jetpt_C", 500,0.,1000.,100, 0, 100);
     hist_n_tracks_cHpt_C = new TH2F("n_tracks_cHpt_C", "n_tracks_cHpt_C", 500,0.,1000.,100, 0, 100);
     hist_n_BCtracks_jetpt_C = new TH2F("n_BCtracks_jetpt_C", "n_BCtracks_jetpt_C", 500,0.,1000.,100, 0, 100);
     hist_n_BCtracks_cHpt_C = new TH2F("n_BCtracks_cHpt_C", "n_BCtracks_cHpt_C", 500,0.,1000.,100, 0, 100);

     hist_nljets = new TH1F("nljets","nljets",1,1,2);
     hist_n_tracks_jetpt_l = new TH2F("n_tracks_jetpt_l", "n_tracks_jetpt_l", 500,0.,1000.,100, 0, 100);
    
   }
}

void DAOD_selector::bookHistosForFlavorLabelStudies()
{ 
  if (fDoFlavorLabelMatrix)
    {
      hist2_jetFlavorMatrix = new TH2F("jetFlavorLabelMatrix","jetFlavorLabelMatrix",6,-0.5,5.5,6,-0.5,5.5);

     //Histograms for labeling studies 
     std::string HistoName[6] = {"1B", "1D0B", "0B0D", "2B", "2D0B", "2B2D"};
     std::string hVariable;
     std::string hNameLab;
     for(uint32_t iX=0; iX<6; iX++){
       for(uint32_t iY=0; iY<6; iY++){
	 hNameLab = HistoName[iX]+"_"+HistoName[iY];
	 hVariable = "Labels_"+hNameLab+"_jetPt";
         //HistopT[iX][iY] = new TH1D(hVariable.c_str(),(HistoName[iX]+" (cone labeling) - "+HistoName[iY]+" (ghost labeling); Jet p_{T} [GeV]; Entries").c_str(),50,0,300);
	 BookHisto(hVariable, (HistoName[iX]+" (cone labeling) - "+HistoName[iY]+" (ghost labeling); Jet p_{T} [GeV]; Entries"), 60,0.,300.);
	 hVariable = "Labels_"+hNameLab+"_jetEta";
	 BookHisto(hVariable, (HistoName[iX]+" (cone labeling) - "+HistoName[iY]+" (ghost labeling); Jet \eta; Entries"), 100,-3.,3.);
	 hVariable = "Labels_"+hNameLab+"_bHPt";
	 BookHisto(hVariable, (HistoName[iX]+" (cone labeling) - "+HistoName[iY]+" (ghost labeling); B-hadron p_{T} [GeV]; Entries"), 60, 0.,300.);
	 hVariable = "Labels_"+hNameLab+"_bHPtFraction";
	 BookHisto(hVariable, (HistoName[iX]+" (cone labeling) - "+HistoName[iY]+" (ghost labeling); B-hadron p_{T}/jet p_{T}; Entries"), 60, 0.,300.);
	 hVariable = "Labels_"+hNameLab+"_bHjetDR";
	 BookHisto(hVariable, (HistoName[iX]+" (cone labeling) - "+HistoName[iY]+" (ghost labeling); B-hadron - jet DR; Entries"), 60, 0.,300.);
	 std::string stringAlgo[5]={"IP2D","IP3D","RNNIP","SV1","JF"};
	 for (unsigned int algoBit=0; algoBit<5; ++algoBit) /// 0=IP2D, 1=IP3D, 2=RNNIP, 3=SV1, 4=JF
	   {
	     hVariable="Labels_"+hNameLab+"_nTrkAlgo_"+stringAlgo[algoBit];
	     BookHisto(hVariable,
		       (HistoName[iX]+" (cone labeling) - "+HistoName[iY]+" (ghost labeling); nTrk used by "+stringAlgo[algoBit]+" ; Entries"),
		       11,-0.5,10.5);
	   }
	 BookHisto("Labels_"+hNameLab+"_SV1nVtx",  (HistoName[iX]+" (cone labeling) - "+HistoName[iY]+" (ghost labeling); SV1 nVertices")       ,4,  -0.5 ,3.5);
	 BookHisto("Labels_"+hNameLab+"_SV1mVtx",  (HistoName[iX]+" (cone labeling) - "+HistoName[iY]+" (ghost labeling); SV1 Vertex mass")     ,50,  0.  ,20.0);
	 BookHisto("Labels_"+hNameLab+"_SV1eFc" ,  (HistoName[iX]+" (cone labeling) - "+HistoName[iY]+" (ghost labeling); SV1 Vtx Energy frac."),50, -0.1  ,1.1);
	 BookHisto("Labels_"+hNameLab+"_SV1sig3d", (HistoName[iX]+" (cone labeling) - "+HistoName[iY]+" (ghost labeling); SV1 Vtx 3d sign")     ,100, 0.  ,1000.0);
	 BookHisto("Labels_"+hNameLab+"_SV1dR",    (HistoName[iX]+" (cone labeling) - "+HistoName[iY]+" (ghost labeling); SV1 Vtx DeltaR")      ,100, 0.  ,1.0);

	 BookHisto("Labels_"+hNameLab+"_JFnVtx",  (HistoName[iX]+" (cone labeling) - "+HistoName[iY]+" (ghost labeling); JF nVertices")       ,10,  -0.5 ,9.5);
	 BookHisto("Labels_"+hNameLab+"_JFmVtx",  (HistoName[iX]+" (cone labeling) - "+HistoName[iY]+" (ghost labeling); JF Vertex mass")     ,50,  0.  ,20.0);
	 BookHisto("Labels_"+hNameLab+"_JFeFc" ,  (HistoName[iX]+" (cone labeling) - "+HistoName[iY]+" (ghost labeling); JF Vtx Energy frac."),50,  -.1  ,1.1);
	 BookHisto("Labels_"+hNameLab+"_JFsig3d", (HistoName[iX]+" (cone labeling) - "+HistoName[iY]+" (ghost labeling); JF Vtx 3d sign")     ,100, 0.  ,1000.0);
	 BookHisto("Labels_"+hNameLab+"_JFdR",    (HistoName[iX]+" (cone labeling) - "+HistoName[iY]+" (ghost labeling); JF Vtx DeltaR")      ,100, 0.  ,1.0);
       }
     }      
    }
  return;
}

void DAOD_selector::bookHistosForDiscriminants()
{
   if(discriminants){

     int nbin=1400;
     hist_ip2d_llr_l = new TH1F("ip2d_llr_l","ip2d_llr_l",nbin, -20., 50.);
     hist_ip2d_llr_B = new TH1F("ip2d_llr_B", "ip2d_llr_B", nbin, -20., 50.);
     hist_ip2d_llr_C = new TH1F("ip2d_llr_C", "ip2d_llr_C", nbin, -20., 50.);
     hist_ip2d_llr_jetpt_B = new TH2F("ip2d_llr_jetpt_B","ip2d_llr_jetpt_B",nbin,-20,50, 500, 0., 1000.);
     hist_ip2d_llr_jetpt_singleB = new TH2F("ip2d_llr_jetpt_singleB","ip2d_llr_jetpt_singleB",nbin,-20,50, 500, 0., 1000.);
     hist_ip2d_llr_jetpt_C = new TH2F("ip2d_llr_jetpt_C","ip2d_llr_jetpt_C",nbin,-20,50, 500, 0., 1000.);
     hist_ip2d_llr_jetpt_singleC = new TH2F("ip2d_llr_jetpt_singleC","ip2d_llr_jetpt_singleC",nbin,-20,50, 500, 0., 1000.);


     hist_ip3d_llr_l = new TH1F("ip3d_llr_l","ip3d_llr_l",nbin, -20., 50.);
     hist_ip3d_llr_B = new TH1F("ip3d_llr_B", "ip3d_llr_B", nbin, -20., 50.);
     hist_ip3d_llr_C = new TH1F("ip3d_llr_C", "ip3d_llr_C", nbin, -20., 50.);
     hist_ip3d_llr_jetpt_B = new TH2F("ip3d_llr_jetpt_B","ip3d_llr_jetpt_B",nbin,-20,50, 500, 0., 1000.);
     hist_ip3d_llr_jetpt_singleB = new TH2F("ip3d_llr_jetpt_singleB","ip3d_llr_jetpt_singleB",nbin,-20,50, 500, 0., 1000.);
     hist_ip3d_llr_jetpt_C = new TH2F("ip3d_llr_jetpt_C","ip3d_llr_jetpt_C",nbin,-20,50, 500, 0., 1000.);
     hist_ip3d_llr_jetpt_singleC = new TH2F("ip3d_llr_jetpt_singleC","ip3d_llr_jetpt_singleC",nbin,-20,50, 500, 0., 1000.);

     hist_rnnip_llr_l = new TH1F("rnnip_llr_l","rnnip_llr_l",nbin, -20., 50.);
     hist_rnnip_llr_B = new TH1F("rnnip_llr_B", "rnnip_llr_B", nbin, -20., 50.);
     hist_rnnip_llr_C = new TH1F("rnnip_llr_C", "rnnip_llr_C", nbin, -20., 50.);
     hist_rnnip_llr_jetpt_B = new TH2F("rnnip_llr_jetpt_B","rnnip_llr_jetpt_B",nbin,-20,50, 500, 0., 1000.);
     hist_rnnip_llr_jetpt_singleB = new TH2F("rnnip_llr_jetpt_singleB","rnnip_llr_jetpt_singleB",nbin,-20,50, 500, 0., 1000.);
     hist_rnnip_llr_jetpt_C = new TH2F("rnnip_llr_jetpt_C","rnnip_llr_jetpt_C",nbin,-20,50, 500, 0., 1000.);
     hist_rnnip_llr_jetpt_singleC = new TH2F("rnnip_llr_jetpt_singleC","rnnip_llr_jetpt_singleC",nbin,-20,50, 500, 0., 1000.);

     hist_sv1_llr_l = new TH1F("sv1_llr_l","sv1_llr_l", nbin, -20., 50.);
     hist_sv1_llr_B = new TH1F("sv1_llr_B", "sv1_llr_B", nbin, -20., 50.);
     hist_sv1_llr_C = new TH1F("sv1_llr_C", "sv1_llr_C", nbin, -20., 50.);
     hist_sv1_llr_jetpt_B = new TH2F("sv1_llr_jetpt_B","sv1_llr_jetpt_B",nbin,-20,50, 500, 0., 1000.);
     hist_sv1_llr_jetpt_singleB = new TH2F("sv1_llr_jetpt_singleB","sv1_llr_jetpt_singleB",nbin,-20,50, 500, 0., 1000.);
     hist_sv1_llr_jetpt_C = new TH2F("sv1_llr_jetpt_C","sv1_llr_jetpt_C",nbin,-20,50, 500, 0., 1000.);
     hist_sv1_llr_jetpt_singleC = new TH2F("sv1_llr_jetpt_singleC","sv1_llr_jetpt_singleC",nbin,-20,50, 500, 0., 1000.);

     hist_jet_sv1_Nvtx_B  = new TH1F("jet_sv1_Nvtx_B","jet_sv1_Nvtx_B",3, 0., 3);
     hist_jet_sv1_ntrkv_B = new TH1F("jet_sv1_ntrkv_B","jet_sv1_ntrkv_B",50, 0., 50.);
     hist_jet_sv1_n2t_B = new TH1F("jet_sv1_n2t_B","jet_sv1_n2t_B",100, 0., 100.);
     hist_jet_sv1_m_B = new TH1F("jet_sv1_m_B","jet_sv1_m_B_B",1000, 0., 10000);
     hist_jet_sv1_efc_B = new TH1F("jet_sv1_efc_B","jet_sv1_efc_B",1000, 0., 1.);
     hist_jet_sv1_sig3d_B = new TH1F("jet_sv1_sig3d_B","jet_sv1_sig3d_B",1000, 0., 1000.);
     hist_jet_sv1_deltaR_B = new TH1F("jet_sv1_deltaR_B","jet_sv1_deltaR_B",1000, 0., 10.);
     hist_jet_sv1_Lxy_B  = new TH1F("jet_sv1_Lxy_B","jet_sv1_Lxy_B",1000, 0., 1000.);
     hist_jet_sv1_L3d_B  = new TH1F("jet_sv1_L3d_B","jet_sv1_L3d_B",1000, 0., 1000.);

     hist_jf_llr_l = new TH1F("jf_llr_l","jf_llr_l",nbin, -20., 50.);
     hist_jf_llr_B = new TH1F("jf_llr_B", "jf_llr_B", nbin, -20., 50.);
     hist_jf_llr_C = new TH1F("jf_llr_C", "jf_llr_C", nbin, -20., 50.);
     hist_jf_llr_jetpt_B = new TH2F("jf_llr_jetpt_B","jf_llr_jetpt_B",nbin,-20,50, 500, 0., 1000.);
     hist_jf_llr_jetpt_singleB = new TH2F("jf_llr_jetpt_singleB","jf_llr_jetpt_singleB",nbin,-20,50, 500, 0., 1000.);
     hist_jf_llr_jetpt_C = new TH2F("jf_llr_jetpt_C","jf_llr_jetpt_C",nbin,-20,50, 500, 0., 1000.);
     hist_jf_llr_jetpt_singleC = new TH2F("jf_llr_jetpt_singleC","jf_llr_jetpt_singleC",nbin,-20,50, 500, 0., 1000.);

     hist_dl1_l = new TH1F("dl1_llr_l","dl1_llr_l",nbin, -20., 50.);
     hist_dl1_C = new TH1F("dl1_llr_C","dl1_llr_C",nbin, -20., 50.);
     hist_dl1_B = new TH1F("dl1_llr_B","dl1_llr_B",nbin, -20., 50.);
     hist_dl1_llr_jetpt_B = new TH2F("dl1_llr_jetpt_B","dl1_llr_jetpt_B",nbin,-20,50, 500, 0., 1000.);
     hist_dl1_llr_jetpt_singleB = new TH2F("dl1_llr_jetpt_singleB","dl1_llr_jetpt_singleB",nbin,-20,50, 500, 0., 1000.);
     hist_dl1_llr_jetpt_C = new TH2F("dl1_llr_jetpt_C","dl1_llr_jetpt_C",nbin,-20,50, 500, 0., 1000.);
     hist_dl1_llr_jetpt_singleC = new TH2F("dl1_llr_jetpt_singleC","dl1_llr_jetpt_singleC",nbin,-20,50, 500, 0., 1000.);
   }

}

void DAOD_selector::bookHistosForShrinkingCone()
{
   if(shrinking_cone){

     hist_n_tracks = new TH1F("n tracks", "n1==1", 100, 0, 50);
     hist_tracks_DR = new TH2F("n tracks-Delta_R", "n1==1", 100, 0., 1., 100., 0., 100.);
     hist_DR_1 = new TH1F("Delta_R", "n1==1", 100, 0., 10.);
//   hist_std_dev_DR_1 = new TH1F("Standard Deviation from Delta_R", "n1==1", 100, 0., 1.);
     jet_DR_pT = new TH2F("jets Delta_R-pT","n1==1", 500, 0., 500., 200, 0., 2.);

   }

}


void DAOD_selector::bookHistosForSelectionAlgos()
{

   if(selection_alg){

     hist_efficiency_B = new TH1F("efficiency","efficiency",120,-0.1,1.1);
     hist_n_trk = new TH1D("n_trk","n_trk",100,0.,100);
     hist_n_child = new TH1D("n_child","n_child",100,0.,100);
     hist_n_match = new TH1D("n_match","n_match",100,0.,100);

     hist_trk_pT_B = new TH1F("trk_pT_B", "trk_pT_B", 3000, 0., 300.);
     hist_trk_eta_B = new TH1F("trk_eta_B","trk_eta_B",260,-2.6,2.6);
     hist_trk_phi_B = new TH1F("trk_phi_B","trk_phi_B",200,-4.,4.);
     hist_trk_Deta_B = new TH1F("trk_Deta_B","trk_Deta_B",200,-1.,1.);
     hist_trk_Dphi_B = new TH1F("trk_Dphi_B","trk_Dphi_B",200,-1.,1.);
     hist_trk_Dphi_Deta_B = new TH2F("trk_Dphi_Deta","trk_Dphi_Deta", 100, -1.,1., 100, -1.,1.);
     hist_trk_DR_B = new TH1F("trk_DR_B","trk_DR_B",200,0.,2.);
     hist_trk_pT_DR_B = new TH2F("trk_pT_DR_B","trk_pT_DR_B",300,0.,150.,200,0.,2.);
     hist_trk_pT_jet_DR_B = new TH2F("trk_pT_jet_DR_B","trk_pT_jet_DR_B",2000,0.,1000,200,0.,2.);
     hist_trk_pdgId_B = new TH1F("trk_pdgID_B","trk_pdgID_B",200000,-100000,100000);
     hist_trk_origin_B = new TH1F("trk_origin_B","trk_origin_B",6,0,6);
     hist_trk_origin_truth_label_B = new TH1F("trk_origin_truth_label_B","trk_origin_truth_label_B",121,-1,120);

     hist_trk_d0_B = new TH1F("trk_d0_B","trk_d0_B",3000,-15.,15.);
     hist_trk_d0_jet_pT_BC_B = new TH2F("trk_d0_jet_pT_BC_B","trk_d0_jet_pT_BC_B",2000,0.,1000,3000,-15.,15.);
     hist_trk_DR_jet_pt_BC = new TH2F("trk_DR_jet_pt_BC","trk_DR_jet_pt_BC", 500,0.,1000.,200,0.,2.);
     hist_trk_sigd0_B = new TH1F("trk_sigd0_B","trk_sigd0_B",1000,0.,10.);
     hist_trk_z0_B = new TH1F("trk_z0_B","trk_z0_B",3000,-15.,15.);
     hist_trk_sigz0_B = new TH1F("trk_sigz0_B","trk_sigz0_B",1000,0.,10.);
     hist_trk_z0sinth_B = new TH1F("trk_z0sinth_B","trk_z0sinth_B",1000,-5.,5.);
     hist_trk_d0sig_B = new TH1F("trk_d0sig_B","trk_d0sig_B",3000,-15.,15.);
     hist_trk_z0sinthsig_B = new TH1F("trk_z0sinthsig_B","trk_z0sinthsig_B",3000,-15.,15.);
     hist_trk_d0sig_origin_B = new TH2F("trk_d0sig_origin_B","trk_d0sig_origin_B",300,-15.,15.,6,0,6);
     hist_trk_d0sig_truthlabel_B = new TH2F("trk_d0sig_truthlabel_B","trk_d0sig_truthlabel_B",300,-15.,15.,121,-1,120);
     hist_trk_z0sinthsig_origin_B = new TH2F("trk_z0sinthsig_origin_B","trk_z0sinthsig_origin_B",300,-15.,15.,6,0,6);
     hist_trk_logpTfrac_origin_B = new TH2F("trk_logpTfrac_origin_B","trk_logpTfrac_origin_B",300,-10.,2,6,0,6);
     hist_trk_logDR_origin_B = new TH2F("trk_logDR_origin_B","trk_logDR_origin_B",200,-10.,2.,6,0,6);
     hist_trk_IBLhits_origin_B = new TH2F("trk_IBLhits_origin_B","trk_IBLhits_origin_B",12,0.,12.,6,0,6);
     hist_trk_NextToIBLhits_origin_B = new TH2F("trk_NextToIBLhits_origin_B","trk_NextToIBLhits_origin_B",12,0.,12.,6,0,6);
     hist_trk_sharedIBLhits_origin_B = new TH2F("trk_sharedIBLhits_origin_B","trk_sharedIBLhits_origin_B",12,0.,12.,6,0,6);
     hist_trk_splitIBLhits_origin_B = new TH2F("trk_splitIBLhits_origin_B","trk_splitIBLhits_origin_B",12,0.,12.,6,0,6);
     hist_trk_nPixhits_origin_B = new TH2F("trk_nPixhits_origin_B","trk_nPixhits_origin_B",12,0.,12.,6,0,6);
     hist_trk_sharedPixhits_origin_B = new TH2F("trk_sharedPixhits_origin_B","trk_sharedPixhits_origin_B",12,0.,12.,6,0,6);
     hist_trk_splitPixhits_origin_B = new TH2F("trk_splitPixhits_origin_B","trk_splitPixhits_origin_B",12,0.,12.,6,0,6);
     hist_trk_nSCThits_origin_B = new TH2F("trk_nSCThits_origin_B","trk_nSCThits_origin_B",24,0.,24.,6,0,6);
     hist_trk_sharedSCThits_origin_B = new TH2F("trk_sharedSCThits_origin_B","trk_sharedSCThits_origin_B",24,0.,24.,6,0,6);

     hist_trk_pT_JF_B = new TH1F("trk_pT_JF_B", "trk_pT_JF_B", 3000, 0., 300.);
     hist_jet_pT_origin_JF_B = new TH2F("trk_jet_pT_origin_JF_B", "trk_jet_pT_origin_JF_B",6,0,6,500, 0., 1000.);
     hist_jet_pT_origin_truth_label_JF_B = new TH2F("trk_jet_pT_origin_truth_label_JF_B", "trk_jet_pT_origin_truth_label_JF_B",121,-1,120,500, 0., 1000.);
     hist_jet_pt_JF_B = new TH1F("trk_jet_pT_JF_B", "trk_jet_pT_JF_B", 500, 0., 1000.);
     hist_trk_eta_JF_B = new TH1F("trk_JF_eta_B","trk_JF_eta_B",260,-2.6,2.6);
     hist_trk_pT_jet_DR_JF_B = new TH2F("trk_pT_jet_DR_JF_B","trk_pT_jet_DR_JF_B",2000,0.,1000,200,0.,2.);
     hist_trk_origin_JF_B = new TH1F("trk_JF_origin_B","trk_JF_origin_B",6,0,6);
     hist_trk_origin_truth_label_JF_B = new TH1F("trk_origin_truth_label_JF_B","trk_origin_truth_label_JF_B",121,-1,120);
     hist_trk_d0_JF_B = new TH1F("trk_d0_JF_B","trk_d0_JF_B",3000,-15.,15.);
     hist_trk_z0sinth_JF_B = new TH1F("trk_z0sinth_JF_B","trk_z0sinth_JF_B",1000,-5.,5.);
     hist_trk_d0sig_JF_B = new TH1F("trk_d0sig_JF_B","trk_d0sig_JF_B",300,-15.,15.);
     hist_trk_z0sinthsig_JF_B = new TH1F("trk_z0sinthsig_JF_B","trk_z0sinthsig_JF_B",300,-15.,15.);
     hist_trk_d0sig_origin_JF_B = new TH2F("trk_d0sig_origin_JF_B","trk_d0sig_origin_JF_B",300,-15.,15.,6,0,6);
     hist_trk_z0sinthsig_origin_JF_B = new TH2F("trk_z0sinthsig_origin_JF_B","trk_z0sinthsig_origin_JF_B",300,-15.,15.,6,0,6);
     hist_trk_logpTfrac_origin_JF_B = new TH2F("trk_logpTfrac_origin_JF_B","trk_logpTfrac_origin_JF_B",300,-10.,2,6,0,6);
     hist_trk_logDR_origin_JF_B = new TH2F("trk_logDR_origin_JF_B","trk_logDR_origin_JF_B",200,-10.,2.,6,0,6);
     hist_trk_IBLhits_origin_JF_B = new TH2F("trk_IBLhits_origin_JF_B","trk_IBLhits_origin_JF_B",12,0.,12.,6,0,6);
     hist_trk_NextToIBLhits_origin_JF_B = new TH2F("trk_NextToIBLhits_origin_JF_B","trk_NextToIBLhits_origin_JF_B",12,0.,12.,6,0,6);
     hist_trk_sharedIBLhits_origin_JF_B = new TH2F("trk_sharedIBLhits_origin_JF_B","trk_sharedIBLhits_origin_JF_B",12,0.,12.,6,0,6);
     hist_trk_splitIBLhits_origin_JF_B = new TH2F("trk_splitIBLhits_origin_JF_B","trk_splitIBLhits_origin_JF_B",12,0.,12.,6,0,6);
     hist_trk_nPixhits_origin_JF_B = new TH2F("trk_nPixhits_origin_JF_B","trk_nPixhits_origin_JF_B",12,0.,12.,6,0,6);
     hist_trk_sharedPixhits_origin_JF_B = new TH2F("trk_sharedPixhits_origin_JF_B","trk_sharedPixhits_origin_JF_B",12,0.,12.,6,0,6);
     hist_trk_splitPixhits_origin_JF_B = new TH2F("trk_splitPixhits_origin_JF_B","trk_splitPixhits_origin_JF_B",12,0.,12.,6,0,6);
     hist_trk_nSCThits_origin_JF_B = new TH2F("trk_nSCThits_origin_JF_B","trk_nSCThits_origin_JF_B",24,0.,24.,6,0,6);
     hist_trk_sharedSCThits_origin_JF_B = new TH2F("trk_sharedSCThits_origin_JF_B","trk_sharedSCThits_origin_JF_B",24,0.,24.,6,0,6);

     hist_trk_pT_SV1_B = new TH1F("trk_pT_SV1_B", "trk_pT_SV1_B", 3000, 0., 300.);
     hist_jet_pT_origin_SV1_B = new TH2F("trk_jet_pT_origin_SV1_B", "trk_jet_pT_origin_SV1_B",6,0,6,500, 0., 1000.);
     hist_jet_pT_origin_truth_label_SV1_B = new TH2F("trk_jet_pT_origin_truth_label_SV1_B", "trk_jet_pT_origin_truth_label_SV1_B",121,-1,120,500, 0., 1000.);
     hist_bH_pT_origin_truth_label_SV1_B = new TH2F("trk_bH_pT_origin_truth_label_SV1_B", "trk_bH_pT_origin_truth_label_SV1_B",121,-1,120,500, 0., 1000.);
     hist_jet_pt_SV1_B = new TH1F("trk_jet_pT_SV1_B", "trk_jet_pT_SV1_B", 500, 0., 1000.);
     hist_trk_eta_SV1_B = new TH1F("trk_SV1_eta_B","trk_SV1_eta_B",260,-2.6,2.6);
     hist_trk_pT_jet_DR_SV1_B = new TH2F("trk_pT_jet_DR_SV1_B","trk_pT_jet_DR_SV1_B",2000,0.,1000,200,0.,2.);
     hist_trk_origin_SV1_B = new TH1F("trk_SV1_origin_B","trk_SV1_origin_B",6,0,6);
     hist_trk_origin_truth_label_SV1_B = new TH1F("trk_origin_truth_label_SV1_B","trk_origin_truth_label_SV1_B",121,-1,120);
     hist_trk_d0_SV1_B = new TH1F("trk_d0_SV1_B","trk_d0_SV1_B",3000,-15.,15.);
     hist_trk_z0sinth_SV1_B = new TH1F("trk_z0sinth_SV1_B","trk_z0sinth_SV1_B",1000,-5.,5.);
     hist_trk_d0sig_SV1_B = new TH1F("trk_d0sig_SV1_B","trk_d0sig_SV1_B",300,-15.,15.);
     hist_trk_z0sinthsig_SV1_B = new TH1F("trk_z0sinthsig_SV1_B","trk_z0sinthsig_SV1_B",300,-15.,15.);
     hist_trk_d0sig_origin_SV1_B = new TH2F("trk_d0sig_origin_SV1_B","trk_d0sig_origin_SV1_B",300,-15.,15.,6,0,6);
     hist_trk_z0sinthsig_origin_SV1_B = new TH2F("trk_z0sinthsig_origin_SV1_B","trk_z0sinthsig_origin_SV1_B",300,-15.,15.,6,0,6);
     hist_trk_logpTfrac_origin_SV1_B = new TH2F("trk_logpTfrac_origin_SV1_B","trk_logpTfrac_origin_SV1_B",300,-10.,2,6,0,6);
     hist_trk_logDR_origin_SV1_B = new TH2F("trk_logDR_origin_SV1_B","trk_logDR_origin_SV1_B",200,-10.,2.,6,0,6);
     hist_trk_IBLhits_origin_SV1_B = new TH2F("trk_IBLhits_origin_SV1_B","trk_IBLhits_origin_SV1_B",12,0.,12.,6,0,6);
     hist_trk_NextToIBLhits_origin_SV1_B = new TH2F("trk_NextToIBLhits_origin_SV1_B","trk_NextToIBLhits_origin_SV1_B",12,0.,12.,6,0,6);
     hist_trk_sharedIBLhits_origin_SV1_B = new TH2F("trk_sharedIBLhits_origin_SV1_B","trk_sharedIBLhits_origin_SV1_B",12,0.,12.,6,0,6);
     hist_trk_splitIBLhits_origin_SV1_B = new TH2F("trk_splitIBLhits_origin_SV1_B","trk_splitIBLhits_origin_SV1_B",12,0.,12.,6,0,6);
     hist_trk_nPixhits_origin_SV1_B = new TH2F("trk_nPixhits_origin_SV1_B","trk_nPixhits_origin_SV1_B",12,0.,12.,6,0,6);
     hist_trk_sharedPixhits_origin_SV1_B = new TH2F("trk_sharedPixhits_origin_SV1_B","trk_sharedPixhits_origin_SV1_B",12,0.,12.,6,0,6);
     hist_trk_splitPixhits_origin_SV1_B = new TH2F("trk_splitPixhits_origin_SV1_B","trk_splitPixhits_origin_SV1_B",12,0.,12.,6,0,6);
     hist_trk_nSCThits_origin_SV1_B = new TH2F("trk_nSCThits_origin_SV1_B","trk_nSCThits_origin_SV1_B",24,0.,24.,6,0,6);
     hist_trk_sharedSCThits_origin_SV1_B = new TH2F("trk_sharedSCThits_origin_SV1_B","trk_sharedSCThits_origin_SV1_B",24,0.,24.,6,0,6);
     hist_trk_chi2_SV1_B = new TH1F("trk_chi2_SV1_B","trk_chi2_SV1_B",400,0,400);
     hist_trk_SV1input_origin_B = new TH1F("trk_SV1input_origin_B","trk_SV1input_origin_B",6,0,6);

     hist_trk_pT_SV0_B = new TH1F("trk_pT_SV0_B", "trk_pT_SV0_B", 3000, 0., 300.);
     hist_trk_eta_SV0_B = new TH1F("trk_SV0_eta_B","trk_SV0_eta_B",260,-2.6,2.6);
     hist_trk_pT_jet_DR_SV0_B = new TH2F("trk_pT_jet_DR_SV0_B","trk_pT_jet_DR_SV0_B",2000,0.,1000,200,0.,2.);
     hist_trk_origin_SV0_B = new TH1F("trk_SV0_origin_B","trk_SV0_origin_B",6,0,6);
     hist_trk_d0_SV0_B = new TH1F("trk_d0_SV0_B","trk_d0_SV0_B",3000,-15.,15.);
     hist_trk_z0sinth_SV0_B = new TH1F("trk_z0sinth_SV0_B","trk_z0sinth_SV0_B",1000,-5.,5.);
     hist_trk_d0sig_SV0_B = new TH1F("trk_d0sig_SV0_B","trk_d0sig_SV0_B",300,-15.,15.);
     hist_trk_z0sinthsig_SV0_B = new TH1F("trk_z0sinthsig_SV0_B","trk_z0sinthsig_SV0_B",300,-15.,15.);
     hist_trk_d0sig_origin_SV0_B = new TH2F("trk_d0sig_origin_SV0_B","trk_d0sig_origin_SV0_B",300,-15.,15.,6,0,6);
     hist_trk_z0sinthsig_origin_SV0_B = new TH2F("trk_z0sinthsig_origin_SV0_B","trk_z0sinthsig_origin_SV0_B",300,-15.,15.,6,0,6);
     hist_trk_logpTfrac_origin_SV0_B = new TH2F("trk_logpTfrac_origin_SV0_B","trk_logpTfrac_origin_SV0_B",300,-10.,2,6,0,6);
     hist_trk_logDR_origin_SV0_B = new TH2F("trk_logDR_origin_SV0_B","trk_logDR_origin_SV0_B",200,-10.,2.,6,0,6);
     hist_trk_IBLhits_origin_SV0_B = new TH2F("trk_IBLhits_origin_SV0_B","trk_IBLhits_origin_SV0_B",12,0.,12.,6,0,6);
     hist_trk_NextToIBLhits_origin_SV0_B = new TH2F("trk_NextToIBLhits_origin_SV0_B","trk_NextToIBLhits_origin_SV0_B",12,0.,12.,6,0,6);
     hist_trk_sharedIBLhits_origin_SV0_B = new TH2F("trk_sharedIBLhits_origin_SV0_B","trk_sharedIBLhits_origin_SV0_B",12,0.,12.,6,0,6);
     hist_trk_splitIBLhits_origin_SV0_B = new TH2F("trk_splitIBLhits_origin_SV0_B","trk_splitIBLhits_origin_SV0_B",12,0.,12.,6,0,6);
     hist_trk_nPixhits_origin_SV0_B = new TH2F("trk_nPixhits_origin_SV0_B","trk_nPixhits_origin_SV0_B",12,0.,12.,6,0,6);
     hist_trk_sharedPixhits_origin_SV0_B = new TH2F("trk_sharedPixhits_origin_SV0_B","trk_sharedPixhits_origin_SV0_B",12,0.,12.,6,0,6);
     hist_trk_splitPixhits_origin_SV0_B = new TH2F("trk_splitPixhits_origin_SV0_B","trk_splitPixhits_origin_SV0_B",12,0.,12.,6,0,6);
     hist_trk_nSCThits_origin_SV0_B = new TH2F("trk_nSCThits_origin_SV0_B","trk_nSCThits_origin_SV0_B",24,0.,24.,6,0,6);
     hist_trk_sharedSCThits_origin_SV0_B = new TH2F("trk_sharedSCThits_origin_SV0_B","trk_sharedSCThits_origin_SV0_B",24,0.,24.,6,0,6);

     hist_trk_pT_IP3D_B = new TH1F("trk_pT_IP3D_B", "trk_pT_IP3D_B", 3000, 0., 300.);
     hist_jet_pT_origin_IP3D_B = new TH2F("trk_jet_pT_origin_IP3D_B", "trk_jet_pT_origin_IP3D_B",6,0,6,500, 0., 1000.);
     hist_jet_pT_origin_truth_label_IP3D_B = new TH2F("trk_jet_pT_origin_truth_label_IP3D_B", "trk_jet_pT_origin_truth_label_IP3D_B",121,-1,120,500, 0., 1000.);
     hist_bH_pT_origin_truth_label_IP3D_B = new TH2F("trk_bH_pT_origin_truth_label_IP3D_B", "trk_bH_pT_origin_truth_label_IP3D_B",121,-1,120,500, 0., 1000.);
     hist_jet_pt_IP3D_B = new TH1F("trk_jet_pT_IP3D_B", "trk_jet_pT_IP3D_B", 500, 0., 1000.);
     hist_trk_eta_IP3D_B = new TH1F("trk_IP3D_eta_B","trk_IP3D_eta_B",260,-2.6,2.6);
     hist_trk_pT_jet_DR_IP3D_B = new TH2F("trk_pT_jet_DR_IP3D_B","trk_pT_jet_DR_IP3D_B",2000,0.,1000,200,0.,2.);
     hist_trk_origin_IP3D_B = new TH1F("trk_IP3D_origin_B","trk_IP3D_origin_B",6,0,6);
     hist_trk_origin_truth_label_IP3D_B = new TH1F("trk_origin_truth_label_IP3D_B","trk_origin_truth_label_IP3D_B",121,-1,120);
     hist_trk_d0_IP3D_B = new TH1F("trk_d0_IP3D_B","trk_d0_IP3D_B",3000,-15.,15.);
     hist_trk_z0sinth_IP3D_B = new TH1F("trk_z0sinth_IP3D_B","trk_z0sinth_IP3D_B",1000,-5.,5.);
     hist_trk_d0sig_IP3D_B = new TH1F("trk_d0sig_IP3D_B","trk_d0sig_IP3D_B",300,-15.,15.);
     hist_trk_z0sinthsig_IP3D_B = new TH1F("trk_z0sinthsig_IP3D_B","trk_z0sinthsig_IP3D_B",300,-15.,15.);
     hist_trk_d0sig_origin_IP3D_B = new TH2F("trk_d0sig_origin_IP3D_B","trk_d0sig_origin_IP3D_B",300,-15.,15.,6,0,6);
     hist_trk_z0sinthsig_origin_IP3D_B = new TH2F("trk_z0sinthsig_origin_IP3D_B","trk_z0sinthsig_origin_IP3D_B",300,-15.,15.,6,0,6);
     hist_trk_logpTfrac_origin_IP3D_B = new TH2F("trk_logpTfrac_origin_IP3D_B","trk_logpTfrac_origin_IP3D_B",300,-10.,2,6,0,6);
     hist_trk_logDR_origin_IP3D_B = new TH2F("trk_logDR_origin_IP3D_B","trk_logDR_origin_IP3D_B",200,-10.,2.,6,0,6);
     hist_trk_IBLhits_origin_IP3D_B = new TH2F("trk_IBLhits_origin_IP3D_B","trk_IBLhits_origin_IP3D_B",12,0.,12.,6,0,6);
     hist_trk_NextToIBLhits_origin_IP3D_B = new TH2F("trk_NextToIBLhits_origin_IP3D_B","trk_NextToIBLhits_origin_IP3D_B",12,0.,12.,6,0,6);
     hist_trk_sharedIBLhits_origin_IP3D_B = new TH2F("trk_sharedIBLhits_origin_IP3D_B","trk_sharedIBLhits_origin_IP3D_B",12,0.,12.,6,0,6);
     hist_trk_splitIBLhits_origin_IP3D_B = new TH2F("trk_splitIBLhits_origin_IP3D_B","trk_splitIBLhits_origin_IP3D_B",12,0.,12.,6,0,6);
     hist_trk_nPixhits_origin_IP3D_B = new TH2F("trk_nPixhits_origin_IP3D_B","trk_nPixhits_origin_IP3D_B",12,0.,12.,6,0,6);
     hist_trk_sharedPixhits_origin_IP3D_B = new TH2F("trk_sharedPixhits_origin_IP3D_B","trk_sharedPixhits_origin_IP3D_B",12,0.,12.,6,0,6);
     hist_trk_splitPixhits_origin_IP3D_B = new TH2F("trk_splitPixhits_origin_IP3D_B","trk_splitPixhits_origin_IP3D_B",12,0.,12.,6,0,6);
     hist_trk_nSCThits_origin_IP3D_B = new TH2F("trk_nSCThits_origin_IP3D_B","trk_nSCThits_origin_IP3D_B",24,0.,24.,6,0,6);
     hist_trk_sharedSCThits_origin_IP3D_B = new TH2F("trk_sharedSCThits_origin_IP3D_B","trk_sharedSCThits_origin_IP3D_B",24,0.,24.,6,0,6);

     hist_trk_pT_IP2D_B = new TH1F("trk_pT_IP2D_B", "trk_pT_IP2D_B", 3000, 0., 300.);
     hist_trk_eta_IP2D_B = new TH1F("trk_IP2D_eta_B","trk_IP2D_eta_B",260,-2.6,2.6);
     hist_trk_pT_jet_DR_IP2D_B = new TH2F("trk_pT_jet_DR_IP2D_B","trk_pT_jet_DR_IP2D_B",2000,0.,1000,200,0.,2.);
     hist_trk_origin_IP2D_B = new TH1F("trk_IP2D_origin_B","trk_IP2D_origin_B",6,0,6);
     hist_trk_d0_IP2D_B = new TH1F("trk_d0_IP2D_B","trk_d0_IP2D_B",3000,-15.,15.);
     hist_trk_z0sinth_IP2D_B = new TH1F("trk_z0sinth_IP2D_B","trk_z0sinth_IP2D_B",1000,-5.,5.);
     hist_trk_d0sig_IP2D_B = new TH1F("trk_d0sig_IP2D_B","trk_d0sig_IP2D_B",300,-15.,15.);
     hist_trk_z0sinthsig_IP2D_B = new TH1F("trk_z0sinthsig_IP2D_B","trk_z0sinthsig_IP2D_B",300,-15.,15.);
     hist_trk_d0sig_origin_IP2D_B = new TH2F("trk_d0sig_origin_IP2D_B","trk_d0sig_origin_IP2D_B",300,-15.,15.,6,0,6);
     hist_trk_z0sinthsig_origin_IP2D_B = new TH2F("trk_z0sinthsig_origin_IP2D_B","trk_z0sinthsig_origin_IP2D_B",300,-15.,15.,6,0,6);
     hist_trk_logpTfrac_origin_IP2D_B = new TH2F("trk_logpTfrac_origin_IP2D_B","trk_logpTfrac_origin_IP2D_B",300,-10.,2,6,0,6);
     hist_trk_logDR_origin_IP2D_B = new TH2F("trk_logDR_origin_IP2D_B","trk_logDR_origin_IP2D_B",200,-10.,2.,6,0,6);
     hist_trk_IBLhits_origin_IP2D_B = new TH2F("trk_IBLhits_origin_IP2D_B","trk_IBLhits_origin_IP2D_B",12,0.,12.,6,0,6);
     hist_trk_NextToIBLhits_origin_IP2D_B = new TH2F("trk_NextToIBLhits_origin_IP2D_B","trk_NextToIBLhits_origin_IP2D_B",12,0.,12.,6,0,6);
     hist_trk_sharedIBLhits_origin_IP2D_B = new TH2F("trk_sharedIBLhits_origin_IP2D_B","trk_sharedIBLhits_origin_IP2D_B",12,0.,12.,6,0,6);
     hist_trk_splitIBLhits_origin_IP2D_B = new TH2F("trk_splitIBLhits_origin_IP2D_B","trk_splitIBLhits_origin_IP2D_B",12,0.,12.,6,0,6);
     hist_trk_nPixhits_origin_IP2D_B = new TH2F("trk_nPixhits_origin_IP2D_B","trk_nPixhits_origin_IP2D_B",12,0.,12.,6,0,6);
     hist_trk_sharedPixhits_origin_IP2D_B = new TH2F("trk_sharedPixhits_origin_IP2D_B","trk_sharedPixhits_origin_IP2D_B",12,0.,12.,6,0,6);
     hist_trk_splitPixhits_origin_IP2D_B = new TH2F("trk_splitPixhits_origin_IP2D_B","trk_splitPixhits_origin_IP2D_B",12,0.,12.,6,0,6);
     hist_trk_nSCThits_origin_IP2D_B = new TH2F("trk_nSCThits_origin_IP2D_B","trk_nSCThits_origin_IP2D_B",24,0.,24.,6,0,6);
     hist_trk_sharedSCThits_origin_IP2D_B = new TH2F("trk_sharedSCThits_origin_IP2D_B","trk_sharedSCThits_origin_IP2D_B",24,0.,24.,6,0,6);

     hist_trk_d0_PUB = new TH1F("trk_d0_PUB","trk_d0_PUB",300,-15.,15.);
     hist_trk_d0_BB = new TH1F("trk_d0_BB","trk_d0_BB",300,-15.,15.);
     hist_trk_d0_CB = new TH1F("trk_d0_CB","trk_d0_CB",300,-15.,15.);
     hist_trk_d0_FRAGB = new TH1F("trk_d0_FRAGB","trk_d0_FRAGB",300,-15.,15.);
     hist_trk_d0_GEANTB = new TH1F("trk_d0_GEANTB","trk_d0_GEANTB",300,-15.,15.);

     hist_child_pT_B = new TH1F("child_pT_B", "child_pT_B", 500,0.,1000.);
     hist_child_eta_B = new TH1F("child_eta_B","child_eta_B",260,-2.6,2.6);
     hist_child_phi_B = new TH1F("child_phi_B","child_phi_B",200,-4.,4.);
     hist_child_Deta_B = new TH1F("child_Deta_B","child_Deta_B",200,-1.,1.);
     hist_child_Dphi_B = new TH1F("child_Dphi_B","child_Dphi_B",200,-1.,1.);
     hist_child_Dphi_Deta_B = new TH2F("child_Dphi_Deta","child_Dphi_Deta", 100, -1.,1., 100, -1.,1.);
     hist_child_DR_B = new TH1F("child_DR_B","child_DR_B",200,0.,2.);
     hist_child_pT_DR_B = new TH2F("child_pT_DR_B","child_pT_DR_B", 500,0.,1000.,200,0.,2.);
     hist_child_pT_jet_DR_B = new TH2F("child_pT_jet_DR_B","child_pT_jet_DR_B", 500,0.,1000.,200,0.,2.);
     hist_child_pdgID_B = new TH1F("child_pdgID_B","child_pdgID_B",200000,-100000,100000);

     hist_child_jetpT_B = new TH1F("child_jetpT_B", "child_jetpT_B",  500,0.,1000.);
     hist_child_DR_jetpt_B = new TH2F("child_DR_jetpt_B","child_DR_jetpt_B", 500,0.,1000.,200,0.,2.);
     hist_child_DR_bHpt_B = new TH2F("child_DR_bHpt_B","child_DR_bHpt_B", 500,0.,1000.,200,0.,2.);
     hist_child_bHpT_B = new TH1F("child_bHpT_B", "child_bHpT_B", 500,0.,1000.);

     hist_child_jetpT_C = new TH1F("child_jetpT_C", "child_jetpT_C",  500,0.,1000.);
     hist_child_DR_jetpt_C = new TH2F("child_DR_jetpt_C","child_DR_jetpt_C", 500,0.,1000.,200,0.,2.);
     hist_child_DR_cHpt_C = new TH2F("child_DR_cHpt_C","child_DR_cHpt_C", 500,0.,1000.,200,0.,2.);
     hist_child_cHpT_C = new TH1F("child_cHpT_C", "child_cHpT_C", 500,0.,1000.);



     hist_child_pi_notD = new TH1F("child_pi_pT_notD_B", "child_pi_pT_notD_B", 300, 0., 150.);
     hist_child_K_notD = new TH1F("child_K_pT_notD_B", "child_K_pT_notD_B", 300, 0., 150.);
     hist_child_pi = new TH1F("child_pi_pT_B", "child_pi_pT_B", 300, 0., 150.);
     hist_child_pi_Lxy_B = new TH1F("child_pi_Lxy_B","child_pi_Lxy_B",1000,0.,1000.);
     hist_child_pi_Lxyz_B = new TH1F("child_pi_Lxyz_B","child_pi_Lxyz_B",1000,0.,1000.);
     hist_child_pi_d0_truth_B = new TH1F("child_pi_d0_truth_B","child_pi_d0_truth_B",300,-15.,15.);
     hist_child_pi_z0_truth_B = new TH1F("child_pi_z0_truth_B","child_pi_z0_truth_B",1000,-1000.,1000.);
     hist_child_K = new TH1F("child_K_pT_B", "child_K_pT_B", 300, 0., 150.);
     hist_child_K_Lxy_B = new TH1F("child_K_Lxy_B","child_K_Lxy_B",1000,0.,1000.);
     hist_child_K_Lxyz_B = new TH1F("child_K_Lxyz_B","child_K_Lxyz_B",1000,0.,1000.);
     hist_child_K_d0_truth_B = new TH1F("child_K_d0_truth_B","child_K_d0_truth_B",300,-15.,15.);
     hist_child_K_z0_truth_B = new TH1F("child_K_z0_truth_B","child_K_z0_truth_B",1000,-1000.,1000.);
     hist_child_mu = new TH1F("child_mu_pT_B", "child_mu_pT_B", 300, 0., 150.);
     hist_child_mu_Lxy_B = new TH1F("child_mu_Lxy_B","child_mu_Lxy_B",1000,0.,1000.);
     hist_child_mu_Lxyz_B = new TH1F("child_mu_Lxyz_B","child_mu_Lxyz_B",1000,0.,1000.);
     hist_child_mu_d0_truth_B = new TH1F("child_mu_d0_truth_B","child_mu_d0_truth_B",300,-15.,15.);
     hist_child_mu_z0_truth_B = new TH1F("child_mu_z0_truth_B","child_mu_z0_truth_B",1000,-1000.,1000.);
     hist_child_p = new TH1F("child_p_pT_B", "child_p_pT_B", 300, 0., 150.);
     hist_child_p_Lxy_B = new TH1F("child_p_Lxy_B","child_p_Lxy_B",1000,0.,1000.);
     hist_child_p_Lxyz_B = new TH1F("child_p_Lxyz_B","child_p_Lxyz_B",1000,0.,1000.);
     hist_child_p_d0_truth_B = new TH1F("child_p_d0_truth_B","child_p_d0_truth_B",300,-15.,15.);
     hist_child_p_z0_truth_B = new TH1F("child_p_z0_truth_B","child_p_z0_truth_B",1000,-1000.,1000.);
     hist_child_e = new TH1F("child_e_pT_B", "child_e_pT_B", 300, 0., 150.);
     hist_child_e_Lxy_B = new TH1F("child_e_Lxy_B","child_e_Lxy_B",1000,0.,1000.);
     hist_child_e_Lxyz_B = new TH1F("child_e_Lxyz_B","child_e_Lxyz_B",1000,0.,1000.);
     hist_child_e_d0_truth_B = new TH1F("child_e_d0_truth_B","child_e_d0_truth_B",300,-15.,15.);
     hist_child_e_z0_truth_B = new TH1F("child_e_z0_truth_B","child_e_z0_truth_B",1000,-1000.,1000.);

     hist_child_Lxy_B = new TH1F("child_Lxy_B","child_Lxy_B",1000,0.,1000.);
     hist_child_Lxyz_B = new TH1F("child_Lxyz_B","child_Lxyz_B",1000,0.,1000.);
//     hist_child_decay_IP = new TH1F("child_decay_IP_B","child_decay_IP_B",300,-15.,15.);
//     hist_child_nodecay_IP = new TH1F("child_nodecay_IP_B","child_nodecay_IP_B",300,-15.,15.);
//     hist_child_linear_IP = new TH1F("child_linear_IP_B","child_linear_IP_B",300,-15.,15.);
     hist_child_d0_truth = new TH1F("child_d0_truth_B","child_d0_truth_B",3000,-15.,15.);
     hist_child_d0 = new TH1F("child_d0_B","child_d0_B",3000,-15.,15.);
     hist_child_z0 = new TH1F("child_z0_B","child_z0_B",3000,-15.,15.);
     hist_child_d0_pT = new TH2F("child_d0_pT_B","child_pT_d0_B",300,-15.,15.,300,0.,150.);
     hist_child_z0sinth_B = new TH1F("child_z0sinth_B","child_z0sinth_B",1000,-5.,5.);
     hist_pT_vs_R0_ratio_B = new TH1F("child_pT_vs_R0_ratio_B","child_pT_vs_R0_ratio_B",300,-1.2,2.4);

     if(origin_selection){

       hist_matched_origin_pT_B = new TH1F("matched_origin_trk_pT_B","matched_origin_trk_pT_B", 500, 0., 1000.);
       hist_matched_origin_eta_B = new TH1F("matched_origin_trk_eta_B","matched_origin_trk_eta_B",260,-2.6,2.6);
       hist_matched_origin_phi_B = new TH1F("matched_origin_trk_phi_B","matched_origin_trk_phi_B",200,-4.,4.);
       hist_matched_origin_Deta_B = new TH1F("matched_origin_trk_Deta_B","matched_origin_trk_Deta_B",200,-1.,1.);
       hist_matched_origin_Dphi_B = new TH1F("matched_origin_trk_Dphi_B","matched_origin_trk_Dphi_B",200,-1.,1.);
       hist_matched_origin_Dphi_Deta_B = new TH2F("matched_origin_trk_Dphi_Deta","matched_origin_trk_Dphi_Deta", 100, -1.,1., 100, -1.,1.);
       hist_matched_origin_DR_B = new TH1F("matched_origin_trk_DR_B","matched_origin_trk_DR_B",200,0.,2.);
       hist_matched_origin_pT_DR_B = new TH2F("matched_origin_trk_pT_DR_B","matched_origin_trk_pT_DR_B",300,0.,150.,200,0.,2.);
       hist_matched_origin_pT_jet_DR_B = new TH2F("matched_origin_trk_pT_jet_DR_B","matched_origin_trk_pT_jet_DR_B",3000, 0., 300.,200,0.,2.);
       hist_matched_origin_pdgId_B = new TH1F("matched_origin_trk_pdgId_B","matched_origin_trk_pdgId_B",200000,-100000,100000);
       hist_matched_origin_origin_B = new TH1F("matched_origin_trk_origin_B","matched_origin_trk_origin_B",7,-2,5);
       hist_matched_origin_d0_B = new TH1F("matched_origin_trk_d0_B","matched_origin_trk_d0_B",300,-15.,15.);
//       hist_matched_origin_Lxy_B = new TH1F("matched_origin_trk_Lxy_B","matched_origin_trk_Lxy_B",300,0.,1000.);
//       hist_matched_origin_Lxyz_B = new TH1F("matched_origin_trk_Lxyz_B","matched_origin_trk_Lxyz_B",300,0.,1000.);
       hist_matched_origin_jetpT_B = new TH1F("matched_origin_trk_jetpT_B","matched_origin_trk_jetpT_B", 500, 0., 1000.);
       hist_trk_DR_jetpt_B = new TH2F("matched_origin_trk_DR_jetpt_B","matched_origin_trk_DR_jetpt_B", 500, 0., 1000.,200,0.,2.);
       hist_trk_DR_bHpt_B = new TH2F("matched_origin_trk_DR_bHpt_B","matched_origin_trk_DR_bHpt_B", 500, 0., 1000.,200,0.,2.);
       hist_matched_origin_bHpT_B = new TH1F("matched_origin_trk_bHpT_B","matched_origin_trk_bHpT_B", 500, 0., 1000.);

       hist_matched_origin_jetpT_C = new TH1F("matched_origin_trk_jetpT_C","matched_origin_trk_jetpT_C", 500, 0., 1000.);
       hist_trk_DR_jetpt_C = new TH2F("matched_origin_trk_DR_jetpt_C","matched_origin_trk_DR_jetpt_C", 500, 0., 1000.,200,0.,2.);
       hist_trk_DR_cHpt_C = new TH2F("matched_origin_trk_DR_cHpt_C","matched_origin_trk_DR_cHpt_C", 500, 0., 1000.,200,0.,2.);
       hist_matched_origin_cHpT_C = new TH1F("matched_origin_trk_cHpT_C","matched_origin_trk_cHpT_C", 500, 0., 1000.);
     }

     if(geometric_selection){

       hist_matched_pT_B = new TH1F("matched_child_pT_B","matched_child_pT_B", 3000, 0., 300.);
       hist_matched_eta_B = new TH1F("matched_child_eta_B","matched_child_eta_B",260,-2.6,2.6);
       hist_matched_phi_B = new TH1F("matched_child_phi_B","matched_child_phi_B",200,-4.,4.);
       hist_matched_Deta_B = new TH1F("matched_child_Deta_B","matched_child_Deta_B",200,-1.,1.);
       hist_matched_Dphi_B = new TH1F("matched_child_Dphi_B","matched_child_Dphi_B",200,-1.,1.);
       hist_matched_Dphi_Deta_B = new TH2F("matched_child_Dphi_Deta","matched_child_Dphi_Deta", 100, -1.,1., 100, -1.,1.);
       hist_matched_DR_B = new TH1F("matched_child_DR_B","matched_child_DR_B",200,0.,2.);
       hist_matched_pT_DR_B = new TH2F("matched_child_pT_DR_B","matched_child_pT_DR_B",300,0.,150.,200,0.,2.);
       hist_matched_pT_jet_DR_B = new TH2F("matched_child_pT_jet_DR_B","matched_child_pT_jet_DR_B",3000, 0., 300.,200,0.,2.);
       hist_matched_pdgId_B = new TH1F("matched_child_pdgId_B","matched_child_pdgId_B",200000,-100000,100000);

       hist_matched_child_pi_notD = new TH1F("matched_child_pi_pT_notD_B", "matched_child_pi_pT_notD_B", 300, 0., 150.);
       hist_matched_child_K_notD = new TH1F("matched_child_K_pT_notD_B", "matched_child_K_pT_notD_B", 300, 0., 150.);
       hist_matched_child_pi = new TH1F("matched_child_pi_pT_B", "matched_child_pi_pT_B", 300, 0., 150.);
       hist_matched_child_pi_Lxy_B = new TH1F("matched_child_pi_Lxy_B","matched_child_pi_Lxy_B",1000,0.,1000.);
       hist_matched_child_pi_Lxyz_B = new TH1F("matched_child_pi_Lxyz_B","matched_child_pi_Lxyz_B",1000,0.,1000.);
       hist_matched_child_pi_d0_truth_B = new TH1F("matched_child_pi_d0_truth_B","matched_child_pi_d0_truth_B",300,-15.,15.);
       hist_matched_child_pi_z0_truth_B = new TH1F("matched_child_pi_z0_truth_B","matched_child_pi_z0_truth_B",1000,-1000.,1000.);
       hist_matched_child_K = new TH1F("matched_child_K_pT_B", "matched_child_K_pT_B", 300, 0., 150.);
       hist_matched_child_K_Lxy_B = new TH1F("matched_child_K_Lxy_B","matched_child_K_Lxy_B",1000,0.,1000.);
       hist_matched_child_K_Lxyz_B = new TH1F("matched_child_K_Lxyz_B","matched_child_K_Lxyz_B",1000,0.,1000.);
       hist_matched_child_K_d0_truth_B = new TH1F("matched_child_K_d0_truth_B","matched_child_K_d0_truth_B",300,-15.,15.);
       hist_matched_child_K_z0_truth_B = new TH1F("matched_child_K_z0_truth_B","matched_child_K_z0_truth_B",1000,-1000.,1000.);
       hist_matched_child_mu = new TH1F("matched_child_mu_pT_B", "matched_child_mu_pT_B", 300, 0., 150.);
       hist_matched_child_mu_Lxy_B = new TH1F("matched_child_mu_Lxy_B","matched_child_mu_Lxy_B",1000,0.,1000.);
       hist_matched_child_mu_Lxyz_B = new TH1F("matched_child_mu_Lxyz_B","matched_child_mu_Lxyz_B",1000,0.,1000.);
       hist_matched_child_mu_d0_truth_B = new TH1F("matched_child_mu_d0_truth_B","matched_child_mu_d0_truth_B",300,-15.,15.);
       hist_matched_child_mu_z0_truth_B = new TH1F("matched_child_mu_z0_truth_B","matched_child_mu_z0_truth_B",1000,-1000.,1000.);
       hist_matched_child_p = new TH1F("matched_child_p_pT_B", "matched_child_p_pT_B", 300, 0., 150.);
       hist_matched_child_p_Lxy_B = new TH1F("matched_child_p_Lxy_B","matched_child_p_Lxy_B",1000,0.,1000.);
       hist_matched_child_p_Lxyz_B = new TH1F("matched_child_p_Lxyz_B","matched_child_p_Lxyz_B",1000,0.,1000.);
       hist_matched_child_p_d0_truth_B = new TH1F("matched_child_p_d0_truth_B","matched_child_p_d0_truth_B",300,-15.,15.);
       hist_matched_child_p_z0_truth_B = new TH1F("matched_child_p_z0_truth_B","matched_child_p_z0_truth_B",1000,-1000.,1000.);
       hist_matched_child_e = new TH1F("matched_child_e_pT_B", "matched_child_e_pT_B", 300, 0., 150.);
       hist_matched_child_e_Lxy_B = new TH1F("matched_child_e_Lxy_B","matched_child_e_Lxy_B",1000,0.,1000.);
       hist_matched_child_e_Lxyz_B = new TH1F("matched_child_e_Lxyz_B","matched_child_e_Lxyz_B",1000,0.,1000.);
       hist_matched_child_e_d0_truth_B = new TH1F("matched_child_e_d0_truth_B","matched_child_e_d0_truth_B",300,-15.,15.);
       hist_matched_child_e_z0_truth_B = new TH1F("matched_child_e_z0_truth_B","matched_child_e_z0_truth_B",1000,-1000.,1000.);


       hist_matched_origin_B = new TH1F("matched_trk_origin_B","matched_trk_origin_B",7,-2,5);
       hist_matched_pT_child_pTfraction_B = new TH2F("matched_pT_child_pTfraction_B","matched_pT_child_vs_DpT/pT_child_B",300,0.,150.,500,-1.,10.);
       hist_matched_DR_trk_B = new TH1F("matched_DR_trk_B","matched_DR_trk_B_B",500,0.,1.);
       hist_matched_DR_trk_pTfraction = new TH2F("matched_DR_trk_pTfraction_B","matched_DR_trk_pTfraction_B",500,0.,1.,500,-1.,10.);
       hist_matched_d0_B = new TH1F("matched_child_d0_B","matched_child_d0_B",3000,-15.,15.);
       hist_matched_Lxy_B = new TH1F("matched_child_Lxy_B","matched_child_Lxy_B",1000,0.,1000.);
       hist_matched_Lxyz_B = new TH1F("matched_child_Lxyz_B","matched_child_Lxyz_B",1000,0.,1000.);
  /*
       hist_nomatched_pT_B = new TH1F("nomatched_child_pT_B", "nomatched_child_pT_B", 500, 0., 150.);
       hist_nomatched_eta_B = new TH1F("nomatched_child_eta_B","nomatched_child_eta_B",500,-2.6,2.6);
       hist_nomatched_phi_B = new TH1F("nomatched_child_phi_B","nomatched_child_phi_B",500,-4.,4.);
       hist_nomatched_Deta_B = new TH1F("nomatched_child_Deta_B","nomatched_child_Deta_B",500,-1.,1.);
       hist_nomatched_Dphi_B = new TH1F("nomatched_child_Dphi_B","nomatched_child_Dphi_B",500,-1.,1.);
       hist_nomatched_Dphi_Deta_B = new TH2F("nomatched_child_Dphi_Deta","nomatched_child_Dphi_Deta", 100, -1.,1., 100, -1.,1.);
       hist_nomatched_DR_B = new TH1F("nomatched_child_DR_B","nomatched_child_DR_B",500,-0.1,2.);
       hist_nomatched_pT_DR_B = new TH2F("nomatched_child_pT_DR_B","nomatched_child_pT_DR_B",80,0.,150.,80,0.,0.6);
       */
       hist_nomatched_pT_jet_DR_B = new TH2F("nomatched_child_pT_jet_DR_B","nomatched_child_pT_jet_DR_B",2000,0.,1000,200,0.,2.);
       hist_nomatchedIN_pT_B = new TH1F("nomatchedIN_child_pT_B", "nomatchedIN_child_pT_B", 500, 0., 150.);
       hist_nomatchedIN_eta_B = new TH1F("nomatchedIN_child_eta_B","nomatchedIN_child_eta_B",500,-2.6,2.6);
       hist_nomatchedIN_phi_B = new TH1F("nomatchedIN_child_phi_B","nomatchedIN_child_phi_B",500,-4.,4.);
       hist_nomatchedIN_DR_B = new TH1F("nomatchedIN_child_DR_B","nomatchedIN_child_DR_B",500,-0.1,2.);
       hist_nomatchedIN_pT_jet_DR_B = new TH2F("nomatchedIN_child_pT_jet_DR_B","nomatchedIN_child_pT_jet_DR_B",3000, 0., 300.,200,0.,2.);
       hist_nomatchedIN_d0_B = new TH1F("nomatchedIN_child_d0_truth_B","nomatchedIN_child_d0_truth_B",3000,-15.,15.);
       hist_nomatchedIN_z0sinth_B = new TH1F("nomatchedIN_child_z0sinth_B","nomatchedIN_child_z0sinth_B",3000,-15.,15.);
       hist_nomatchedOUT_pT_B = new TH1F("nomatchedOUT_child_pT_B", "nomatchedOUT_child_pT_B", 500, 0., 150.);
       hist_nomatchedOUT_eta_B = new TH1F("nomatchedOUT_child_eta_B","nomatchedOUT_child_eta_B",500,-2.6,2.6);
       hist_nomatchedOUT_phi_B = new TH1F("nomatchedOUT_child_phi_B","nomatchedOUT_child_phi_B",500,-4.,4.);
       hist_nomatchedOUT_DR_B = new TH1F("nomatchedOUT_child_DR_B","nomatchedOUT_child_DR_B",500,-0.1,2.);
       hist_nomatchedOUT_pT_jet_DR_B = new TH2F("nomatchedOUT_child_pT_jet_DR_B","nomatchedOUT_child_pT_jet_DR_B",2000,0.,1000,200,0.,2.);
       hist_nomatchedOUT_d0_B = new TH1F("nomatchedOUT_child_d0_truth_B","nomatchedOUT_child_d0_truth_B",3000,-15.,15.);
       hist_nomatchedOUT_z0sinth_B = new TH1F("nomatchedOUT_child_z0sinth_B","nomatchedOUT_child_z0sinth_B",3000,-15.,15.);
//       hist_nomatched_pdgId_B = new TH1F("nomatched_child_pdgId_B","nomatched_child_pdgId_B",2000,-1000,1000);

  /*
       hist_single_matched_pT_B = new TH1F("single_matched_pT_B","single_matched_pT_B", 500, 0., 150.);
       hist_single_matched_eta_B = new TH1F("single_matched_child_eta_B","single_matched_child_eta_B",500,-2.6,2.6);
       hist_single_matched_phi_B = new TH1F("single_matched_child_phi_B","single_matched_child_phi_B",500,-4.,4.);
       hist_single_matched_Deta_B = new TH1F("single_matched_child_Deta_B","single_matched_child_Deta_B",500,-1.,1.);
       hist_single_matched_Dphi_B = new TH1F("single_matched_child_Dphi_B","single_matched_child_Dphi_B",500,-1.,1.);
       hist_single_matched_Dphi_Deta_B = new TH2F("single_matched_child_Dphi_Deta","single_matched_child_Dphi_Deta", 100, -1.,1., 100, -1.1,1.);
       hist_single_matched_DR_B = new TH1F("single_matched_DR_B","single_matched_DR_B",500,-0.1,2.);
       hist_single_matched_pT_DR_B = new TH2F("single_matched_pT_DR_B","single_matched_pT_DR_B",80,0.,150.,80,0.,0.6);
       hist_single_matched_pT_jet_DR_B = new TH2F("single_matched_child_pT_jet_DR_B","single_matched_child_pT_jet_DR_B",100,0.,500.,100,0.,2.);
       hist_single_matched_pdgId_B = new TH1F("single_matched_child_pdgId_B","single_matched_child_pdgId_B",2000,-1000,1000);
       hist_single_matched_origin_B = new TH1F("single_matched_trk_origin_B","single_matched_trk_origin_B",7,-2,5);
       hist_single_matched_pT_child_pTfraction_B = new TH2F("single_matched_pT_child_pTfraction_B","single_matched_pT_child_vs_DpT/pT_child_B",500,0.,150.,500,-1.1,10.);
       hist_single_matched_DR_trk_B = new TH1F("single_matched_DR_trk_B","single_matched_DR_trk_B", 500, 0., 1.);
       hist_single_matched_DR_trk_pTfraction = new TH2F("single_matched_DR_trk_pTfraction_B","single_matched_DR_trk_pTfraction_B",500,0.,1.,500,-1.1,10.);
       hist_single_matched_d0_B = new TH1F("single_matched_child_d0_B","single_matched_child_d0_B",300,-15.,15.);
  */
     }
   }

}
void DAOD_selector::initFlagsAndCuts()
{
  std::cout<<"\n In DAOD_selector::initFlagsAndCuts"<<std::endl;


  m_cut=1.;  // the working point !!! maybe
  m_fc=0.08,m_fcRNNIP=0.07; // fraction of c jets in the non-b jet sample
  m_pTfraction_cut=1.,m_DRcut=0.1,m_pTfraction_nocut=1e6,m_DRnocut=1e6; // cuts for geometrical matching
  m_track_cut=10;


  //D_phi=0.,D_eta=0.,DR=0.,px=0.,py=0.,pz=0.,Dx_1=0.,Dy_1=0.,Dz_1=0.,Dx_2=0.,Dy_2=0.,Dz_2=0,Lxy=0,Lxyz=0,Dxy_1=0,x0=0,y0=0,Dx_3=0.,Dy_3=0.,Dxy_3=0.,rand_n=0.,R0=0,d0=0;
  c=2.99792458e8;//,nx=0,ny=0;


  m_pt_max_shrCone=500.;
  m_pt_min_shrCone=0.;
  m_Delta_pt_shrCone=(m_pt_max_shrCone-m_pt_min_shrCone)/bin_1;



  return;
}
void DAOD_selector::initCounters()
{
  std::cout<<"\n In DAOD_selector::initCounters"<<std::endl;

  m_Ntot=0; // number of events read/processed
  cnt_13=0,cnt_103=0,cnt_113=0,cnt_15=0,cnt_105=0,cnt_115=0;
  m_noB=0,m_bb=0,m_b=0,m_bc_overlap=0,m_sc=0,m_sc2=0,m_sc3=0,m_match=0,m_nomatch=0,m_match_overlap=0,m_match_notoverlap=0,n_trk_pT_cut=0,n_trk_PU_pT_cut=0,n_trk_FRAG_pT_cut=0,n_trk_GEANT_pT_cut=0,n_trk_B=0,n_trk_C=0;

  m_qc=0,m_qj=0,q=0,a=0,b=0,sc=0,sgn=0,sgn_d0=0,sgn_z0=0;
  m_den=0; // some denominator

  // Number of jets, jets with at least 1 B hadron, al least 1 C hadron, otherwise light (in general and after the jet selection [variables with _2] )
  m_njets=0,m_njets_2=0;
  m_nBjets=0,m_nCjets=0,m_nljets=0;
  m_nBjets_2=0,m_nCjets_2=0,m_nljets_2=0;
  // Number of jets with c and b hadrons
  m_nJetBCoverlap = 0;   // ov_1=0,ov_2=0;
  m_nJetBCoverlap_postJetSel = 0;
  m_nBcheck=0,m_nCcheck=0,m_nlcheck=0,ov_check=0;

  // count tracks that have a mismatch between geometry-based matching and origin based matching to the jet
  m_GeomNOr_PU=0,m_GeomNOr_F=0,m_GeomNOr_G=0;

  //cut flow on jets
     m_njets_2_passPtMin=0;
     m_njets_2_passPtMax=0;
     m_njets_2_passEtaRange=0;
     m_njets_2_passBadMedium=0;
     m_njets_2_passOR=0;
     m_njets_2_passORmu=0;
     m_njets_2_passJVT=0;
     m_njets_2_passIsol=0;

  return;
}

void DAOD_selector::OverlapRemoval(std::vector<int>& isJet, std::vector<int>& isJet_OR)
{
//  isJet_OR=isJet;
  double D_eta = 0.,D_phi = 0.;
  double DeltaR=0;
  int idx_l=0,idx_k=0.;
  unsigned n_good=0;
    for(unsigned l=0;l<isJet.size();l++){
      n_good=0;
      idx_l=isJet[l];
      for(unsigned k=l+1;k<isJet.size();k++){
        idx_k=isJet[k];
        D_eta=jet_eta[idx_l]-jet_eta[idx_k];
        if(abs(jet_phi[idx_l]-jet_phi[idx_k])>M_PI){
          D_phi=2*M_PI-abs(jet_phi[idx_l]-jet_phi[idx_k]);
        }
        else{
          D_phi=jet_phi[idx_l]-jet_phi[idx_k];
        }
        DeltaR=sqrt(D_eta*D_eta+D_phi*D_phi);
        if(DeltaR>jet_Isol_cut){
//          std::cout<<"Overlap\t"<<isJet.size()<<" "<<isJet_OR.size()<<" "<<l<<"\n";
//          isJet_OR.erase(isJet_OR.begin()+l);
//          continue;
//          isJet_OR.push_back(l);
          n_good=n_good+1;
        }
      }
      if(n_good==isJet.size()-l-1){
        isJet_OR.push_back(l);
      }
    }
//    std::cout<<isJet.size()<<"\t"<<isJet_OR.size()<<"\n";
}

void DAOD_selector::getGhostExtJetFlavourLabel(std::vector<int>& isJet, std::vector<int>& isBcheck, std::vector<int>& isCcheck, std::vector<int>& islcheck, bool diag_trms){
  isBcheck.clear();
  isCcheck.clear();
  islcheck.clear();

  int ic=-1;
  int ig=-1;
  int jetIndex=0;
  for(std::vector<int>::iterator it = isJet.begin(); it != isJet.end(); ++it){
    jetIndex = *it;

    ig = getGhosFlavLabel(jetIndex);
    if (diag_trms)
      {
	ic= getConeFlavLabel(jetIndex);
	if      (ig==0 && ic==0) isBcheck.push_back(*it);
	else if (ig==1 && ic==1) isCcheck.push_back(*it);
	else if (ig==2 && ic==2) islcheck.push_back(*it);
      }
    else
      {
	if      (ig==0) isBcheck.push_back(*it);
	else if (ig==1) isCcheck.push_back(*it);
	else if (ig==2) islcheck.push_back(*it);
      }
  }
  
  /*
  bool jet_labelled = false;
  if(diag_trms==false){
    for(std::vector<int>::iterator it = isJet.begin(); it != isJet.end(); ++it){
      jet_labelled=false;

      if(jet_ghostBHadCount[*it]==1){
        isBcheck.push_back(*it);
        jet_labelled=true;
      }
      if(jet_ghostCHadCount[*it]==1 && jet_ghostBHadCount[*it]==0){
        isCcheck.push_back(*it);
        jet_labelled=true;
      }
      if(jet_ghostCHadCount[*it]==0 && jet_ghostBHadCount[*it]==0){
        islcheck.push_back(*it);
        jet_labelled=true;
      }
      if(jet_labelled)  continue;
    }
  }
  if(diag_trms==true){
    for(std::vector<int>::iterator it = isJet.begin(); it != isJet.end(); ++it){
      jet_labelled=false;

      if(jet_ghostBHadCount[*it]==1 && jet_DoubleHadLabel[*it]==5){
        isBcheck.push_back(*it);
        jet_labelled=true;
      }
      if(jet_ghostCHadCount[*it]==1 && jet_ghostBHadCount[*it]==0 && jet_DoubleHadLabel[*it]==4){
        isCcheck.push_back(*it);
        jet_labelled=true;
      }
      //SS { 
      //if(jet_ghostCHadCount[*it]==0 && jet_DoubleHadLabel[*it]==0){
      if(jet_ghostBHadCount[*it]==0 && jet_ghostCHadCount[*it]==0 && (jet_DoubleHadLabel[*it]==0 || jet_DoubleHadLabel[*it]==15)){
	//SS }
        islcheck.push_back(*it);
        jet_labelled=true;
      }
      // SS the line below is irrelevant 
      if(jet_labelled)  continue;
    }
  }
  */

}
void DAOD_selector::getGhostJetFlavourLabel(std::vector<int>& isJet, std::vector<int>& isBcheck, std::vector<int>& isCcheck, std::vector<int>& islcheck, bool diag_trms){
  isBcheck.clear();
  isCcheck.clear();
  islcheck.clear();


  int ic=-1;
  int ig=-1;
  int jetIndex=0;
  for(std::vector<int>::iterator it = isJet.begin(); it != isJet.end(); ++it){
    jetIndex = *it;

    ig = getGhosFlavLabel(jetIndex);
    if (diag_trms)
      {
	ic= getConeFlavLabel(jetIndex);
	if      ((ig==0 || ig==3 || ig==5) && (ic==0 || ic==3 || ic==5)) isBcheck.push_back(*it);
	else if ((ig==1 || ig==4         ) && (ic==1 || ic==4         )) isCcheck.push_back(*it);
	else if ( ig==2                    &&  ic==2                   ) islcheck.push_back(*it);
      }
    else
      {
	if      (ig==0 || ig==3 || ig==5) isBcheck.push_back(*it);
	else if (ig==1 || ig==4         ) isCcheck.push_back(*it);
	else if (ig==2                  ) islcheck.push_back(*it);
      }
  }
  /*
  bool jet_labelled = false;
  if(diag_trms==false){
    for(std::vector<int>::iterator it = isJet.begin(); it != isJet.end(); ++it){
      jet_labelled=false;

      if(jet_ghostBHadCount[*it]>0){
        isBcheck.push_back(*it);
        jet_labelled=true;
      }
      if(jet_ghostCHadCount[*it]>0 && jet_ghostBHadCount[*it]==0){
        isCcheck.push_back(*it);
        jet_labelled=true;
      }
      if(jet_ghostCHadCount[*it]==0 && jet_ghostBHadCount[*it]==0){
        islcheck.push_back(*it);
        jet_labelled=true;
      }
      if(jet_labelled)  continue;
    }
  }
  if(diag_trms==true){
    for(std::vector<int>::iterator it = isJet.begin(); it != isJet.end(); ++it){
      jet_labelled=false;

      if(jet_ghostBHadCount[*it]>0 && jet_LabDr_HadF[*it]==5){
        isBcheck.push_back(*it);
        jet_labelled=true;
      }
      if(jet_ghostCHadCount[*it]>0 && jet_ghostBHadCount[*it]==0 && jet_LabDr_HadF[*it]==4){
        isCcheck.push_back(*it);
        jet_labelled=true;
      }
      //SS { 
      //if(jet_ghostCHadCount[*it]==0 && jet_DoubleHadLabel[*it]==0){
      if(jet_ghostBHadCount[*it]==0 && jet_ghostCHadCount[*it]==0 && jet_LabDr_HadF[*it]==0){
	//SS }
        islcheck.push_back(*it);
        jet_labelled=true;
      }
      // SS the line below is irrelevant 
      if(jet_labelled)  continue;
    }
  }
  */

}

void DAOD_selector::getHadronConeFlavourLabel(std::vector<int>& isJet, std::vector<int>& isBcheck, std::vector<int>& isCcheck, std::vector<int>& islcheck)
{
  isBcheck.clear();
  isCcheck.clear();
  islcheck.clear();


  int ic=-1;
  int ig=-1;
  int jetIndex=0;
  for(std::vector<int>::iterator it = isJet.begin(); it != isJet.end(); ++it){
    jetIndex = *it;

    ic= getConeFlavLabel(jetIndex);
    if      (ic==0 || ic==3 || ic==5) isBcheck.push_back(*it);
    else if (ic==1 || ic==4         ) isCcheck.push_back(*it);
    else if (ic==2                  ) islcheck.push_back(*it);
  }

  /*
  bool jet_labelled = false;

  for(std::vector<int>::iterator it = isJet.begin(); it != isJet.end(); ++it){
    jet_labelled=false;

    if(jet_LabDr_HadF[*it]==5){
      isBcheck.push_back(*it);
      jet_labelled=true;
    }
    if(jet_LabDr_HadF[*it]==4){
      isCcheck.push_back(*it);
      jet_labelled=true;
    }
    if(jet_LabDr_HadF[*it]==0 || jet_LabDr_HadF[*it]==15){
      islcheck.push_back(*it);
      jet_labelled=true;
    }
    if(jet_labelled)  continue;

  }
  */
}
void DAOD_selector::getHadronConeExtFlavourLabel(std::vector<int>& isJet, std::vector<int>& isBcheck, std::vector<int>& isCcheck, std::vector<int>& islcheck)
{
  isBcheck.clear();
  isCcheck.clear();
  islcheck.clear();

  int ic=-1;
  int ig=-1;
  int jetIndex=0;
  for(std::vector<int>::iterator it = isJet.begin(); it != isJet.end(); ++it){
    jetIndex = *it;

    ic= getConeFlavLabel(jetIndex);
    if      (ic==0) isBcheck.push_back(*it);
    else if (ic==1) isCcheck.push_back(*it);
    else if (ic==2) islcheck.push_back(*it);
  }


  /*
  bool jet_labelled = false;

  for(std::vector<int>::iterator it = isJet.begin(); it != isJet.end(); ++it){
    jet_labelled=false;

    if(jet_DoubleHadLabel[*it]==5){
      isBcheck.push_back(*it);
      jet_labelled=true;
    }
    if(jet_DoubleHadLabel[*it]==4){
      isCcheck.push_back(*it);
      jet_labelled=true;
    }
    if(jet_DoubleHadLabel[*it]==0 || jet_DoubleHadLabel[*it]==15){
      islcheck.push_back(*it);
      jet_labelled=true;
    }
    if(jet_labelled)  continue;

  }
  */
}
   // jet labelling (b-, c-, l-) based on truth (with pt, Deta cuts), Clusive samples
void DAOD_selector::getTrueJetFlavourLabel(std::vector<int>& isJet, std::vector<int>& isBcheck, std::vector<int>& isCcheck, std::vector<int>& islcheck)
{

   isBcheck.clear();
   isCcheck.clear();
   islcheck.clear();

   bool jet_labelled = false;
   int n=0;
   double D_eta = 0.,D_phi = 0.;
   double DeltaR=0.,pt=0.;
   for(std::vector<int>::iterator it = isJet.begin(); it != isJet.end(); ++it){
     jet_labelled=false;

     if(jet_nBHadr[*it]>0){
       if(!decay_mode.compare("leptonic") || !decay_mode.compare("hadronic")){

         for(unsigned k=0;k<jet_bH_child_E[*it].size();k++){
           if(abs(jet_bH_child_pdg_id[*it].at(k))==11 || abs(jet_bH_child_pdg_id[*it].at(k))==13)
            n++;
         }

/*
         for(unsigned k=0;k<jet_trk_pt[*it].size();k++){
           if(jet_trk_orig[*it].at(k)==0 || jet_trk_orig[*it].at(k)==1){
             if(abs(jet_trk_pdg_id[*it].at(k))==11 || abs(jet_trk_pdg_id[*it].at(k))==13)
               n++;
           }
         }
*/
         if(!decay_mode.compare("leptonic"))
           if(n==0)  continue;//with n==0 we select only "leptonic" b-jets (b hadrons decaying leptonically)
         if(!decay_mode.compare("hadronic"))
           if(n>0)   continue;//with n>0 we select only "hadronic" b-jets (b hadrons decaying only hadronically)
       }

       for(unsigned i=0;i<jet_bH_pt[*it].size();i++){
         D_eta=jet_bH_eta[*it].at(i)-jet_eta[*it];
         if(abs(jet_bH_phi[*it].at(i)-jet_phi[*it])>M_PI){
           D_phi=2*M_PI-abs(jet_bH_phi[*it].at(i)-jet_phi[*it]);
         }
         else{
           D_phi=jet_bH_phi[*it].at(i)-jet_phi[*it];
         }

         DeltaR=sqrt(D_eta*D_eta+D_phi*D_phi);
         pt=jet_bH_pt[*it].at(i);
      	 if(DeltaR < m_DR_bcH_truth_cut && pt > m_pT_bcH_truth_cut){
      	     jet_labelled = true;
      	     isBcheck.push_back(*it);//// vector of jets with >=1 b hadron and labelled as b-jets
      	     i=jet_bH_pt[*it].size(); /// ==>>> no need to look for another b-hadron matching the "labelling requirements"
    	   }
       } // end loop over b-hadrons
     }
     if (jet_labelled) continue;
     //if (isBcheck.size() > 0) continue; /// ==>>> this jet is labelled as b-jet; it cannot be labelled as C/l anymore

     if(jet_nCHadr[*it]>0){
       for(unsigned i=0;i<jet_cH_pt[*it].size();i++){
         D_eta=jet_cH_eta[*it].at(i)-jet_eta[*it];
         if(abs(jet_cH_phi[*it].at(i)-jet_phi[*it])>M_PI){
           D_phi=2*M_PI-abs(jet_cH_phi[*it].at(i)-jet_phi[*it]);
         }
         else {
           D_phi=jet_cH_phi[*it].at(i)-jet_phi[*it];
         }

         DeltaR=sqrt(D_eta*D_eta+D_phi*D_phi);
         pt=jet_cH_pt[*it].at(i);
      	 if(DeltaR < m_DR_bcH_truth_cut && pt > m_pT_bcH_truth_cut){
      	     jet_labelled = true;
      	     isCcheck.push_back(*it);  //// vector of jets with non b-hadrons and >=1 c hadron-->  labelled as c-jet
      	     i=jet_cH_pt[*it].size(); /// ==>>> no need to look for another c-hadron matching the "labelling requirements"
      	  }
       } // end loop over c-hadrons
     }
     if (jet_labelled) continue;
     //if (isCcheck.size() > 0) continue; /// ==>>> this jet is labelled as c-jet; it cannot be labelled as l anymore

     // if we get here the jet failed b-labelling and c-labellig => it's a l-jet, according to truth labelling
     islcheck.push_back(*it);
   }


   return;
}
void DAOD_selector::getFlavorLabelMatrix(std::vector<int>& isJet)
{
  for (unsigned int i1=0; i1<6; ++i1)
    for (unsigned int i2=0; i2<6; ++i2)
      m_jetFlavorMatrix[i1][i2].clear();


  bool jet_labelled = false;

  int iflc = -1;
  int iflg = -1;
  for(std::vector<int>::iterator it = isJet.begin(); it != isJet.end(); ++it){
    jet_labelled=false;

    iflc = getConeFlavLabel(*it);
    iflg = getGhosFlavLabel(*it);
    if (fDoFlavorLabelMatrix) hist2_jetFlavorMatrix->Fill(double(iflc), double(iflg));
    //std::cout<<" iflc, g " << iflc<<" "<<iflg<<" double hadr "<<jet_DoubleHadLabel[*it]<<" ... nB/C hdr "<< jet_ghostBHadCount[*it]<<" "<<jet_ghostCHadCount[*it]<<std::endl;
    if (iflc > -1 && iflg > -1)
      {
	m_jetFlavorMatrix[iflc][iflg].push_back(*it);
      }
    else
      {
	std::cout<<" Flavor Label not assigned for this jet,  DoubleHadLabel, nGhost B/C hadrons   "
	  << jet_DoubleHadLabel[*it] <<" "<<jet_ghostBHadCount[*it]<<"/"<<jet_ghostCHadCount[*it]<<std::endl;
      }
  }
  
  return;
}
int DAOD_selector::getConeFlavLabel(int jetIndex)
{
  int iflav = -1;
  // if (jet_DoubleHadLabel[jetIndex]==5) return 0;
  // if (jet_DoubleHadLabel[jetIndex]==4) return 1;
  // if (jet_DoubleHadLabel[jetIndex]==0 || jet_DoubleHadLabel[jetIndex]==15) return 2;
  // if (jet_DoubleHadLabel[jetIndex]==55) return 3;
  // if (jet_DoubleHadLabel[jetIndex]==44) return 4;
  // if (jet_DoubleHadLabel[jetIndex]==54) return 5;

  if (jet_DoubleHadLabel[jetIndex]==5) return 0;
  if (jet_DoubleHadLabel[jetIndex]==4) return 1;
  if (jet_DoubleHadLabel[jetIndex]==0 || jet_DoubleHadLabel[jetIndex]==15) return 2; //includes tau 
  if (jet_DoubleHadLabel[jetIndex]==55) return 3;
  if (jet_DoubleHadLabel[jetIndex]==44) return 4;
  else
    return 5; // does not include taus 
  return iflav;
}
int DAOD_selector::getGhosFlavLabel(int jetIndex)
{
  int iflav = -1;
  // if (jet_ghostBHadCount[jetIndex]>1  && jet_ghostCHadCount[jetIndex]>1 ) return 5; 
  // if (jet_ghostBHadCount[jetIndex]==1 ) return 0; 
  // if (jet_ghostBHadCount[jetIndex]==0 && jet_ghostCHadCount[jetIndex]==1) return 1; 
  // if (jet_ghostBHadCount[jetIndex]==0 && jet_ghostCHadCount[jetIndex]==0) return 2; 
  // if (jet_ghostBHadCount[jetIndex]>1  ) return 3; 
  // if (jet_ghostBHadCount[jetIndex]==0 && jet_ghostCHadCount[jetIndex]>1 ) return 4; 
  if (jet_ghostBHadCount[jetIndex]==1 ) return 0; 
  if (jet_ghostBHadCount[jetIndex]==0 && jet_ghostCHadCount[jetIndex]==1) return 1; 
  if (jet_ghostBHadCount[jetIndex]==0 && jet_ghostCHadCount[jetIndex]==0) return 2; // includes tau 
  if (jet_ghostBHadCount[jetIndex]==2 ) return 3; 
  if (jet_ghostBHadCount[jetIndex]==0 && jet_ghostCHadCount[jetIndex]==2) return 4; 
  else
    return 5; 
  return iflav;
}
void DAOD_selector::getJetFeaturesInFlavorLabelMatrix()
{
  // here use m_jetFlavorMatrix[icone][ighost] to fill histograms
  std::string HistoName[6] = {"1B", "1D0B", "0B0D", "2B", "2D0B", "2B2D"};
  std::string hNameLab;
  std::string hVariable;
  std::vector<int> jVec;
  int jetIndex=-1;
  for(uint32_t iX=0; iX<6; iX++){
    for(uint32_t iY=0; iY<6; iY++){
      hNameLab = HistoName[iX]+"_"+HistoName[iY];
      jVec = m_jetFlavorMatrix[iX][iY];
      if (jVec.size()>0)
	{
	  for (unsigned int ijet=0; ijet<jVec.size(); ++ijet)
	    {
	      jetIndex = jVec[ijet];
	      if (jetIndex < 0)
		{
		  std::cout<<"jetIndex not found ... skipping"<<std::endl;
		  continue;
		}

	      // jet pT
	      // name of the histogram => hVariable
	      hVariable = "Labels_"+hNameLab+"_jetPt";
	      //HistopT[iX][iY]->Fill( jet_pt[jetIndex] );
	      FillHisto(hVariable,                    0.001*jet_pt[jetIndex]);
	      FillHisto("Labels_"+hNameLab+"_jetEta",       jet_eta[jetIndex]);

	      
	      //HistopT[iX][iY]->Fill( jet_pt[jetIndex] );
	      double ptH = 0.;
	      double jetEta = jet_eta[jetIndex];
	      double jetPhi = jet_phi[jetIndex];
	      double dR=0.;double dPhi=0.;double dEta=0.;
	      if (jet_nBHadr[jetIndex]>0)
		{
		  for (int i=0; i<jet_nBHadr[jetIndex]; ++i)
		    {
		      ptH = jet_bH_pt[jetIndex].at(i);
		      FillHisto("Labels_"+hNameLab+"_bHPt", 0.001*ptH);
		      FillHisto("Labels_"+hNameLab+"_bHPtFraction", ptH/jet_pt[jetIndex]);
		      dEta = jetEta - jet_bH_eta[jetIndex].at(i);
		      dPhi = jetPhi - jet_bH_phi[jetIndex].at(i);
		      if (fabs(dPhi)>M_PI) dPhi=2*M_PI-abs(dPhi);
		      dR = sqrt( dEta*dEta + dPhi*dPhi );
		      FillHisto("Labels_"+hNameLab+"_bHjetDR", dR);
		    }
		}
	      // track variables
	      int ntIPxD=0;
	      int ntSV1=0;
	      int ntJF=0;
	      int ntInAlgo[5]={0,0,0,0,0};// count tracks used by the vavious algorithms
	      std::string stringAlgo[5]={"IP2D","IP3D","RNNIP","SV1","JF"};
	      for (unsigned int algoBit=0; algoBit<5; ++algoBit) /// 0=IP2D, 1=IP3D, 2=RNNIP, 3=SV1, 4=JF
		{
		  for (unsigned int i=0; i<jet_trk_pt[jetIndex].size(); ++i)
		    {
		      int algoBits=jet_trk_algo[jetIndex].at(i);
		      if (CHECK_BIT(algoBits,algoBit))
			{
			  ntInAlgo[algoBit]+=1;
			}
		    }
		  FillHisto("Labels_"+hNameLab+"_nTrkAlgo_"+stringAlgo[algoBit], ntInAlgo[algoBit]);
		}
	      /// algo results 
	      /// vertices 
	      //jet_sv1_Nvtx;
	      FillHisto("Labels_"+hNameLab+"_SV1nVtx", jet_sv1_Nvtx[jetIndex]);
	      FillHisto("Labels_"+hNameLab+"_SV1mVtx", 0.001*jet_sv1_m[jetIndex]);
	      FillHisto("Labels_"+hNameLab+"_SV1eFc",  jet_sv1_efc[jetIndex]);
	      FillHisto("Labels_"+hNameLab+"_SV1sig3d",jet_sv1_sig3d[jetIndex]);
	      FillHisto("Labels_"+hNameLab+"_SV1dR",   jet_sv1_deltaR[jetIndex]);
	      //jet_jf_nvtx;
	      FillHisto("Labels_"+hNameLab+"_JFnVtx", jet_jf_nvtx[jetIndex]);
	      FillHisto("Labels_"+hNameLab+"_JFmVtx", 0.001*jet_jf_m[jetIndex]);
	      FillHisto("Labels_"+hNameLab+"_JFeFc",  jet_jf_efc[jetIndex]);
	      FillHisto("Labels_"+hNameLab+"_JFsig3d",jet_jf_sig3d[jetIndex]);
	      FillHisto("Labels_"+hNameLab+"_JFdR",   jet_jf_dR[jetIndex]);
	    }
	}
    }
  }      



  return;
}

void DAOD_selector::BookHisto (std::string name, std::string nametitle,  
			       Int_t nx, Double_t xlow, Double_t xup)
{
  if (nametitle=="") nametitle=name;
  TH1D* hT1y = new TH1D(name.c_str(), nametitle.c_str(), nx, xlow, xup);
  fHisto1DMap[name]=hT1y;
  return;
}
void DAOD_selector::BookHisto2(std::string name, std::string nametitle,  
			       Int_t nx, Double_t xlow, Double_t xup, Int_t ny, Double_t ylow, Double_t yup)
{
  if (nametitle=="") nametitle=name;
  TH2D* hT1y = new TH2D(name.c_str(), nametitle.c_str(), nx, xlow, xup, ny, ylow, yup);
  fHisto2DMap[name]=hT1y;
  return;
}
void DAOD_selector::FillHisto (std::string hname, Double_t bin, Double_t weight)
{
  if (fHisto1DMap.find(hname)==fHisto1DMap.end()) {
    std::cout << "---> warning from HistoSvc::FillHisto() : histo " << hname
	      << " does not exist. (xbin=" << bin << " weight=" << weight << ")"
	      << std::endl;
    return;
  }
  if  (fHisto1DMap[hname]) {fHisto1DMap[hname]->Fill(bin, weight); }
  else
    {
      std::cout << "---> warning from HistoSvc::FillHisto() : histo " << hname
		<< " found in the map with a NULL pointer "
		<< std::endl;
    }
  return;
}
void DAOD_selector::FillHisto2(std::string hname, Double_t xbin, Double_t ybin, Double_t weight)
{
  if (fHisto2DMap.find(hname)==fHisto2DMap.end()) {
    std::cout << "---> warning from HistoSvc::FillHisto() : histo2D " << hname
	      << " does not exist. (xbin=" << xbin << " ybin=" << ybin <<" weight=" << weight << ")"
	      << std::endl;
    return;
  }
  if  (fHisto2DMap[hname]) {fHisto2DMap[hname]->Fill(xbin, ybin, weight); }
  else
    {
      std::cout << "---> warning from HistoSvc::FillHisto() : histo " << hname
		<< " found in the map with a NULL pointer "
		<< std::endl;
    }  
  return;
}
void DAOD_selector::Normalize (std::string hname, Double_t fac)
{
   if (fHisto1DMap.find(hname)!=fHisto1DMap.end()) {
     if  (fHisto1DMap[hname]) {
       fHisto1DMap[hname]->Scale(fac);
     }
     else{
       std::cout << "---> warning from HistoSvc::FillHisto() : histo " << hname
		 << " found in the map with a NULL pointer "
		 << std::endl;
     }
     //    std::cout << "---> warning from HistoSvc::FillHisto() : histo2D " << hname
     //           << ", to be normalized, does not exist. " << xbin << " ybin=" << ybin <<" weight=" << weight << ")"
     //           << std::endl;
    return;
   }
   else{
     if (fHisto2DMap.find(hname)!=fHisto2DMap.end()) {
       if  (fHisto2DMap[hname]) {
	 fHisto2DMap[hname]->Scale(fac);
       }
       else{
	 std::cout << "---> warning from HistoSvc::NormalizeHisto() : histo2D " << hname
		   << " found in the map with a NULL pointer "
		   << std::endl;
       }
       return;
     }
     else
       {
	 std::cout << "---> warning from HistoSvc::NormalizeHisto() : histo " << hname
		   << " NOT found in the map of 1D histos nor in the map of 2D histos"
		   << std::endl;
	 return;
       }
   }
  return;
} 
void DAOD_selector::saveHistosInMaps()
{
  for (std::map<std::string, TH1D*>::const_iterator it=fHisto1DMap.begin(); it!=fHisto1DMap.end(); ++it )
    (*it).second->Write();
  for (std::map<std::string, TH2D*>::const_iterator it=fHisto2DMap.begin(); it!=fHisto2DMap.end(); ++it )
    (*it).second->Write();

}
