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


void DAOD_selector::setFlags(bool lxplus_flag, bool debug_flag, bool selections_flag, bool discriminants_flag, bool shrinking_cone_flag, bool selection_alg_flag, bool origin_selection_flag, bool geometric_selection_flag, bool cut_flag, bool retagT_flag)
{
   std::cout<<"\n In DAOD_selector::setFlags"<<std::endl;

   lxplus=lxplus_flag;
   debug=debug_flag;
   selections=selections_flag;
   discriminants=discriminants_flag;
   shrinking_cone=shrinking_cone_flag;
   selection_alg=selection_alg_flag;
   origin_selection=origin_selection_flag;
   geometric_selection=geometric_selection_flag;
   cut=cut_flag;
   retagT=retagT_flag;
   return;
}

void DAOD_selector::setCuts(float m_jet_pT_infcut, float m_jet_pT_supcut, float m_jet_eta_cut, float m_jet_JVT_cut, float m_DR_bcH_cut, float m_pT_bcH_cut, float m_trk_pT_cut, float m_trk_eta_cut, float m_trk_d0_cut)
{
   std::cout<<"\n In DAOD_selector::setCuts"<<std::endl;
   jet_pT_infcut=m_jet_pT_infcut;
   jet_pT_supcut=m_jet_pT_supcut;
   jet_eta_cut=m_jet_eta_cut;
   jet_JVT_cut=m_jet_JVT_cut;
   m_DR_bcH_truth_cut=m_DR_bcH_cut;
   m_pT_bcH_truth_cut=m_pT_bcH_cut;
   trk_pT_cut=m_trk_pT_cut;
   trk_eta_cut=m_trk_eta_cut;
   trk_d0_cut=m_trk_d0_cut;
   return;
}

void DAOD_selector::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).


   std::cout<<"\n In DOAD_selector::Begin ... "<<std::endl;

   initFlagsAndCuts();
   initCounters();

   //tmp_pTfraction=0.,tmp_DR=0.,tmp_min_pTfraction=1.,tmp_min_DR=1.;

   m1=0,m2=0,m1_ex=0,m2_ex=0,mm1_ex=0,mm2_ex=0,m1_ov=0,m2_ov=0,mm=0;



   JF_ntrk=0,SV1_ntrk=0,SV0_ntrk=0,IP2D_ntrk=0,IP3D_ntrk=0;

   //   std::cout<<"\n";
   openOutputFile();

   bookHistosForSelections();

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
   std::vector<int> isJet,isBcheck,isCcheck,islcheck;

   if (m_Ntot%1000==0) std::cout<<".... in Process ... events processed so far "<<m_Ntot<<"/"<<fReader.GetCurrentEntry()<<"  out of "<<fReader.GetEntries()<<" Run / event nb = "<<*runnb<<"/"<<*eventnb<<std::endl;

   m_Ntot++;

   double D_phi=0.,D_eta=0.,DR=0.,px=0.,py=0.,pz=0.,Dx_1=0.,Dy_1=0.,Dz_1=0.,Dx_2=0.,Dy_2=0.,Dz_2=0,Lxy=0,Lxyz=0,Dxy_1=0.,x0=0.,y0=0.,Dx_3=0.,Dy_3=0.,Dxy_3=0.,rand_n=0.,R0=0.,d0=0.;

   double D_phi_trk=0.,D_eta_trk=0.,DR_trk=0.,DpT_trk=0.;
   double tmp_pTfraction=0.,tmp_DR=0.,tmp_min_pTfraction=1.,tmp_min_DR=1.;
   int match=0,max_size=0;

   isJet.clear();


   // here loop over jets
   std::vector<int> isB_1,isB_2,isC_1,isC_2;
   for(int i=0;i<*njets;i++) {

     if(jet_nBHadr[i]>0){
       nBjets++;
       isB_1.push_back(i);
     }

     if(jet_nCHadr[i]>0){
       nCjets++;
       isC_1.push_back(i);
     }

     if(jet_nBHadr[i]==0 && jet_nCHadr[i]==0){
       nljets++;
     }


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


// jet labelling (b-, c-, l-) based on truth (with pt, Deta cuts), exclusive samples
   getTrueJetFlavourLabel(isJet, isBcheck, isCcheck, islcheck);

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

     if(nBcheck>0){
       for(std::vector<int>::iterator it = isBcheck.begin(); it != isBcheck.end(); ++it){
         hist_pt_inB->Fill(jet_pt[*it]*0.001);
         hist_eta_inB->Fill(jet_eta[*it]);
         hist_phi_inB->Fill(jet_phi[*it]);
         hist_E_inB->Fill(jet_E[*it]);
         for(unsigned i=0;i<jet_trk_orig[*it].size();i++){
           hist_Bjet_origin->Fill(jet_trk_orig[*it].at(i));
         }
       }
     }

     if(nCcheck>0){
       for(std::vector<int>::iterator it = isBcheck.begin(); it != isBcheck.end(); ++it){
         hist_pt_inC->Fill(jet_pt[*it]*0.001);
         hist_eta_inC->Fill(jet_eta[*it]);
         hist_phi_inC->Fill(jet_phi[*it]);
         hist_E_inC->Fill(jet_E[*it]);
         for(unsigned i=0;i<jet_trk_orig[*it].size();i++){
           hist_Cjet_origin->Fill(jet_trk_orig[*it].at(i));
         }
       }
     }

     if(nlcheck>0){
       for(std::vector<int>::iterator it = isBcheck.begin(); it != isBcheck.end(); ++it){
         hist_pt_l->Fill(jet_pt[*it]*0.001);
         hist_eta_l->Fill(jet_eta[*it]);
         hist_phi_l->Fill(jet_phi[*it]);
         hist_E_l->Fill(jet_E[*it]);
         for(unsigned i=0;i<jet_trk_orig[*it].size();i++){
           hist_ljet_origin->Fill(jet_trk_orig[*it].at(i));
         }
       }
     }
/*
     if(n1==1 && n2==0 && n3==0){
       int i=is1B.at(0);
       m_b++;
       hist_pt_1->Fill(jet_pt[i]*0.001);
       hist_eta_1->Fill(jet_eta[i]);
       hist_phi_1->Fill(jet_phi[i]);
       hist_E_1->Fill(jet_E[i]);

     }

     if(n1==2 && n2==0 && n3==0){
       m_bb++;
       for(std::vector<int>::iterator it = is1B.begin(); it != is1B.end(); ++it){
         hist_pt_2->Fill(jet_pt[*it]*0.001);
         hist_eta_2->Fill(jet_eta[*it]);
         hist_phi_2->Fill(jet_phi[*it]);
         hist_E_2->Fill(jet_E[*it]);
       }
     }

     //select by: n1==0 && n2==1 && n3==0
     if(n1==0 && n2==1 && n3==0){
       for(std::vector<int>::iterator it = is2B.begin(); it != is2B.end(); ++it){
         hist_pt_2b->Fill(jet_pt[*it]*0.001);
         hist_eta_2b->Fill(jet_eta[*it]);
         hist_phi_2b->Fill(jet_phi[*it]);
         hist_E_2b->Fill(jet_E[*it]);
       }
     }

     //select by: n1==3 && n2==0 && n3==0
     if(n1==3 && n2==0 && n3==0){
       for(std::vector<int>::iterator it = is1B.begin(); it != is1B.end(); ++it){
         hist_pt_3a->Fill(jet_pt[*it]*0.001);
         hist_eta_3a->Fill(jet_eta[*it]);
         hist_phi_3a->Fill(jet_phi[*it]);
         hist_E_3a->Fill(jet_E[*it]);
       }
     }
     //select by: n1==0 && n2==0 && n3==1
     if(n1==0 && n2==0 && n3==1){
       for(std::vector<int>::iterator it = is3B.begin(); it != is3B.end(); ++it){
         hist_pt_3b->Fill(jet_pt[*it]*0.001);
         hist_eta_3b->Fill(jet_eta[*it]);
         hist_phi_3b->Fill(jet_phi[*it]);
         hist_E_3b->Fill(jet_E[*it]);
        }
     }

     //select by: n1==4 && n2==0 && n3==0
     if(n1==4 && n2==0 && n3==0){
       for(std::vector<int>::iterator it = is1B.begin(); it != is1B.end(); ++it){
         hist_pt_4->Fill(jet_pt[*it]*0.001);
         hist_eta_4->Fill(jet_eta[*it]);
         hist_phi_4->Fill(jet_phi[*it]);
         hist_E_4->Fill(jet_E[*it]);
       }
     }
*/
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

//DISCRIMINANTS
   if(discriminants){
     double DL1=0,RNNIP=0;
     if(nBcheck>0){
       for(std::vector<int>::iterator it = isBcheck.begin(); it != isBcheck.end(); ++it){
         if(jet_ip2d_pb[*it]!=-99){
           hist_ip2d_llr_exB->Fill(jet_ip2d_llr[*it]); //llr is computed as log(pb/pu)
         }
         if(jet_ip3d_pb[*it]!=-99){
           hist_ip3d_llr_exB->Fill(jet_ip3d_llr[*it]); //llr is computed as log(pb/pu)
         }
         if(jet_rnnip_pb[*it]!=-99){
           RNNIP=log(jet_rnnip_pb[*it]/(m_fcRNNIP*jet_rnnip_pc[*it]+(1.-m_fcRNNIP)*jet_rnnip_pu[*it]));
           hist_rnnip_llr_exB->Fill(RNNIP); //llr is computed as log(pb/pu)
         }
         if(sv1_llr[*it]!=-99){
           hist_sv1_llr_exB->Fill(sv1_llr[*it]); //llr is computed as log(pb/pu)
         }
         if(jet_jf_llr[*it]!=-99){
           hist_jf_llr_exB->Fill(jet_jf_llr[*it]); //llr is computed as log(pb/pu)
         }
         if(jet_dl1_pb[*it]!=-99){
           DL1=log(jet_dl1_pb[*it]/(m_fc*jet_dl1_pc[*it]+(1-m_fc)*jet_dl1_pu[*it]));
           hist_dl1_exB->Fill(DL1);
         }
       }
     }


  //inclusive plots
     if(nCcheck>0){
       for(std::vector<int>::iterator it = isCcheck.begin(); it != isCcheck.end(); ++it){
         if(jet_ip2d_pb[*it]!=-99){
           hist_ip2d_llr_exC->Fill(jet_ip2d_llr[*it]); //llr is computed as log(pb/pu)
         }
         if(jet_ip3d_pb[*it]!=-99){
           hist_ip3d_llr_exC->Fill(jet_ip3d_llr[*it]); //llr is computed as log(pb/pu)
         }
         if(jet_rnnip_pb[*it]!=-99){
           RNNIP=log(jet_rnnip_pb[*it]/(m_fcRNNIP*jet_rnnip_pc[*it]+(1.-m_fcRNNIP)*jet_rnnip_pu[*it]));
           hist_rnnip_llr_exC->Fill(RNNIP); //llr is computed as log(pb/pu)
         }
         if(sv1_llr[*it]!=-99){
           hist_sv1_llr_exC->Fill(sv1_llr[*it]); //llr is computed as log(pb/pu)
         }
         if(jet_jf_llr[*it]!=-99){
           hist_jf_llr_exC->Fill(jet_jf_llr[*it]); //llr is computed as log(pb/pu)
         }
         if(jet_dl1_pb[*it]!=-99){
           DL1=log(jet_dl1_pb[*it]/(m_fc*jet_dl1_pc[*it]+(1-m_fc)*jet_dl1_pu[*it]));
           hist_dl1_exC->Fill(DL1);
         }
       }
     }

     if(nlcheck>0){
       for(std::vector<int>::iterator it = islcheck.begin(); it != islcheck.end(); ++it){
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
       }
     }



   }


   unsigned size_jet=0,size_child=0;

   if(selection_alg){
     if(nBcheck>0){
       for(std::vector<int>::iterator it = isBcheck.begin(); it != isBcheck.end(); ++it){//isBcheck isB

         std::vector<int> matching_trk,matching_child,origin_trk;
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

	         int den=0;
           for(unsigned i=0;i<size_jet;i++){
             if(abs(jet_trk_eta[*it].at(i))<trk_eta_cut && jet_trk_pt[*it].at(i)>trk_pT_cut && abs(jet_trk_d0[*it].at(i))<trk_d0_cut ){// && abs(jet_trk_z0[*it].at(i)*sin(jet_trk_theta[*it].at(i)))<1.5

               n_trk_pT_cut++;

               D_eta=jet_trk_eta[*it][i]-jet_eta[*it];
               if(abs(jet_trk_phi[*it].at(i)-jet_phi[*it])>M_PI){
                 D_phi=2*M_PI-abs(jet_trk_phi[*it].at(i)-jet_phi[*it]);
               }
               if(abs(jet_trk_phi[*it].at(i)-jet_phi[*it])<M_PI){
                 D_phi=jet_trk_phi[*it].at(i)-jet_phi[*it];
               }
               DR=sqrt(D_eta*D_eta+D_phi*D_phi);
               hist_trk_pT_inB->Fill(1e-3*jet_trk_pt[*it].at(i));
               hist_trk_eta_inB->Fill(jet_trk_eta[*it].at(i));
               hist_trk_phi_inB->Fill(jet_trk_phi[*it].at(i));
               hist_trk_Deta_inB->Fill(D_eta);
               hist_trk_Dphi_inB->Fill(D_phi);
               hist_trk_Dphi_Deta_inB->Fill(D_phi,D_eta);
               hist_trk_DR_inB->Fill(DR);
               hist_trk_pT_DR_inB->Fill(1e-3*jet_trk_pt[*it].at(i),DR);
               hist_trk_pT_jet_DR_inB->Fill(1e-3*jet_pt[*it],DR);
               hist_trk_pdgId_inB->Fill(jet_trk_pdg_id[*it].at(i));
               hist_trk_origin_inB->Fill(jet_trk_orig[*it].at(i));

               A=sin(jet_phi[*it]-jet_trk_phi[*it].at(i))*jet_trk_d0[*it].at(i);
               sgn=A/abs(A);
               d0=sgn*abs(jet_trk_d0[*it].at(i));
               hist_trk_d0_inB->Fill(d0);
               hist_trk_z0sinth_inB->Fill(sgn*jet_trk_z0[*it].at(i)*sin(jet_trk_theta[*it].at(i)));
               hist_trk_d0sig_inB->Fill(d0/jet_trk_d0sig[*it].at(i));
               hist_trk_z0sinthsig_inB->Fill(sgn*jet_trk_z0[*it].at(i)*sin(jet_trk_theta[*it].at(i))/jet_trk_d0sig[*it].at(i));
               hist_trk_d0sig_origin_inB->Fill(d0/jet_trk_d0sig[*it].at(i),jet_trk_orig[*it].at(i));
               hist_trk_z0sinthsig_origin_inB->Fill(sgn*jet_trk_z0[*it].at(i)*sin(jet_trk_theta[*it].at(i)),jet_trk_orig[*it].at(i));
               hist_trk_logpTfrac_origin_inB->Fill(log(jet_trk_pt[*it].at(i)/jet_pt[*it]),jet_trk_orig[*it].at(i));
               hist_trk_logDR_origin_inB->Fill(log(DR),jet_trk_orig[*it].at(i));
               hist_trk_IBLhits_origin_inB->Fill(jet_trk_nInnHits[*it].at(i),jet_trk_orig[*it].at(i));
               hist_trk_NextToIBLhits_origin_inB->Fill(jet_trk_nNextToInnHits[*it].at(i),jet_trk_orig[*it].at(i));
               hist_trk_sharedIBLhits_origin_inB->Fill(jet_trk_nsharedBLHits[*it].at(i),jet_trk_orig[*it].at(i));
               hist_trk_splitIBLhits_origin_inB->Fill(jet_trk_nsplitBLHits[*it].at(i),jet_trk_orig[*it].at(i));//jet_trk_nsplitPixHits
               hist_trk_nPixhits_origin_inB->Fill(jet_trk_nPixHits[*it].at(i),jet_trk_orig[*it].at(i));
               hist_trk_sharedPixhits_origin_inB->Fill(jet_trk_nsharedPixHits[*it].at(i),jet_trk_orig[*it].at(i));
               hist_trk_splitPixhits_origin_inB->Fill(jet_trk_nsplitPixHits[*it].at(i),jet_trk_orig[*it].at(i));
               hist_trk_nSCThits_origin_inB->Fill(jet_trk_nSCTHits[*it].at(i),jet_trk_orig[*it].at(i));
               hist_trk_sharedSCThits_origin_inB->Fill(jet_trk_nsharedSCTHits[*it].at(i),jet_trk_orig[*it].at(i));

//               std::cout<<jet_trk_algo[*it].at(i)<<"\n";
               std::vector<int> bin_alg(5);
               int x=jet_trk_algo[*it].at(i);

               int idx=0,p=0;
               while(idx<5){
                 p=x%2;
                 bin_alg.at(idx)=p;
                 if(idx==4) {
                   JF_ntrk=JF_ntrk+p;
                   if(p==1){
                     hist_trk_d0sig_JF_inB->Fill(d0/jet_trk_d0sig[*it].at(i));
                     hist_trk_z0sinthsig_JF_inB->Fill(sgn*jet_trk_z0[*it].at(i)*sin(jet_trk_theta[*it].at(i))/jet_trk_d0sig[*it].at(i));
                     hist_trk_d0sig_origin_JF_inB->Fill(d0/jet_trk_d0sig[*it].at(i),jet_trk_orig[*it].at(i));
                     hist_trk_z0sinthsig_origin_JF_inB->Fill(sgn*jet_trk_z0[*it].at(i)*sin(jet_trk_theta[*it].at(i))/jet_trk_d0sig[*it].at(i),jet_trk_orig[*it].at(i));
                     hist_trk_logpTfrac_origin_JF_inB->Fill(log(jet_trk_pt[*it].at(i)/jet_pt[*it]),jet_trk_orig[*it].at(i));
                     hist_trk_logDR_origin_JF_inB->Fill(log(DR),jet_trk_orig[*it].at(i));
                     hist_trk_IBLhits_origin_JF_inB->Fill(jet_trk_nInnHits[*it].at(i),jet_trk_orig[*it].at(i));
                     hist_trk_NextToIBLhits_origin_JF_inB->Fill(jet_trk_nNextToInnHits[*it].at(i),jet_trk_orig[*it].at(i));
                     hist_trk_sharedIBLhits_origin_JF_inB->Fill(jet_trk_nsharedBLHits[*it].at(i),jet_trk_orig[*it].at(i));
                     hist_trk_splitIBLhits_origin_JF_inB->Fill(jet_trk_nsplitBLHits[*it].at(i),jet_trk_orig[*it].at(i));//jet_trk_nsplitPixHits
                     hist_trk_nPixhits_origin_JF_inB->Fill(jet_trk_nPixHits[*it].at(i),jet_trk_orig[*it].at(i));
                     hist_trk_sharedPixhits_origin_JF_inB->Fill(jet_trk_nsharedPixHits[*it].at(i),jet_trk_orig[*it].at(i));
                     hist_trk_splitPixhits_origin_JF_inB->Fill(jet_trk_nsplitPixHits[*it].at(i),jet_trk_orig[*it].at(i));
                     hist_trk_nSCThits_origin_JF_inB->Fill(jet_trk_nSCTHits[*it].at(i),jet_trk_orig[*it].at(i));
                     hist_trk_sharedSCThits_origin_JF_inB->Fill(jet_trk_nsharedSCTHits[*it].at(i),jet_trk_orig[*it].at(i));
                   }
                 }
                 if(idx==3) {
                   SV1_ntrk=SV1_ntrk+p;
                   if(p==1){
                     hist_trk_d0sig_SV1_inB->Fill(d0/jet_trk_d0sig[*it].at(i));
                     hist_trk_z0sinthsig_SV1_inB->Fill(sgn*jet_trk_z0[*it].at(i)*sin(jet_trk_theta[*it].at(i))/jet_trk_d0sig[*it].at(i));
                     hist_trk_d0sig_origin_SV1_inB->Fill(d0/jet_trk_d0sig[*it].at(i),jet_trk_orig[*it].at(i));
                     hist_trk_z0sinthsig_origin_SV1_inB->Fill(sgn*jet_trk_z0[*it].at(i)*sin(jet_trk_theta[*it].at(i))/jet_trk_d0sig[*it].at(i),jet_trk_orig[*it].at(i));
                     hist_trk_logpTfrac_origin_SV1_inB->Fill(log(jet_trk_pt[*it].at(i)/jet_pt[*it]),jet_trk_orig[*it].at(i));
                     hist_trk_logDR_origin_SV1_inB->Fill(log(DR),jet_trk_orig[*it].at(i));
                     hist_trk_IBLhits_origin_SV1_inB->Fill(jet_trk_nInnHits[*it].at(i),jet_trk_orig[*it].at(i));
                     hist_trk_NextToIBLhits_origin_SV1_inB->Fill(jet_trk_nNextToInnHits[*it].at(i),jet_trk_orig[*it].at(i));
                     hist_trk_sharedIBLhits_origin_SV1_inB->Fill(jet_trk_nsharedBLHits[*it].at(i),jet_trk_orig[*it].at(i));
                     hist_trk_splitIBLhits_origin_SV1_inB->Fill(jet_trk_nsplitBLHits[*it].at(i),jet_trk_orig[*it].at(i));//jet_trk_nsplitPixHits
                     hist_trk_nPixhits_origin_SV1_inB->Fill(jet_trk_nPixHits[*it].at(i),jet_trk_orig[*it].at(i));
                     hist_trk_sharedPixhits_origin_SV1_inB->Fill(jet_trk_nsharedPixHits[*it].at(i),jet_trk_orig[*it].at(i));
                     hist_trk_splitPixhits_origin_SV1_inB->Fill(jet_trk_nsplitPixHits[*it].at(i),jet_trk_orig[*it].at(i));
                     hist_trk_nSCThits_origin_SV1_inB->Fill(jet_trk_nSCTHits[*it].at(i),jet_trk_orig[*it].at(i));
                     hist_trk_sharedSCThits_origin_SV1_inB->Fill(jet_trk_nsharedSCTHits[*it].at(i),jet_trk_orig[*it].at(i));
                   }
                 }
                 if(idx==2) {
                   SV0_ntrk=SV0_ntrk+p;
                   if(p==1){
                     hist_trk_d0sig_SV0_inB->Fill(d0/jet_trk_d0sig[*it].at(i));
                     hist_trk_z0sinthsig_SV0_inB->Fill(sgn*jet_trk_z0[*it].at(i)*sin(jet_trk_theta[*it].at(i))/jet_trk_d0sig[*it].at(i));
                     hist_trk_d0sig_origin_SV0_inB->Fill(d0/jet_trk_d0sig[*it].at(i),jet_trk_orig[*it].at(i));
                     hist_trk_z0sinthsig_origin_SV0_inB->Fill(sgn*jet_trk_z0[*it].at(i)*sin(jet_trk_theta[*it].at(i))/jet_trk_d0sig[*it].at(i),jet_trk_orig[*it].at(i));
                     hist_trk_logpTfrac_origin_SV0_inB->Fill(log(jet_trk_pt[*it].at(i)/jet_pt[*it]),jet_trk_orig[*it].at(i));
                     hist_trk_logDR_origin_SV0_inB->Fill(log(DR),jet_trk_orig[*it].at(i));
                     hist_trk_IBLhits_origin_SV0_inB->Fill(jet_trk_nInnHits[*it].at(i),jet_trk_orig[*it].at(i));
                     hist_trk_NextToIBLhits_origin_SV0_inB->Fill(jet_trk_nNextToInnHits[*it].at(i),jet_trk_orig[*it].at(i));
                     hist_trk_sharedIBLhits_origin_SV0_inB->Fill(jet_trk_nsharedBLHits[*it].at(i),jet_trk_orig[*it].at(i));
                     hist_trk_splitIBLhits_origin_SV0_inB->Fill(jet_trk_nsplitBLHits[*it].at(i),jet_trk_orig[*it].at(i));//jet_trk_nsplitPixHits
                     hist_trk_nPixhits_origin_SV0_inB->Fill(jet_trk_nPixHits[*it].at(i),jet_trk_orig[*it].at(i));
                     hist_trk_sharedPixhits_origin_SV0_inB->Fill(jet_trk_nsharedPixHits[*it].at(i),jet_trk_orig[*it].at(i));
                     hist_trk_splitPixhits_origin_SV0_inB->Fill(jet_trk_nsplitPixHits[*it].at(i),jet_trk_orig[*it].at(i));
                     hist_trk_nSCThits_origin_SV0_inB->Fill(jet_trk_nSCTHits[*it].at(i),jet_trk_orig[*it].at(i));
                     hist_trk_sharedSCThits_origin_SV0_inB->Fill(jet_trk_nsharedSCTHits[*it].at(i),jet_trk_orig[*it].at(i));
                   }
                 }
                 if(idx==1) {
                   IP3D_ntrk=IP3D_ntrk+p;
                   if(p==1){
                     hist_trk_d0sig_IP3D_inB->Fill(d0/jet_trk_d0sig[*it].at(i));
                     hist_trk_z0sinthsig_IP3D_inB->Fill(sgn*jet_trk_z0[*it].at(i)*sin(jet_trk_theta[*it].at(i))/jet_trk_d0sig[*it].at(i));
                     hist_trk_d0sig_origin_IP3D_inB->Fill(d0/jet_trk_d0sig[*it].at(i),jet_trk_orig[*it].at(i));
                     hist_trk_z0sinthsig_origin_IP3D_inB->Fill(sgn*jet_trk_z0[*it].at(i)*sin(jet_trk_theta[*it].at(i))/jet_trk_d0sig[*it].at(i),jet_trk_orig[*it].at(i));
                     hist_trk_logpTfrac_origin_IP3D_inB->Fill(log(jet_trk_pt[*it].at(i)/jet_pt[*it]),jet_trk_orig[*it].at(i));
                     hist_trk_logDR_origin_IP3D_inB->Fill(log(DR),jet_trk_orig[*it].at(i));
                     hist_trk_IBLhits_origin_IP3D_inB->Fill(jet_trk_nInnHits[*it].at(i),jet_trk_orig[*it].at(i));
                     hist_trk_NextToIBLhits_origin_IP3D_inB->Fill(jet_trk_nNextToInnHits[*it].at(i),jet_trk_orig[*it].at(i));
                     hist_trk_sharedIBLhits_origin_IP3D_inB->Fill(jet_trk_nsharedBLHits[*it].at(i),jet_trk_orig[*it].at(i));
                     hist_trk_splitIBLhits_origin_IP3D_inB->Fill(jet_trk_nsplitBLHits[*it].at(i),jet_trk_orig[*it].at(i));//jet_trk_nsplitPixHits
                     hist_trk_nPixhits_origin_IP3D_inB->Fill(jet_trk_nPixHits[*it].at(i),jet_trk_orig[*it].at(i));
                     hist_trk_sharedPixhits_origin_IP3D_inB->Fill(jet_trk_nsharedPixHits[*it].at(i),jet_trk_orig[*it].at(i));
                     hist_trk_splitPixhits_origin_IP3D_inB->Fill(jet_trk_nsplitPixHits[*it].at(i),jet_trk_orig[*it].at(i));
                     hist_trk_nSCThits_origin_IP3D_inB->Fill(jet_trk_nSCTHits[*it].at(i),jet_trk_orig[*it].at(i));
                     hist_trk_sharedSCThits_origin_IP3D_inB->Fill(jet_trk_nsharedSCTHits[*it].at(i),jet_trk_orig[*it].at(i));
                   }
                 }
                 if(idx==0) {
                   IP2D_ntrk=IP2D_ntrk+p;
                   if(p==1){
                     hist_trk_d0sig_IP2D_inB->Fill(d0/jet_trk_d0sig[*it].at(i));
                     hist_trk_z0sinthsig_IP2D_inB->Fill(sgn*jet_trk_z0[*it].at(i)*sin(jet_trk_theta[*it].at(i))/jet_trk_d0sig[*it].at(i));
                     hist_trk_d0sig_origin_IP2D_inB->Fill(d0/jet_trk_d0sig[*it].at(i),jet_trk_orig[*it].at(i));
                     hist_trk_z0sinthsig_origin_IP2D_inB->Fill(sgn*jet_trk_z0[*it].at(i)*sin(jet_trk_theta[*it].at(i))/jet_trk_d0sig[*it].at(i),jet_trk_orig[*it].at(i));
                     hist_trk_logpTfrac_origin_IP2D_inB->Fill(log(jet_trk_pt[*it].at(i)/jet_pt[*it]),jet_trk_orig[*it].at(i));
                     hist_trk_logDR_origin_IP2D_inB->Fill(log(DR),jet_trk_orig[*it].at(i));
                     hist_trk_IBLhits_origin_IP2D_inB->Fill(jet_trk_nInnHits[*it].at(i),jet_trk_orig[*it].at(i));
                     hist_trk_NextToIBLhits_origin_IP2D_inB->Fill(jet_trk_nNextToInnHits[*it].at(i),jet_trk_orig[*it].at(i));
                     hist_trk_sharedIBLhits_origin_IP2D_inB->Fill(jet_trk_nsharedBLHits[*it].at(i),jet_trk_orig[*it].at(i));
                     hist_trk_splitIBLhits_origin_IP2D_inB->Fill(jet_trk_nsplitBLHits[*it].at(i),jet_trk_orig[*it].at(i));//jet_trk_nsplitPixHits
                     hist_trk_nPixhits_origin_IP2D_inB->Fill(jet_trk_nPixHits[*it].at(i),jet_trk_orig[*it].at(i));
                     hist_trk_sharedPixhits_origin_IP2D_inB->Fill(jet_trk_nsharedPixHits[*it].at(i),jet_trk_orig[*it].at(i));
                     hist_trk_splitPixhits_origin_IP2D_inB->Fill(jet_trk_nsplitPixHits[*it].at(i),jet_trk_orig[*it].at(i));
                     hist_trk_nSCThits_origin_IP2D_inB->Fill(jet_trk_nSCTHits[*it].at(i),jet_trk_orig[*it].at(i));
                     hist_trk_sharedSCThits_origin_IP2D_inB->Fill(jet_trk_nsharedSCTHits[*it].at(i),jet_trk_orig[*it].at(i));
                   }
                 }
                 x=x/2;
                 idx++;
               }
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

               if(jet_trk_orig[*it].at(i)==-1){
                 hist_trk_d0_PUinB->Fill(d0);
                 n_trk_PU_pT_cut++;
               }
               if(jet_trk_orig[*it].at(i)==0){
                 hist_trk_d0_BinB->Fill(d0);
                 n_trk_B++;
               }
               if(jet_trk_orig[*it].at(i)==1){
                 hist_trk_d0_CinB->Fill(d0);
                 n_trk_C++;
               }
               if(jet_trk_orig[*it].at(i)==2){
                 hist_trk_d0_FRAGinB->Fill(d0);
                 n_trk_FRAG_pT_cut++;
               }
               if(jet_trk_orig[*it].at(i)==3){
                 hist_trk_d0_GEANTinB->Fill(d0);
                 n_trk_GEANT_pT_cut++;
               }

               if(origin_selection){
                 if(jet_trk_orig[*it].at(i)==0 || jet_trk_orig[*it].at(i)==1){
                   origin_trk.push_back(i);

                   hist_matched_origin_pT_inB->Fill(1e-3*jet_trk_pt[*it].at(i));
                   hist_matched_origin_eta_inB->Fill(jet_trk_eta[*it].at(i));
                   hist_matched_origin_phi_inB->Fill(jet_trk_phi[*it].at(i));
                   hist_matched_origin_Deta_inB->Fill(D_eta);
                   hist_matched_origin_Dphi_inB->Fill(D_phi);
                   hist_matched_origin_Dphi_Deta_inB->Fill(D_phi,D_eta);
                   hist_matched_origin_DR_inB->Fill(DR);
                   hist_matched_origin_pT_DR_inB->Fill(1e-3*jet_trk_pt[*it].at(i),DR);
                   hist_matched_origin_pT_jet_DR_inB->Fill(1e-3*jet_pt[*it],DR);
                   hist_matched_origin_pdgId_inB->Fill(jet_trk_pdg_id[*it].at(i));
                   hist_matched_origin_origin_inB->Fill(jet_trk_orig[*it].at(i));
                   hist_matched_origin_d0_inB->Fill(d0);
//                   hist_matched_origin_Lxy_inB->Fill();
//                   hist_matched_origin_Lxyz_inB->Fill();

                 }

               }

             }
           }

           std::vector<double> child_Pt,child_Eta,child_Phi;
           std::vector<int> child_idx;

           double child_IP[size_child],child_Lxy[size_child],child_Lxyz[size_child];
           den=0;//n of child_400 per jet

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

//             float z0sinth=(jet_bH_child_z0[*it].at(j)-(*PVz))*sin(jet_bH_child_theta[*it].at(j));


             if(abs(child_Eta.at(j))<trk_eta_cut && child_Pt.at(j)>trk_pT_cut && abs(jet_bH_child_d0[*it].at(j))<trk_d0_cut ){//CHILD SELECTION CRITERIA && abs(z0sinth)<1.5
               if(jet_bH_child_prod_x[*it].size()!=size_child || jet_bH_child_decay_x[*it].size()!=size_child){
                 std::cout<<"WARNING\n";
               }
               den++;

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

                 hist_child_Lxy_inB->Fill(Lxy);
                 hist_child_Lxyz_inB->Fill(Lxyz);
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

                 hist_child_Lxy_inB->Fill(Lxy);
                 hist_child_Lxyz_inB->Fill(Lxyz);
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

                 hist_pT_vs_R0_ratio_inB->Fill(child_Pt.at(j)/abs(R0));

//                 ncorr+=jet_bH_child_charge[*it].at(j)*d0/abs(d0);
//                 std::cout<<jet_bH_child_charge[*it].at(j)*R0/abs(jet_bH_child_charge[*it].at(j)*R0)<<"\n";
//                 std::cout<<jet_bH_child_charge[*it].at(j)<<"\t"<<d0/abs(d0)<<"\n";
               }//DECAY CHILD

               d0=jet_bH_child_d0[*it].at(j);
               hist_child_d0->Fill(d0);
               hist_child_d0_pT->Fill(d0,1e-3*sqrt(px*px+py*py));
               child_IP[j]=d0;

               A=sin(jet_phi[*it]-child_Phi.at(j))*jet_bH_child_d0[*it].at(j);
               sgn=A/abs(A);

               hist_child_d0_truth->Fill(sgn*abs(jet_bH_child_d0[*it].at(j)));
               if(jet_bH_child_theta[*it].at(j)==-99)
               std::cout<<"W\n";
               hist_child_z0sinth_inB->Fill((jet_bH_child_z0[*it].at(j)-(*PVz))*sin(jet_bH_child_theta[*it].at(j)));//*sqrt(px*px+py*py)/sqrt(px*px+py*py+pz*pz)

               D_eta=child_Eta[j]-jet_eta[*it];
               if(abs(child_Phi[j]-jet_phi[*it])>M_PI){
                 D_phi=2*M_PI-abs(child_Phi[j]-jet_phi[*it]);
               }
               if(abs(child_Phi[j]-jet_phi[*it])<M_PI){
                 D_phi=child_Phi[j]-jet_phi[*it];
               }
               DR=sqrt(D_eta*D_eta+D_phi*D_phi);
               hist_child_pT_inB->Fill(1e-3*child_Pt[j]);
               hist_child_eta_inB->Fill(child_Eta[j]);
               hist_child_phi_inB->Fill(child_Phi[j]);
               hist_child_Deta_inB->Fill(D_eta);
               hist_child_Dphi_inB->Fill(D_phi);
               hist_child_Dphi_Deta_inB->Fill(D_phi,D_eta);
               hist_child_DR_inB->Fill(DR);
               hist_child_pT_DR_inB->Fill(1e-3*child_Pt[j],DR);
               hist_child_pT_jet_DR_inB->Fill(1e-3*jet_pt[*it],DR);
               hist_child_pdgID_inB->Fill(jet_bH_child_pdg_id[*it].at(j));

               if(abs(jet_bH_child_pdg_id[*it].at(j))==211){
                 hist_child_pi->Fill(1e-3*child_Pt[j]);
                 hist_child_pi_Lxy_inB->Fill(Lxy);
                 hist_child_pi_Lxyz_inB->Fill(Lxyz);//CHECK
                 hist_child_pi_d0_truth_inB->Fill(sgn*abs(jet_bH_child_d0[*it].at(j)));
                 hist_child_pi_z0_truth_inB->Fill(jet_bH_child_z0[*it].at(j));
               }
               if(abs(jet_bH_child_pdg_id[*it].at(j))==321){
                 hist_child_K->Fill(1e-3*child_Pt[j]);
                 hist_child_K_Lxy_inB->Fill(Lxy);
                 hist_child_K_Lxyz_inB->Fill(Lxyz);
                 hist_child_K_d0_truth_inB->Fill(sgn*abs(jet_bH_child_d0[*it].at(j)));
                 hist_child_K_z0_truth_inB->Fill(jet_bH_child_z0[*it].at(j));
               }
               if(abs(jet_bH_child_pdg_id[*it].at(j))==13){
                 hist_child_mu->Fill(1e-3*child_Pt[j]);
                 hist_child_mu_Lxy_inB->Fill(Lxy);
                 hist_child_mu_Lxyz_inB->Fill(Lxyz);
                 hist_child_mu_d0_truth_inB->Fill(sgn*abs(jet_bH_child_d0[*it].at(j)));
                 hist_child_mu_z0_truth_inB->Fill(jet_bH_child_z0[*it].at(j));
               }
               if(abs(jet_bH_child_pdg_id[*it].at(j))==2212){
                 hist_child_p->Fill(1e-3*child_Pt[j]);
                 hist_child_p_Lxy_inB->Fill(Lxy);
                 hist_child_p_Lxyz_inB->Fill(Lxyz);
                 hist_child_p_d0_truth_inB->Fill(sgn*abs(jet_bH_child_d0[*it].at(j)));
                 hist_child_p_z0_truth_inB->Fill(jet_bH_child_z0[*it].at(j));
               }
               if(abs(jet_bH_child_pdg_id[*it].at(j))==11){
                 hist_child_e->Fill(1e-3*child_Pt[j]);
                 hist_child_e_Lxy_inB->Fill(Lxy);
                 hist_child_e_Lxyz_inB->Fill(Lxyz);
                 hist_child_e_d0_truth_inB->Fill(sgn*abs(jet_bH_child_d0[*it].at(j)));
                 hist_child_e_z0_truth_inB->Fill(jet_bH_child_z0[*it].at(j));
               }
             }
           }



           if(geometric_selection){

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
//               float sintheta=child_Pt.at(j)/sqrt(child_Pt.at(j)*child_Pt.at(j)+jet_bH_child_pz[*it].at(j)*jet_bH_child_pz[*it].at(j));
               if(abs(child_Eta.at(j))<trk_eta_cut && child_Pt.at(j)>trk_pT_cut && abs(jet_bH_child_d0[*it].at(j))<trk_d0_cut){// && abs(jet_bH_child_z0[*it].at(j)*sintheta)<1.5
                 if(jet_bH_child_charge[*it].at(j)==0){
                   std::cout<<"CHARGE 0\n";
                 }
                 for(unsigned i=0;i<size_jet;i++){
                   if(abs(jet_trk_eta[*it].at(i))<trk_eta_cut && jet_trk_pt[*it].at(i)>trk_pT_cut && abs(jet_trk_d0[*it].at(i))<trk_d0_cut){// && abs(jet_trk_z0[*it].at(i)*sin(jet_trk_theta[*it].at(i)))<1.5
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
                   std::cout<<"EXCLUSED CHILDS:\t"<<m_Ntot<<"\t"<<tmp_DpT<<"\t"<<tmp_DR<<"\n";
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
                 hist_matched_pT_inB->Fill(1e-3*child_Pt.at(b));
                 hist_matched_eta_inB->Fill(child_Eta.at(b));
                 hist_matched_phi_inB->Fill(child_Phi.at(b));
                 hist_matched_Deta_inB->Fill(D_eta);
                 hist_matched_Dphi_inB->Fill(D_phi);
                 hist_matched_Dphi_Deta_inB->Fill(D_phi,D_eta);
                 hist_matched_DR_inB->Fill(DR);
                 hist_matched_pT_DR_inB->Fill(1e-3*child_Pt.at(b),DR);
                 hist_matched_pT_jet_DR_inB->Fill(1e-3*jet_pt[*it],DR);
                 hist_matched_pdgId_inB->Fill(jet_bH_child_pdg_id[*it].at(b));

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
                   hist_matched_child_pi_Lxy_inB->Fill(Lxy);
                   hist_matched_child_pi_Lxyz_inB->Fill(Lxyz);
                   hist_matched_child_pi_d0_truth_inB->Fill(sgn*abs(jet_bH_child_d0[*it].at(b)));
                   hist_matched_child_pi_z0_truth_inB->Fill(jet_bH_child_z0[*it].at(b));
                 }
                 if(abs(jet_bH_child_pdg_id[*it].at(b))==321){
                   hist_matched_child_K->Fill(1e-3*child_Pt.at(b));
                   hist_matched_child_K_Lxy_inB->Fill(Lxy);
                   hist_matched_child_K_Lxyz_inB->Fill(Lxyz);
                   hist_matched_child_K_d0_truth_inB->Fill(sgn*abs(jet_bH_child_d0[*it].at(b)));
                   hist_matched_child_K_z0_truth_inB->Fill(jet_bH_child_z0[*it].at(b));
                 }
                 if(abs(jet_bH_child_pdg_id[*it].at(b))==13){
                   hist_matched_child_mu->Fill(1e-3*child_Pt.at(b));
                   hist_matched_child_mu_Lxy_inB->Fill(Lxy);
                   hist_matched_child_mu_Lxyz_inB->Fill(Lxyz);
                   hist_matched_child_mu_d0_truth_inB->Fill(sgn*abs(jet_bH_child_d0[*it].at(b)));
                   hist_matched_child_mu_z0_truth_inB->Fill(jet_bH_child_z0[*it].at(b));
                 }
                 if(abs(jet_bH_child_pdg_id[*it].at(b))==2212){
                   hist_matched_child_p->Fill(1e-3*child_Pt.at(b));
                   hist_matched_child_p_Lxy_inB->Fill(Lxy);
                   hist_matched_child_p_Lxyz_inB->Fill(Lxyz);
                   hist_matched_child_p_d0_truth_inB->Fill(sgn*abs(jet_bH_child_d0[*it].at(b)));
                   hist_matched_child_p_z0_truth_inB->Fill(jet_bH_child_z0[*it].at(b));
                 }
                 if(abs(jet_bH_child_pdg_id[*it].at(b))==11){
                   hist_matched_child_e->Fill(1e-3*child_Pt.at(b));
                   hist_matched_child_e_Lxy_inB->Fill(Lxy);
                   hist_matched_child_e_Lxyz_inB->Fill(Lxyz);
                   hist_matched_child_e_d0_truth_inB->Fill(sgn*abs(jet_bH_child_d0[*it].at(b)));
                   hist_matched_child_e_z0_truth_inB->Fill(jet_bH_child_z0[*it].at(b));
                 }

                 hist_matched_origin_inB->Fill(jet_trk_orig[*it].at(a));
                 hist_matched_d0_inB->Fill(child_IP[b]);
                 hist_matched_Lxy_inB->Fill(child_Lxy[b]);
                 hist_matched_Lxyz_inB->Fill(child_Lxyz[b]);

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
                 hist_matched_pT_child_pTfraction_inB->Fill(1e-3*child_Pt.at(b),-1.*DpT_trk/(1e-3*child_Pt.at(b)));
                 hist_matched_DR_trk_inB->Fill(DR_trk);
                 hist_matched_DR_trk_pTfraction->Fill(DR_trk,-1.*DpT_trk/(1e-3*child_Pt.at(b)));

                 if(sc==1){
                   m_sc+=1;
                   /*
                   hist_single_matched_pT_inB->Fill(1e-3*child_Pt.at(b));
                   hist_single_matched_eta_inB->Fill(child_Eta.at(b));
                   hist_single_matched_phi_inB->Fill(child_Phi.at(b));
                   hist_single_matched_Deta_inB->Fill(D_eta);
                   hist_single_matched_Dphi_inB->Fill(D_phi);
                   hist_single_matched_Dphi_Deta_inB->Fill(D_phi,D_eta);
                   hist_single_matched_DR_inB->Fill(DR);//RISPETTO AL JET
                   hist_single_matched_pT_DR_inB->Fill(1e-3*child_Pt.at(b),DR);
                   hist_single_matched_pT_jet_DR_inB->Fill(1e-3*jet_pt[*it],DR);
                   hist_single_matched_pdgId_inB->Fill(jet_bH_child_pdg_id[*it].at(b));
                   hist_single_matched_origin_inB->Fill(jet_trk_orig[*it].at(a));

                   hist_single_matched_pT_child_pTfraction_inB->Fill(1e-3*child_Pt.at(b),-1.*DpT_trk/(1e-3*child_Pt.at(b)));
                   hist_single_matched_DR_trk_inB->Fill(DR_trk);
                   hist_single_matched_DR_trk_pTfraction->Fill(DR_trk,-1.*DpT_trk/(1e-3*child_Pt.at(b)));

                   hist_single_matched_d0_inB->Fill(child_IP[b]);
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
  /*
             int j=0;
             for(unsigned l=0;l<child_idx.size();l++){
               j=child_idx.at(l);
               //CINEMATICA RISPETTO AL JET
               if(abs(child_Eta.at(j))<trk_eta_cut && child_Pt.at(j)>trk_pT_cut){
                 D_eta=child_Eta.at(j)-jet_eta[*it];
                 if(abs(child_Phi.at(j)-jet_phi[*it])>M_PI){
                   D_phi=2*M_PI-abs(child_Phi.at(j)-jet_phi[*it]);
                 }
                 if(abs(child_Phi.at(j)-jet_phi[*it])<M_PI){
                   D_phi=child_Phi.at(j)-jet_phi[*it];
                 }
                 DR=sqrt(D_eta*D_eta+D_phi*D_phi);
                 hist_nomatched_pT_inB->Fill(1e-3*child_Pt.at(j));
                 hist_nomatched_eta_inB->Fill(child_Eta.at(j));
                 hist_nomatched_phi_inB->Fill(child_Phi.at(j));
                 hist_nomatched_Deta_inB->Fill(D_eta);
                 hist_nomatched_Dphi_inB->Fill(D_phi);
                 hist_nomatched_Dphi_Deta_inB->Fill(D_phi,D_eta);
                 hist_nomatched_DR_inB->Fill(DR);
                 hist_nomatched_pT_DR_inB->Fill(1e-3*child_Pt.at(j),DR);
                 hist_nomatched_pT_jet_DR_inB->Fill(1e-3*jet_pt[*it],DR);
                 hist_nomatched_pdgId_inB->Fill(jet_bH_child_pdg_id[*it].at(j));
               }
             }
  */

    //std::cout<<"event:\t"<<m_Ntot<< "\tn of matched tracks:\t"<< match << "\tn of child tracks(denominator):\t" << den << "\tratio:\t" <<(float) match/den<<"\t";
             m_match+=match;
             m_den+=den;
             hist_n_match->Fill(match);
             hist_efficiency_inB->Fill((float) match/den);
             match=0;

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


       }// end of if(isBcheck)

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
   std::cout<<"\n In DAOD_selector::Terminate"<<std::endl;
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.

   if(selections){

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

     hist_pt_inB->Write();
     hist_eta_inB->Write();
     hist_phi_inB->Write();
     hist_E_inB->Write();
     hist_Bjet_origin->Write();

     hist_pt_inC->Write();
     hist_eta_inC->Write();
     hist_phi_inC->Write();
     hist_E_inC->Write();
     hist_Cjet_origin->Write();

     hist_pt_l->Write();
     hist_eta_l->Write();
     hist_phi_l->Write();
     hist_E_l->Write();
     hist_ljet_origin->Write();
   }


   if(discriminants){

//     hist_jet_IP2_inB->Write();
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

     hist_ip2d_llr_exB->Write();
     hist_ip2d_llr_exC->Write();
     hist_ip2d_llr_l->Write();
     hist_ip3d_llr_exB->Write();
     hist_ip3d_llr_exC->Write();
     hist_ip3d_llr_l->Write();
     hist_rnnip_llr_exB->Write();
     hist_rnnip_llr_exC->Write();
     hist_rnnip_llr_l->Write();
     hist_sv1_llr_exB->Write();
     hist_sv1_llr_exC->Write();
     hist_sv1_llr_l->Write();
     hist_jf_llr_exB->Write();
     hist_jf_llr_exC->Write();
     hist_jf_llr_l->Write();
     hist_dl1_exB->Write();
     hist_dl1_exC->Write();
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

      hist_trk_pT_inB->Write();
      hist_trk_eta_inB->Write();
      hist_trk_phi_inB->Write();
      hist_trk_Deta_inB->Write();
      hist_trk_Dphi_inB->Write();
      hist_trk_Dphi_Deta_inB->Write();
      hist_trk_DR_inB->Write();
      hist_trk_pT_DR_inB->Write();
      hist_trk_pT_jet_DR_inB->Write();
      hist_trk_pdgId_inB->Write();
      hist_trk_origin_inB->Write();

      hist_trk_d0_inB->Write();
      hist_trk_z0sinth_inB->Write();
      hist_trk_d0sig_inB->Write();
      hist_trk_z0sinthsig_inB->Write();
      hist_trk_d0sig_origin_inB->Write();
      hist_trk_z0sinthsig_origin_inB->Write();
      hist_trk_logpTfrac_origin_inB->Write();
      hist_trk_logDR_origin_inB->Write();
      hist_trk_IBLhits_origin_inB->Write();
      hist_trk_NextToIBLhits_origin_inB->Write();
      hist_trk_sharedIBLhits_origin_inB->Write();
      hist_trk_splitIBLhits_origin_inB->Write();
      hist_trk_nPixhits_origin_inB->Write();
      hist_trk_sharedPixhits_origin_inB->Write();
      hist_trk_splitPixhits_origin_inB->Write();
      hist_trk_nSCThits_origin_inB->Write();
      hist_trk_sharedSCThits_origin_inB->Write();

      hist_trk_d0sig_JF_inB->Write();
      hist_trk_z0sinthsig_JF_inB->Write();
      hist_trk_d0sig_origin_JF_inB->Write();
      hist_trk_z0sinthsig_origin_JF_inB->Write();
      hist_trk_logpTfrac_origin_JF_inB->Write();
      hist_trk_logDR_origin_JF_inB->Write();
      hist_trk_IBLhits_origin_JF_inB->Write();
      hist_trk_NextToIBLhits_origin_JF_inB->Write();
      hist_trk_sharedIBLhits_origin_JF_inB->Write();
      hist_trk_splitIBLhits_origin_JF_inB->Write();
      hist_trk_nPixhits_origin_JF_inB->Write();
      hist_trk_sharedPixhits_origin_JF_inB->Write();
      hist_trk_splitPixhits_origin_JF_inB->Write();
      hist_trk_nSCThits_origin_JF_inB->Write();
      hist_trk_sharedSCThits_origin_JF_inB->Write();

      hist_trk_d0sig_SV1_inB->Write();
      hist_trk_z0sinthsig_SV1_inB->Write();
      hist_trk_d0sig_origin_SV1_inB->Write();
      hist_trk_z0sinthsig_origin_SV1_inB->Write();
      hist_trk_logpTfrac_origin_SV1_inB->Write();
      hist_trk_logDR_origin_SV1_inB->Write();
      hist_trk_IBLhits_origin_SV1_inB->Write();
      hist_trk_NextToIBLhits_origin_SV1_inB->Write();
      hist_trk_sharedIBLhits_origin_SV1_inB->Write();
      hist_trk_splitIBLhits_origin_SV1_inB->Write();
      hist_trk_nPixhits_origin_SV1_inB->Write();
      hist_trk_sharedPixhits_origin_SV1_inB->Write();
      hist_trk_splitPixhits_origin_SV1_inB->Write();
      hist_trk_nSCThits_origin_SV1_inB->Write();
      hist_trk_sharedSCThits_origin_SV1_inB->Write();

      hist_trk_d0sig_SV0_inB->Write();
      hist_trk_z0sinthsig_SV0_inB->Write();
      hist_trk_d0sig_origin_SV0_inB->Write();
      hist_trk_z0sinthsig_origin_SV0_inB->Write();
      hist_trk_logpTfrac_origin_SV0_inB->Write();
      hist_trk_logDR_origin_SV0_inB->Write();
      hist_trk_IBLhits_origin_SV0_inB->Write();
      hist_trk_NextToIBLhits_origin_SV0_inB->Write();
      hist_trk_sharedIBLhits_origin_SV0_inB->Write();
      hist_trk_splitIBLhits_origin_SV0_inB->Write();
      hist_trk_nPixhits_origin_SV0_inB->Write();
      hist_trk_sharedPixhits_origin_SV0_inB->Write();
      hist_trk_splitPixhits_origin_SV0_inB->Write();
      hist_trk_nSCThits_origin_SV0_inB->Write();
      hist_trk_sharedSCThits_origin_SV0_inB->Write();

      hist_trk_d0sig_IP3D_inB->Write();
      hist_trk_z0sinthsig_IP3D_inB->Write();
      hist_trk_d0sig_origin_IP3D_inB->Write();
      hist_trk_z0sinthsig_origin_IP3D_inB->Write();
      hist_trk_logpTfrac_origin_IP3D_inB->Write();
      hist_trk_logDR_origin_IP3D_inB->Write();
      hist_trk_IBLhits_origin_IP3D_inB->Write();
      hist_trk_NextToIBLhits_origin_IP3D_inB->Write();
      hist_trk_sharedIBLhits_origin_IP3D_inB->Write();
      hist_trk_splitIBLhits_origin_IP3D_inB->Write();
      hist_trk_nPixhits_origin_IP3D_inB->Write();
      hist_trk_sharedPixhits_origin_IP3D_inB->Write();
      hist_trk_splitPixhits_origin_IP3D_inB->Write();
      hist_trk_nSCThits_origin_IP3D_inB->Write();
      hist_trk_sharedSCThits_origin_IP3D_inB->Write();

      hist_trk_d0sig_IP2D_inB->Write();
      hist_trk_z0sinthsig_IP2D_inB->Write();
      hist_trk_d0sig_origin_IP2D_inB->Write();
      hist_trk_z0sinthsig_origin_IP2D_inB->Write();
      hist_trk_logpTfrac_origin_IP2D_inB->Write();
      hist_trk_logDR_origin_IP2D_inB->Write();
      hist_trk_IBLhits_origin_IP2D_inB->Write();
      hist_trk_NextToIBLhits_origin_IP2D_inB->Write();
      hist_trk_sharedIBLhits_origin_IP2D_inB->Write();
      hist_trk_splitIBLhits_origin_IP2D_inB->Write();
      hist_trk_nPixhits_origin_IP2D_inB->Write();
      hist_trk_sharedPixhits_origin_IP2D_inB->Write();
      hist_trk_splitPixhits_origin_IP2D_inB->Write();
      hist_trk_nSCThits_origin_IP2D_inB->Write();
      hist_trk_sharedSCThits_origin_IP2D_inB->Write();

      hist_trk_d0_PUinB->Write();
      hist_trk_d0_BinB->Write();
      hist_trk_d0_CinB->Write();
      hist_trk_d0_FRAGinB->Write();
      hist_trk_d0_GEANTinB->Write();

      hist_child_pT_inB->Write();
      hist_child_eta_inB->Write();
      hist_child_phi_inB->Write();
      hist_child_Deta_inB->Write();
      hist_child_Dphi_inB->Write();
      hist_child_Dphi_Deta_inB->Write();
      hist_child_DR_inB->Write();
      hist_child_pT_DR_inB->Write();
      hist_child_pT_jet_DR_inB->Write();
      hist_child_pdgID_inB->Write();

      hist_child_pi_notD->Write();
      hist_child_K_notD->Write();
      hist_child_pi->Write();
      hist_child_pi_Lxy_inB->Write();
      hist_child_pi_Lxyz_inB->Write();
      hist_child_pi_d0_truth_inB->Write();
      hist_child_pi_z0_truth_inB->Write();
      hist_child_K->Write();
      hist_child_K_Lxy_inB->Write();
      hist_child_K_Lxyz_inB->Write();
      hist_child_K_d0_truth_inB->Write();
      hist_child_K_z0_truth_inB->Write();
      hist_child_mu->Write();
      hist_child_mu_Lxy_inB->Write();
      hist_child_mu_Lxyz_inB->Write();
      hist_child_mu_d0_truth_inB->Write();
      hist_child_mu_z0_truth_inB->Write();
      hist_child_p->Write();
      hist_child_p_Lxy_inB->Write();
      hist_child_p_Lxyz_inB->Write();
      hist_child_p_d0_truth_inB->Write();
      hist_child_p_z0_truth_inB->Write();
      hist_child_e->Write();
      hist_child_e_Lxy_inB->Write();
      hist_child_e_Lxyz_inB->Write();
      hist_child_e_d0_truth_inB->Write();
      hist_child_e_z0_truth_inB->Write();

      hist_child_Lxy_inB->Write();
      hist_child_Lxyz_inB->Write();
//      hist_child_decay_IP->Write();
//      hist_child_nodecay_IP->Write();
//      hist_child_linear_IP->Write();
      hist_child_d0_truth->Write();
      hist_child_d0->Write();
      hist_child_d0_pT->Write();
      hist_child_z0sinth_inB->Write();
      hist_pT_vs_R0_ratio_inB->Write();

      if(origin_selection){

        hist_matched_origin_pT_inB->Write();
        hist_matched_origin_eta_inB->Write();
        hist_matched_origin_phi_inB->Write();
        hist_matched_origin_Deta_inB->Write();
        hist_matched_origin_Dphi_inB->Write();
        hist_matched_origin_Dphi_Deta_inB->Write();
        hist_matched_origin_DR_inB->Write();
        hist_matched_origin_pT_DR_inB->Write();
        hist_matched_origin_pT_jet_DR_inB->Write();
        hist_matched_origin_pdgId_inB->Write();
        hist_matched_origin_origin_inB->Write();
//        hist_matched_origin_Lxy_inB->Write();
//        hist_matched_origin_Lxyz_inB->Write();
        hist_matched_origin_d0_inB->Write();

      }

      if(geometric_selection){

        hist_efficiency_inB->Write();
        hist_n_trk->Write();
        hist_n_child->Write();
        hist_n_match->Write();

        hist_matched_pT_inB->Write();
        hist_matched_eta_inB->Write();
        hist_matched_phi_inB->Write();
        hist_matched_Deta_inB->Write();
        hist_matched_Dphi_inB->Write();
        hist_matched_Dphi_Deta_inB->Write();
        hist_matched_DR_inB->Write();
        hist_matched_pT_DR_inB->Write();
        hist_matched_pT_jet_DR_inB->Write();
        hist_matched_pdgId_inB->Write();

        hist_matched_child_pi_notD->Write();
        hist_matched_child_K_notD->Write();
        hist_matched_child_pi->Write();
        hist_matched_child_pi_Lxy_inB->Write();
        hist_matched_child_pi_Lxyz_inB->Write();
        hist_matched_child_pi_d0_truth_inB->Write();
        hist_matched_child_pi_z0_truth_inB->Write();
        hist_matched_child_K->Write();
        hist_matched_child_K_Lxy_inB->Write();
        hist_matched_child_K_Lxyz_inB->Write();
        hist_matched_child_K_d0_truth_inB->Write();
        hist_matched_child_K_z0_truth_inB->Write();
        hist_matched_child_mu->Write();
        hist_matched_child_mu_Lxy_inB->Write();
        hist_matched_child_mu_Lxyz_inB->Write();
        hist_matched_child_mu_d0_truth_inB->Write();
        hist_matched_child_mu_z0_truth_inB->Write();
        hist_matched_child_p->Write();
        hist_matched_child_p_Lxy_inB->Write();
        hist_matched_child_p_Lxyz_inB->Write();
        hist_matched_child_p_d0_truth_inB->Write();
        hist_matched_child_p_z0_truth_inB->Write();
        hist_matched_child_e->Write();
        hist_matched_child_e_Lxy_inB->Write();
        hist_matched_child_e_Lxyz_inB->Write();
        hist_matched_child_e_d0_truth_inB->Write();
        hist_matched_child_e_z0_truth_inB->Write();

        hist_matched_origin_inB->Write();
        hist_matched_pT_child_pTfraction_inB->Write();
        hist_matched_DR_trk_inB->Write();
        hist_matched_DR_trk_pTfraction->Write();
        hist_matched_Lxy_inB->Write();
        hist_matched_Lxyz_inB->Write();
        hist_matched_d0_inB->Write();
  /*
        hist_nomatched_pT_inB->Write();
        hist_nomatched_eta_inB->Write();
        hist_nomatched_phi_inB->Write();
        hist_nomatched_Deta_inB->Write();
        hist_nomatched_Dphi_inB->Write();
        hist_nomatched_Dphi_Deta_inB->Write();
        hist_nomatched_DR_inB->Write();
        hist_nomatched_pT_DR_inB->Write();
        hist_nomatched_pT_jet_DR_inB->Write();
        hist_nomatched_pdgId_inB->Write();
  */
  /*
        hist_single_matched_pT_inB->Write();
        hist_single_matched_eta_inB->Write();
        hist_single_matched_phi_inB->Write();
        hist_single_matched_Deta_inB->Write();
        hist_single_matched_Dphi_inB->Write();
        hist_single_matched_Dphi_Deta_inB->Write();
        hist_single_matched_DR_inB->Write();
        hist_single_matched_pT_DR_inB->Write();
        hist_single_matched_pT_jet_DR_inB->Write();
        hist_single_matched_pdgId_inB->Write();
        hist_single_matched_origin_inB->Write();
        hist_single_matched_pT_child_pTfraction_inB->Write();
        hist_single_matched_DR_trk_inB->Write();
        hist_single_matched_DR_trk_pTfraction->Write();
        hist_single_matched_d0_inB->Write();
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

    if(selection_alg){
      std::cout<<"\nTRUTH LEVEL COMPOSITION\n# of jets: "<< m_njets << " of which:\n" << "b jets: " << m_nBjets << ", c jets: " << m_nCjets << ", light jets: "<< m_nljets << "\nB-C overlap: " << m_nJetBCoverlap<<" ("<< (float) 100*m_nJetBCoverlap/m_nBjets<<" %)\n";

      std::cout<<"\nB-C JET CUTS\n"<<1e-3*jet_pT_infcut<<" < Jet pT [GeV] < "<<1e-3*jet_pT_supcut<<"\nJet JVT > "<<jet_JVT_cut<<" (for jet pT in [20, 60] GeV, |eta| < 2.4)\nJet |eta| < "<<jet_eta_cut<<", OverlapRemoval (e+-/mu+-), NOT isBad\n";

      std::cout<<" CutFlow: passing PtMin           "<<m_njets_2_passPtMin<<std::endl;
      std::cout<<" CutFlow: passing PtMax           "<<m_njets_2_passPtMax<<std::endl;
      std::cout<<" CutFlow: passing Eta Range       "<<m_njets_2_passEtaRange<<std::endl;
      std::cout<<" CutFlow: passing BadMedium       "<<m_njets_2_passBadMedium<<std::endl;
      std::cout<<" CutFlow: passing OR              "<<m_njets_2_passOR<<std::endl;
      std::cout<<" CutFlow: passing ORmu            "<<m_njets_2_passORmu<<std::endl;
      std::cout<<" CutFlow: passing JVT             "<<m_njets_2_passJVT<<std::endl;


      std::cout<<"\n# of jets: "<< m_njets_2 << " of which:\n" << "b jets: " << m_nBjets_2 << ", c jets: " << m_nCjets_2 << ", light jets: "<< m_nljets_2 << "\nB-C overlap: " <<m_nJetBCoverlap_postJetSel<<" ("<< (float) 100*m_nJetBCoverlap_postJetSel/m_nBjets_2<<" %)\n";


      std::cout<<"\nB-C HADRONS CUTS\n"<<"b-c Hadr pT > "<<1e-3*m_pT_bcH_truth_cut<<" GeV\n"<<"b-c Hadr DR < "<<m_DR_bcH_truth_cut<<"\n";
      std::cout<<"\nb jets: " << m_nBcheck << ", c jets: "<< m_nCcheck << ", light jets "<<m_nlcheck<<"\nB-C overlap: "<<ov_check<<" ("<< (float) 100*ov_check/m_nBcheck<<" %)\n";
      std::cout<<"\nTRACK CUTS\ntrk pT > "<<1e-3*trk_pT_cut<<" GeV\n"<<"trk |eta| < "<<trk_eta_cut<<"\n"<<"trk |d0| < "<<trk_d0_cut<<" mm"<<"\n";
      std::cout<<"\n(JF_ntrk, SV1_ntrk, SV0_ntrk, IP2D_ntrk, IP3D_ntrk)\n"<<JF_ntrk<<", \t"<<SV1_ntrk<<", \t"<<SV0_ntrk<<", \t"<<IP2D_ntrk<<", \t"<<IP3D_ntrk<<"\n";
      if(origin_selection){
        std::cout<<"\nORIGIN SELECTION\n";
        std::cout<< std::fixed << std::setprecision(3) << "\n";
        std::cout<< "number of tracks BC:\t" << n_trk_B+n_trk_C << "\t(" << (float) 100*(n_trk_B+n_trk_C)/n_trk_pT_cut << " %)\n";
        std::cout<< "number of tracks PU:\t" << n_trk_PU_pT_cut << "\t(" << (float) 100*n_trk_PU_pT_cut/n_trk_pT_cut << " %)\n";
        std::cout<< "number of tracks FRAG:\t" << n_trk_FRAG_pT_cut << "\t(" << (float) 100*n_trk_FRAG_pT_cut/n_trk_pT_cut << " %)\n";
        std::cout<< "number of tracks GEANT:\t" << n_trk_GEANT_pT_cut << "\t(" << (float) 100*n_trk_GEANT_pT_cut/n_trk_pT_cut << " %)\n";
        std::cout<< "\nmatches:\t\t"<< n_trk_B+n_trk_C <<"\t("<<mm+mm1_ex<<")"<<"\n";
//        std::cout<< "no matches:\t\t"<< m_den-(n_trk_B+n_trk_C) <<"\n";
        std::cout<< "\naverage efficiency:\t" << (float) (n_trk_B+n_trk_C)/m_den << "\n";
      }

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
        std::cout<< "overlap/origin_matches:\t\t" << (float) m_match_overlap/(n_trk_B+n_trk_C) <<"\n";
        std::cout<< "overlap/geometric_matches:\t" << (float) m_match_overlap/m_match <<"\n";
        std::cout<<"origin - not-geometry: "<<mm1_ex<<"\n";
        std::cout<<"geometry - not-origin: "<<mm2_ex<<"\n";
        std::cout<<"geometry - not-origin from PU: "<<m_GeomNOr_PU<<"\t("<<(float) 100*m_GeomNOr_PU/mm2_ex<<" %)"<<"\n";
        std::cout<<"geometry - not-origin from Frag: "<<m_GeomNOr_F<<"\t("<<(float) 100*m_GeomNOr_F/mm2_ex<<" %)"<<"\n";
        std::cout<<"geometry - not-origin from Geant: "<<m_GeomNOr_G<<"\t("<<(float) 100*m_GeomNOr_G/mm2_ex<<" %)"<<"\n";
      }
    }

   std::cout<<"Out of DAOD_selector::Terminate"<<std::endl;
   return;
}
void DAOD_selector::openOutputFile(std::string fileNameStringID)
{
   std::cout<<"\n In DAOD_selector::openOutputFile"<<std::endl;
   if(lxplus && !debug){
     if(retagT){
       std::cout<<"RETAG TRUE\n";
       if(cut){
         std::cout<<"CUT\n";
         file = new TFile(("../output_files/lxplus_output_doRetagT_cut"+fileNameStringID+".root").c_str(),"RECREATE");
       }
       if(!cut){
         std::cout<<"NO CUT\n";
         file = new TFile(("../output_files/lxplus_output_doRetagT_nocut"+fileNameStringID+".root").c_str(),"RECREATE");
       }
     }
     if(!retagT){
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
     file = new TFile(("debug_output"+fileNameStringID+".root").c_str(),"RECREATE");
   }
   if(!lxplus && !debug){
     if(retagT){
       std::cout<<"RETAG TRUE\n";
       if(cut){
         std::cout<<"CUT\n";
         file = new TFile(("../output_files/output_doRetagT_cut"+fileNameStringID+".root").c_str(),"RECREATE");
       }
       if(!cut){
         std::cout<<"NO CUT\n";
         file = new TFile(("../output_files/output_doRetagT_nocut"+fileNameStringID+".root").c_str(),"RECREATE");
       }
     }
     if(!retagT){
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
   if(selections){

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

     hist_pt_inB = new TH1F("pT_inB", "inB", 100, 0., 1000.);
     hist_eta_inB = new TH1F("eta_inB", "inB", 100, -5., 5.);
     hist_phi_inB = new TH1F("phi_inB", "inB", 100, -4.,4.);
     hist_E_inB = new TH1F("E_inB", "inB", 100, 0., 10e5);
     hist_Bjet_origin = new TH1F("Bjet_origin","Bjet_origin",5,-1,4);

     hist_pt_inC = new TH1F("pT_inC", "inC", 100, 0., 1000.);
     hist_eta_inC = new TH1F("eta_inC", "inC", 100, -5., 5.);
     hist_phi_inC = new TH1F("phi_inC", "inC", 100, -4.,4.);
     hist_E_inC = new TH1F("E_inC", "inC", 100, 0., 10e5);
     hist_Cjet_origin = new TH1F("Cjet_origin","Cjet_origin",5,-1,4);

     hist_pt_l = new TH1F("pT_l", "l", 100, 0., 1000.);
     hist_eta_l = new TH1F("eta_l", "l", 100, -5., 5.);
     hist_phi_l = new TH1F("phi_l", "l", 100, -4.,4.);
     hist_E_l = new TH1F("E_l", "l", 100, 0., 10e5);
     hist_ljet_origin = new TH1F("ljet_origin","ljet_origin",5,-1,4);
   }
}

void DAOD_selector::bookHistosForDiscriminants()
{

   if(discriminants){

     hist_ip2d_llr_l = new TH1F("ip2d_llr_l","ip2d_llr_l",280., -20., 50.);
     hist_ip2d_llr_exB = new TH1F("ip2d_llr_exB", "ip2d_llr_exB", 280., -20., 50.);
     hist_ip2d_llr_exC = new TH1F("ip2d_llr_exC", "ip2d_llr_exC", 280., -20., 50.);

     hist_ip3d_llr_l = new TH1F("ip3d_llr_l","ip3d_llr_l",280., -20., 50.);
     hist_ip3d_llr_exB = new TH1F("ip3d_llr_exB", "ip3d_llr_exB", 280., -20., 50.);
     hist_ip3d_llr_exC = new TH1F("ip3d_llr_exC", "ip3d_llr_exC", 280., -20., 50.);

     hist_rnnip_llr_l = new TH1F("rnnip_llr_l","rnnip_llr_l",280., -20., 50.);
     hist_rnnip_llr_exB = new TH1F("rnnip_llr_exB", "rnnip_llr_exB", 280., -20., 50.);
     hist_rnnip_llr_exC = new TH1F("rnnip_llr_exC", "rnnip_llr_exC", 280., -20., 50.);

     hist_sv1_llr_l = new TH1F("sv1_llr_l","sv1_llr_l",280., -20., 50.);
     hist_sv1_llr_exB = new TH1F("sv1_llr_exB", "sv1_llr_exB", 280., -20., 50.);
     hist_sv1_llr_exC = new TH1F("sv1_llr_exC", "sv1_llr_exC", 280., -20., 50.);

     hist_jf_llr_l = new TH1F("jf_llr_l","jf_llr_l",280., -20., 50.);
     hist_jf_llr_exB = new TH1F("jf_llr_exB", "jf_llr_exB", 280., -20., 50.);
     hist_jf_llr_exC = new TH1F("jf_llr_exC", "jf_llr_exC", 280., -20., 50.);

     hist_dl1_l = new TH1F("DL1_light","DL1_light",260,-8.,18.);
     hist_dl1_exC = new TH1F("DL1_exC","DL1_exC",260,-8.,18.);
     hist_dl1_exB = new TH1F("DL1_exB","DL1_exB",260,-8.,18.);

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

     hist_efficiency_inB = new TH1F("efficiency","efficiency",120,-0.1,1.1);
     hist_n_trk = new TH1D("n_trk","n_trk",100,0.,100);
     hist_n_child = new TH1D("n_child","n_child",100,0.,100);
     hist_n_match = new TH1D("n_match","n_match",100,0.,100);

     hist_trk_pT_inB = new TH1F("trk_pT_inB", "trk_pT_inB", 300, 0., 150.);
     hist_trk_eta_inB = new TH1F("trk_eta_inB","trk_eta_inB",260,-2.6,2.6);
     hist_trk_phi_inB = new TH1F("trk_phi_inB","trk_phi_inB",200,-4.,4.);
     hist_trk_Deta_inB = new TH1F("trk_Deta_inB","trk_Deta_inB",200,-1.,1.);
     hist_trk_Dphi_inB = new TH1F("trk_Dphi_inB","trk_Dphi_inB",200,-1.,1.);
     hist_trk_Dphi_Deta_inB = new TH2F("trk_Dphi_Deta","trk_Dphi_Deta", 100, -1.,1., 100, -1.,1.);
     hist_trk_DR_inB = new TH1F("trk_DR_inB","trk_DR_inB",200,0.,2.);
     hist_trk_pT_DR_inB = new TH2F("trk_pT_DR_inB","trk_pT_DR_inB",300,0.,150.,200,0.,2.);
     hist_trk_pT_jet_DR_inB = new TH2F("trk_pT_jet_DR_inB","trk_pT_jet_DR_inB",1000,0.,500.,200,0.,2.);
     hist_trk_pdgId_inB = new TH1F("trk_pdgID_inB","trk_pdgID_inB",200000,-100000,100000);
     hist_trk_origin_inB = new TH1F("trk_origin_inB","trk_origin_inB",5,-1,4);

     hist_trk_d0_inB = new TH1F("trk_d0_inB","trk_d0_inB",300,-15.,15.);
     hist_trk_z0sinth_inB = new TH1F("trk_z0sinth_inB","trk_z0sinth_inB",300,-15.,15.);
     hist_trk_d0sig_inB = new TH1F("trk_d0sig_inB","trk_d0sig_inB",300,-15.,15.);
     hist_trk_z0sinthsig_inB = new TH1F("trk_z0sinthsig_inB","trk_z0sinthsig_inB",300,-15.,15.);
     hist_trk_d0sig_origin_inB = new TH2F("trk_d0sig_origin_inB","trk_d0sig_origin_inB",300,-15.,15.,5,-1,4);
     hist_trk_z0sinthsig_origin_inB = new TH2F("trk_z0sinthsig_origin_inB","trk_z0sinthsig_origin_inB",300,-15.,15.,5,-1,4);
     hist_trk_logpTfrac_origin_inB = new TH2F("trk_logpTfrac_origin_inB","trk_logpTfrac_origin_inB",300,-10.,2,5,-1,4);
     hist_trk_logDR_origin_inB = new TH2F("trk_logDR_origin_inB","trk_logDR_origin_inB",200,-10.,2.,5,-1,4);
     hist_trk_IBLhits_origin_inB = new TH2F("trk_IBLhits_origin_inB","trk_IBLhits_origin_inB",12,0.,12.,5,-1,4);
     hist_trk_NextToIBLhits_origin_inB = new TH2F("trk_NextToIBLhits_origin_inB","trk_NextToIBLhits_origin_inB",12,0.,12.,5,-1,4);
     hist_trk_sharedIBLhits_origin_inB = new TH2F("trk_sharedIBLhits_origin_inB","trk_sharedIBLhits_origin_inB",12,0.,12.,5,-1,4);
     hist_trk_splitIBLhits_origin_inB = new TH2F("trk_splitIBLhits_origin_inB","trk_splitIBLhits_origin_inB",12,0.,12.,5,-1,4);
     hist_trk_nPixhits_origin_inB = new TH2F("trk_nPixhits_origin_inB","trk_nPixhits_origin_inB",12,0.,12.,5,-1,4);
     hist_trk_sharedPixhits_origin_inB = new TH2F("trk_sharedPixhits_origin_inB","trk_sharedPixhits_origin_inB",12,0.,12.,5,-1,4);
     hist_trk_splitPixhits_origin_inB = new TH2F("trk_splitPixhits_origin_inB","trk_splitPixhits_origin_inB",12,0.,12.,5,-1,4);
     hist_trk_nSCThits_origin_inB = new TH2F("trk_nSCThits_origin_inB","trk_nSCThits_origin_inB",24,0.,24.,5,-1,4);
     hist_trk_sharedSCThits_origin_inB = new TH2F("trk_sharedSCThits_origin_inB","trk_sharedSCThits_origin_inB",24,0.,24.,5,-1,4);

     hist_trk_d0sig_JF_inB = new TH1F("trk_d0sig_JF_inB","trk_d0sig_JF_inB",300,-15.,15.);
     hist_trk_z0sinthsig_JF_inB = new TH1F("trk_z0sinthsig_JF_inB","trk_z0sinthsig_JF_inB",300,-15.,15.);
     hist_trk_d0sig_origin_JF_inB = new TH2F("trk_d0sig_origin_JF_inB","trk_d0sig_origin_JF_inB",300,-15.,15.,5,-1,4);
     hist_trk_z0sinthsig_origin_JF_inB = new TH2F("trk_z0sinthsig_origin_JF_inB","trk_z0sinthsig_origin_JF_inB",300,-15.,15.,5,-1,4);
     hist_trk_logpTfrac_origin_JF_inB = new TH2F("trk_logpTfrac_origin_JF_inB","trk_logpTfrac_origin_JF_inB",300,-10.,2,5,-1,4);
     hist_trk_logDR_origin_JF_inB = new TH2F("trk_logDR_origin_JF_inB","trk_logDR_origin_JF_inB",200,-10.,2.,5,-1,4);
     hist_trk_IBLhits_origin_JF_inB = new TH2F("trk_IBLhits_origin_JF_inB","trk_IBLhits_origin_JF_inB",12,0.,12.,5,-1,4);
     hist_trk_NextToIBLhits_origin_JF_inB = new TH2F("trk_NextToIBLhits_origin_JF_inB","trk_NextToIBLhits_origin_JF_inB",12,0.,12.,5,-1,4);
     hist_trk_sharedIBLhits_origin_JF_inB = new TH2F("trk_sharedIBLhits_origin_JF_inB","trk_sharedIBLhits_origin_JF_inB",12,0.,12.,5,-1,4);
     hist_trk_splitIBLhits_origin_JF_inB = new TH2F("trk_splitIBLhits_origin_JF_inB","trk_splitIBLhits_origin_JF_inB",12,0.,12.,5,-1,4);
     hist_trk_nPixhits_origin_JF_inB = new TH2F("trk_nPixhits_origin_JF_inB","trk_nPixhits_origin_JF_inB",12,0.,12.,5,-1,4);
     hist_trk_sharedPixhits_origin_JF_inB = new TH2F("trk_sharedPixhits_origin_JF_inB","trk_sharedPixhits_origin_JF_inB",12,0.,12.,5,-1,4);
     hist_trk_splitPixhits_origin_JF_inB = new TH2F("trk_splitPixhits_origin_JF_inB","trk_splitPixhits_origin_JF_inB",12,0.,12.,5,-1,4);
     hist_trk_nSCThits_origin_JF_inB = new TH2F("trk_nSCThits_origin_JF_inB","trk_nSCThits_origin_JF_inB",24,0.,24.,5,-1,4);
     hist_trk_sharedSCThits_origin_JF_inB = new TH2F("trk_sharedSCThits_origin_JF_inB","trk_sharedSCThits_origin_JF_inB",24,0.,24.,5,-1,4);

     hist_trk_d0sig_SV1_inB = new TH1F("trk_d0sig_SV1_inB","trk_d0sig_SV1_inB",300,-15.,15.);
     hist_trk_z0sinthsig_SV1_inB = new TH1F("trk_z0sinthsig_SV1_inB","trk_z0sinthsig_SV1_inB",300,-15.,15.);
     hist_trk_d0sig_origin_SV1_inB = new TH2F("trk_d0sig_origin_SV1_inB","trk_d0sig_origin_SV1_inB",300,-15.,15.,5,-1,4);
     hist_trk_z0sinthsig_origin_SV1_inB = new TH2F("trk_z0sinthsig_origin_SV1_inB","trk_z0sinthsig_origin_SV1_inB",300,-15.,15.,5,-1,4);
     hist_trk_logpTfrac_origin_SV1_inB = new TH2F("trk_logpTfrac_origin_SV1_inB","trk_logpTfrac_origin_SV1_inB",300,-10.,2,5,-1,4);
     hist_trk_logDR_origin_SV1_inB = new TH2F("trk_logDR_origin_SV1_inB","trk_logDR_origin_SV1_inB",200,-10.,2.,5,-1,4);
     hist_trk_IBLhits_origin_SV1_inB = new TH2F("trk_IBLhits_origin_SV1_inB","trk_IBLhits_origin_SV1_inB",12,0.,12.,5,-1,4);
     hist_trk_NextToIBLhits_origin_SV1_inB = new TH2F("trk_NextToIBLhits_origin_SV1_inB","trk_NextToIBLhits_origin_SV1_inB",12,0.,12.,5,-1,4);
     hist_trk_sharedIBLhits_origin_SV1_inB = new TH2F("trk_sharedIBLhits_origin_SV1_inB","trk_sharedIBLhits_origin_SV1_inB",12,0.,12.,5,-1,4);
     hist_trk_splitIBLhits_origin_SV1_inB = new TH2F("trk_splitIBLhits_origin_SV1_inB","trk_splitIBLhits_origin_SV1_inB",12,0.,12.,5,-1,4);
     hist_trk_nPixhits_origin_SV1_inB = new TH2F("trk_nPixhits_origin_SV1_inB","trk_nPixhits_origin_SV1_inB",12,0.,12.,5,-1,4);
     hist_trk_sharedPixhits_origin_SV1_inB = new TH2F("trk_sharedPixhits_origin_SV1_inB","trk_sharedPixhits_origin_SV1_inB",12,0.,12.,5,-1,4);
     hist_trk_splitPixhits_origin_SV1_inB = new TH2F("trk_splitPixhits_origin_SV1_inB","trk_splitPixhits_origin_SV1_inB",12,0.,12.,5,-1,4);
     hist_trk_nSCThits_origin_SV1_inB = new TH2F("trk_nSCThits_origin_SV1_inB","trk_nSCThits_origin_SV1_inB",24,0.,24.,5,-1,4);
     hist_trk_sharedSCThits_origin_SV1_inB = new TH2F("trk_sharedSCThits_origin_SV1_inB","trk_sharedSCThits_origin_SV1_inB",24,0.,24.,5,-1,4);

     hist_trk_d0sig_SV0_inB = new TH1F("trk_d0sig_SV0_inB","trk_d0sig_SV0_inB",300,-15.,15.);
     hist_trk_z0sinthsig_SV0_inB = new TH1F("trk_z0sinthsig_SV0_inB","trk_z0sinthsig_SV0_inB",300,-15.,15.);
     hist_trk_d0sig_origin_SV0_inB = new TH2F("trk_d0sig_origin_SV0_inB","trk_d0sig_origin_SV0_inB",300,-15.,15.,5,-1,4);
     hist_trk_z0sinthsig_origin_SV0_inB = new TH2F("trk_z0sinthsig_origin_SV0_inB","trk_z0sinthsig_origin_SV0_inB",300,-15.,15.,5,-1,4);
     hist_trk_logpTfrac_origin_SV0_inB = new TH2F("trk_logpTfrac_origin_SV0_inB","trk_logpTfrac_origin_SV0_inB",300,-10.,2,5,-1,4);
     hist_trk_logDR_origin_SV0_inB = new TH2F("trk_logDR_origin_SV0_inB","trk_logDR_origin_SV0_inB",200,-10.,2.,5,-1,4);
     hist_trk_IBLhits_origin_SV0_inB = new TH2F("trk_IBLhits_origin_SV0_inB","trk_IBLhits_origin_SV0_inB",12,0.,12.,5,-1,4);
     hist_trk_NextToIBLhits_origin_SV0_inB = new TH2F("trk_NextToIBLhits_origin_SV0_inB","trk_NextToIBLhits_origin_SV0_inB",12,0.,12.,5,-1,4);
     hist_trk_sharedIBLhits_origin_SV0_inB = new TH2F("trk_sharedIBLhits_origin_SV0_inB","trk_sharedIBLhits_origin_SV0_inB",12,0.,12.,5,-1,4);
     hist_trk_splitIBLhits_origin_SV0_inB = new TH2F("trk_splitIBLhits_origin_SV0_inB","trk_splitIBLhits_origin_SV0_inB",12,0.,12.,5,-1,4);
     hist_trk_nPixhits_origin_SV0_inB = new TH2F("trk_nPixhits_origin_SV0_inB","trk_nPixhits_origin_SV0_inB",12,0.,12.,5,-1,4);
     hist_trk_sharedPixhits_origin_SV0_inB = new TH2F("trk_sharedPixhits_origin_SV0_inB","trk_sharedPixhits_origin_SV0_inB",12,0.,12.,5,-1,4);
     hist_trk_splitPixhits_origin_SV0_inB = new TH2F("trk_splitPixhits_origin_SV0_inB","trk_splitPixhits_origin_SV0_inB",12,0.,12.,5,-1,4);
     hist_trk_nSCThits_origin_SV0_inB = new TH2F("trk_nSCThits_origin_SV0_inB","trk_nSCThits_origin_SV0_inB",24,0.,24.,5,-1,4);
     hist_trk_sharedSCThits_origin_SV0_inB = new TH2F("trk_sharedSCThits_origin_SV0_inB","trk_sharedSCThits_origin_SV0_inB",24,0.,24.,5,-1,4);

     hist_trk_d0sig_IP3D_inB = new TH1F("trk_d0sig_IP3D_inB","trk_d0sig_IP3D_inB",300,-15.,15.);
     hist_trk_z0sinthsig_IP3D_inB = new TH1F("trk_z0sinthsig_IP3D_inB","trk_z0sinthsig_IP3D_inB",300,-15.,15.);
     hist_trk_d0sig_origin_IP3D_inB = new TH2F("trk_d0sig_origin_IP3D_inB","trk_d0sig_origin_IP3D_inB",300,-15.,15.,5,-1,4);
     hist_trk_z0sinthsig_origin_IP3D_inB = new TH2F("trk_z0sinthsig_origin_IP3D_inB","trk_z0sinthsig_origin_IP3D_inB",300,-15.,15.,5,-1,4);
     hist_trk_logpTfrac_origin_IP3D_inB = new TH2F("trk_logpTfrac_origin_IP3D_inB","trk_logpTfrac_origin_IP3D_inB",300,-10.,2,5,-1,4);
     hist_trk_logDR_origin_IP3D_inB = new TH2F("trk_logDR_origin_IP3D_inB","trk_logDR_origin_IP3D_inB",200,-10.,2.,5,-1,4);
     hist_trk_IBLhits_origin_IP3D_inB = new TH2F("trk_IBLhits_origin_IP3D_inB","trk_IBLhits_origin_IP3D_inB",12,0.,12.,5,-1,4);
     hist_trk_NextToIBLhits_origin_IP3D_inB = new TH2F("trk_NextToIBLhits_origin_IP3D_inB","trk_NextToIBLhits_origin_IP3D_inB",12,0.,12.,5,-1,4);
     hist_trk_sharedIBLhits_origin_IP3D_inB = new TH2F("trk_sharedIBLhits_origin_IP3D_inB","trk_sharedIBLhits_origin_IP3D_inB",12,0.,12.,5,-1,4);
     hist_trk_splitIBLhits_origin_IP3D_inB = new TH2F("trk_splitIBLhits_origin_IP3D_inB","trk_splitIBLhits_origin_IP3D_inB",12,0.,12.,5,-1,4);
     hist_trk_nPixhits_origin_IP3D_inB = new TH2F("trk_nPixhits_origin_IP3D_inB","trk_nPixhits_origin_IP3D_inB",12,0.,12.,5,-1,4);
     hist_trk_sharedPixhits_origin_IP3D_inB = new TH2F("trk_sharedPixhits_origin_IP3D_inB","trk_sharedPixhits_origin_IP3D_inB",12,0.,12.,5,-1,4);
     hist_trk_splitPixhits_origin_IP3D_inB = new TH2F("trk_splitPixhits_origin_IP3D_inB","trk_splitPixhits_origin_IP3D_inB",12,0.,12.,5,-1,4);
     hist_trk_nSCThits_origin_IP3D_inB = new TH2F("trk_nSCThits_origin_IP3D_inB","trk_nSCThits_origin_IP3D_inB",24,0.,24.,5,-1,4);
     hist_trk_sharedSCThits_origin_IP3D_inB = new TH2F("trk_sharedSCThits_origin_IP3D_inB","trk_sharedSCThits_origin_IP3D_inB",24,0.,24.,5,-1,4);

     hist_trk_d0sig_IP2D_inB = new TH1F("trk_d0sig_IP2D_inB","trk_d0sig_IP2D_inB",300,-15.,15.);
     hist_trk_z0sinthsig_IP2D_inB = new TH1F("trk_z0sinthsig_IP2D_inB","trk_z0sinthsig_IP2D_inB",300,-15.,15.);
     hist_trk_d0sig_origin_IP2D_inB = new TH2F("trk_d0sig_origin_IP2D_inB","trk_d0sig_origin_IP2D_inB",300,-15.,15.,5,-1,4);
     hist_trk_z0sinthsig_origin_IP2D_inB = new TH2F("trk_z0sinthsig_origin_IP2D_inB","trk_z0sinthsig_origin_IP2D_inB",300,-15.,15.,5,-1,4);
     hist_trk_logpTfrac_origin_IP2D_inB = new TH2F("trk_logpTfrac_origin_IP2D_inB","trk_logpTfrac_origin_IP2D_inB",300,-10.,2,5,-1,4);
     hist_trk_logDR_origin_IP2D_inB = new TH2F("trk_logDR_origin_IP2D_inB","trk_logDR_origin_IP2D_inB",200,-10.,2.,5,-1,4);
     hist_trk_IBLhits_origin_IP2D_inB = new TH2F("trk_IBLhits_origin_IP2D_inB","trk_IBLhits_origin_IP2D_inB",12,0.,12.,5,-1,4);
     hist_trk_NextToIBLhits_origin_IP2D_inB = new TH2F("trk_NextToIBLhits_origin_IP2D_inB","trk_NextToIBLhits_origin_IP2D_inB",12,0.,12.,5,-1,4);
     hist_trk_sharedIBLhits_origin_IP2D_inB = new TH2F("trk_sharedIBLhits_origin_IP2D_inB","trk_sharedIBLhits_origin_IP2D_inB",12,0.,12.,5,-1,4);
     hist_trk_splitIBLhits_origin_IP2D_inB = new TH2F("trk_splitIBLhits_origin_IP2D_inB","trk_splitIBLhits_origin_IP2D_inB",12,0.,12.,5,-1,4);
     hist_trk_nPixhits_origin_IP2D_inB = new TH2F("trk_nPixhits_origin_IP2D_inB","trk_nPixhits_origin_IP2D_inB",12,0.,12.,5,-1,4);
     hist_trk_sharedPixhits_origin_IP2D_inB = new TH2F("trk_sharedPixhits_origin_IP2D_inB","trk_sharedPixhits_origin_IP2D_inB",12,0.,12.,5,-1,4);
     hist_trk_splitPixhits_origin_IP2D_inB = new TH2F("trk_splitPixhits_origin_IP2D_inB","trk_splitPixhits_origin_IP2D_inB",12,0.,12.,5,-1,4);
     hist_trk_nSCThits_origin_IP2D_inB = new TH2F("trk_nSCThits_origin_IP2D_inB","trk_nSCThits_origin_IP2D_inB",24,0.,24.,5,-1,4);
     hist_trk_sharedSCThits_origin_IP2D_inB = new TH2F("trk_sharedSCThits_origin_IP2D_inB","trk_sharedSCThits_origin_IP2D_inB",24,0.,24.,5,-1,4);

     hist_trk_d0_PUinB = new TH1F("trk_d0_PUinB","trk_d0_PUinB",300,-15.,15.);
     hist_trk_d0_BinB = new TH1F("trk_d0_BinB","trk_d0_BinB",300,-15.,15.);
     hist_trk_d0_CinB = new TH1F("trk_d0_CinB","trk_d0_CinB",300,-15.,15.);
     hist_trk_d0_FRAGinB = new TH1F("trk_d0_FRAGinB","trk_d0_FRAGinB",300,-15.,15.);
     hist_trk_d0_GEANTinB = new TH1F("trk_d0_GEANTinB","trk_d0_GEANTinB",300,-15.,15.);

     hist_child_pT_inB = new TH1F("child_pT_inB", "child_pT_inB", 300, 0., 150.);
     hist_child_eta_inB = new TH1F("child_eta_inB","child_eta_inB",260,-2.6,2.6);
     hist_child_phi_inB = new TH1F("child_phi_inB","child_phi_inB",200,-4.,4.);
     hist_child_Deta_inB = new TH1F("child_Deta_inB","child_Deta_inB",200,-1.,1.);
     hist_child_Dphi_inB = new TH1F("child_Dphi_inB","child_Dphi_inB",200,-1.,1.);
     hist_child_Dphi_Deta_inB = new TH2F("child_Dphi_Deta","child_Dphi_Deta", 100, -1.,1., 100, -1.,1.);
     hist_child_DR_inB = new TH1F("child_DR_inB","child_DR_inB",200,0.,2.);
     hist_child_pT_DR_inB = new TH2F("child_pT_DR_inB","child_pT_DR_inB",300,0.,150.,200,0.,2.);
     hist_child_pT_jet_DR_inB = new TH2F("child_pT_jet_DR_inB","child_pT_jet_DR_inB",1000,0.,500.,200,0.,2.);
     hist_child_pdgID_inB = new TH1F("child_pdgID_inB","child_pdgID_inB",200000,-100000,100000);

     hist_child_pi_notD = new TH1F("child_pi_pT_notD_inB", "child_pi_pT_notD_inB", 300, 0., 150.);
     hist_child_K_notD = new TH1F("child_K_pT_notD_inB", "child_K_pT_notD_inB", 300, 0., 150.);
     hist_child_pi = new TH1F("child_pi_pT_inB", "child_pi_pT_inB", 300, 0., 150.);
     hist_child_pi_Lxy_inB = new TH1F("child_pi_Lxy_inB","child_pi_Lxy_inB",1000,0.,1000.);
     hist_child_pi_Lxyz_inB = new TH1F("child_pi_Lxyz_inB","child_pi_Lxyz_inB",1000,0.,1000.);
     hist_child_pi_d0_truth_inB = new TH1F("child_pi_d0_truth_inB","child_pi_d0_truth_inB",300,-15.,15.);
     hist_child_pi_z0_truth_inB = new TH1F("child_pi_z0_truth_inB","child_pi_z0_truth_inB",1000,-1000.,1000.);
     hist_child_K = new TH1F("child_K_pT_inB", "child_K_pT_inB", 300, 0., 150.);
     hist_child_K_Lxy_inB = new TH1F("child_K_Lxy_inB","child_K_Lxy_inB",1000,0.,1000.);
     hist_child_K_Lxyz_inB = new TH1F("child_K_Lxyz_inB","child_K_Lxyz_inB",1000,0.,1000.);
     hist_child_K_d0_truth_inB = new TH1F("child_K_d0_truth_inB","child_K_d0_truth_inB",300,-15.,15.);
     hist_child_K_z0_truth_inB = new TH1F("child_K_z0_truth_inB","child_K_z0_truth_inB",1000,-1000.,1000.);
     hist_child_mu = new TH1F("child_mu_pT_inB", "child_mu_pT_inB", 300, 0., 150.);
     hist_child_mu_Lxy_inB = new TH1F("child_mu_Lxy_inB","child_mu_Lxy_inB",1000,0.,1000.);
     hist_child_mu_Lxyz_inB = new TH1F("child_mu_Lxyz_inB","child_mu_Lxyz_inB",1000,0.,1000.);
     hist_child_mu_d0_truth_inB = new TH1F("child_mu_d0_truth_inB","child_mu_d0_truth_inB",300,-15.,15.);
     hist_child_mu_z0_truth_inB = new TH1F("child_mu_z0_truth_inB","child_mu_z0_truth_inB",1000,-1000.,1000.);
     hist_child_p = new TH1F("child_p_pT_inB", "child_p_pT_inB", 300, 0., 150.);
     hist_child_p_Lxy_inB = new TH1F("child_p_Lxy_inB","child_p_Lxy_inB",1000,0.,1000.);
     hist_child_p_Lxyz_inB = new TH1F("child_p_Lxyz_inB","child_p_Lxyz_inB",1000,0.,1000.);
     hist_child_p_d0_truth_inB = new TH1F("child_p_d0_truth_inB","child_p_d0_truth_inB",300,-15.,15.);
     hist_child_p_z0_truth_inB = new TH1F("child_p_z0_truth_inB","child_p_z0_truth_inB",1000,-1000.,1000.);
     hist_child_e = new TH1F("child_e_pT_inB", "child_e_pT_inB", 300, 0., 150.);
     hist_child_e_Lxy_inB = new TH1F("child_e_Lxy_inB","child_e_Lxy_inB",1000,0.,1000.);
     hist_child_e_Lxyz_inB = new TH1F("child_e_Lxyz_inB","child_e_Lxyz_inB",1000,0.,1000.);
     hist_child_e_d0_truth_inB = new TH1F("child_e_d0_truth_inB","child_e_d0_truth_inB",300,-15.,15.);
     hist_child_e_z0_truth_inB = new TH1F("child_e_z0_truth_inB","child_e_z0_truth_inB",1000,-1000.,1000.);

     hist_child_Lxy_inB = new TH1F("child_Lxy_inB","child_Lxy_inB",1000,0.,1000.);
     hist_child_Lxyz_inB = new TH1F("child_Lxyz_inB","child_Lxyz_inB",1000,0.,1000.);
//     hist_child_decay_IP = new TH1F("child_decay_IP_inB","child_decay_IP_inB",300,-15.,15.);
//     hist_child_nodecay_IP = new TH1F("child_nodecay_IP_inB","child_nodecay_IP_inB",300,-15.,15.);
//     hist_child_linear_IP = new TH1F("child_linear_IP_inB","child_linear_IP_inB",300,-15.,15.);
     hist_child_d0_truth = new TH1F("child_d0_truth_inB","child_d0_truth_inB",300,-15.,15.);
     hist_child_d0 = new TH1F("child_d0_inB","child_d0_inB",300,-15.,15.);
     hist_child_d0_pT = new TH2F("child_d0_pT_inB","child_pT_d0_inB",300,-15.,15.,300,0.,150.);
     hist_child_z0sinth_inB = new TH1F("child_z0sinth_inB","child_z0sinth_inB",300,-150.,150.);
     hist_pT_vs_R0_ratio_inB = new TH1F("child_pT_vs_R0_ratio_inB","child_pT_vs_R0_ratio_inB",300,-1.2,2.4);

     if(origin_selection){

       hist_matched_origin_pT_inB = new TH1F("matched_origin_trk_pT_inB","matched_origin_trk_pT_inB", 300, 0., 150.);
       hist_matched_origin_eta_inB = new TH1F("matched_origin_trk_eta_inB","matched_origin_trk_eta_inB",260,-2.6,2.6);
       hist_matched_origin_phi_inB = new TH1F("matched_origin_trk_phi_inB","matched_origin_trk_phi_inB",200,-4.,4.);
       hist_matched_origin_Deta_inB = new TH1F("matched_origin_trk_Deta_inB","matched_origin_trk_Deta_inB",200,-1.,1.);
       hist_matched_origin_Dphi_inB = new TH1F("matched_origin_trk_Dphi_inB","matched_origin_trk_Dphi_inB",200,-1.,1.);
       hist_matched_origin_Dphi_Deta_inB = new TH2F("matched_origin_trk_Dphi_Deta","matched_origin_trk_Dphi_Deta", 100, -1.,1., 100, -1.,1.);
       hist_matched_origin_DR_inB = new TH1F("matched_origin_trk_DR_inB","matched_origin_trk_DR_inB",200,0.,2.);
       hist_matched_origin_pT_DR_inB = new TH2F("matched_origin_trk_pT_DR_inB","matched_origin_trk_pT_DR_inB",300,0.,150.,200,0.,2.);
       hist_matched_origin_pT_jet_DR_inB = new TH2F("matched_origin_trk_pT_jet_DR_inB","matched_origin_trk_pT_jet_DR_inB",1000,0.,500.,200,0.,2.);
       hist_matched_origin_pdgId_inB = new TH1F("matched_origin_trk_pdgId_inB","matched_origin_trk_pdgId_inB",200000,-100000,100000);
       hist_matched_origin_origin_inB = new TH1F("matched_origin_trk_origin_inB","matched_origin_trk_origin_inB",7,-2,5);
       hist_matched_origin_d0_inB = new TH1F("matched_origin_trk_d0_inB","matched_origin_trk_d0_inB",300,-15.,15.);
//       hist_matched_origin_Lxy_inB = new TH1F("matched_origin_trk_Lxy_inB","matched_origin_trk_Lxy_inB",300,0.,1000.);
//       hist_matched_origin_Lxyz_inB = new TH1F("matched_origin_trk_Lxyz_inB","matched_origin_trk_Lxyz_inB",300,0.,1000.);
     }

     if(geometric_selection){

       hist_matched_pT_inB = new TH1F("matched_child_pT_inB","matched_child_pT_inB", 300, 0., 150.);
       hist_matched_eta_inB = new TH1F("matched_child_eta_inB","matched_child_eta_inB",260,-2.6,2.6);
       hist_matched_phi_inB = new TH1F("matched_child_phi_inB","matched_child_phi_inB",200,-4.,4.);
       hist_matched_Deta_inB = new TH1F("matched_child_Deta_inB","matched_child_Deta_inB",200,-1.,1.);
       hist_matched_Dphi_inB = new TH1F("matched_child_Dphi_inB","matched_child_Dphi_inB",200,-1.,1.);
       hist_matched_Dphi_Deta_inB = new TH2F("matched_child_Dphi_Deta","matched_child_Dphi_Deta", 100, -1.,1., 100, -1.,1.);
       hist_matched_DR_inB = new TH1F("matched_child_DR_inB","matched_child_DR_inB",200,0.,2.);
       hist_matched_pT_DR_inB = new TH2F("matched_child_pT_DR_inB","matched_child_pT_DR_inB",300,0.,150.,200,0.,2.);
       hist_matched_pT_jet_DR_inB = new TH2F("matched_child_pT_jet_DR_inB","matched_child_pT_jet_DR_inB",1000,0.,500.,200,0.,2.);
       hist_matched_pdgId_inB = new TH1F("matched_child_pdgId_inB","matched_child_pdgId_inB",200000,-100000,100000);

       hist_matched_child_pi_notD = new TH1F("matched_child_pi_pT_notD_inB", "matched_child_pi_pT_notD_inB", 300, 0., 150.);
       hist_matched_child_K_notD = new TH1F("matched_child_K_pT_notD_inB", "matched_child_K_pT_notD_inB", 300, 0., 150.);
       hist_matched_child_pi = new TH1F("matched_child_pi_pT_inB", "matched_child_pi_pT_inB", 300, 0., 150.);
       hist_matched_child_pi_Lxy_inB = new TH1F("matched_child_pi_Lxy_inB","matched_child_pi_Lxy_inB",1000,0.,1000.);
       hist_matched_child_pi_Lxyz_inB = new TH1F("matched_child_pi_Lxyz_inB","matched_child_pi_Lxyz_inB",1000,0.,1000.);
       hist_matched_child_pi_d0_truth_inB = new TH1F("matched_child_pi_d0_truth_inB","matched_child_pi_d0_truth_inB",300,-15.,15.);
       hist_matched_child_pi_z0_truth_inB = new TH1F("matched_child_pi_z0_truth_inB","matched_child_pi_z0_truth_inB",1000,-1000.,1000.);
       hist_matched_child_K = new TH1F("matched_child_K_pT_inB", "matched_child_K_pT_inB", 300, 0., 150.);
       hist_matched_child_K_Lxy_inB = new TH1F("matched_child_K_Lxy_inB","matched_child_K_Lxy_inB",1000,0.,1000.);
       hist_matched_child_K_Lxyz_inB = new TH1F("matched_child_K_Lxyz_inB","matched_child_K_Lxyz_inB",1000,0.,1000.);
       hist_matched_child_K_d0_truth_inB = new TH1F("matched_child_K_d0_truth_inB","matched_child_K_d0_truth_inB",300,-15.,15.);
       hist_matched_child_K_z0_truth_inB = new TH1F("matched_child_K_z0_truth_inB","matched_child_K_z0_truth_inB",1000,-1000.,1000.);
       hist_matched_child_mu = new TH1F("matched_child_mu_pT_inB", "matched_child_mu_pT_inB", 300, 0., 150.);
       hist_matched_child_mu_Lxy_inB = new TH1F("matched_child_mu_Lxy_inB","matched_child_mu_Lxy_inB",1000,0.,1000.);
       hist_matched_child_mu_Lxyz_inB = new TH1F("matched_child_mu_Lxyz_inB","matched_child_mu_Lxyz_inB",1000,0.,1000.);
       hist_matched_child_mu_d0_truth_inB = new TH1F("matched_child_mu_d0_truth_inB","matched_child_mu_d0_truth_inB",300,-15.,15.);
       hist_matched_child_mu_z0_truth_inB = new TH1F("matched_child_mu_z0_truth_inB","matched_child_mu_z0_truth_inB",1000,-1000.,1000.);
       hist_matched_child_p = new TH1F("matched_child_p_pT_inB", "matched_child_p_pT_inB", 300, 0., 150.);
       hist_matched_child_p_Lxy_inB = new TH1F("matched_child_p_Lxy_inB","matched_child_p_Lxy_inB",1000,0.,1000.);
       hist_matched_child_p_Lxyz_inB = new TH1F("matched_child_p_Lxyz_inB","matched_child_p_Lxyz_inB",1000,0.,1000.);
       hist_matched_child_p_d0_truth_inB = new TH1F("matched_child_p_d0_truth_inB","matched_child_p_d0_truth_inB",300,-15.,15.);
       hist_matched_child_p_z0_truth_inB = new TH1F("matched_child_p_z0_truth_inB","matched_child_p_z0_truth_inB",1000,-1000.,1000.);
       hist_matched_child_e = new TH1F("matched_child_e_pT_inB", "matched_child_e_pT_inB", 300, 0., 150.);
       hist_matched_child_e_Lxy_inB = new TH1F("matched_child_e_Lxy_inB","matched_child_e_Lxy_inB",1000,0.,1000.);
       hist_matched_child_e_Lxyz_inB = new TH1F("matched_child_e_Lxyz_inB","matched_child_e_Lxyz_inB",1000,0.,1000.);
       hist_matched_child_e_d0_truth_inB = new TH1F("matched_child_e_d0_truth_inB","matched_child_e_d0_truth_inB",300,-15.,15.);
       hist_matched_child_e_z0_truth_inB = new TH1F("matched_child_e_z0_truth_inB","matched_child_e_z0_truth_inB",1000,-1000.,1000.);


       hist_matched_origin_inB = new TH1F("matched_trk_origin_inB","matched_trk_origin_inB",7,-2,5);
       hist_matched_pT_child_pTfraction_inB = new TH2F("matched_pT_child_pTfraction_inB","matched_pT_child_vs_DpT/pT_child_inB",300,0.,150.,500,-1.,10.);
       hist_matched_DR_trk_inB = new TH1F("matched_DR_trk_inB","matched_DR_trk_inB_inB",500,0.,1.);
       hist_matched_DR_trk_pTfraction = new TH2F("matched_DR_trk_pTfraction_inB","matched_DR_trk_pTfraction_inB",500,0.,1.,500,-1.,10.);
       hist_matched_d0_inB = new TH1F("matched_child_d0_inB","matched_child_d0_inB",300,-15.,15.);
       hist_matched_Lxy_inB = new TH1F("matched_child_Lxy_inB","matched_child_Lxy_inB",1000,0.,1000.);
       hist_matched_Lxyz_inB = new TH1F("matched_child_Lxyz_inB","matched_child_Lxyz_inB",1000,0.,1000.);
  /*
       hist_nomatched_pT_inB = new TH1F("nomatched_child_pT_inB", "nomatched_child_pT_inB", 500, 0., 150.);
       hist_nomatched_eta_inB = new TH1F("nomatched_child_eta_inB","nomatched_child_eta_inB",500,-2.6,2.6);
       hist_nomatched_phi_inB = new TH1F("nomatched_child_phi_inB","nomatched_child_phi_inB",500,-4.,4.);
       hist_nomatched_Deta_inB = new TH1F("nomatched_child_Deta_inB","nomatched_child_Deta_inB",500,-1.,1.);
       hist_nomatched_Dphi_inB = new TH1F("nomatched_child_Dphi_inB","nomatched_child_Dphi_inB",500,-1.,1.);
       hist_nomatched_Dphi_Deta_inB = new TH2F("nomatched_child_Dphi_Deta","nomatched_child_Dphi_Deta", 100, -1.,1., 100, -1.,1.);
       hist_nomatched_DR_inB = new TH1F("nomatched_child_DR_inB","nomatched_child_DR_inB",500,-0.1,2.);
       hist_nomatched_pT_DR_inB = new TH2F("nomatched_child_pT_DR_inB","nomatched_child_pT_DR_inB",80,0.,150.,80,0.,0.6);
       hist_nomatched_pT_jet_DR_inB = new TH2F("nomatched_child_pT_jet_DR_inB","nomatched_child_pT_jet_DR_inB",100,0.,500.,100,0.,2.);
       hist_nomatched_pdgId_inB = new TH1F("nomatched_child_pdgId_inB","nomatched_child_pdgId_inB",2000,-1000,1000);
  */
  /*
       hist_single_matched_pT_inB = new TH1F("single_matched_pT_inB","single_matched_pT_inB", 500, 0., 150.);
       hist_single_matched_eta_inB = new TH1F("single_matched_child_eta_inB","single_matched_child_eta_inB",500,-2.6,2.6);
       hist_single_matched_phi_inB = new TH1F("single_matched_child_phi_inB","single_matched_child_phi_inB",500,-4.,4.);
       hist_single_matched_Deta_inB = new TH1F("single_matched_child_Deta_inB","single_matched_child_Deta_inB",500,-1.,1.);
       hist_single_matched_Dphi_inB = new TH1F("single_matched_child_Dphi_inB","single_matched_child_Dphi_inB",500,-1.,1.);
       hist_single_matched_Dphi_Deta_inB = new TH2F("single_matched_child_Dphi_Deta","single_matched_child_Dphi_Deta", 100, -1.,1., 100, -1.1,1.);
       hist_single_matched_DR_inB = new TH1F("single_matched_DR_inB","single_matched_DR_inB",500,-0.1,2.);
       hist_single_matched_pT_DR_inB = new TH2F("single_matched_pT_DR_inB","single_matched_pT_DR_inB",80,0.,150.,80,0.,0.6);
       hist_single_matched_pT_jet_DR_inB = new TH2F("single_matched_child_pT_jet_DR_inB","single_matched_child_pT_jet_DR_inB",100,0.,500.,100,0.,2.);
       hist_single_matched_pdgId_inB = new TH1F("single_matched_child_pdgId_inB","single_matched_child_pdgId_inB",2000,-1000,1000);
       hist_single_matched_origin_inB = new TH1F("single_matched_trk_origin_inB","single_matched_trk_origin_inB",7,-2,5);
       hist_single_matched_pT_child_pTfraction_inB = new TH2F("single_matched_pT_child_pTfraction_inB","single_matched_pT_child_vs_DpT/pT_child_inB",500,0.,150.,500,-1.1,10.);
       hist_single_matched_DR_trk_inB = new TH1F("single_matched_DR_trk_inB","single_matched_DR_trk_inB", 500, 0., 1.);
       hist_single_matched_DR_trk_pTfraction = new TH2F("single_matched_DR_trk_pTfraction_inB","single_matched_DR_trk_pTfraction_inB",500,0.,1.,500,-1.1,10.);
       hist_single_matched_d0_inB = new TH1F("single_matched_child_d0_inB","single_matched_child_d0_inB",300,-15.,15.);
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

  m_noB=0,m_bb=0,m_b=0,m_bc_overlap=0,m_sc=0,m_sc2=0,m_sc3=0,m_match=0,m_nomatch=0,m_match_overlap=0,m_match_notoverlap=0,n_trk_pT_cut=0,n_trk_PU_pT_cut=0,n_trk_FRAG_pT_cut=0,n_trk_GEANT_pT_cut=0,n_trk_B=0,n_trk_C=0;

  m_qc=0,m_qj=0,q=0,a=0,b=0,sc=0,sgn=0;
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

  return;
}

   // jet labelling (b-, c-, l-) based on truth (with pt, Deta cuts), exclusive samples
void DAOD_selector::getTrueJetFlavourLabel(std::vector<int>& isJet, std::vector<int>& isBcheck, std::vector<int>& isCcheck, std::vector<int>& islcheck)
{

  isBcheck.clear();
  isCcheck.clear();
  islcheck.clear();

  bool jet_labelled = false;
  double D_eta = 0.,D_phi = 0.;
  double DeltaR=0.,pt=0.;
   for(std::vector<int>::iterator it = isJet.begin(); it != isJet.end(); ++it){
     jet_labelled=false;
     if(jet_nBHadr[*it]>0){
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
