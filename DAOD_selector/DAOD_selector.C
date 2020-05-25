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


void DAOD_selector::setFlags(bool lxplus_flag, bool debug_flag, bool selections_flag, bool discriminants_flag, bool shrinking_cone_flag, bool selection_alg_flag, bool cut_flag, bool retagT_flag)
{
   lxplus=lxplus_flag;
   debug=debug_flag;
   selections=selections_flag;
   discriminants=discriminants_flag;
   shrinking_cone=shrinking_cone_flag;
   selection_alg=selection_alg_flag;
   cut=cut_flag;
   retagT=retagT_flag;
   return;
}

void DAOD_selector::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

   m_cut=1.,m_fc=0.08;
   m_N=0,m_Ntot=0,m_b2d=0,m_b3d=0,m_c2d=0,m_c3d=0,m_noB=0,m_bb=0,m_b=0,m_bc_overlap=0,m_nbjets=0,m_nl=0,m_sc=0,m_sc2=0,m_sc3=0,m_match=0,m_nomatch=0,m_truth_match=0,m_trk_BC=0;
   m_qc=0,m_qj=0,q=0,a=0,b=0,sc=0,sgn=0;
   D_phi=0.,D_eta=0.,DR=0.,px=0.,py=0.,Dx_1=0.,Dy_1=0.,Dz_1=0.,Dx_2=0.,Dy_2=0.,Dz_2=0,Dxy_1=0,x0=0,y0=0,Dx_3=0.,Dy_3=0.,Dxy_3=0.,rand_n=0.,R0=0,d0=0,c=2.99792458e8;//,nx=0,ny=0;
   D_phi_trk=0.,D_eta_trk=0.,DR_trk=0.,DpT_trk=0.;
   match=0,max_size=0;
   tmp_pTfraction=0.,tmp_DR=0.,tmp_min_pTfraction=1.,tmp_min_DR=1.,m_pTfraction_cut=1.,m_DRcut=0.1,m_pTfraction_nocut=1e6,m_DRnocut=1e6;
   size_jet=0,size_child=0;
   den=0,m_den=0;

   m_track_cut=10;

   pt_max=500., pt_min=0.;
   Delta=(pt_max-pt_min)/bin_1;

   std::cout<<"\n";

   if(lxplus && !debug){
     if(retagT){
       std::cout<<"RETAG TRUE\n";
       if(cut){
         std::cout<<"CUT\n";
         file = new TFile("../output_files/lxplus_output_doRetagT_cut.root","RECREATE");
       }
       if(!cut){
         std::cout<<"NO CUT\n";
         file = new TFile("../output_files/lxplus_output_doRetagT_nocut.root","RECREATE");
       }
     }
     if(!retagT){
       std::cout<<"RETAG FALSE\n";
       if(cut){
         std::cout<<"CUT\n";
         file = new TFile("../output_files/lxplus_output_doRetagF_cut.root","RECREATE");
       }
       if(!cut){
         std::cout<<"NO CUT\n";
         file = new TFile("../output_files/lxplus_output_doRetagF_nocut.root","RECREATE");
       }
     }
   }
   if(debug){
     std::cout<<"DEBUG MODE\n";
     file = new TFile("debug.root","RECREATE");
   }
   if(!lxplus && !debug){
     if(retagT){
       std::cout<<"RETAG TRUE\n";
       if(cut){
         std::cout<<"CUT\n";
         file = new TFile("../output_files/output_13022127_000030_doRetagT_cut.root","RECREATE");
       }
       if(!cut){
         std::cout<<"NO CUT\n";
         file = new TFile("../output_files/output_13022127_000030_doRetagT_nocut.root","RECREATE");
       }
     }
     if(!retagT){
       std::cout<<"RETAG FALSE\n";
       if(cut){
         std::cout<<"CUT\n";
         file = new TFile("../output_files/output_13022127_000030_doRetagF_cut.root","RECREATE");
       }
       if(!cut){
         std::cout<<"NO CUT\n";
         file = new TFile("../output_files/output_13022127_000030_doRetagF_nocut.root","RECREATE");
       }
     }
   }



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

   }


   if(shrinking_cone){

     hist_n_tracks = new TH1F("n tracks", "n1==1", 100, 0, 50);
     hist_tracks_DR = new TH2F("n tracks-Delta_R", "n1==1", 100, 0., 1., 100., 0., 100.);
     hist_DR_1 = new TH1F("Delta_R", "n1==1", 100, 0., 10.);
//   hist_std_dev_DR_1 = new TH1F("Standard Deviation from Delta_R", "n1==1", 100, 0., 1.);
     jet_DR_pT = new TH2F("jets Delta_R-pT","n1==1", 100, 0., 500., 100, 0., 1.);

   }


   if(selection_alg){

     hist_efficiency_inB = new TH1F("efficiency","efficiency",50,-0.1,1.1);
     hist_n_trk = new TH1D("n_trk","n_trk",100,0.,50);
     hist_n_child = new TH1D("n_child","n_child",100,0.,50);
     hist_n_match = new TH1D("n_match","n_match",100,0.,50);

     hist_trk_pT_inB = new TH1F("trk_pT_400_inB", "trk_pT_400_inB", 500, 0., 150.);
     hist_trk_eta_inB = new TH1F("trk_eta_400_inB","trk_eta_400_inB",500,-2.6,2.6);
     hist_trk_phi_inB = new TH1F("trk_phi_400_inB","trk_phi_400_inB",500,-4.,4.);
     hist_trk_Deta_inB = new TH1F("trk_Deta_400_inB","trk_Deta_400_inB",500,-1.,1.);
     hist_trk_Dphi_inB = new TH1F("trk_Dphi_400_inB","trk_Dphi_400_inB",500,-1.,1.);
     hist_trk_Dphi_Deta_inB = new TH2F("trk_Dphi_Deta_400","trk_Dphi_Deta_400", 100, -1.,1., 100, -1.,1.);
     hist_trk_DR_inB = new TH1F("trk_DR_400_inB","trk_DR_400_inB",500,-0.1,2.);
     hist_trk_pT_DR_inB = new TH2F("trk_pT_DR_400_inB","trk_pT_DR_400_inB",80,0.,150.,80,0.,0.6);
     hist_trk_origin_inB = new TH1F("trk_origin_inB","trk_origin_inB",7,-2,5);
     hist_trk_d0_inB = new TH1F("trk_400_d0_inB","trk_400_d0_inB",300,-15.,15.);
     hist_trk_d0_PUinB = new TH1F("trk_400_d0_PUinB","trk_400_d0_PUinB",300,-15.,15.);
     hist_trk_d0_BinB = new TH1F("trk_400_d0_BinB","trk_400_d0_BinB",300,-15.,15.);
     hist_trk_d0_CinB = new TH1F("trk_400_d0_CinB","trk_400_d0_CinB",300,-15.,15.);
     hist_trk_d0_FRAGinB = new TH1F("trk_400_d0_FRAGinB","trk_400_d0_FRAGinB",300,-15.,15.);
     hist_trk_d0_GEANTinB = new TH1F("trk_400_d0_GEANTinB","trk_400_d0_GEANTinB",300,-15.,15.);

     hist_child_pT_inB = new TH1F("child_pT_400_inB", "child_pT_400_inB", 500, 0., 150.);
     hist_child_eta_inB = new TH1F("child_eta_400_inB","child_eta_400_inB",500,-2.6,2.6);
     hist_child_phi_inB = new TH1F("child_phi_400_inB","child_phi_400_inB",500,-4.,4.);
     hist_child_Deta_inB = new TH1F("child_Deta_400_inB","child_Deta_400_inB",500,-1.,1.);
     hist_child_Dphi_inB = new TH1F("child_Dphi_400_inB","child_Dphi_400_inB",500,-1.,1.);
     hist_child_Dphi_Deta_inB = new TH2F("child_Dphi_Deta_400","child_Dphi_Deta_400", 100, -1.,1., 100, -1.,1.);
     hist_child_DR_inB = new TH1F("child_DR_400_inB","child_DR_400_inB",500,-0.1,2.);
     hist_child_pT_DR_inB = new TH2F("child_pT_DR_400_inB","child_pT_DR_400_inB",80,0.,150.,80,0.,0.6);
     hist_child_pT_jet_DR_inB = new TH2F("child_pT_jet_DR_400_inB","child_pT_jet_DR_400_inB",100,0.,500.,100,0.,2.);
     hist_child_Lxy_inB = new TH1F("child_400_Lxy_inB","child_400_Lxy_inB",300,0.,1000.);
     hist_child_Lxyz_inB = new TH1F("child_400_Lxyz_inB","child_400_Lxyz_inB",300,0.,1000.);
     hist_child_decay_IP = new TH1F("child_400_decay_IP_inB","child_400_decay_IP_inB",300,-15.,15.);
     hist_child_nodecay_IP = new TH1F("child_400_nodecay_IP_inB","child_400_nodecay_IP_inB",300,-15.,15.);
//     hist_child_linear_IP = new TH1F("child_400_linear_IP_inB","child_400_linear_IP_inB",300,-15.,15.);
     hist_child_d0_truth = new TH1F("child_400_d0_truth_inB","child_400_d0_truth_inB",300,-15.,15.);
     hist_child_d0 = new TH1F("child_400_d0_inB","child_400_d0_inB",300,-15.,15.);
     hist_child_d0_pT = new TH2F("child_400_d0_pT_inB","child_400_pT_d0_inB",300,-15.,15.,300,0.,150.);
     hist_pT_vs_R0_ratio_inB = new TH1F("child_pT_vs_R0_ratio_inB","child_pT_vs_R0_ratio_inB",300,-1.2,2.4);

     hist_matched_pT_inB = new TH1F("matched_child_pT_400_inB","matched_child_pT_400_inB", 500, 0., 150.);
     hist_matched_eta_inB = new TH1F("matched_child_eta_400_inB","matched_child_eta_400_inB",500,-2.6,2.6);
     hist_matched_phi_inB = new TH1F("matched_child_phi_400_inB","matched_child_phi_400_inB",500,-4.,4.);
     hist_matched_Deta_inB = new TH1F("matched_child_Deta_400_inB","matched_child_Deta_400_inB",500,-1.,1.);
     hist_matched_Dphi_inB = new TH1F("matched_child_Dphi_400_inB","matched_child_Dphi_400_inB",500,-1.,1.);
     hist_matched_Dphi_Deta_inB = new TH2F("matched_child_Dphi_Deta_400","matched_child_Dphi_Deta_400", 100, -1.,1., 100, -1.,1.);
     hist_matched_DR_inB = new TH1F("matched_child_DR_400_inB","matched_child_DR_400_inB",500,-0.1,2.);
     hist_matched_pT_DR_inB = new TH2F("matched_child_pT_DR_400_inB","matched_child_pT_DR_400_inB",80,0.,150.,80,0.,0.6);
     hist_matched_pT_jet_DR_inB = new TH2F("matched_child_pT_jet_DR_400_inB","matched_child_pT_jet_DR_400_inB",100,0.,500.,100,0.,2.);
     hist_matched_pdgId_inB = new TH1F("matched_child_pdgId_400_inB","matched_child_pdgId_400_inB",2000,-1000,1000);
     hist_matched_origin_inB = new TH1F("matched_trk_origin_inB","matched_trk_origin_inB",7,-2,5);
     hist_matched_pT_child_pTfraction_inB = new TH2F("matched_pT_child_pTfraction_400_inB","matched_pT_child_vs_DpT/pT_child_400_inB",500,0.,150.,500,-1.1,10.);
     hist_matched_DR_trk_inB = new TH1F("matched_DR_trk_400_inB","matched_DR_trk_inB_400_inB",500,0.,1.);
     hist_matched_DR_trk_pTfraction = new TH2F("matched_DR_trk_pTfraction_400_inB","matched_DR_trk_pTfraction_400_inB",500,0.,1.,500,-1.,10.);
     hist_matched_d0_inB = new TH1F("matched_child_400_d0_inB","matched_child_400_d0_inB",300,-15.,15.);
     hist_matched_Lxy_inB = new TH1F("matched_child_400_Lxy_inB","matched_child_400_Lxy_inB",300,0.,1000.);
     hist_matched_Lxyz_inB = new TH1F("matched_child_400_Lxyz_inB","matched_child_400_Lxyz_inB",300,0.,1000.);

     hist_nomatched_pT_inB = new TH1F("nomatched_child_pT_400_inB", "nomatched_child_pT_400_inB", 500, 0., 150.);
     hist_nomatched_eta_inB = new TH1F("nomatched_child_eta_400_inB","nomatched_child_eta_400_inB",500,-2.6,2.6);
     hist_nomatched_phi_inB = new TH1F("nomatched_child_phi_400_inB","nomatched_child_phi_400_inB",500,-4.,4.);
     hist_nomatched_Deta_inB = new TH1F("nomatched_child_Deta_400_inB","nomatched_child_Deta_400_inB",500,-1.,1.);
     hist_nomatched_Dphi_inB = new TH1F("nomatched_child_Dphi_400_inB","nomatched_child_Dphi_400_inB",500,-1.,1.);
     hist_nomatched_Dphi_Deta_inB = new TH2F("nomatched_child_Dphi_Deta_400","nomatched_child_Dphi_Deta_400", 100, -1.,1., 100, -1.,1.);
     hist_nomatched_DR_inB = new TH1F("nomatched_child_DR_400_inB","nomatched_child_DR_400_inB",500,-0.1,2.);
     hist_nomatched_pT_DR_inB = new TH2F("nomatched_child_pT_DR_400_inB","nomatched_child_pT_DR_400_inB",80,0.,150.,80,0.,0.6);
     hist_nomatched_pT_jet_DR_inB = new TH2F("nomatched_child_pT_jet_DR_400_inB","nomatched_child_pT_jet_DR_400_inB",100,0.,500.,100,0.,2.);
     hist_nomatched_pdgId_inB = new TH1F("nomatched_child_pdgId_400_inB","nomatched_child_pdgId_400_inB",2000,-1000,1000);

     hist_single_matched_pT_inB = new TH1F("single_matched_pT_400_inB","single_matched_pT_400_inB", 500, 0., 150.);
     hist_single_matched_eta_inB = new TH1F("single_matched_child_eta_400_inB","single_matched_child_eta_400_inB",500,-2.6,2.6);
     hist_single_matched_phi_inB = new TH1F("single_matched_child_phi_400_inB","single_matched_child_phi_400_inB",500,-4.,4.);
     hist_single_matched_Deta_inB = new TH1F("single_matched_child_Deta_400_inB","single_matched_child_Deta_400_inB",500,-1.,1.);
     hist_single_matched_Dphi_inB = new TH1F("single_matched_child_Dphi_400_inB","single_matched_child_Dphi_400_inB",500,-1.,1.);
     hist_single_matched_Dphi_Deta_inB = new TH2F("single_matched_child_Dphi_Deta_400","single_matched_child_Dphi_Deta_400", 100, -1.,1., 100, -1.1,1.);
     hist_single_matched_DR_inB = new TH1F("single_matched_DR_400_inB","single_matched_DR_400_inB",500,-0.1,2.);
     hist_single_matched_pT_DR_inB = new TH2F("single_matched_pT_DR_400_inB","single_matched_pT_DR_400_inB",80,0.,150.,80,0.,0.6);
     hist_single_matched_pT_jet_DR_inB = new TH2F("single_matched_child_pT_jet_DR_400_inB","single_matched_child_pT_jet_DR_400_inB",100,0.,500.,100,0.,2.);
     hist_single_matched_pdgId_inB = new TH1F("single_matched_child_pdgId_400_inB","single_matched_child_pdgId_400_inB",2000,-1000,1000);
     hist_single_matched_origin_inB = new TH1F("single_matched_trk_origin_inB","single_matched_trk_origin_inB",7,-2,5);
     hist_single_matched_pT_child_pTfraction_inB = new TH2F("single_matched_pT_child_pTfraction_400_inB","single_matched_pT_child_vs_DpT/pT_child_400_inB",500,0.,150.,500,-1.1,10.);
     hist_single_matched_DR_trk_inB = new TH1F("single_matched_DR_trk_400_inB","single_matched_DR_trk_400_inB", 500, 0., 1.);
     hist_single_matched_DR_trk_pTfraction = new TH2F("single_matched_DR_trk_pTfraction_400_inB","single_matched_DR_trk_pTfraction_400_inB",500,0.,1.,500,-1.1,10.);
     hist_single_matched_d0_inB = new TH1F("single_matched_child_400_d0_inB","single_matched_child_400_d0_inB",300,-15.,15.);

   }


   TString option = GetOption();
}

void DAOD_selector::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();
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

   int nB=0,nBjets=0;
   int n1=0,n2=0,n3=0; //<- the file has max=3 b-tagged particles in a single jet
   std::vector<int> isB,is1B,is2B,is3B;

   int nl=0;
   std::vector<int> isl;

   m_Ntot++;

   for(int i=0;i<*njets;i++) {
     if(jet_nBHadr[i]==0 && jet_nCHadr[i]==0){
       nl++;  //m_nl++;
       isl.push_back(i);
     }
   }

   for(unsigned i=0;i<jet_nBHadr.GetSize();i++) {

     nB+=jet_nBHadr[i];
     if(jet_nBHadr[i]>0){
       isB.push_back(i);
       nBjets++;  m_nbjets++;
       if(jet_nCHadr[i]>0)  m_bc_overlap++;
     }
     if(jet_nBHadr[i]==1)  {is1B.push_back(i); n1++;}
     if(jet_nBHadr[i]==2)  {is2B.push_back(i); n2++;}
     if(jet_nBHadr[i]==3)  {is3B.push_back(i); n3++;}

   }

   if(nB==0)  {m_noB++;}

   if (n1+2*n2+3*n3!=nB) {
        std::cout << "Warning: n1+n2+n3!=nB\t" << n1<<"\t"<<n2<<"\t"<<n3<<"\t"<<nB<<"\n";
    }

   int nC=0,nCjets=0;
   int nC1=0,nC2=0,nC3=0; //<- the file has max=3 b-tagged particles in a single jet
   std::vector<int> isC,is1C,is2C,is3C;

   for(unsigned i=0;i<jet_nCHadr.GetSize();i++) {

     if(jet_nCHadr[i]>0){
       isC.push_back(i);
       nCjets++;
     }
     if(jet_nCHadr[i]==1)  {is1C.push_back(i); nC1++;}
     if(jet_nCHadr[i]==2)  {is2C.push_back(i); nC2++;}
     if(jet_nCHadr[i]==3)  {is3C.push_back(i); nC3++;}
   }


   if(selections){


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
         m_N++;
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




     if(nBjets>0){
       for(std::vector<int>::iterator it = isB.begin(); it != isB.end(); ++it){

         hist_pt_inB->Fill(jet_pt[*it]*0.001);
         hist_eta_inB->Fill(jet_eta[*it]);
         hist_phi_inB->Fill(jet_phi[*it]);
         hist_E_inB->Fill(jet_E[*it]);
       }
     }



   }



   if(shrinking_cone){
 //select by: n1==1 && n2==0 && n3==0
     if(n1==1 && n2==0 && n3==0){

       int i=is1B.at(0);
       if(jet_pt[i]*0.001<pt_max){

         std::vector<float> eta=jet_trk_eta[i];
         std::vector<float> phi=jet_trk_phi[i];
         if(eta.size()!=phi.size()) std::cout<< "ERROR"<< "\n";
         int size=eta.size();//size is the number of tracks, on which we make a cut
 //        if(size>=m_track_cut){
         std::vector<float> R(size);
         for(int j=0;j<size;j++){
           R.at(j)=sqrt(eta.at(j)*eta.at(j)+phi.at(j)*phi.at(j));
         }
         double R_m;
 //          R_m=std::accumulate(R.begin(), R.end(), 0.0)/R.size();
         R_m=sqrt(jet_eta[i]*jet_eta[i]+jet_phi[i]*jet_phi[i]);

         float D_phi=0;
         float DR_Max=0,tmp_M=0,sq_sum=0,std_dev=0;
         for(int j=0;j<size;j++){
 //            tmp_M=abs(R.at(j)-R_m);//<---------------------------------------WRONG!!!!!!!!!!!!!
           if(abs(jet_trk_phi[i].at(j)-jet_phi[i])>M_PI){
             D_phi=2*M_PI-abs(jet_trk_phi[i].at(j)-jet_phi[i]);
           }
           if(abs(jet_trk_phi[i].at(j)-jet_phi[i])<M_PI){
             D_phi=jet_trk_phi[i].at(j)-jet_phi[i];
           }
             //if(jet_trk_phi[i].at(j)>0 && jet_phi[i]<0) std::cout<<D_phi<<"\n";
           tmp_M=sqrt((eta.at(j)-jet_eta[i])*(eta.at(j)-jet_eta[i])+D_phi*D_phi);
 //            sq_sum+=tmp_M*tmp_M;
           if(tmp_M>DR_Max){
             DR_Max=tmp_M;
           }
         }

 //          std_dev=sqrt(sq_sum/(size-1));
 //          hist_std_dev_DR_1->Fill(std_dev);
         hist_n_tracks->Fill(size);
 //        g->SetPoint(g->GetN(), jet_pt[i]*0.001, DR_Max);
         jet_DR_pT->Fill(jet_pt[i]*0.001, DR_Max);
         hist_tracks_DR->Fill(DR_Max, size);
         hist_DR_1->Fill(DR_Max);

         int quot=(int) jet_pt[i]*0.001/Delta;
         bin_v.at(quot).push_back(DR_Max);

 //        }
       }
     }
   }



   if(selection_alg){
     if(nBjets>0){
       for(std::vector<int>::iterator it = isB.begin(); it != isB.end(); ++it){

  //         int b=(jet_bH_pdgId[*it][i]/100)%10;
         size_jet=jet_trk_pt[*it].size();

         if(size_jet!=jet_trk_eta[*it].size() || size_jet!=jet_trk_phi[*it].size()){
           std::cout<<"WARNING\n";
         }
         size_child=jet_bH_child_px[*it].size();
         if(size_child!=jet_bH_child_py[*it].size() || size_child!=jet_bH_child_pz[*it].size() || size_child!=jet_bH_child_E[*it].size()){
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

           for(unsigned i=0;i<size_jet;i++){
             if(abs(jet_trk_eta[*it].at(i))<2.5 && jet_trk_pt[*it].at(i)>400.){

               A=sin(jet_phi[*it]-jet_trk_phi[*it].at(i))*jet_trk_d0[*it].at(i);
               sgn=A/abs(A);
               d0=sgn*abs(jet_trk_d0[*it].at(i));
               hist_trk_d0_inB->Fill(d0);
               if(jet_trk_orig[*it].at(i)==-1){
                 hist_trk_d0_PUinB->Fill(d0);
               }
               if(jet_trk_orig[*it].at(i)==0){
                 hist_trk_d0_BinB->Fill(d0);
                 m_trk_BC++;
               }
               if(jet_trk_orig[*it].at(i)==1){
                 hist_trk_d0_CinB->Fill(d0);
                 m_trk_BC++;
               }
               if(jet_trk_orig[*it].at(i)==2){
                 hist_trk_d0_FRAGinB->Fill(d0);
               }
               if(jet_trk_orig[*it].at(i)==3){
                 hist_trk_d0_GEANTinB->Fill(d0);
               }

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
               hist_trk_origin_inB->Fill(jet_trk_orig[*it].at(i));
             }
           }

           std::vector<double> child_Pt,child_Eta,child_Phi;
           std::vector<int> child_idx;

           double child_IP[size_child],child_Lxy[size_child],child_Lxyz[size_child];

           int bool_matrix[size_jet][size_child];
           for(unsigned l=0;l<size_jet;l++){
             for(unsigned k=0;k<size_child;k++){
                bool_matrix[l][k]=0;
             }
           }
           den=0;

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

             if(abs(child_Eta.at(j))<2.5 && child_Pt.at(j)>400.){//CHILD SELECTION CRITERIA
               if(jet_bH_child_prod_x[*it].size()!=size_child || jet_bH_child_decay_x[*it].size()!=size_child){
                 std::cout<<"WARNING\n";
               }

               den++;
               child_idx.push_back(j);

               Dx_1=(*PVx)-jet_bH_child_prod_x[*it].at(j);//x_PV - x_SV
               Dy_1=(*PVx)-jet_bH_child_prod_y[*it].at(j);
               Dz_1=(*PVx)-jet_bH_child_prod_z[*it].at(j);
               Dxy_1=sqrt(Dx_1*Dx_1+Dy_1*Dy_1);
               hist_child_Lxy_inB->Fill(Dxy_1);
               hist_child_Lxyz_inB->Fill(sqrt(Dxy_1*Dxy_1+Dz_1*Dz_1));
               child_Lxy[j]=Dxy_1;
               child_Lxyz[j]=sqrt(Dxy_1*Dxy_1+Dz_1*Dz_1);
               px=jet_bH_child_px[*it].at(j);
               py=jet_bH_child_py[*it].at(j);

               if(jet_bH_child_decay_x[*it].at(j)==-999){//NO DECAY CHILD
                 if(jet_bH_child_charge[*it].at(j)>1.){
                   c=-1;
                 }
                 if(jet_bH_child_charge[*it].at(j)<1.){
                   c=+1;
                 }

//                 TRandom3 *rand = new TRandom3(0);
//                 rand_n=rand->Gaus(0.59,0.05943);
                 rand_n=0.6;
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

                 hist_child_nodecay_IP->Fill(d0);

               }

               if(jet_bH_child_decay_x[*it].at(j)!=-999){//DECAY CHILD
//                 std::cout<<jet_bH_child_decay_x[*it].at(j)<<"\n";
                 Dx_2=jet_bH_child_decay_x[*it].at(j)-jet_bH_child_prod_x[*it].at(j);//x_DV - x_SV
                 Dy_2=jet_bH_child_decay_y[*it].at(j)-jet_bH_child_prod_y[*it].at(j);
                 Dz_2=jet_bH_child_decay_z[*it].at(j)-jet_bH_child_prod_z[*it].at(j);
//                 Dxy_2=sqrt((Dx_1*Dx_1)+(Dy_1*Dy_1));

//                 vx=1e3*c*px/jet_bH_child_E[*it].at(j);
//                 vy=1e3*c*py/jet_bH_child_E[*it].at(j);
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

                 hist_child_decay_IP->Fill(d0);

                 hist_pT_vs_R0_ratio_inB->Fill(child_Pt.at(j)/abs(R0));

//                 ncorr+=jet_bH_child_charge[*it].at(j)*d0/abs(d0);
//                 std::cout<<jet_bH_child_charge[*it].at(j)*R0/abs(jet_bH_child_charge[*it].at(j)*R0)<<"\n";
//                 std::cout<<jet_bH_child_charge[*it].at(j)<<"\t"<<d0/abs(d0)<<"\n";
               }

//               d0=(px*(-Dy_1)-py*(-Dx_1))/sqrt(px*px+py*py);//LINEAR D0
//               t1=sin(jet_phi[*it]-child_Phi.at(j))*d0;
//               t2=t1/abs(t1);
//               hist_child_linear_IP->Fill(t2*abs(d0));

               A=sin(jet_phi[*it]-child_Phi.at(j))*jet_bH_child_d0[*it].at(j);
               sgn=A/abs(A);

               hist_child_d0_truth->Fill(sgn*abs(jet_bH_child_d0[*it].at(j)));
               hist_child_d0->Fill(d0);
               hist_child_d0_pT->Fill(d0,1e-3*sqrt(px*px+py*py));
               child_IP[j]=d0;

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

               for(unsigned i=0;i<size_jet;i++){
                 if((jet_trk_parent_pdgid[*it].at(i)==jet_bH_child_parent_pdg_id[*it].at(j) && jet_trk_pdg_id[*it].at(i)==jet_bH_child_pdg_id[*it].at(j)) && (abs(jet_trk_parent_pdgid[*it].at(i))!=99 && abs(jet_trk_pdg_id[*it].at(i))!=99)){
//                 if(jet_trk_pdg_id[*it].at(i)==jet_bH_child_pdg_id[*it].at(j)){
                   bool_matrix[i][j]=1;
                 }
               }

             }
           }
           if(child_Pt.size()!=size_child || child_Eta.size()!=size_child || child_Phi.size()!=size_child){
             std::cout<<"WARNING\n";
           }



  /*
  std::cout<<"\nCORRESPONDENCE MATRIX\t"<<size_jet<<"x"<<size_child<<"\tevent:"<<m_Ntot;
           for(int i=0;i<size_jet;i++){
             std::cout<<"\n";
             for(int j=0;j<size_child;j++){
               std::cout << bool_matrix[i][j]<<"\t";
             }
           }
           std::cout<<"\n";*/

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
               if(jet_trk_orig[*it].at(a)==0 || jet_trk_orig[*it].at(a)==1){
                 m_truth_match++;
               }

               for(unsigned i=0;i<child_idx.size();i++){
                 if(child_idx.at(i)==b){
                   child_idx.erase(child_idx.begin()+i);
                 }
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
               }
               if(sc==2){
                 m_sc2++;
               }
               if(sc>=3){
                 m_sc3++;
               }
             }

  //std::cout<<"\n"<<"q="<<q<<"\t(a,b)=("<<a+1<<","<<b+1<<")"<<"\n";

             for(unsigned j=0;j<size_child;j++){
               bool_matrix[a][j]=0;
             }
             for(unsigned i=0;i<size_jet;i++){
               bool_matrix[i][b]=0;
             }

  /*
  std::cout<<"ELIMINATION STARTS:\n";
             for(int i=0;i<size_jet;i++){
               std::cout<<"\n";
               for(int j=0;j<size_child;j++){
                 std::cout << bool_matrix[i][j]<<"\t";
               }
             }
           std::cout<<"\n";*/
           }

           m_nomatch+=child_idx.size();
           int j=0;
           for(unsigned l=0;l<child_idx.size();l++){
             j=child_idx.at(l);
             //CINEMATICA RISPETTO AL JET
             if(abs(child_Eta.at(j))<2.5 && child_Pt.at(j)>400.){
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

  //std::cout<<"event:\t"<<m_Ntot<< "\tn of matched tracks:\t"<< match << "\tn of child tracks(denominator):\t" << den << "\tratio:\t" <<(float) match/den<<"\t";
           m_match+=match;
           m_den+=den;
           hist_n_match->Fill(match);
           hist_efficiency_inB->Fill((float) match/den);
           match=0;


         }
       }
     }
   }

//inspect tracks

//   for(TTreeReaderArray<std::vector<float, std::allocator<float>>>::iterator it=jet_trk_eta.begin(); it!=jet_trk_eta.end(); ++it){
//      std::vector v=*it;

//      double m;
//      v.size()==0 ? m=-999. : m=std::accumulate(v.begin(), v.end(), 0.0)/v.size();
//       std::cout << m << "\t" << v.size() << "\n";
//      cnt_1++;
//  }

/*
   std::vector<int> count(*njets);
   int n=0;
   for(TTreeReaderArray<std::vector<int, std::allocator<int>>>::iterator it=jet_trk_pdg_id.begin(); it!=jet_trk_pdg_id.end(); ++it){
      std::vector v=*it;
      int k=0;
      count.at(n)=0;
      for(std::vector<int>::iterator itv=v.begin(); itv!=v.end(); ++itv){
          if(*itv==11){
              k++;
              count[n]=k;
          }
//          std::cout<< *itv << "\n";
      }
      n++;
//    std::cout<< "\n";
   }


//    std::cout << count.size() << "\t" << n << "\n";
   int i=0;
   for(std::vector<int>::iterator it=count.begin(); it!=count.end(); ++it) {
     std::cout << *it << "\t" << jet_nBHadr[i] << "\n";
     if(jet_nBHadr[i]-*it){
       std::cout << "WARNING" << "\n";
     }
     i++;
   }

  std::cout << "\n";


//   if (n1>3 || n2>1 || n3>1) std::cout << "Warning" << "\t"<< n1<<"\t"<<n2<<"\t"<<n3<<"\t"<<nB<<"\n";

//  for(int i=0;i<*njets;i++){
//    std::cout << jet_ip2d_pb[i] << "\n";




//   std::cout << "\n";

   for(TTreeReaderArray<std::vector<float, std::allocator<float>>>::iterator it=jet_trk_vtx_Z.begin(); it!=jet_trk_vtx_Z.end(); ++it){
      std::vector v=*it;
      double m=0.,sum=0.;
      int cnt=0;
      if(v.size()==0)   m=-999.;
      if(v.size()!=0) {
          for(std::vector<float>::iterator itv=v.begin(); itv!=v.end(); ++itv){
              if(*itv!=-999){
                  sum+=*itv;
                  cnt++;
              }
           }
      m=sum/cnt;
      cnt=0;
      }
      std::cout << m << "\t" << v.size() << "\n";
//        std::cout << v.size() << "\n";
      cnt_2++;
   }


//       std::cout << "\n";
//       std::cout << cnt_1 << "\t" << "\n";
//       std::cout << "\n";

*/
   //secondary vtx inspection: jet_sv1_Nvtx is a boolean vector (0,1), with jet_sv1_Nvtx.GetSize()=*njet
   //idea: get the pdgId && jet_sv1_Nvtx==1 distribution


   return kTRUE;
}

void DAOD_selector::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.

}

void DAOD_selector::Terminate()
{
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
      hist_trk_origin_inB->Write();
      hist_trk_d0_inB->Write();
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
      hist_child_Lxy_inB->Write();
      hist_child_Lxyz_inB->Write();
      hist_child_decay_IP->Write();
      hist_child_nodecay_IP->Write();
//      hist_child_linear_IP->Write();
      hist_child_d0_truth->Write();
      hist_child_d0->Write();
      hist_child_d0_pT->Write();
      hist_pT_vs_R0_ratio_inB->Write();

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
      hist_matched_origin_inB->Write();
      hist_matched_pT_child_pTfraction_inB->Write();
      hist_matched_DR_trk_inB->Write();
      hist_matched_DR_trk_pTfraction->Write();
      hist_matched_Lxy_inB->Write();
      hist_matched_Lxyz_inB->Write();
      hist_matched_d0_inB->Write();

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
    }


    file->Close();

    if(selections){
      std::cout<< "\nSELECTIONS\n\n";
      std::cout<< "fraction of events with one single b:\t" << (double) m_b/m_Ntot << "\n";
      std::cout<< "fraction of events with two single b:\t" << (double) m_bb/m_Ntot << "\n";
      std::cout<< "fraction of events without b:\t\t" << (double) m_noB/m_Ntot << "\n";
      std::cout<< "total b-c overlap:\t" << (double) m_bc_overlap/m_nbjets << "\n";
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
          std::cout << std::fixed << std::setprecision(5) << R_bin[i] << "\t+- " << std_dev_R[i] << "\t" << "with\t" << sup << "\tpoints with pT in\t" << "["  << std::fixed << std::setprecision(1) << i*Delta << ","  << (i+1)*Delta << "]\tGeV" << "\n";

          x[i]=i*Delta+Delta*0.5;
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
      if(cut && selection_alg)  std::cout<<"\nSELECTION ALGORITHM - CUT\n";
      if(!cut && selection_alg) std::cout<<"\nSELECTION ALGORITHM - NO CUT\n";
      std::cout<< std::fixed << std::setprecision(3) << "\n";
      std::cout<< "number of childs_400:\t" << m_den << "\n";
      std::cout<< "matches:\t\t"<< m_match <<"\n";
      std::cout<< "no matches:\t\t"<< m_nomatch <<"\n";
      std::cout<< "average efficiency:\t" << (float) m_match/m_den << "\n";
      std::cout<< "single matches:\t\t" << m_sc <<"\tfraction:\t"<<(float) m_sc/m_match << "\n";
      std::cout<< "double matches:\t\t" << m_sc2 <<"\tfraction:\t"<<(float) m_sc2/m_match << "\n";
      std::cout<< "triple+ matches:\t" << m_sc3 <<"\tfraction:\t"<<(float) m_sc3/m_match << "\n";
      std::cout<< "algorithm average eff.:\t\t" << (float) m_truth_match/m_match <<"\n";
      std::cout<< "relative matching freq.:\t" << (float) m_truth_match/m_trk_BC <<"\n";
      std::cout<< "efficiency from origin:\t" << (float) m_trk_BC/m_den <<"\n";
    }

    std::cout<<"\n";

}
