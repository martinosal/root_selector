//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon May 25 17:19:26 2020 by ROOT version 6.18/04
// from TTree bTag_AntiKt4EMPFlowJets_BTagging201903/bTagAntiKt4EMPFlowJets_BTagging201903
// found on file: flav_Akt4EMPf_BTagging201903_IP_doRetagF.root
//////////////////////////////////////////////////////////

#ifndef DAOD_selector_h
#define DAOD_selector_h

#define bin_1 50
#define tracksize 10

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <TH2.h>
#include <TH2F.h>
#include <TH1F.h>
#include <TLorentzVector.h>
#include <TMultiGraph.h>
#include <TGraphErrors.h>
#include <TGraph.h>
#include <TRandom.h>
#include <random>
// Headers needed by this particular selector
#include <math.h>
#include <numeric>
#include <vector>
#include <iostream>
#include <cmath>



class DAOD_selector : public TSelector {
public :
   TTreeReader     fReader;  //!the tree reader
   TTree          *fChain = 0;   //!pointer to the analyzed TTree or TChain

   // Readers to access the data (delete the ones you do not need).
   TTreeReaderValue<Int_t> runnb = {fReader, "runnb"};
   TTreeReaderValue<Int_t> eventnb = {fReader, "eventnb"};
   TTreeReaderValue<Int_t> mcchan = {fReader, "mcchan"};
   TTreeReaderValue<Float_t> mcwg = {fReader, "mcwg"};
   TTreeReaderValue<Float_t> avgmu = {fReader, "avgmu"};
   TTreeReaderValue<Int_t> actmu = {fReader, "actmu"};
   TTreeReaderValue<Float_t> PVx = {fReader, "PVx"};
   TTreeReaderValue<Float_t> PVy = {fReader, "PVy"};
   TTreeReaderValue<Float_t> PVz = {fReader, "PVz"};
   TTreeReaderValue<Float_t> truth_PVx = {fReader, "truth_PVx"};
   TTreeReaderValue<Float_t> truth_PVy = {fReader, "truth_PVy"};
   TTreeReaderValue<Float_t> truth_PVz = {fReader, "truth_PVz"};
   TTreeReaderValue<Int_t> njets = {fReader, "njets"};
   TTreeReaderArray<float> jet_pt = {fReader, "jet_pt"};
   TTreeReaderArray<float> jet_eta = {fReader, "jet_eta"};
   TTreeReaderArray<float> jet_phi = {fReader, "jet_phi"};
   TTreeReaderArray<float> jet_E = {fReader, "jet_E"};
   TTreeReaderArray<float> jet_pt_orig = {fReader, "jet_pt_orig"};
   TTreeReaderArray<float> jet_eta_orig = {fReader, "jet_eta_orig"};
   TTreeReaderArray<float> jet_phi_orig = {fReader, "jet_phi_orig"};
   TTreeReaderArray<float> jet_E_orig = {fReader, "jet_E_orig"};
   TTreeReaderArray<int> jet_LabDr_HadF = {fReader, "jet_LabDr_HadF"};
   TTreeReaderArray<int> jet_DoubleHadLabel = {fReader, "jet_DoubleHadLabel"};
   TTreeReaderArray<float> jet_JVT = {fReader, "jet_JVT"};
   TTreeReaderArray<float> jet_m = {fReader, "jet_m"};
   TTreeReaderArray<float> jet_nConst = {fReader, "jet_nConst"};
   TTreeReaderArray<float> jet_dRiso = {fReader, "jet_dRiso"};
   TTreeReaderArray<int> jet_truthMatch = {fReader, "jet_truthMatch"};
   TTreeReaderArray<int> jet_isPU = {fReader, "jet_isPU"};
   TTreeReaderArray<int> jet_aliveAfterOR = {fReader, "jet_aliveAfterOR"};
   TTreeReaderArray<int> jet_aliveAfterORmu = {fReader, "jet_aliveAfterORmu"};
   TTreeReaderArray<int> jet_isBadMedium = {fReader, "jet_isBadMedium"};
   TTreeReaderArray<float> jet_truthPt = {fReader, "jet_truthPt"};
   TTreeReaderArray<float> jet_dRminToB = {fReader, "jet_dRminToB"};
   TTreeReaderArray<float> jet_dRminToC = {fReader, "jet_dRminToC"};
   TTreeReaderArray<float> jet_dRminToT = {fReader, "jet_dRminToT"};
   TTreeReaderArray<int> jet_nBHadr = {fReader, "jet_nBHadr"};
   TTreeReaderArray<int> jet_nCHadr = {fReader, "jet_nCHadr"};
   TTreeReaderArray<vector<int>> jet_bH_pdgId = {fReader, "jet_bH_pdgId"};
   TTreeReaderArray<vector<int>> jet_bH_parent_pdgId = {fReader, "jet_bH_parent_pdgId"};
   TTreeReaderArray<vector<float>> jet_bH_pt = {fReader, "jet_bH_pt"};
   TTreeReaderArray<vector<float>> jet_bH_eta = {fReader, "jet_bH_eta"};
   TTreeReaderArray<vector<float>> jet_bH_phi = {fReader, "jet_bH_phi"};
   TTreeReaderArray<vector<float>> jet_bH_E = {fReader, "jet_bH_E"};
   TTreeReaderArray<vector<float>> jet_bH_charge = {fReader, "jet_bH_charge"};
   TTreeReaderArray<vector<float>> jet_bH_Lxy = {fReader, "jet_bH_Lxy"};
   TTreeReaderArray<vector<float>> jet_bH_x = {fReader, "jet_bH_x"};
   TTreeReaderArray<vector<float>> jet_bH_y = {fReader, "jet_bH_y"};
   TTreeReaderArray<vector<float>> jet_bH_z = {fReader, "jet_bH_z"};
   TTreeReaderArray<vector<float>> jet_bH_dRjet = {fReader, "jet_bH_dRjet"};
   TTreeReaderArray<vector<int>> jet_cH_pdgId = {fReader, "jet_cH_pdgId"};
   TTreeReaderArray<vector<int>> jet_cH_parent_pdgId = {fReader, "jet_cH_parent_pdgId"};
   TTreeReaderArray<vector<float>> jet_cH_pt = {fReader, "jet_cH_pt"};
   TTreeReaderArray<vector<float>> jet_cH_eta = {fReader, "jet_cH_eta"};
   TTreeReaderArray<vector<float>> jet_cH_phi = {fReader, "jet_cH_phi"};
   TTreeReaderArray<vector<float>> jet_cH_E = {fReader, "jet_cH_E"};
   TTreeReaderArray<vector<float>> jet_cH_charge = {fReader, "jet_cH_charge"};
   TTreeReaderArray<vector<float>> jet_cH_Lxy = {fReader, "jet_cH_Lxy"};
   TTreeReaderArray<vector<float>> jet_cH_x = {fReader, "jet_cH_x"};
   TTreeReaderArray<vector<float>> jet_cH_y = {fReader, "jet_cH_y"};
   TTreeReaderArray<vector<float>> jet_cH_z = {fReader, "jet_cH_z"};
   TTreeReaderArray<vector<float>> jet_cH_dRjet = {fReader, "jet_cH_dRjet"};
   TTreeReaderArray<vector<int>> jet_trk_orig = {fReader, "jet_trk_orig"};
   TTreeReaderArray<vector<float>> jet_bH_prod_x = {fReader, "jet_bH_prod_x"};
   TTreeReaderArray<vector<float>> jet_bH_prod_y = {fReader, "jet_bH_prod_y"};
   TTreeReaderArray<vector<float>> jet_bH_prod_z = {fReader, "jet_bH_prod_z"};
   TTreeReaderArray<vector<float>> jet_bH_PtTrk = {fReader, "jet_bH_PtTrk"};
   TTreeReaderArray<vector<float>> jet_bH_MTrk = {fReader, "jet_bH_MTrk"};
   TTreeReaderArray<vector<int>> jet_bH_nBtracks = {fReader, "jet_bH_nBtracks"};
   TTreeReaderArray<vector<int>> jet_bH_nCtracks = {fReader, "jet_bH_nCtracks"};
   TTreeReaderArray<vector<int>> jet_bH_nBtracks_400 = {fReader, "jet_bH_nBtracks_400"};
   TTreeReaderArray<vector<int>> jet_bH_nCtracks_400 = {fReader, "jet_bH_nCtracks_400"};
   TTreeReaderArray<vector<int>> jet_bH_child_hadron_idx = {fReader, "jet_bH_child_hadron_idx"};
   TTreeReaderArray<vector<int>> jet_bH_child_pdg_id = {fReader, "jet_bH_child_pdg_id"};
   TTreeReaderArray<vector<int>> jet_bH_child_parent_pdg_id = {fReader, "jet_bH_child_parent_pdg_id"};
   TTreeReaderArray<vector<int>> jet_bH_child_barcode = {fReader, "jet_bH_child_barcode"};
   TTreeReaderArray<vector<float>> jet_bH_child_charge = {fReader, "jet_bH_child_charge"};
   TTreeReaderArray<vector<float>> jet_bH_child_px = {fReader, "jet_bH_child_px"};
   TTreeReaderArray<vector<float>> jet_bH_child_py = {fReader, "jet_bH_child_py"};
   TTreeReaderArray<vector<float>> jet_bH_child_pz = {fReader, "jet_bH_child_pz"};
   TTreeReaderArray<vector<float>> jet_bH_child_E = {fReader, "jet_bH_child_E"};
   TTreeReaderArray<vector<float>> jet_bH_child_prod_x = {fReader, "jet_bH_child_prod_x"};
   TTreeReaderArray<vector<float>> jet_bH_child_prod_y = {fReader, "jet_bH_child_prod_y"};
   TTreeReaderArray<vector<float>> jet_bH_child_prod_z = {fReader, "jet_bH_child_prod_z"};
   TTreeReaderArray<vector<float>> jet_bH_child_decay_x = {fReader, "jet_bH_child_decay_x"};
   TTreeReaderArray<vector<float>> jet_bH_child_decay_y = {fReader, "jet_bH_child_decay_y"};
   TTreeReaderArray<vector<float>> jet_bH_child_decay_z = {fReader, "jet_bH_child_decay_z"};
   TTreeReaderArray<vector<float>> jet_bH_child_d0 = {fReader, "jet_bH_child_d0"};
//   TTreeReaderArray<vector<float>> jet_bH_child_z0 = {fReader, "jet_bH_child_z0"};
   TTreeReaderArray<vector<float>> jet_cH_prod_x = {fReader, "jet_cH_prod_x"};
   TTreeReaderArray<vector<float>> jet_cH_prod_y = {fReader, "jet_cH_prod_y"};
   TTreeReaderArray<vector<float>> jet_cH_prod_z = {fReader, "jet_cH_prod_z"};
   TTreeReaderArray<vector<float>> jet_cH_PtTrk = {fReader, "jet_cH_PtTrk"};
   TTreeReaderArray<vector<float>> jet_cH_MTrk = {fReader, "jet_cH_MTrk"};
   TTreeReaderArray<vector<int>> jet_cH_nCtracks = {fReader, "jet_cH_nCtracks"};
   TTreeReaderArray<vector<int>> jet_cH_nCtracks_400 = {fReader, "jet_cH_nCtracks_400"};
   TTreeReaderArray<vector<int>> jet_cH_child_hadron_idx = {fReader, "jet_cH_child_hadron_idx"};
   TTreeReaderArray<vector<int>> jet_cH_child_pdg_id = {fReader, "jet_cH_child_pdg_id"};
   TTreeReaderArray<vector<int>> jet_cH_child_parent_pdg_id = {fReader, "jet_cH_child_parent_pdg_id"};
   TTreeReaderArray<vector<int>> jet_cH_child_barcode = {fReader, "jet_cH_child_barcode"};
   TTreeReaderArray<vector<float>> jet_cH_child_charge = {fReader, "jet_cH_child_charge"};
   TTreeReaderArray<vector<float>> jet_cH_child_px = {fReader, "jet_cH_child_px"};
   TTreeReaderArray<vector<float>> jet_cH_child_py = {fReader, "jet_cH_child_py"};
   TTreeReaderArray<vector<float>> jet_cH_child_pz = {fReader, "jet_cH_child_pz"};
   TTreeReaderArray<vector<float>> jet_cH_child_E = {fReader, "jet_cH_child_E"};
   TTreeReaderArray<vector<float>> jet_cH_child_prod_x = {fReader, "jet_cH_child_prod_x"};
   TTreeReaderArray<vector<float>> jet_cH_child_prod_y = {fReader, "jet_cH_child_prod_y"};
   TTreeReaderArray<vector<float>> jet_cH_child_prod_z = {fReader, "jet_cH_child_prod_z"};
   TTreeReaderArray<vector<float>> jet_cH_child_decay_x = {fReader, "jet_cH_child_decay_x"};
   TTreeReaderArray<vector<float>> jet_cH_child_decay_y = {fReader, "jet_cH_child_decay_y"};
   TTreeReaderArray<vector<float>> jet_cH_child_decay_z = {fReader, "jet_cH_child_decay_z"};
   TTreeReaderArray<int> jet_nGhostBHadrFromParent = {fReader, "jet_nGhostBHadrFromParent"};
   TTreeReaderArray<int> jet_nGhostCHadrFromParent = {fReader, "jet_nGhostCHadrFromParent"};
   TTreeReaderArray<int> jet_nGhostCHadrFromParentNotFromB = {fReader, "jet_nGhostCHadrFromParentNotFromB"};
   TTreeReaderArray<int> jet_nGhostTauFromParent = {fReader, "jet_nGhostTauFromParent"};
   TTreeReaderArray<int> jet_nGhostHBosoFromParent = {fReader, "jet_nGhostHBosoFromParent"};
   TTreeReaderArray<int> jet_nGhostBHadr = {fReader, "jet_nGhostBHadr"};
   TTreeReaderArray<int> jet_nGhostCHadr = {fReader, "jet_nGhostCHadr"};
   TTreeReaderArray<int> jet_nGhostCHadrNotFromB = {fReader, "jet_nGhostCHadrNotFromB"};
   TTreeReaderArray<int> jet_nGhostTau = {fReader, "jet_nGhostTau"};
   TTreeReaderArray<int> jet_nGhostHBoso = {fReader, "jet_nGhostHBoso"};
   TTreeReaderArray<float> bH1FromParent_pt = {fReader, "bH1FromParent_pt"};
   TTreeReaderArray<float> bH1FromParent_eta = {fReader, "bH1FromParent_eta"};
   TTreeReaderArray<float> bH1FromParent_phi = {fReader, "bH1FromParent_phi"};
   TTreeReaderArray<float> bH1FromParent_Lxy = {fReader, "bH1FromParent_Lxy"};
   TTreeReaderArray<float> bH1FromParent_x = {fReader, "bH1FromParent_x"};
   TTreeReaderArray<float> bH1FromParent_y = {fReader, "bH1FromParent_y"};
   TTreeReaderArray<float> bH1FromParent_z = {fReader, "bH1FromParent_z"};
   TTreeReaderArray<float> bH1FromParent_dRjet = {fReader, "bH1FromParent_dRjet"};
   TTreeReaderArray<float> bH2FromParent_pt = {fReader, "bH2FromParent_pt"};
   TTreeReaderArray<float> bH2FromParent_eta = {fReader, "bH2FromParent_eta"};
   TTreeReaderArray<float> bH2FromParent_phi = {fReader, "bH2FromParent_phi"};
   TTreeReaderArray<float> bH2FromParent_Lxy = {fReader, "bH2FromParent_Lxy"};
   TTreeReaderArray<float> bH2FromParent_x = {fReader, "bH2FromParent_x"};
   TTreeReaderArray<float> bH2FromParent_y = {fReader, "bH2FromParent_y"};
   TTreeReaderArray<float> bH2FromParent_z = {fReader, "bH2FromParent_z"};
   TTreeReaderArray<float> bH2FromParent_dRjet = {fReader, "bH2FromParent_dRjet"};
   TTreeReaderArray<float> bH1_pt = {fReader, "bH1_pt"};
   TTreeReaderArray<float> bH1_eta = {fReader, "bH1_eta"};
   TTreeReaderArray<float> bH1_phi = {fReader, "bH1_phi"};
   TTreeReaderArray<float> bH1_Lxy = {fReader, "bH1_Lxy"};
   TTreeReaderArray<float> bH1_x = {fReader, "bH1_x"};
   TTreeReaderArray<float> bH1_y = {fReader, "bH1_y"};
   TTreeReaderArray<float> bH1_z = {fReader, "bH1_z"};
   TTreeReaderArray<float> bH1_dRjet = {fReader, "bH1_dRjet"};
   TTreeReaderArray<float> bH2_pt = {fReader, "bH2_pt"};
   TTreeReaderArray<float> bH2_eta = {fReader, "bH2_eta"};
   TTreeReaderArray<float> bH2_phi = {fReader, "bH2_phi"};
   TTreeReaderArray<float> bH2_Lxy = {fReader, "bH2_Lxy"};
   TTreeReaderArray<float> bH2_x = {fReader, "bH2_x"};
   TTreeReaderArray<float> bH2_y = {fReader, "bH2_y"};
   TTreeReaderArray<float> bH2_z = {fReader, "bH2_z"};
   TTreeReaderArray<float> bH2_dRjet = {fReader, "bH2_dRjet"};
   TTreeReaderArray<vector<float>> jet_trk_pt = {fReader, "jet_trk_pt"};
   TTreeReaderArray<vector<float>> jet_trk_eta = {fReader, "jet_trk_eta"};
   TTreeReaderArray<vector<float>> jet_trk_theta = {fReader, "jet_trk_theta"};
   TTreeReaderArray<vector<float>> jet_trk_phi = {fReader, "jet_trk_phi"};
   TTreeReaderArray<vector<float>> jet_trk_qoverp = {fReader, "jet_trk_qoverp"};
   TTreeReaderArray<vector<float>> jet_trk_charge = {fReader, "jet_trk_charge"};
   TTreeReaderArray<vector<float>> jet_trk_chi2 = {fReader, "jet_trk_chi2"};
   TTreeReaderArray<vector<float>> jet_trk_ndf = {fReader, "jet_trk_ndf"};
   TTreeReaderArray<vector<int>> jet_trk_nNextToInnHits = {fReader, "jet_trk_nNextToInnHits"};
   TTreeReaderArray<vector<int>> jet_trk_nInnHits = {fReader, "jet_trk_nInnHits"};
   TTreeReaderArray<vector<int>> jet_trk_nBLHits = {fReader, "jet_trk_nBLHits"};
   TTreeReaderArray<vector<int>> jet_trk_nsharedBLHits = {fReader, "jet_trk_nsharedBLHits"};
   TTreeReaderArray<vector<int>> jet_trk_nsplitBLHits = {fReader, "jet_trk_nsplitBLHits"};
   TTreeReaderArray<vector<int>> jet_trk_nPixHits = {fReader, "jet_trk_nPixHits"};
   TTreeReaderArray<vector<int>> jet_trk_nPixHoles = {fReader, "jet_trk_nPixHoles"};
   TTreeReaderArray<vector<int>> jet_trk_nsharedPixHits = {fReader, "jet_trk_nsharedPixHits"};
   TTreeReaderArray<vector<int>> jet_trk_nsplitPixHits = {fReader, "jet_trk_nsplitPixHits"};
   TTreeReaderArray<vector<int>> jet_trk_nSCTHits = {fReader, "jet_trk_nSCTHits"};
   TTreeReaderArray<vector<int>> jet_trk_nSCTHoles = {fReader, "jet_trk_nSCTHoles"};
   TTreeReaderArray<vector<int>> jet_trk_nsharedSCTHits = {fReader, "jet_trk_nsharedSCTHits"};
   TTreeReaderArray<vector<int>> jet_trk_expectBLayerHit = {fReader, "jet_trk_expectBLayerHit"};
   TTreeReaderArray<vector<int>> jet_trk_expectInnermostPixelLayerHit = {fReader, "jet_trk_expectInnermostPixelLayerHit"};
   TTreeReaderArray<vector<int>> jet_trk_expectNextToInnermostPixelLayerHit = {fReader, "jet_trk_expectNextToInnermostPixelLayerHit"};
   TTreeReaderArray<vector<float>> jet_trk_d0 = {fReader, "jet_trk_d0"};
   TTreeReaderArray<vector<float>> jet_trk_z0 = {fReader, "jet_trk_z0"};
   TTreeReaderArray<vector<float>> jet_trk_ip3d_d0 = {fReader, "jet_trk_ip3d_d0"};
   TTreeReaderArray<vector<float>> jet_trk_ip3d_d02D = {fReader, "jet_trk_ip3d_d02D"};
   TTreeReaderArray<vector<float>> jet_trk_ip3d_z0 = {fReader, "jet_trk_ip3d_z0"};
   TTreeReaderArray<vector<float>> jet_trk_ip3d_d0sig = {fReader, "jet_trk_ip3d_d0sig"};
   TTreeReaderArray<vector<float>> jet_trk_ip3d_z0sig = {fReader, "jet_trk_ip3d_z0sig"};
   TTreeReaderArray<vector<int>> jet_trk_algo = {fReader, "jet_trk_algo"};
   TTreeReaderArray<vector<float>> jet_trk_vtx_X = {fReader, "jet_trk_vtx_X"};
   TTreeReaderArray<vector<float>> jet_trk_vtx_Y = {fReader, "jet_trk_vtx_Y"};
   TTreeReaderArray<vector<float>> jet_trk_vtx_Z = {fReader, "jet_trk_vtx_Z"};
   TTreeReaderArray<vector<int>> jet_trk_pdg_id = {fReader, "jet_trk_pdg_id"};
   TTreeReaderArray<vector<int>> jet_trk_barcode = {fReader, "jet_trk_barcode"};
   TTreeReaderArray<vector<int>> jet_trk_parent_pdgid = {fReader, "jet_trk_parent_pdgid"};
   TTreeReaderArray<int> jet_btag_ntrk = {fReader, "jet_btag_ntrk"};
   TTreeReaderArray<int> jet_ip3d_ntrk = {fReader, "jet_ip3d_ntrk"};
   TTreeReaderArray<vector<int>> jet_trk_ip3d_grade = {fReader, "jet_trk_ip3d_grade"};
   TTreeReaderArray<vector<float>> jet_trk_ip3d_llr = {fReader, "jet_trk_ip3d_llr"};


   DAOD_selector(TTree * /*tree*/ =0) { }
   virtual ~DAOD_selector() { }
   virtual Int_t   Version() const { return 2; }
   virtual void    setFlags(bool, bool, bool, bool, bool, bool, bool, bool, bool, bool);
   virtual void    Begin(TTree *tree);
   virtual void    SlaveBegin(TTree *tree);
   virtual void    Init(TTree *tree);
   virtual Bool_t  Notify();
   virtual Bool_t  Process(Long64_t entry);
   virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
   virtual void    SetOption(const char *option) { fOption = option; }
   virtual void    SetObject(TObject *obj) { fObject = obj; }
   virtual void    SetInputList(TList *input) { fInput = input; }
   virtual TList  *GetOutputList() const { return fOutput; }
   virtual void    SlaveTerminate();
   virtual void    Terminate();

   ClassDef(DAOD_selector,0);

   private:

     bool selections,discriminants,shrinking_cone,selection_alg,origin_selection,geometric_selection,cut,retagT,debug,lxplus;

     double m_cut,m_fc;
     int m_N,m_Ntot,m_b2d,m_b3d,m_c2d,m_c3d,m_noB,m_bb,m_b,m_bc_overlap,m_nbjets,m_nl,m_sc,m_sc2,m_sc3,m_match,m_nomatch,m_match_overlap,m_trk_400,m_trk_B,m_trk_C,m_trk_PU_400,m_trk_FRAG_400,m_trk_GEANT_400;
     int m_qc,m_qj,q,a,b,sc,sgn;
     double D_phi,D_eta,DR,px,py,Dx_1,Dy_1,Dz_1,Dx_2,Dy_2,Dz_2,Dxy_1,x0,y0,Dx_3,Dy_3,Dxy_3,rand_n,R0,d0,c,A,gamma;//,nx=0,ny=0;
     double D_phi_trk,D_eta_trk,DR_trk,DpT_trk;
     int match,max_size;
     double tmp_pTfraction,tmp_DR,tmp_min_pTfraction,tmp_min_DR,m_pTfraction_cut,m_DRcut,m_pTfraction_nocut,m_DRnocut;
     unsigned size_jet,size_child;
     int den,m_den;

     unsigned m_track_cut;

     float pt_max, pt_min;
     float Delta;

     std::vector< std::vector<float> > bin_v = std::vector< std::vector<float> >(bin_1);
     TLorentzVector jet;

     TFile *file;
  //   TGraph *g = new TGraph ();
  //   TGraphErrors *g_E = new TGraphErrors ();
     TH2F *jet_DR_pT;

     TH1F *hist_pt_1;
     TH1F *hist_eta_1;
     TH1F *hist_phi_1;
     TH1F *hist_E_1;
     TH1F *hist_pt_2;
     TH1F *hist_eta_2;
     TH1F *hist_phi_2;
     TH1F *hist_E_2;
     TH1F *hist_pt_2b;
     TH1F *hist_eta_2b;
     TH1F *hist_phi_2b;
     TH1F *hist_E_2b;
     TH1F *hist_pt_3a;
     TH1F *hist_eta_3a;
     TH1F *hist_phi_3a;
     TH1F *hist_E_3a;
     TH1F *hist_pt_3b;
     TH1F *hist_eta_3b;
     TH1F *hist_phi_3b;
     TH1F *hist_E_3b;
     TH1F *hist_pt_4;
     TH1F *hist_eta_4;
     TH1F *hist_phi_4;
     TH1F *hist_E_4;

     TH1F *hist_pt_inB;
     TH1F *hist_eta_inB;
     TH1F *hist_phi_inB;
     TH1F *hist_E_inB;

     TH1F *hist_n_tracks;
     TH2F *hist_tracks_DR;
     TH1F *hist_DR_1;
  //   TH1F *hist_std_dev_DR_1;


     TH1F *hist_trk_pT_inB;
     TH1F *hist_trk_Deta_inB;
     TH1F *hist_trk_eta_inB;
     TH1F *hist_trk_Dphi_inB;
     TH1F *hist_trk_phi_inB;
     TH2F *hist_trk_Dphi_Deta_inB;
     TH1F *hist_trk_DR_inB;
     TH2F *hist_trk_pT_DR_inB;
     TH2F *hist_trk_pT_jet_DR_inB;     
     TH1F *hist_trk_origin_inB;
     TH1F *hist_trk_d0_inB;
     TH1F *hist_trk_d0_PUinB;
     TH1F *hist_trk_d0_BinB;
     TH1F *hist_trk_d0_CinB;
     TH1F *hist_trk_d0_FRAGinB;
     TH1F *hist_trk_d0_GEANTinB;

     TH1F *hist_child_pT_inB;
     TH1F *hist_child_Deta_inB;
     TH1F *hist_child_eta_inB;
     TH1F *hist_child_Dphi_inB;
     TH1F *hist_child_phi_inB;
     TH2F *hist_child_Dphi_Deta_inB;
     TH1F *hist_child_DR_inB;
     TH2F *hist_child_pT_DR_inB;
     TH2F *hist_child_pT_jet_DR_inB;
     TH1F *hist_child_Lxy_inB;
     TH1F *hist_child_Lxyz_inB;
     TH1F *hist_child_decay_IP;
     TH1F *hist_child_nodecay_IP;
  //   TH1F *hist_child_linear_IP;
     TH1F *hist_child_d0_truth;
     TH1F *hist_child_d0;
     TH2F *hist_child_d0_pT;
     TH1F *hist_pT_vs_R0_ratio_inB;
  //   TH1F *hist_child_linearIP;

     TH1F *hist_efficiency_inB;
     TH1D *hist_n_child;
     TH1D *hist_n_trk;
     TH1D *hist_n_match;

     TH1F *hist_matched_origin_pT_inB;
     TH1F *hist_matched_origin_eta_inB;
     TH1F *hist_matched_origin_phi_inB;
     TH1F *hist_matched_origin_Deta_inB;
     TH1F *hist_matched_origin_Dphi_inB;
     TH2F *hist_matched_origin_Dphi_Deta_inB;
     TH1F *hist_matched_origin_DR_inB;
     TH2F *hist_matched_origin_pT_DR_inB;
     TH2F *hist_matched_origin_pT_jet_DR_inB;
     TH1F *hist_matched_origin_pdgId_inB;
     TH1F *hist_matched_origin_origin_inB;
     TH1F *hist_matched_origin_d0_inB;
//     TH1F *hist_matched_origin_Lxy_inB;
//     TH1F *hist_matched_origin_Lxyz_inB;

     TH1F *hist_matched_pT_inB;
     TH1F *hist_matched_eta_inB;
     TH1F *hist_matched_phi_inB;
     TH1F *hist_matched_Deta_inB;
     TH1F *hist_matched_Dphi_inB;
     TH2F *hist_matched_Dphi_Deta_inB;
     TH1F *hist_matched_DR_inB;
     TH2F *hist_matched_pT_DR_inB;
     TH2F *hist_matched_pT_jet_DR_inB;
     TH1F *hist_matched_pdgId_inB;
     TH1F *hist_matched_origin_inB;
     TH2F *hist_matched_pT_child_pTfraction_inB;
     TH1F *hist_matched_DR_trk_inB;
     TH2F *hist_matched_DR_trk_pTfraction;
     TH1F *hist_matched_d0_inB;
     TH1F *hist_matched_Lxy_inB;
     TH1F *hist_matched_Lxyz_inB;

     TH1F *hist_nomatched_pT_inB;
     TH1F *hist_nomatched_eta_inB;
     TH1F *hist_nomatched_phi_inB;
     TH1F *hist_nomatched_Deta_inB;
     TH1F *hist_nomatched_Dphi_inB;
     TH2F *hist_nomatched_Dphi_Deta_inB;
     TH1F *hist_nomatched_DR_inB;
     TH2F *hist_nomatched_pT_DR_inB;
     TH2F *hist_nomatched_pT_jet_DR_inB;
     TH1F *hist_nomatched_pdgId_inB;

     TH1F *hist_single_matched_pT_inB;
     TH1F *hist_single_matched_eta_inB;
     TH1F *hist_single_matched_phi_inB;
     TH1F *hist_single_matched_Deta_inB;
     TH1F *hist_single_matched_Dphi_inB;
     TH2F *hist_single_matched_Dphi_Deta_inB;
     TH1F *hist_single_matched_DR_inB;
     TH2F *hist_single_matched_pT_DR_inB;
     TH2F *hist_single_matched_pT_jet_DR_inB;
     TH1F *hist_single_matched_pdgId_inB;
     TH1F *hist_single_matched_origin_inB;
     TH2F *hist_single_matched_pT_child_pTfraction_inB;
     TH1F *hist_single_matched_DR_trk_inB;
     TH2F *hist_single_matched_DR_trk_pTfraction;
     TH1F *hist_single_matched_d0_inB;




};

#endif

#ifdef DAOD_selector_cxx
void DAOD_selector::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the reader is initialized.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   fReader.SetTree(tree);
}

Bool_t DAOD_selector::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}


#endif // #ifdef DAOD_selector_cxx
