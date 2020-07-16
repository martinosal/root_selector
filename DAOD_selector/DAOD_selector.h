//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon May 25 17:19:26 2020 by ROOT version 6.18/04
// from TTree bTag_AntiKt4EMPFlowJets_BTagging201903/bTagAntiKt4EMPFlowJets_BTagging201903
// found on file: flav_Akt4EMPf_BTagging201903_IP_doRetagF.root
//////////////////////////////////////////////////////////

#ifndef DAOD_selector_h
#define DAOD_selector_h

#define bin_1 50
#define tracksize 11

#include "functions.h"

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
   TTreeReaderArray<double> jet_dl1_pb = {fReader, "jet_dl1_pb"};//<-TaggerScoreBranches.cxx
   TTreeReaderArray<double> jet_dl1_pc = {fReader, "jet_dl1_pc"};//<-TaggerScoreBranches.cxx
   TTreeReaderArray<double> jet_dl1_pu = {fReader, "jet_dl1_pu"};//<-TaggerScoreBranches.cxx
   TTreeReaderArray<double> jet_mv2c10 = {fReader, "jet_mv2c10"};//<-TaggerScoreBranches.cxx
   TTreeReaderArray<double> jet_mv2c10mu = {fReader, "jet_mv2c10mu"};//<-TaggerScoreBranches.cxx
   TTreeReaderArray<double> jet_mv2c10rnn = {fReader, "jet_mv2c10rnn"};//<-TaggerScoreBranches.cxx
   TTreeReaderArray<double> jet_mv2c100 = {fReader, "jet_mv2c100"};//<-TaggerScoreBranches.cxx
   TTreeReaderArray<double> jet_mv2cl100 = {fReader, "jet_mv2cl100"};//<-TaggerScoreBranches.cxx
   TTreeReaderArray<float> jet_ip2d_pb = {fReader, "jet_ip2d_pb"};//<-ImpactParameterBranches.cxx
   TTreeReaderArray<float> jet_ip2d_pc = {fReader, "jet_ip2d_pc"};//<-ImpactParameterBranches.cxx
   TTreeReaderArray<float> jet_ip2d_pu = {fReader, "jet_ip2d_pu"};//<-ImpactParameterBranches.cxx
   TTreeReaderArray<float> jet_ip2d_llr = {fReader, "jet_ip2d_llr"};//<-ImpactParameterBranches.cxx
   TTreeReaderArray<float> jet_ip3d_pb = {fReader, "jet_ip3d_pb"};//<-ImpactParameterBranches.cxx
   TTreeReaderArray<float> jet_ip3d_pc = {fReader, "jet_ip3d_pc"};//<-ImpactParameterBranches.cxx
   TTreeReaderArray<float> jet_ip3d_pu = {fReader, "jet_ip3d_pu"};//<-ImpactParameterBranches.cxx
   TTreeReaderArray<float> jet_ip3d_llr = {fReader, "jet_ip3d_llr"};//<-ImpactParameterBranches.cxx
   TTreeReaderArray<float> jet_ip2 = {fReader, "jet_ip2"};//<-ImpactParameterBranches.cxx
   TTreeReaderArray<float> jet_ip2_c = {fReader, "jet_ip2_c"};//<-ImpactParameterBranches.cxx
   TTreeReaderArray<float> jet_ip2_cu = {fReader, "jet_ip2_cu"};//<-ImpactParameterBranches.cxx
   TTreeReaderArray<float> jet_ip2_nan = {fReader, "jet_ip2_nan"};//<-ImpactParameterBranches.cxx
   TTreeReaderArray<float> jet_ip2_c_nan = {fReader, "jet_ip2_c_nan"};//<-ImpactParameterBranches.cxx
   TTreeReaderArray<float> jet_ip2_cu_nan = {fReader, "jet_ip2_cu_nan"};//<-ImpactParameterBranches.cxx
   TTreeReaderArray<float> jet_ip3 = {fReader, "jet_ip3"};//<-ImpactParameterBranches.cxx
   TTreeReaderArray<float> jet_ip3_c = {fReader, "jet_ip3_c"};//<-ImpactParameterBranches.cxx
   TTreeReaderArray<float> jet_ip3_cu = {fReader, "jet_ip3_cu"};//<-ImpactParameterBranches.cxx
   TTreeReaderArray<float> jet_ip3_nan = {fReader, "jet_ip3_nan"};//<-ImpactParameterBranches.cxx
   TTreeReaderArray<float> jet_ip3_c_nan = {fReader, "jet_ip3_c_nan"};//<-ImpactParameterBranches.cxx
   TTreeReaderArray<float> jet_ip3_cu_nan = {fReader, "jet_ip3_cu_nan"};//<-ImpactParameterBranches.cxx
   TTreeReaderArray<float> jet_rnnip_pb = {fReader, "jet_rnnip_pb"};//<-ImpactParameterBranches.cxx
   TTreeReaderArray<float> jet_rnnip_pc = {fReader, "jet_rnnip_pc"};//<-ImpactParameterBranches.cxx
   TTreeReaderArray<float> jet_rnnip_pu = {fReader, "jet_rnnip_pu"};//<-ImpactParameterBranches.cxx
   TTreeReaderArray<float> jet_rnnip_ptau = {fReader, "jet_rnnip_ptau"};//<-ImpactParameterBranches.cxx
   TTreeReaderValue<Float_t> PV_jf_x = {fReader, "PV_jf_x"};//<-JetFitterBranches.cxx
   TTreeReaderValue<Float_t> PV_jf_y = {fReader, "PV_jf_y"};//<-JetFitterBranches.cxx
   TTreeReaderValue<Float_t> PV_jf_z = {fReader, "PV_jf_z"};//<-JetFitterBranches.cxx
   TTreeReaderArray<vector<int>> jet_trk_jf_Vertex = {fReader, "jet_trk_jf_Vertex"};//<-JetFitterBranches.cxx
   TTreeReaderArray<float> jet_jf_pb = {fReader, "jet_jf_pb"};//<-JetFitterBranches.cxx
   TTreeReaderArray<float> jet_jf_pc = {fReader, "jet_jf_pc"};//<-JetFitterBranches.cxx
   TTreeReaderArray<float> jet_jf_pu = {fReader, "jet_jf_pu"};//<-JetFitterBranches.cxx
   TTreeReaderArray<float> jet_jf_llr = {fReader, "jet_jf_llr"};//<-JetFitterBranches.cxx
   TTreeReaderArray<float> jet_jf_m = {fReader, "jet_jf_m"};//<-JetFitterBranches.cxx
   TTreeReaderArray<float> jet_jf_mUncorr = {fReader, "jet_jf_mUncorr"};//<-JetFitterBranches.cxx
   TTreeReaderArray<float> jet_jf_efc = {fReader, "jet_jf_efc"};//<-JetFitterBranches.cxx
   TTreeReaderArray<float> jet_jf_deta = {fReader, "jet_jf_deta"};//<-JetFitterBranches.cxx
   TTreeReaderArray<float> jet_jf_dphi = {fReader, "jet_jf_dphi"};//<-JetFitterBranches.cxx
   TTreeReaderArray<float> jet_jf_dR = {fReader, "jet_jf_dR"};//<-JetFitterBranches.cxx
   TTreeReaderArray<float> jet_jf_dRFlightDir = {fReader, "jet_jf_dRFlightDir"};//<-JetFitterBranches.cxx
   TTreeReaderArray<float> jet_jf_ntrkAtVx = {fReader, "jet_jf_ntrkAtVx"};//<-JetFitterBranches.cxx
   TTreeReaderArray<float> jet_jf_nvtx = {fReader, "jet_jf_nvtx"};//<-JetFitterBranches.cxx
   TTreeReaderArray<float> jet_jf_sig3d = {fReader, "jet_jf_sig3d"};//<-JetFitterBranches.cxx
   TTreeReaderArray<float> jet_jf_nvtx1t = {fReader, "jet_jf_nvtx1t"};//<-JetFitterBranches.cxx
   TTreeReaderArray<float> jet_jf_n2t = {fReader, "jet_jf_n2t"};//<-JetFitterBranches.cxx
   TTreeReaderArray<float> jet_jf_VTXsize = {fReader, "jet_jf_VTXsize"};//<-JetFitterBranches.cxx
   TTreeReaderArray<vector<float>> jet_jf_vtx_chi2 = {fReader, "jet_jf_vtx_chi2"};//<-JetFitterBranches.cxx
   TTreeReaderArray<vector<float>> jet_jf_vtx_ndf = {fReader, "jet_jf_vtx_ndf"};//<-JetFitterBranches.cxx
   TTreeReaderArray<vector<int>> jet_jf_vtx_ntrk = {fReader, "jet_jf_vtx_ntrk"};//<-JetFitterBranches.cxx
   TTreeReaderArray<vector<float>> jet_jf_vtx_L3D = {fReader, "jet_jf_vtx_L3D"};//<-JetFitterBranches.cxx
   TTreeReaderArray<vector<float>> jet_jf_vtx_sig3D = {fReader, "jet_jf_vtx_sig3D"};//<-JetFitterBranches.cxx
   TTreeReaderArray<float> jet_jf_phi = {fReader, "jet_jf_phi"};//<-JetFitterBranches.cxx
   TTreeReaderArray<float> jet_jf_theta = {fReader, "jet_jf_theta"};//<-JetFitterBranches.cxx
   TTreeReaderArray<vector<float>> jet_jf_vtx_sigTrans = {fReader, "jet_jf_vtx_sigTrans"};//<-JetFitterBranches.cxx
   TTreeReaderArray<vector<float>> jet_jf_vtx_x = {fReader, "jet_jf_vtx_x"};//<-JetFitterBranches.cxx
   TTreeReaderArray<vector<float>> jet_jf_vtx_x_err = {fReader, "jet_jf_vtx_x_err"};//<-JetFitterBranches.cxx
   TTreeReaderArray<vector<float>> jet_jf_vtx_y = {fReader, "jet_jf_vtx_y"};//<-JetFitterBranches.cxx
   TTreeReaderArray<vector<float>> jet_jf_vtx_y_err = {fReader, "jet_jf_vtx_y_err"};//<-JetFitterBranches.cxx
   TTreeReaderArray<vector<float>> jet_jf_vtx_z = {fReader, "jet_jf_vtx_z"};//<-JetFitterBranches.cxx
   TTreeReaderArray<vector<float>> jet_jf_vtx_z_err = {fReader, "jet_jf_vtx_z_err"};//<-JetFitterBranches.cxx
   TTreeReaderArray<float> jet_jf_theta_err = {fReader, "jet_jf_theta_err"};//<-JetFitterBranches.cxx
   TTreeReaderArray<float> jet_jf_phi_err = {fReader, "jet_jf_phi_err"};//<-JetFitterBranches.cxx
   TTreeReaderArray<float> nTrk_vtx1 = {fReader, "nTrk_vtx1"};//<-JetFitterBranches.cxx
   TTreeReaderArray<float> mass_first_vtx = {fReader, "mass_first_vtx"};//<-JetFitterBranches.cxx
   TTreeReaderArray<float> e_first_vtx = {fReader, "e_first_vtx"};//<-JetFitterBranches.cxx
   TTreeReaderArray<float> e_frac_vtx1 = {fReader, "e_frac_vtx1"};//<-JetFitterBranches.cxx
   TTreeReaderArray<float> closestVtx_L3D = {fReader, "closestVtx_L3D"};//<-JetFitterBranches.cxx
   TTreeReaderArray<float> JF_Lxy1 = {fReader, "JF_Lxy1"};//<-JetFitterBranches.cxx
   TTreeReaderArray<float> vtx1_MaxTrkRapidity = {fReader, "vtx1_MaxTrkRapidity"};//<-JetFitterBranches.cxx
   TTreeReaderArray<float> vtx1_AvgTrkRapidity = {fReader, "vtx1_AvgTrkRapidity"};//<-JetFitterBranches.cxx
   TTreeReaderArray<float> vtx1_MinTrkRapidity = {fReader, "vtx1_MinTrkRapidity"};//<-JetFitterBranches.cxx
   TTreeReaderArray<float> nTrk_vtx2 = {fReader, "nTrk_vtx2"};//<-JetFitterBranches.cxx
   TTreeReaderArray<float> mass_second_vtx = {fReader, "mass_second_vtx"};//<-JetFitterBranches.cxx//SV mass?
   TTreeReaderArray<float> e_second_vtx = {fReader, "e_second_vtx"};//<-JetFitterBranches.cxx
   TTreeReaderArray<float> e_frac_vtx2 = {fReader, "e_frac_vtx2"};//<-JetFitterBranches.cxx//SV energy fraction?
   TTreeReaderArray<float> second_closestVtx_L3D = {fReader, "second_closestVtx_L3D"};//<-JetFitterBranches.cxx
   TTreeReaderArray<float> JF_Lxy2 = {fReader, "JF_Lxy2"};//<-JetFitterBranches.cxx
   TTreeReaderArray<float> vtx2_MaxTrkRapidity = {fReader, "vtx2_MaxTrkRapidity"};//<-JetFitterBranches.cxx
   TTreeReaderArray<float> vtx2_AvgTrkRapidity = {fReader, "vtx2_AvgTrkRapidity"};//<-JetFitterBranches.cxx
   TTreeReaderArray<float> vtx2_MinTrkRapidity = {fReader, "vtx2_MinTrkRapidity"};//<-JetFitterBranches.cxx
   TTreeReaderArray<float> MaxTrkRapidity = {fReader, "MaxTrkRapidity"};//<-JetFitterBranches.cxx
   TTreeReaderArray<float> MinTrkRapidity = {fReader, "MinTrkRapidity"};//<-JetFitterBranches.cxx
   TTreeReaderArray<float> AvgTrkRapidity = {fReader, "AvgTrkRapidity"};//<-JetFitterBranches.cxx
   TTreeReaderArray<int> jet_sv1_Nvtx = {fReader, "jet_sv1_Nvtx"};//<-SVBranches.cxx
   TTreeReaderArray<float> jet_sv1_ntrkv = {fReader, "jet_sv1_ntrkv"};//<-SVBranches.cxx
   TTreeReaderArray<float> jet_sv1_n2t = {fReader, "jet_sv1_n2t"};//<-SVBranches.cxx//number of two-track vertices
   TTreeReaderArray<float> jet_sv1_m = {fReader, "jet_sv1_m"};//<-SVBranches.cxx//SV mass
   TTreeReaderArray<float> jet_sv1_efc = {fReader, "jet_sv1_efc"};//<-SVBranches.cxx//SV energy fraction
   TTreeReaderArray<float> jet_sv1_sig3d = {fReader, "jet_sv1_sig3d"};//<-SVBranches.cxx
   TTreeReaderArray<float> sv1_llr = {fReader, "sv1_llr"};//<-SVBranches.cxx
   TTreeReaderArray<float> jet_sv1_normdist = {fReader, "jet_sv1_normdist"};//<-SVBranches.cxx
   TTreeReaderArray<float> jet_sv1_deltaR = {fReader, "jet_sv1_deltaR"};//<-SVBranches.cxx
   TTreeReaderArray<float> jet_sv1_Lxy = {fReader, "jet_sv1_Lxy"};//<-SVBranches.cxx
   TTreeReaderArray<float> jet_sv1_L3d = {fReader, "jet_sv1_L3d"};//<-SVBranches.cxx
   TTreeReaderArray<vector<float>> jet_sv1_vtx_x = {fReader, "jet_sv1_vtx_x"};//<-SVBranches.cxx
   TTreeReaderArray<vector<float>> jet_sv1_vtx_y = {fReader, "jet_sv1_vtx_y"};//<-SVBranches.cxx
   TTreeReaderArray<vector<float>> jet_sv1_vtx_z = {fReader, "jet_sv1_vtx_z"};//<-SVBranches.cxx
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
   TTreeReaderArray<vector<float>> jet_bH_child_theta = {fReader, "jet_bH_child_theta"};
   TTreeReaderArray<vector<float>> jet_bH_child_prod_x = {fReader, "jet_bH_child_prod_x"};
   TTreeReaderArray<vector<float>> jet_bH_child_prod_y = {fReader, "jet_bH_child_prod_y"};
   TTreeReaderArray<vector<float>> jet_bH_child_prod_z = {fReader, "jet_bH_child_prod_z"};
   TTreeReaderArray<vector<float>> jet_bH_child_decay_x = {fReader, "jet_bH_child_decay_x"};
   TTreeReaderArray<vector<float>> jet_bH_child_decay_y = {fReader, "jet_bH_child_decay_y"};
   TTreeReaderArray<vector<float>> jet_bH_child_decay_z = {fReader, "jet_bH_child_decay_z"};
   TTreeReaderArray<vector<float>> jet_bH_child_d0 = {fReader, "jet_bH_child_d0"};
   TTreeReaderArray<vector<float>> jet_bH_child_z0 = {fReader, "jet_bH_child_z0"};
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
   TTreeReaderArray<vector<float>> jet_trk_d0sig = {fReader, "jet_trk_d0sig"};
   TTreeReaderArray<vector<float>> jet_trk_z0 = {fReader, "jet_trk_z0"};
   TTreeReaderArray<vector<float>> jet_trk_z0sig = {fReader, "jet_trk_z0sig"};
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
   virtual void    setCuts(float, float, float, float, float, float, float, float, float);
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

   void initFlagsAndCuts();
   void initCounters();
   void openOutputFile(std::string fileNameStringID="");

   void bookHistosForSelections();
   void bookHistosForDiscriminants();
   void bookHistosForShrinkingCone();
   void bookHistosForSelectionAlgos();

   void getTrueJetFlavourLabel(std::vector<int>& isJet, std::vector<int>& isJetB, std::vector<int>& isJetC, std::vector<int>& isJetl);



   ClassDef(DAOD_selector,0);

   private:

     // Flags selectong the running mode
     bool selections,discriminants,shrinking_cone,selection_alg,origin_selection,geometric_selection,cut,retagT,debug,lxplus;

     // selection cuts
     float jet_pT_infcut,jet_pT_supcut,jet_eta_cut,jet_JVT_cut,m_pT_bcH_truth_cut,m_DR_bcH_truth_cut;
     float trk_pT_cut,trk_eta_cut,trk_d0_cut;

     // Service variables
     double m_pt_max_shrCone, m_pt_min_shrCone, m_Delta_pt_shrCone;

     //float pt_bH,DeltaR_bH,pt_cH,DeltaR_cH;
     double m_cut,m_fc,m_fcRNNIP;
     int m_Ntot,m_noB,m_bb,m_b,m_bc_overlap,m_sc,m_sc2,m_sc3,m_match,m_nomatch,m_match_overlap,m_match_notoverlap,n_trk_pT_cut,n_trk_B,n_trk_C,n_trk_PU_pT_cut,n_trk_FRAG_pT_cut,n_trk_GEANT_pT_cut;
     int m_qc,m_qj,q,a,b,sc,sgn,sgn_d0,sgn_z0;
     //double D_phi,D_eta,DR,px,py,pz,Dx_1,Dy_1,Dz_1,Dx_2,Dy_2,Dz_2,Lxy,Lxyz,Dxy_1,x0,y0,Dx_3,Dy_3,Dxy_3,rand_n,R0,d0;
     double c,A,gamma;//,nx=0,ny=0;
     //double D_phi_trk,D_eta_trk,DR_trk,DpT_trk;
     int mm,m1,m2,m1_ex,m2_ex,mm1_ex,mm2_ex,m1_ov,m2_ov,m_GeomNOr_PU,m_GeomNOr_F,m_GeomNOr_G;
     //double tmp_pTfraction,tmp_DR,tmp_min_pTfraction,tmp_min_DR;
     double m_pTfraction_cut,m_DRcut,m_pTfraction_nocut,m_DRnocut;
     //unsigned size_jet,size_child;
     int m_den;
     // counters for jets
     int m_njets,m_njets_2,m_nBjets,m_nCjets,m_nljets,m_nBjets_2,m_nCjets_2,m_nljets_2;
       //cut flow on jets
     int m_njets_2_passPtMin;
     int m_njets_2_passPtMax;
     int m_njets_2_passEtaRange;
     int m_njets_2_passBadMedium;
     int m_njets_2_passOR;
     int m_njets_2_passORmu;
     int m_njets_2_passJVT;

     int m_nBcheck,m_nCcheck,m_nlcheck;
     int m_nJetBCoverlap,m_nJetBCoverlap_postJetSel,ov_check;
     int JF_ntrk,SV1_ntrk,SV0_ntrk,IP2D_ntrk,IP3D_ntrk;

     unsigned m_track_cut;

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

     TH1F *hist_pt_B;
     TH1F *hist_eta_B;
     TH1F *hist_phi_B;
     TH1F *hist_E_B;
     TH1F *hist_Bjet_origin;

     TH1F *hist_pt_C;
     TH1F *hist_eta_C;
     TH1F *hist_phi_C;
     TH1F *hist_E_C;
     TH1F *hist_Cjet_origin;

     TH1F *hist_pt_l;
     TH1F *hist_eta_l;
     TH1F *hist_phi_l;
     TH1F *hist_E_l;
     TH1F *hist_ljet_origin;

     TH1F *hist_n_tracks;
     TH2F *hist_tracks_DR;
     TH1F *hist_DR_1;
  //   TH1F *hist_std_dev_DR_1;


     TH1F *hist_jet_IP2_B;

     TH1F *hist_ip2d_pb;
     TH1F *hist_ip2d_pc;
     TH1F *hist_ip2d_pu;
     TH1F *hist_ip2d_llr;
     TH1F *hist_ip3d_pb;
     TH1F *hist_ip3d_pc;
     TH1F *hist_ip3d_pu;
     TH1F *hist_ip3d_llr;
     TH1F *hist_dl1_pb;
     TH1F *hist_dl1_pc;
     TH1F *hist_dl1_pu;

     TH1F *hist_ip2d_llr_l;
     TH1F *hist_ip3d_llr_l;
     TH1F *hist_rnnip_llr_l;
     TH1F *hist_sv1_llr_l;
     TH1F *hist_jf_llr_l;
     TH1F *hist_ip2d_llr_B;
     TH1F *hist_ip3d_llr_B;
     TH1F *hist_rnnip_llr_B;
     TH1F *hist_sv1_llr_B;
     TH1F *hist_jf_llr_B;
     TH1F *hist_ip2d_llr_C;
     TH1F *hist_ip3d_llr_C;
     TH1F *hist_rnnip_llr_C;
     TH1F *hist_sv1_llr_C;
     TH1F *hist_jf_llr_C;

     TH1F *hist_dl1_l;
     TH1F *hist_dl1_C;
     TH1F *hist_dl1_B;

     TH1F *hist_trk_pT_B;
     TH1F *hist_trk_Deta_B;
     TH1F *hist_trk_eta_B;
     TH1F *hist_trk_Dphi_B;
     TH1F *hist_trk_phi_B;
     TH2F *hist_trk_Dphi_Deta_B;
     TH1F *hist_trk_DR_B;
     TH2F *hist_trk_pT_DR_B;
     TH2F *hist_trk_pT_jet_DR_B;
     TH1F *hist_trk_pdgId_B;
     TH1F *hist_trk_origin_B;

     TH1F *hist_trk_d0_B;
     TH1F *hist_trk_sigd0_B;
     TH1F *hist_trk_z0_B;
     TH1F *hist_trk_sigz0_B;
     TH1F *hist_trk_z0sinth_B;
     TH1F *hist_trk_d0sig_B;
     TH1F *hist_trk_z0sinthsig_B;
     TH2F *hist_trk_d0sig_origin_B;
     TH2F *hist_trk_z0sinthsig_origin_B;
     TH2F *hist_trk_logpTfrac_origin_B;
     TH2F *hist_trk_logDR_origin_B;
     TH2F *hist_trk_IBLhits_origin_B;
     TH2F *hist_trk_NextToIBLhits_origin_B;
     TH2F *hist_trk_sharedIBLhits_origin_B;
     TH2F *hist_trk_splitIBLhits_origin_B;
     TH2F *hist_trk_nPixhits_origin_B;
     TH2F *hist_trk_sharedPixhits_origin_B;
     TH2F *hist_trk_splitPixhits_origin_B;
     TH2F *hist_trk_nSCThits_origin_B;
     TH2F *hist_trk_sharedSCThits_origin_B;

     TH1F *hist_trk_origin_JF_B;
     TH1F *hist_trk_d0sig_JF_B;
     TH1F *hist_trk_z0sinthsig_JF_B;
     TH2F *hist_trk_d0sig_origin_JF_B;
     TH2F *hist_trk_z0sinthsig_origin_JF_B;
     TH2F *hist_trk_logpTfrac_origin_JF_B;
     TH2F *hist_trk_logDR_origin_JF_B;
     TH2F *hist_trk_IBLhits_origin_JF_B;
     TH2F *hist_trk_NextToIBLhits_origin_JF_B;
     TH2F *hist_trk_sharedIBLhits_origin_JF_B;
     TH2F *hist_trk_splitIBLhits_origin_JF_B;
     TH2F *hist_trk_nPixhits_origin_JF_B;
     TH2F *hist_trk_sharedPixhits_origin_JF_B;
     TH2F *hist_trk_splitPixhits_origin_JF_B;
     TH2F *hist_trk_nSCThits_origin_JF_B;
     TH2F *hist_trk_sharedSCThits_origin_JF_B;

     TH1F *hist_trk_origin_SV1_B;
     TH1F *hist_trk_d0sig_SV1_B;
     TH1F *hist_trk_z0sinthsig_SV1_B;
     TH2F *hist_trk_d0sig_origin_SV1_B;
     TH2F *hist_trk_z0sinthsig_origin_SV1_B;
     TH2F *hist_trk_logpTfrac_origin_SV1_B;
     TH2F *hist_trk_logDR_origin_SV1_B;
     TH2F *hist_trk_IBLhits_origin_SV1_B;
     TH2F *hist_trk_NextToIBLhits_origin_SV1_B;
     TH2F *hist_trk_sharedIBLhits_origin_SV1_B;
     TH2F *hist_trk_splitIBLhits_origin_SV1_B;
     TH2F *hist_trk_nPixhits_origin_SV1_B;
     TH2F *hist_trk_sharedPixhits_origin_SV1_B;
     TH2F *hist_trk_splitPixhits_origin_SV1_B;
     TH2F *hist_trk_nSCThits_origin_SV1_B;
     TH2F *hist_trk_sharedSCThits_origin_SV1_B;

     TH1F *hist_trk_origin_SV0_B;
     TH1F *hist_trk_d0sig_SV0_B;
     TH1F *hist_trk_z0sinthsig_SV0_B;
     TH2F *hist_trk_d0sig_origin_SV0_B;
     TH2F *hist_trk_z0sinthsig_origin_SV0_B;
     TH2F *hist_trk_logpTfrac_origin_SV0_B;
     TH2F *hist_trk_logDR_origin_SV0_B;
     TH2F *hist_trk_IBLhits_origin_SV0_B;
     TH2F *hist_trk_NextToIBLhits_origin_SV0_B;
     TH2F *hist_trk_sharedIBLhits_origin_SV0_B;
     TH2F *hist_trk_splitIBLhits_origin_SV0_B;
     TH2F *hist_trk_nPixhits_origin_SV0_B;
     TH2F *hist_trk_sharedPixhits_origin_SV0_B;
     TH2F *hist_trk_splitPixhits_origin_SV0_B;
     TH2F *hist_trk_nSCThits_origin_SV0_B;
     TH2F *hist_trk_sharedSCThits_origin_SV0_B;

     TH1F *hist_trk_origin_IP3D_B;
     TH1F *hist_trk_d0sig_IP3D_B;
     TH1F *hist_trk_z0sinthsig_IP3D_B;
     TH2F *hist_trk_d0sig_origin_IP3D_B;
     TH2F *hist_trk_z0sinthsig_origin_IP3D_B;
     TH2F *hist_trk_logpTfrac_origin_IP3D_B;
     TH2F *hist_trk_logDR_origin_IP3D_B;
     TH2F *hist_trk_IBLhits_origin_IP3D_B;
     TH2F *hist_trk_NextToIBLhits_origin_IP3D_B;
     TH2F *hist_trk_sharedIBLhits_origin_IP3D_B;
     TH2F *hist_trk_splitIBLhits_origin_IP3D_B;
     TH2F *hist_trk_nPixhits_origin_IP3D_B;
     TH2F *hist_trk_sharedPixhits_origin_IP3D_B;
     TH2F *hist_trk_splitPixhits_origin_IP3D_B;
     TH2F *hist_trk_nSCThits_origin_IP3D_B;
     TH2F *hist_trk_sharedSCThits_origin_IP3D_B;

     TH1F *hist_trk_origin_IP2D_B;
     TH1F *hist_trk_d0sig_IP2D_B;
     TH1F *hist_trk_z0sinthsig_IP2D_B;
     TH2F *hist_trk_d0sig_origin_IP2D_B;
     TH2F *hist_trk_z0sinthsig_origin_IP2D_B;
     TH2F *hist_trk_logpTfrac_origin_IP2D_B;
     TH2F *hist_trk_logDR_origin_IP2D_B;
     TH2F *hist_trk_IBLhits_origin_IP2D_B;
     TH2F *hist_trk_NextToIBLhits_origin_IP2D_B;
     TH2F *hist_trk_sharedIBLhits_origin_IP2D_B;
     TH2F *hist_trk_splitIBLhits_origin_IP2D_B;
     TH2F *hist_trk_nPixhits_origin_IP2D_B;
     TH2F *hist_trk_sharedPixhits_origin_IP2D_B;
     TH2F *hist_trk_splitPixhits_origin_IP2D_B;
     TH2F *hist_trk_nSCThits_origin_IP2D_B;
     TH2F *hist_trk_sharedSCThits_origin_IP2D_B;

     TH1F *hist_trk_d0_PUB;
     TH1F *hist_trk_d0_BB;
     TH1F *hist_trk_d0_CB;
     TH1F *hist_trk_d0_FRAGB;
     TH1F *hist_trk_d0_GEANTB;

     TH1F *hist_child_pT_B;
     TH1F *hist_child_Deta_B;
     TH1F *hist_child_eta_B;
     TH1F *hist_child_Dphi_B;
     TH1F *hist_child_phi_B;
     TH2F *hist_child_Dphi_Deta_B;
     TH1F *hist_child_DR_B;
     TH2F *hist_child_pT_DR_B;
     TH2F *hist_child_pT_jet_DR_B;
     TH1F *hist_child_pdgID_B;

     TH1F *hist_child_pi_notD;
     TH1F *hist_child_K_notD;
     TH1F *hist_child_pi;
     TH1F *hist_child_pi_Lxy_B;
     TH1F *hist_child_pi_Lxyz_B;
     TH1F *hist_child_pi_d0_truth_B;
     TH1F *hist_child_pi_z0_truth_B;
     TH1F *hist_child_K;
     TH1F *hist_child_K_Lxy_B;
     TH1F *hist_child_K_Lxyz_B;
     TH1F *hist_child_K_d0_truth_B;
     TH1F *hist_child_K_z0_truth_B;
     TH1F *hist_child_mu;
     TH1F *hist_child_mu_Lxy_B;
     TH1F *hist_child_mu_Lxyz_B;
     TH1F *hist_child_mu_d0_truth_B;
     TH1F *hist_child_mu_z0_truth_B;
     TH1F *hist_child_p;
     TH1F *hist_child_p_Lxy_B;
     TH1F *hist_child_p_Lxyz_B;
     TH1F *hist_child_p_d0_truth_B;
     TH1F *hist_child_p_z0_truth_B;
     TH1F *hist_child_e;
     TH1F *hist_child_e_Lxy_B;
     TH1F *hist_child_e_Lxyz_B;
     TH1F *hist_child_e_d0_truth_B;
     TH1F *hist_child_e_z0_truth_B;

     TH1F *hist_child_Lxy_B;
     TH1F *hist_child_Lxyz_B;
     TH1F *hist_child_decay_IP;
     TH1F *hist_child_nodecay_IP;
  //   TH1F *hist_child_linear_IP;
     TH1F *hist_child_d0_truth;
     TH1F *hist_child_d0;
     TH2F *hist_child_d0_pT;
     TH1F *hist_child_z0sinth_B;
     TH1F *hist_pT_vs_R0_ratio_B;
  //   TH1F *hist_child_linearIP;

     TH1F *hist_efficiency_B;
     TH1D *hist_n_child;
     TH1D *hist_n_trk;
     TH1D *hist_n_match;

     TH1F *hist_matched_origin_pT_B;
     TH1F *hist_matched_origin_eta_B;
     TH1F *hist_matched_origin_phi_B;
     TH1F *hist_matched_origin_Deta_B;
     TH1F *hist_matched_origin_Dphi_B;
     TH2F *hist_matched_origin_Dphi_Deta_B;
     TH1F *hist_matched_origin_DR_B;
     TH2F *hist_matched_origin_pT_DR_B;
     TH2F *hist_matched_origin_pT_jet_DR_B;
     TH1F *hist_matched_origin_pdgId_B;
     TH1F *hist_matched_origin_origin_B;
     TH1F *hist_matched_origin_d0_B;
//     TH1F *hist_matched_origin_Lxy_B;
//     TH1F *hist_matched_origin_Lxyz_B;

     TH1F *hist_matched_pT_B;
     TH1F *hist_matched_eta_B;
     TH1F *hist_matched_phi_B;
     TH1F *hist_matched_Deta_B;
     TH1F *hist_matched_Dphi_B;
     TH2F *hist_matched_Dphi_Deta_B;
     TH1F *hist_matched_DR_B;
     TH2F *hist_matched_pT_DR_B;
     TH2F *hist_matched_pT_jet_DR_B;
     TH1F *hist_matched_pdgId_B;

     TH1F *hist_matched_child_pi_notD;
     TH1F *hist_matched_child_K_notD;
     TH1F *hist_matched_child_pi;
     TH1F *hist_matched_child_pi_Lxy_B;
     TH1F *hist_matched_child_pi_Lxyz_B;
     TH1F *hist_matched_child_pi_d0_truth_B;
     TH1F *hist_matched_child_pi_z0_truth_B;
     TH1F *hist_matched_child_K;
     TH1F *hist_matched_child_K_Lxy_B;
     TH1F *hist_matched_child_K_Lxyz_B;
     TH1F *hist_matched_child_K_d0_truth_B;
     TH1F *hist_matched_child_K_z0_truth_B;
     TH1F *hist_matched_child_mu;
     TH1F *hist_matched_child_mu_Lxy_B;
     TH1F *hist_matched_child_mu_Lxyz_B;
     TH1F *hist_matched_child_mu_d0_truth_B;
     TH1F *hist_matched_child_mu_z0_truth_B;
     TH1F *hist_matched_child_p;
     TH1F *hist_matched_child_p_Lxy_B;
     TH1F *hist_matched_child_p_Lxyz_B;
     TH1F *hist_matched_child_p_d0_truth_B;
     TH1F *hist_matched_child_p_z0_truth_B;
     TH1F *hist_matched_child_e;
     TH1F *hist_matched_child_e_Lxy_B;
     TH1F *hist_matched_child_e_Lxyz_B;
     TH1F *hist_matched_child_e_d0_truth_B;
     TH1F *hist_matched_child_e_z0_truth_B;

     TH1F *hist_matched_origin_B;
     TH2F *hist_matched_pT_child_pTfraction_B;
     TH1F *hist_matched_DR_trk_B;
     TH2F *hist_matched_DR_trk_pTfraction;
     TH1F *hist_matched_d0_B;
     TH1F *hist_matched_Lxy_B;
     TH1F *hist_matched_Lxyz_B;

     TH1F *hist_nomatched_pT_B;
     TH1F *hist_nomatched_eta_B;
     TH1F *hist_nomatched_phi_B;
     TH1F *hist_nomatched_Deta_B;
     TH1F *hist_nomatched_Dphi_B;
     TH2F *hist_nomatched_Dphi_Deta_B;
     TH1F *hist_nomatched_DR_B;
     TH2F *hist_nomatched_pT_DR_B;
     TH2F *hist_nomatched_pT_jet_DR_B;
     TH1F *hist_nomatched_pdgId_B;

     TH1F *hist_single_matched_pT_B;
     TH1F *hist_single_matched_eta_B;
     TH1F *hist_single_matched_phi_B;
     TH1F *hist_single_matched_Deta_B;
     TH1F *hist_single_matched_Dphi_B;
     TH2F *hist_single_matched_Dphi_Deta_B;
     TH1F *hist_single_matched_DR_B;
     TH2F *hist_single_matched_pT_DR_B;
     TH2F *hist_single_matched_pT_jet_DR_B;
     TH1F *hist_single_matched_pdgId_B;
     TH1F *hist_single_matched_origin_B;
     TH2F *hist_single_matched_pT_child_pTfraction_B;
     TH1F *hist_single_matched_DR_trk_B;
     TH2F *hist_single_matched_DR_trk_pTfraction;
     TH1F *hist_single_matched_d0_B;




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
