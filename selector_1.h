//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sun Mar 15 17:00:23 2020 by ROOT version 6.18/04
// from TTree bTag_AntiKt4EMTopoJets/bTagAntiKt4EMTopoJets
// found on file: flav_Akt4EMTo_13022127_taggers.root
//////////////////////////////////////////////////////////

#ifndef selector_1_h
#define selector_1_h
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
//#include "matrix.h"

//#include </usr/include/eigen3/Eigen/Dense>

class selector_1 : public TSelector {
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
   TTreeReaderArray<float> jet_pt = {fReader, "jet_pt"};//<-JetPropertiesBranches.cxx
   TTreeReaderArray<float> jet_eta = {fReader, "jet_eta"};//<-JetPropertiesBranches.cxx
   TTreeReaderArray<float> jet_phi = {fReader, "jet_phi"};//<-JetPropertiesBranches.cxx
   TTreeReaderArray<float> jet_E = {fReader, "jet_E"};//<-JetPropertiesBranches.cxx
   TTreeReaderArray<float> jet_pt_orig = {fReader, "jet_pt_orig"};//<-JetPropertiesBranches.cxx
   TTreeReaderArray<float> jet_eta_orig = {fReader, "jet_eta_orig"};//<-JetPropertiesBranches.cxx
   TTreeReaderArray<float> jet_phi_orig = {fReader, "jet_phi_orig"};//<-JetPropertiesBranches.cxx
   TTreeReaderArray<float> jet_E_orig = {fReader, "jet_E_orig"};//<-JetPropertiesBranches.cxx
   TTreeReaderArray<int> jet_LabDr_HadF = {fReader, "jet_LabDr_HadF"};//<-JetPropertiesBranches.cxx
   TTreeReaderArray<int> jet_DoubleHadLabel = {fReader, "jet_DoubleHadLabel"};//<-JetPropertiesBranches.cxx
   TTreeReaderArray<float> jet_JVT = {fReader, "jet_JVT"};//<-JetPropertiesBranches.cxx
   TTreeReaderArray<float> jet_m = {fReader, "jet_m"};//<-JetPropertiesBranches.cxx
   TTreeReaderArray<float> jet_nConst = {fReader, "jet_nConst"};//<-JetPropertiesBranches.cxx
   TTreeReaderArray<float> jet_dRiso = {fReader, "jet_dRiso"};//<-JetPropertiesBranches.cxx
   TTreeReaderArray<int> jet_truthMatch = {fReader, "jet_truthMatch"};//<-JetPropertiesBranches.cxx
   TTreeReaderArray<int> jet_isPU = {fReader, "jet_isPU"};//<-JetPropertiesBranches.cxx
   TTreeReaderArray<int> jet_aliveAfterOR = {fReader, "jet_aliveAfterOR"};//<-JetPropertiesBranches.cxx
   TTreeReaderArray<int> jet_aliveAfterORmu = {fReader, "jet_aliveAfterORmu"};//<-JetPropertiesBranches.cxx
   TTreeReaderArray<int> jet_isBadMedium = {fReader, "jet_isBadMedium"};//<-JetPropertiesBranches.cxx
   TTreeReaderArray<float> jet_truthPt = {fReader, "jet_truthPt"};//<-JetPropertiesBranches.cxx
   TTreeReaderArray<float> jet_dRminToB = {fReader, "jet_dRminToB"};//<-JetPropertiesBranches.cxx
   TTreeReaderArray<float> jet_dRminToC = {fReader, "jet_dRminToC"};//<-JetPropertiesBranches.cxx
   TTreeReaderArray<float> jet_dRminToT = {fReader, "jet_dRminToT"};//<-JetPropertiesBranches.cxx
   TTreeReaderArray<double> jet_dl1_pb = {fReader, "jet_dl1_pb"};//<-TaggerScoreBranches.cxx
   TTreeReaderArray<double> jet_dl1_pc = {fReader, "jet_dl1_pc"};//<-TaggerScoreBranches.cxx
   TTreeReaderArray<double> jet_dl1_pu = {fReader, "jet_dl1_pu"};//<-TaggerScoreBranches.cxx
   TTreeReaderArray<double> jet_dl1mu_pb = {fReader, "jet_dl1mu_pb"};//<-TaggerScoreBranches.cxx
   TTreeReaderArray<double> jet_dl1mu_pc = {fReader, "jet_dl1mu_pc"};//<-TaggerScoreBranches.cxx
   TTreeReaderArray<double> jet_dl1mu_pu = {fReader, "jet_dl1mu_pu"};//<-TaggerScoreBranches.cxx
   TTreeReaderArray<double> jet_dl1rnn_pc = {fReader, "jet_dl1rnn_pc"};//<-TaggerScoreBranches.cxx
   TTreeReaderArray<double> jet_dl1rnn_pb = {fReader, "jet_dl1rnn_pb"};//<-TaggerScoreBranches.cxx
   TTreeReaderArray<double> jet_dl1rnn_pu = {fReader, "jet_dl1rnn_pu"};//<-TaggerScoreBranches.cxx
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
   TTreeReaderArray<double> jet_mu_smt = {fReader, "jet_mu_smt"};
   TTreeReaderArray<float> jet_mu_dR = {fReader, "jet_mu_dR"};
   TTreeReaderArray<float> jet_mu_pTrel = {fReader, "jet_mu_pTrel"};
   TTreeReaderArray<float> jet_mu_qOverPratio = {fReader, "jet_mu_qOverPratio"};
   TTreeReaderArray<float> jet_mu_mombalsignif = {fReader, "jet_mu_mombalsignif"};
   TTreeReaderArray<float> jet_mu_scatneighsignif = {fReader, "jet_mu_scatneighsignif"};
   TTreeReaderArray<float> jet_mu_VtxTyp = {fReader, "jet_mu_VtxTyp"};
   TTreeReaderArray<float> jet_mu_pt = {fReader, "jet_mu_pt"};
   TTreeReaderArray<float> jet_mu_eta = {fReader, "jet_mu_eta"};
   TTreeReaderArray<float> jet_mu_phi = {fReader, "jet_mu_phi"};
   TTreeReaderArray<float> jet_mu_d0 = {fReader, "jet_mu_d0"};
   TTreeReaderArray<float> jet_mu_z0 = {fReader, "jet_mu_z0"};
   TTreeReaderArray<float> jet_mu_parent_pdgid = {fReader, "jet_mu_parent_pdgid"};
   TTreeReaderArray<float> jet_mu_ID_qOverP_var = {fReader, "jet_mu_ID_qOverP_var"};
   TTreeReaderArray<float> jet_mu_muonType = {fReader, "jet_mu_muonType"};
   TTreeReaderArray<int> jet_nBHadr = {fReader, "jet_nBHadr"};//<-BHadronBranches.cxx/REGULAR
   TTreeReaderArray<int> jet_nCHadr = {fReader, "jet_nCHadr"};//<-BHadronBranches.cxx/REGULAR
   TTreeReaderArray<vector<int>> jet_bH_pdgId = {fReader, "jet_bH_pdgId"};//<-BHadronBranches.cxx/REGULAR
   TTreeReaderArray<vector<int>> jet_bH_parent_pdgId = {fReader, "jet_bH_parent_pdgId"};//<-BHadronBranches.cxx/REGULAR
   TTreeReaderArray<vector<float>> jet_bH_pt = {fReader, "jet_bH_pt"};//<-BHadronBranches.cxx/REGULAR
   TTreeReaderArray<vector<float>> jet_bH_eta = {fReader, "jet_bH_eta"};//<-BHadronBranches.cxx/REGULAR
   TTreeReaderArray<vector<float>> jet_bH_phi = {fReader, "jet_bH_phi"};//<-BHadronBranches.cxx/REGULAR
   TTreeReaderArray<vector<float>> jet_bH_E = {fReader, "jet_bH_E"};//<-BHadronBranches.cxx/REGULAR
   TTreeReaderArray<vector<float>> jet_bH_charge = {fReader, "jet_bH_charge"};//<-BHadronBranches.cxx/REGULAR
   TTreeReaderArray<vector<float>> jet_bH_Lxy = {fReader, "jet_bH_Lxy"};//<-BHadronBranches.cxx/REGULAR
   TTreeReaderArray<vector<float>> jet_bH_x = {fReader, "jet_bH_x"};//<-BHadronBranches.cxx/REGULAR
   TTreeReaderArray<vector<float>> jet_bH_y = {fReader, "jet_bH_y"};//<-BHadronBranches.cxx/REGULAR
   TTreeReaderArray<vector<float>> jet_bH_z = {fReader, "jet_bH_z"};//<-BHadronBranches.cxx/REGULAR
   TTreeReaderArray<vector<float>> jet_bH_dRjet = {fReader, "jet_bH_dRjet"};//<-BHadronBranches.cxx/REGULAR
   TTreeReaderArray<vector<int>> jet_cH_pdgId = {fReader, "jet_cH_pdgId"};//<-BHadronBranches.cxx/REGULAR
   TTreeReaderArray<vector<int>> jet_cH_parent_pdgId = {fReader, "jet_cH_parent_pdgId"};//<-BHadronBranches.cxx/REGULAR
   TTreeReaderArray<vector<float>> jet_cH_pt = {fReader, "jet_cH_pt"};//<-BHadronBranches.cxx/REGULAR
   TTreeReaderArray<vector<float>> jet_cH_eta = {fReader, "jet_cH_eta"};//<-BHadronBranches.cxx/REGULAR
   TTreeReaderArray<vector<float>> jet_cH_phi = {fReader, "jet_cH_phi"};//<-BHadronBranches.cxx/REGULAR
   TTreeReaderArray<vector<float>> jet_cH_E = {fReader, "jet_cH_E"};//<-BHadronBranches.cxx/REGULAR
   TTreeReaderArray<vector<float>> jet_cH_charge = {fReader, "jet_cH_charge"};//<-BHadronBranches.cxx/REGULAR
   TTreeReaderArray<vector<float>> jet_cH_Lxy = {fReader, "jet_cH_Lxy"};//<-BHadronBranches.cxx/REGULAR
   TTreeReaderArray<vector<float>> jet_cH_x = {fReader, "jet_cH_x"};//<-BHadronBranches.cxx/REGULAR
   TTreeReaderArray<vector<float>> jet_cH_y = {fReader, "jet_cH_y"};//<-BHadronBranches.cxx/REGULAR
   TTreeReaderArray<vector<float>> jet_cH_z = {fReader, "jet_cH_z"};//<-BHadronBranches.cxx/REGULAR
   TTreeReaderArray<vector<float>> jet_cH_dRjet = {fReader, "jet_cH_dRjet"};//<-BHadronBranches.cxx/REGULAR
   TTreeReaderArray<vector<int>> jet_trk_orig = {fReader, "jet_trk_orig"};//<-BHadronBranches.cxx/REGULAR
   TTreeReaderArray<vector<float>> jet_bH_prod_x = {fReader, "jet_bH_prod_x"};//<-BHadronBranches.cxx/REGULAR
   TTreeReaderArray<vector<float>> jet_bH_prod_y = {fReader, "jet_bH_prod_y"};//<-BHadronBranches.cxx/REGULAR
   TTreeReaderArray<vector<float>> jet_bH_prod_z = {fReader, "jet_bH_prod_z"};//<-BHadronBranches.cxx/REGULAR
   TTreeReaderArray<vector<float>> jet_bH_PtTrk = {fReader, "jet_bH_PtTrk"};//<-BHadronBranches.cxx/REGULAR
   TTreeReaderArray<vector<float>> jet_bH_MTrk = {fReader, "jet_bH_MTrk"};//<-BHadronBranches.cxx/REGULAR
   TTreeReaderArray<vector<int>> jet_bH_nBtracks = {fReader, "jet_bH_nBtracks"};//<-BHadronBranches.cxx/REGULAR
   TTreeReaderArray<vector<int>> jet_bH_nCtracks = {fReader, "jet_bH_nCtracks"};//<-BHadronBranches.cxx/REGULAR
   TTreeReaderArray<vector<int>> jet_bH_nBtracks_400 = {fReader, "jet_bH_nBtracks_400"};//<-BHadronBranches.cxx/REGULAR
   TTreeReaderArray<vector<int>> jet_bH_nCtracks_400 = {fReader, "jet_bH_nCtracks_400"};//<-BHadronBranches.cxx/REGULAR
   TTreeReaderArray<vector<int>> jet_bH_child_hadron_idx = {fReader, "jet_bH_child_hadron_idx"};//<-BHadronBranches.cxx/REGULAR
   TTreeReaderArray<vector<int>> jet_bH_child_pdg_id = {fReader, "jet_bH_child_pdg_id"};//<-BHadronBranches.cxx/REGULAR
   TTreeReaderArray<vector<int>> jet_bH_child_parent_pdg_id = {fReader, "jet_bH_child_parent_pdg_id"};//<-BHadronBranches.cxx/REGULAR
   TTreeReaderArray<vector<int>> jet_bH_child_barcode = {fReader, "jet_bH_child_barcode"};//<-BHadronBranches.cxx/REGULAR
   TTreeReaderArray<vector<float>> jet_bH_child_charge = {fReader, "jet_bH_child_charge"};//<-BHadronBranches.cxx/REGULAR
   TTreeReaderArray<vector<float>> jet_bH_child_px = {fReader, "jet_bH_child_px"};//<-BHadronBranches.cxx/REGULAR
   TTreeReaderArray<vector<float>> jet_bH_child_py = {fReader, "jet_bH_child_py"};//<-BHadronBranches.cxx/REGULAR
   TTreeReaderArray<vector<float>> jet_bH_child_pz = {fReader, "jet_bH_child_pz"};//<-BHadronBranches.cxx/REGULAR
   TTreeReaderArray<vector<float>> jet_bH_child_E = {fReader, "jet_bH_child_E"};//<-BHadronBranches.cxx/REGULAR
   TTreeReaderArray<vector<float>> jet_bH_child_prod_x = {fReader, "jet_bH_child_prod_x"};//<-BHadronBranches.cxx/REGULAR
   TTreeReaderArray<vector<float>> jet_bH_child_prod_y = {fReader, "jet_bH_child_prod_y"};//<-BHadronBranches.cxx/REGULAR
   TTreeReaderArray<vector<float>> jet_bH_child_prod_z = {fReader, "jet_bH_child_prod_z"};//<-BHadronBranches.cxx/REGULAR
   TTreeReaderArray<vector<float>> jet_bH_child_decay_x = {fReader, "jet_bH_child_decay_x"};//<-BHadronBranches.cxx/REGULAR
   TTreeReaderArray<vector<float>> jet_bH_child_decay_y = {fReader, "jet_bH_child_decay_y"};//<-BHadronBranches.cxx/REGULAR
   TTreeReaderArray<vector<float>> jet_bH_child_decay_z = {fReader, "jet_bH_child_decay_z"};//<-BHadronBranches.cxx/REGULAR
   TTreeReaderArray<vector<float>> jet_cH_prod_x = {fReader, "jet_cH_prod_x"};//<-BHadronBranches.cxx/REGULAR
   TTreeReaderArray<vector<float>> jet_cH_prod_y = {fReader, "jet_cH_prod_y"};//<-BHadronBranches.cxx/REGULAR
   TTreeReaderArray<vector<float>> jet_cH_prod_z = {fReader, "jet_cH_prod_z"};//<-BHadronBranches.cxx/REGULAR
   TTreeReaderArray<vector<float>> jet_cH_PtTrk = {fReader, "jet_cH_PtTrk"};//<-BHadronBranches.cxx/REGULAR
   TTreeReaderArray<vector<float>> jet_cH_MTrk = {fReader, "jet_cH_MTrk"};//<-BHadronBranches.cxx/REGULAR
   TTreeReaderArray<vector<int>> jet_cH_nCtracks = {fReader, "jet_cH_nCtracks"};//<-BHadronBranches.cxx/REGULAR
   TTreeReaderArray<vector<int>> jet_cH_nCtracks_400 = {fReader, "jet_cH_nCtracks_400"};//<-BHadronBranches.cxx/REGULAR
   TTreeReaderArray<vector<int>> jet_cH_child_hadron_idx = {fReader, "jet_cH_child_hadron_idx"};//<-BHadronBranches.cxx/REGULAR
   TTreeReaderArray<vector<int>> jet_cH_child_pdg_id = {fReader, "jet_cH_child_pdg_id"};//<-BHadronBranches.cxx/REGULAR
   TTreeReaderArray<vector<int>> jet_cH_child_parent_pdg_id = {fReader, "jet_cH_child_parent_pdg_id"};//<-BHadronBranches.cxx/REGULAR
   TTreeReaderArray<vector<int>> jet_cH_child_barcode = {fReader, "jet_cH_child_barcode"};//<-BHadronBranches.cxx/REGULAR
   TTreeReaderArray<vector<float>> jet_cH_child_charge = {fReader, "jet_cH_child_charge"};//<-BHadronBranches.cxx/REGULAR
   TTreeReaderArray<vector<float>> jet_cH_child_px = {fReader, "jet_cH_child_px"};//<-BHadronBranches.cxx/REGULAR
   TTreeReaderArray<vector<float>> jet_cH_child_py = {fReader, "jet_cH_child_py"};//<-BHadronBranches.cxx/REGULAR
   TTreeReaderArray<vector<float>> jet_cH_child_pz = {fReader, "jet_cH_child_pz"};//<-BHadronBranches.cxx/REGULAR
   TTreeReaderArray<vector<float>> jet_cH_child_E = {fReader, "jet_cH_child_E"};//<-BHadronBranches.cxx/REGULAR
   TTreeReaderArray<vector<float>> jet_cH_child_prod_x = {fReader, "jet_cH_child_prod_x"};//<-BHadronBranches.cxx/REGULAR
   TTreeReaderArray<vector<float>> jet_cH_child_prod_y = {fReader, "jet_cH_child_prod_y"};//<-BHadronBranches.cxx/REGULAR
   TTreeReaderArray<vector<float>> jet_cH_child_prod_z = {fReader, "jet_cH_child_prod_z"};//<-BHadronBranches.cxx/REGULAR
   TTreeReaderArray<vector<float>> jet_cH_child_decay_x = {fReader, "jet_cH_child_decay_x"};//<-BHadronBranches.cxx/REGULAR
   TTreeReaderArray<vector<float>> jet_cH_child_decay_y = {fReader, "jet_cH_child_decay_y"};//<-BHadronBranches.cxx/REGULAR
   TTreeReaderArray<vector<float>> jet_cH_child_decay_z = {fReader, "jet_cH_child_decay_z"};//<-BHadronBranches.cxx/REGULAR
   TTreeReaderArray<int> jet_nGhostBHadrFromParent = {fReader, "jet_nGhostBHadrFromParent"};//<-BHadronBranches.cxx/DOUBLE
   TTreeReaderArray<int> jet_nGhostCHadrFromParent = {fReader, "jet_nGhostCHadrFromParent"};//<-BHadronBranches.cxx/DOUBLE
   TTreeReaderArray<int> jet_nGhostCHadrFromParentNotFromB = {fReader, "jet_nGhostCHadrFromParentNotFromB"};//<-BHadronBranches.cxx/DOUBLE
   TTreeReaderArray<int> jet_nGhostTauFromParent = {fReader, "jet_nGhostTauFromParent"};//<-BHadronBranches.cxx/DOUBLE
   TTreeReaderArray<int> jet_nGhostHBosoFromParent = {fReader, "jet_nGhostHBosoFromParent"};//<-BHadronBranches.cxx/DOUBLE
   TTreeReaderArray<int> jet_nGhostBHadr = {fReader, "jet_nGhostBHadr"};//<-BHadronBranches.cxx/DOUBLE
   TTreeReaderArray<int> jet_nGhostCHadr = {fReader, "jet_nGhostCHadr"};//<-BHadronBranches.cxx/DOUBLE
   TTreeReaderArray<int> jet_nGhostCHadrNotFromB = {fReader, "jet_nGhostCHadrNotFromB"};//<-BHadronBranches.cxx/DOUBLE
   TTreeReaderArray<int> jet_nGhostTau = {fReader, "jet_nGhostTau"};//<-BHadronBranches.cxx/DOUBLE
   TTreeReaderArray<int> jet_nGhostHBoso = {fReader, "jet_nGhostHBoso"};//<-BHadronBranches.cxx/DOUBLE
   TTreeReaderArray<float> bH1FromParent_pt = {fReader, "bH1FromParent_pt"};//<-BHadronBranches.cxx/DOUBLE
   TTreeReaderArray<float> bH1FromParent_eta = {fReader, "bH1FromParent_eta"};//<-BHadronBranches.cxx/DOUBLE
   TTreeReaderArray<float> bH1FromParent_phi = {fReader, "bH1FromParent_phi"};//<-BHadronBranches.cxx/DOUBLE
   TTreeReaderArray<float> bH1FromParent_Lxy = {fReader, "bH1FromParent_Lxy"};//<-BHadronBranches.cxx/DOUBLE
   TTreeReaderArray<float> bH1FromParent_x = {fReader, "bH1FromParent_x"};//<-BHadronBranches.cxx/DOUBLE
   TTreeReaderArray<float> bH1FromParent_y = {fReader, "bH1FromParent_y"};//<-BHadronBranches.cxx/DOUBLE
   TTreeReaderArray<float> bH1FromParent_z = {fReader, "bH1FromParent_z"};//<-BHadronBranches.cxx/DOUBLE
   TTreeReaderArray<float> bH1FromParent_dRjet = {fReader, "bH1FromParent_dRjet"};//<-BHadronBranches.cxx/DOUBLE
   TTreeReaderArray<float> bH2FromParent_pt = {fReader, "bH2FromParent_pt"};//<-BHadronBranches.cxx/DOUBLE
   TTreeReaderArray<float> bH2FromParent_eta = {fReader, "bH2FromParent_eta"};//<-BHadronBranches.cxx/DOUBLE
   TTreeReaderArray<float> bH2FromParent_phi = {fReader, "bH2FromParent_phi"};//<-BHadronBranches.cxx/DOUBLE
   TTreeReaderArray<float> bH2FromParent_Lxy = {fReader, "bH2FromParent_Lxy"};//<-BHadronBranches.cxx/DOUBLE
   TTreeReaderArray<float> bH2FromParent_x = {fReader, "bH2FromParent_x"};//<-BHadronBranches.cxx/DOUBLE
   TTreeReaderArray<float> bH2FromParent_y = {fReader, "bH2FromParent_y"};//<-BHadronBranches.cxx/DOUBLE
   TTreeReaderArray<float> bH2FromParent_z = {fReader, "bH2FromParent_z"};//<-BHadronBranches.cxx/DOUBLE
   TTreeReaderArray<float> bH2FromParent_dRjet = {fReader, "bH2FromParent_dRjet"};//<-BHadronBranches.cxx/DOUBLE
   TTreeReaderArray<float> bH1_pt = {fReader, "bH1_pt"};//<-BHadronBranches.cxx/DOUBLE
   TTreeReaderArray<float> bH1_eta = {fReader, "bH1_eta"};//<-BHadronBranches.cxx/DOUBLE
   TTreeReaderArray<float> bH1_phi = {fReader, "bH1_phi"};//<-BHadronBranches.cxx/DOUBLE
   TTreeReaderArray<float> bH1_Lxy = {fReader, "bH1_Lxy"};//<-BHadronBranches.cxx/DOUBLE
   TTreeReaderArray<float> bH1_x = {fReader, "bH1_x"};//<-BHadronBranches.cxx/DOUBLE
   TTreeReaderArray<float> bH1_y = {fReader, "bH1_y"};//<-BHadronBranches.cxx/DOUBLE
   TTreeReaderArray<float> bH1_z = {fReader, "bH1_z"};//<-BHadronBranches.cxx/DOUBLE
   TTreeReaderArray<float> bH1_dRjet = {fReader, "bH1_dRjet"};//<-BHadronBranches.cxx/DOUBLE
   TTreeReaderArray<float> bH2_pt = {fReader, "bH2_pt"};//<-BHadronBranches.cxx/DOUBLE
   TTreeReaderArray<float> bH2_eta = {fReader, "bH2_eta"};//<-BHadronBranches.cxx/DOUBLE
   TTreeReaderArray<float> bH2_phi = {fReader, "bH2_phi"};//<-BHadronBranches.cxx/DOUBLE
   TTreeReaderArray<float> bH2_Lxy = {fReader, "bH2_Lxy"};//<-BHadronBranches.cxx/DOUBLE
   TTreeReaderArray<float> bH2_x = {fReader, "bH2_x"};//<-BHadronBranches.cxx/DOUBLE
   TTreeReaderArray<float> bH2_y = {fReader, "bH2_y"};//<-BHadronBranches.cxx/DOUBLE
   TTreeReaderArray<float> bH2_z = {fReader, "bH2_z"};//<-BHadronBranches.cxx/DOUBLE
   TTreeReaderArray<float> bH2_dRjet = {fReader, "bH2_dRjet"};//<-BHadronBranches.cxx/DOUBLE
   TTreeReaderArray<vector<float>> jet_trk_pt = {fReader, "jet_trk_pt"};//<-TrackBranches.cxx
   TTreeReaderArray<vector<float>> jet_trk_eta = {fReader, "jet_trk_eta"};//<-TrackBranches.cxx
   TTreeReaderArray<vector<float>> jet_trk_theta = {fReader, "jet_trk_theta"};//<-TrackBranches.cxx
   TTreeReaderArray<vector<float>> jet_trk_phi = {fReader, "jet_trk_phi"};//<-TrackBranches.cxx
   TTreeReaderArray<vector<float>> jet_trk_qoverp = {fReader, "jet_trk_qoverp"};//<-TrackBranches.cxx
   TTreeReaderArray<vector<float>> jet_trk_charge = {fReader, "jet_trk_charge"};//<-TrackBranches.cxx
   TTreeReaderArray<vector<float>> jet_trk_chi2 = {fReader, "jet_trk_chi2"};//<-TrackBranches.cxx
   TTreeReaderArray<vector<float>> jet_trk_ndf = {fReader, "jet_trk_ndf"};//<-TrackBranches.cxx
   TTreeReaderArray<vector<int>> jet_trk_nNextToInnHits = {fReader, "jet_trk_nNextToInnHits"};//<-TrackBranches.cxx
   TTreeReaderArray<vector<int>> jet_trk_nInnHits = {fReader, "jet_trk_nInnHits"};//<-TrackBranches.cxx
   TTreeReaderArray<vector<int>> jet_trk_nBLHits = {fReader, "jet_trk_nBLHits"};//<-TrackBranches.cxx
   TTreeReaderArray<vector<int>> jet_trk_nsharedBLHits = {fReader, "jet_trk_nsharedBLHits"};//<-TrackBranches.cxx
   TTreeReaderArray<vector<int>> jet_trk_nsplitBLHits = {fReader, "jet_trk_nsplitBLHits"};//<-TrackBranches.cxx
   TTreeReaderArray<vector<int>> jet_trk_nPixHits = {fReader, "jet_trk_nPixHits"};//<-TrackBranches.cxx
   TTreeReaderArray<vector<int>> jet_trk_nPixHoles = {fReader, "jet_trk_nPixHoles"};//<-TrackBranches.cxx
   TTreeReaderArray<vector<int>> jet_trk_nsharedPixHits = {fReader, "jet_trk_nsharedPixHits"};//<-TrackBranches.cxx
   TTreeReaderArray<vector<int>> jet_trk_nsplitPixHits = {fReader, "jet_trk_nsplitPixHits"};//<-TrackBranches.cxx
   TTreeReaderArray<vector<int>> jet_trk_nSCTHits = {fReader, "jet_trk_nSCTHits"};//<-TrackBranches.cxx
   TTreeReaderArray<vector<int>> jet_trk_nSCTHoles = {fReader, "jet_trk_nSCTHoles"};//<-TrackBranches.cxx
   TTreeReaderArray<vector<int>> jet_trk_nsharedSCTHits = {fReader, "jet_trk_nsharedSCTHits"};//<-TrackBranches.cxx
   TTreeReaderArray<vector<int>> jet_trk_expectBLayerHit = {fReader, "jet_trk_expectBLayerHit"};//<-TrackBranches.cxx
   TTreeReaderArray<vector<int>> jet_trk_expectInnermostPixelLayerHit = {fReader, "jet_trk_expectInnermostPixelLayerHit"};//<-TrackBranches.cxx
   TTreeReaderArray<vector<int>> jet_trk_expectNextToInnermostPixelLayerHit = {fReader, "jet_trk_expectNextToInnermostPixelLayerHit"};//<-TrackBranches.cxx
   TTreeReaderArray<vector<float>> jet_trk_d0 = {fReader, "jet_trk_d0"};//<-TrackBranches.cxx
   TTreeReaderArray<vector<float>> jet_trk_z0 = {fReader, "jet_trk_z0"};//<-TrackBranches.cxx
   TTreeReaderArray<vector<float>> jet_trk_ip3d_d0 = {fReader, "jet_trk_ip3d_d0"};//<-TrackBranches.cxx
   TTreeReaderArray<vector<float>> jet_trk_ip3d_d02D = {fReader, "jet_trk_ip3d_d02D"};//<-TrackBranches.cxx
   TTreeReaderArray<vector<float>> jet_trk_ip3d_z0 = {fReader, "jet_trk_ip3d_z0"};//<-TrackBranches.cxx
   TTreeReaderArray<vector<float>> jet_trk_ip3d_d0sig = {fReader, "jet_trk_ip3d_d0sig"};//<-TrackBranches.cxx
   TTreeReaderArray<vector<float>> jet_trk_ip3d_z0sig = {fReader, "jet_trk_ip3d_z0sig"};//<-TrackBranches.cxx
   TTreeReaderArray<vector<int>> jet_trk_algo = {fReader, "jet_trk_algo"};//<-TrackBranches.cxx
   TTreeReaderArray<vector<float>> jet_trk_vtx_X = {fReader, "jet_trk_vtx_X"};//<-TrackBranches.cxx
   TTreeReaderArray<vector<float>> jet_trk_vtx_Y = {fReader, "jet_trk_vtx_Y"};//<-TrackBranches.cxx
   TTreeReaderArray<vector<float>> jet_trk_vtx_Z = {fReader, "jet_trk_vtx_Z"};//<-TrackBranches.cxx
   TTreeReaderArray<vector<int>> jet_trk_pdg_id = {fReader, "jet_trk_pdg_id"};//<-TrackBranches.cxx
   TTreeReaderArray<vector<int>> jet_trk_barcode = {fReader, "jet_trk_barcode"};//<-TrackBranches.cxx
   TTreeReaderArray<vector<int>> jet_trk_parent_pdgid = {fReader, "jet_trk_parent_pdgid"};//<-TrackBranches.cxx
   TTreeReaderArray<int> jet_btag_ntrk = {fReader, "jet_btag_ntrk"};//<-TrackBranches.cxx
   TTreeReaderArray<int> jet_ip3d_ntrk = {fReader, "jet_ip3d_ntrk"};//<-TrackBranches.cxx
   TTreeReaderArray<vector<int>> jet_trk_ip3d_grade = {fReader, "jet_trk_ip3d_grade"};//<-TrackBranches.cxx
   TTreeReaderArray<vector<float>> jet_trk_ip3d_llr = {fReader, "jet_trk_ip3d_llr"};//<-TrackBranches.cxx


   selector_1(TTree * /*tree*/ =0) { }
   virtual ~selector_1() { }
   virtual Int_t   Version() const { return 2; }
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

   ClassDef(selector_1,0);

private:

   bool selections=1,discriminants=0,shrinking_cone=0,selection_alg=1,cut=1,retagT=0,lxplus=0,debug=1;

   double m_cut=1.,m_fc=0.08;
   int m_N=0,m_Ntot=0,m_b2d=0,m_b3d=0,m_c2d=0,m_c3d=0,m_noB=0,m_bb=0,m_b=0,m_bc_overlap=0,m_nbjets=0,m_nl=0,m_sc=0,m_sc2=0,m_sc3=0,m_match=0,m_nomatch=0,m_truth_match=0,m_truth_tot=0;
   int m_qc=0,m_qj=0,q=0,a=0,b=0,sc=0,sgn=0;
   double D_phi=0.,D_eta=0.,DR=0.,px=0.,py=0.,Dx_1=0.,Dy_1=0.,Dz_1=0.,Dx_2=0.,Dy_2=0.,Dz_2=0,Dxy_1=0,x0=0,y0=0,vx=0,vy=0,v_pa=0,Dx_3=0.,Dy_3=0.,Dxy_3=0.,rand_n=0.,R0=0,d0=0,c=2.99792458e8;//,nx=0,ny=0;
   double D_phi_trk=0.,D_eta_trk=0.,DR_trk=0.,DpT_trk=0.;
   int match=0,max_size=0;
   double tmp_pTfraction=0.,tmp_DR=0.,tmp_min_pTfraction=1.,tmp_min_DR=1.,m_pTfraction_cut=1.,m_DRcut=0.1,m_pTfraction_nocut=1e6,m_DRnocut=1e6;
   unsigned size_jet=0,size_child=0;
   int den=0,m_den=0;

   unsigned m_track_cut=10;

   float pt_max=500., pt_min=0.;
   float Delta=(pt_max-pt_min)/bin_1;

   std::vector< std::vector<float> > bin_v = std::vector< std::vector<float> >(bin_1);

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
   TH1F *hist_dl1rnn_pb;
   TH1F *hist_dl1rnn_pc;
   TH1F *hist_dl1rnn_pu;

   TH1F *hist_ip2d_llr_l;
   TH1F *hist_ip3d_llr_l;
   TH1F *hist_ip2d_llr_inB;
   TH1F *hist_ip3d_llr_inB;
   TH1F *hist_ip2d_llr_inC;
   TH1F *hist_ip3d_llr_inC;
   TH1F *hist_ip2d_llr_exB;
   TH1F *hist_ip3d_llr_exB;
   TH1F *hist_ip2d_llr_exC;
   TH1F *hist_ip3d_llr_exC;

   TH1F *hist_dl1_l;
   TH1F *hist_dl1_inC;
   TH1F *hist_dl1_inB;
   TH1F *hist_dl1_exC;
   TH1F *hist_dl1_exB;

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
   TH1F *hist_trk_origin_inB;

   TH1F *hist_child_pT_inB;
   TH1F *hist_child_Deta_inB;
   TH1F *hist_child_eta_inB;
   TH1F *hist_child_Dphi_inB;
   TH1F *hist_child_phi_inB;
   TH2F *hist_child_Dphi_Deta_inB;
   TH1F *hist_child_DR_inB;
   TH2F *hist_child_pT_DR_inB;
   TH2F *hist_child_pT_jet_DR_inB;
   TH1F *hist_child_Lxyz_inB;
   TH1F *hist_child_decay_IP;
   TH1F *hist_child_nodecay_IP;
   TH1F *hist_pT_vs_R0_ratio_inB;
//   TH1F *hist_child_linearIP;

   TH1F *hist_efficiency_inB;
   TH1D *hist_n_child;
   TH1D *hist_n_trk;
   TH1D *hist_n_match;

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



};

#endif

#ifdef selector_1_cxx
void selector_1::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the reader is initialized.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   fReader.SetTree(tree);
}

Bool_t selector_1::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}


#endif // #ifdef selector_1_cxx
