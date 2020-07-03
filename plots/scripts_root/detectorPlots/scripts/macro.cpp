#include "plots.cpp"
#define N_alg 4
void macro(TFile* fData=_file0){

 std::vector<std::string> alg{"JF","IP2D","SV0","SV1"};


 IBL_plot("trk_IBLhits_origin",3,fData);
 nextToIBL_plot("trk_NextToIBLhits_origin",4,fData);
 sharedIBL_plot("trk_sharedIBLhits_origin", 3,fData);
 splitIBLhits("trk_splitIBLhits_origin", 3,fData);
 nPixhits("trk_nPixhits_origin", 10,fData);
 sharedPixhits("trk_sharedPixhits_origin", 6,fData);
 splitPixhits("trk_splitPixhits_origin", 7,fData);
 nSCThits("trk_nSCThits_origin", 19,fData);
 sharedSCThits("trk_sharedSCThits_origin", 7,fData);
 for(unsigned i=0;i<N_alg;i++)
 {
   IBL_plot("trk_IBLhits_origin",alg.at(i),3,fData);
   nextToIBL_plot("trk_NextToIBLhits_origin",alg.at(i),4,fData);
   sharedIBL_plot("trk_sharedIBLhits_origin",alg.at(i), 3,fData);
   splitIBLhits("trk_splitIBLhits_origin",alg.at(i), 3,fData);
   nPixhits("trk_nPixhits_origin",alg.at(i), 10,fData);
   sharedPixhits("trk_sharedPixhits_origin",alg.at(i), 6,fData);
   splitPixhits("trk_splitPixhits_origin",alg.at(i), 7,fData);
   nSCThits("trk_nSCThits_origin",alg.at(i), 19,fData);
   sharedSCThits("trk_sharedSCThits_origin",alg.at(i), 7,fData);
 }

 return;
}


/*
trk_IBLhits_origin_inB 3
trk_NextToIBLhits_origin_inB 4
trk_sharedIBLhits_origin_inB 3
trk_splitIBLhits_origin_inB 3
trk_nPixhits_origin_inB 10
trk_sharedPixhits_origin_inB 6
trk_splitPixhits_origin_inB 7
trk_nSCThits_origin_inB 19
trk_sharedSCThits_origin_inB 7
*/
