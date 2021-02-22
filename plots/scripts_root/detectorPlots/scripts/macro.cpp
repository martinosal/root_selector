#include "plots.cpp"
#include "plots_2.cpp"
#define N_alg 3
void macro(TFile* fData=_file0){

 std::vector<std::string> alg{"JF","IP3D","SV1"};
 std::vector<std::string> det{"d0sig","z0sinthsig","logpTfrac","logDR","IBLhits","NextToIBLhits","sharedIBLhits","splitIBLhits","nPixhits","sharedPixhits","splitPixhits","nSCThits","sharedSCThits"};
 std::vector<int> n_p{0,0,0,0,3,4,3,3,10,6,7,19,7};

 for(unsigned i=0;i<13;i++)//was 13
 {
   plot((det.at(i)).c_str(),n_p.at(i),fData);
   for(unsigned j=0;j<N_alg;j++)
   {
     plot((det.at(i)).c_str(),alg.at(j),n_p.at(i),fData);
   }
 }


 sig("sigd0",fData);
 sig("sigz0",fData);

 for(unsigned j=0;j<N_alg;j++)
 {
   orig(alg.at(j),fData);
 }
 return;
}
