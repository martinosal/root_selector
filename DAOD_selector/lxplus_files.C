#include <sstream>
#include <iomanip>

void lxplus_files(TChain *f, std::string jcoll)
{
  std::string fileLocation = "/eos/user/m/mcentonz/File/FTPF/FTPF_taggers_output/QT/";

  std::string fileName;
  if      (jcoll=="bTag_AntiKt4EMPFlowJets_BTagging201903")             fileName="flav_Akt4EMPf_BTagging201903_derorigin_";
  else if (jcoll=="bTag_AntiKtVR30Rmax4Rmin02TrackJets_BTagging201903") fileName="flav_AktVR30Rmax4Rmin02Tr_BTagging201903_derorigin_";
  else if (jcoll=="bTag_AntiKtVR30Rmax4Rmin02TrackGhostTagJets")        fileName="flav_AktVR30Rmax4Rmin02TrackGhostTagJets_derorigin_";
  else
    {
      std::cout<<"Unknown jetCollection "<<  jcoll <<" ... unable to build the chain; exit"<<std::endl;
      return;
    }

  for(int i=1;i<=14;i++){
    std::stringstream ss;
    ss << std::setw(3) << std::setfill('0') << i;
    std::string s = ss.str();

    f->Add((fileLocation+fileName+s+".root").c_str());
  }
  /*
  f->Add("/eos/user/m/mcentonz/File/FTPF/FTPF_taggers_output/QT/flav_Akt4EMPf_BTagging201903_derorigin_001.root");
  f->Add("/eos/user/m/mcentonz/File/FTPF/FTPF_taggers_output/QT/flav_Akt4EMPf_BTagging201903_derorigin_002.root");
  f->Add("/eos/user/m/mcentonz/File/FTPF/FTPF_taggers_output/QT/flav_Akt4EMPf_BTagging201903_derorigin_003.root");
  f->Add("/eos/user/m/mcentonz/File/FTPF/FTPF_taggers_output/QT/flav_Akt4EMPf_BTagging201903_derorigin_004.root");
  f->Add("/eos/user/m/mcentonz/File/FTPF/FTPF_taggers_output/QT/flav_Akt4EMPf_BTagging201903_derorigin_005.root");
  f->Add("/eos/user/m/mcentonz/File/FTPF/FTPF_taggers_output/QT/flav_Akt4EMPf_BTagging201903_derorigin_006.root");
  f->Add("/eos/user/m/mcentonz/File/FTPF/FTPF_taggers_output/QT/flav_Akt4EMPf_BTagging201903_derorigin_007.root");
  f->Add("/eos/user/m/mcentonz/File/FTPF/FTPF_taggers_output/QT/flav_Akt4EMPf_BTagging201903_derorigin_008.root");
  f->Add("/eos/user/m/mcentonz/File/FTPF/FTPF_taggers_output/QT/flav_Akt4EMPf_BTagging201903_derorigin_009.root");
  f->Add("/eos/user/m/mcentonz/File/FTPF/FTPF_taggers_output/QT/flav_Akt4EMPf_BTagging201903_derorigin_010.root");
  f->Add("/eos/user/m/mcentonz/File/FTPF/FTPF_taggers_output/QT/flav_Akt4EMPf_BTagging201903_derorigin_011.root");
  f->Add("/eos/user/m/mcentonz/File/FTPF/FTPF_taggers_output/QT/flav_Akt4EMPf_BTagging201903_derorigin_012.root");
  f->Add("/eos/user/m/mcentonz/File/FTPF/FTPF_taggers_output/QT/flav_Akt4EMPf_BTagging201903_derorigin_013.root");
  f->Add("/eos/user/m/mcentonz/File/FTPF/FTPF_taggers_output/QT/flav_Akt4EMPf_BTagging201903_derorigin_014.root");

  f->Add("/eos/user/m/mcentonz/File/FTPF/FTPF_taggers_output/QT/flav_AktVR30Rmax4Rmin02TrackGhostTagJets_derorigin_001.root");
  f->Add("/eos/user/m/mcentonz/File/FTPF/FTPF_taggers_output/QT/flav_AktVR30Rmax4Rmin02TrackGhostTagJets_derorigin_002.root");
  f->Add("/eos/user/m/mcentonz/File/FTPF/FTPF_taggers_output/QT/flav_AktVR30Rmax4Rmin02TrackGhostTagJets_derorigin_003.root");
  f->Add("/eos/user/m/mcentonz/File/FTPF/FTPF_taggers_output/QT/flav_AktVR30Rmax4Rmin02TrackGhostTagJets_derorigin_004.root");
  f->Add("/eos/user/m/mcentonz/File/FTPF/FTPF_taggers_output/QT/flav_AktVR30Rmax4Rmin02TrackGhostTagJets_derorigin_005.root");
  f->Add("/eos/user/m/mcentonz/File/FTPF/FTPF_taggers_output/QT/flav_AktVR30Rmax4Rmin02TrackGhostTagJets_derorigin_006.root");
  f->Add("/eos/user/m/mcentonz/File/FTPF/FTPF_taggers_output/QT/flav_AktVR30Rmax4Rmin02TrackGhostTagJets_derorigin_007.root");
  f->Add("/eos/user/m/mcentonz/File/FTPF/FTPF_taggers_output/QT/flav_AktVR30Rmax4Rmin02TrackGhostTagJets_derorigin_008.root");
  f->Add("/eos/user/m/mcentonz/File/FTPF/FTPF_taggers_output/QT/flav_AktVR30Rmax4Rmin02TrackGhostTagJets_derorigin_009.root");
  f->Add("/eos/user/m/mcentonz/File/FTPF/FTPF_taggers_output/QT/flav_AktVR30Rmax4Rmin02TrackGhostTagJets_derorigin_010.root");
  f->Add("/eos/user/m/mcentonz/File/FTPF/FTPF_taggers_output/QT/flav_AktVR30Rmax4Rmin02TrackGhostTagJets_derorigin_011.root");
  f->Add("/eos/user/m/mcentonz/File/FTPF/FTPF_taggers_output/QT/flav_AktVR30Rmax4Rmin02TrackGhostTagJets_derorigin_012.root");
  f->Add("/eos/user/m/mcentonz/File/FTPF/FTPF_taggers_output/QT/flav_AktVR30Rmax4Rmin02TrackGhostTagJets_derorigin_013.root");
  f->Add("/eos/user/m/mcentonz/File/FTPF/FTPF_taggers_output/QT/flav_AktVR30Rmax4Rmin02TrackGhostTagJets_derorigin_014.root");

  f->Add("/eos/user/m/mcentonz/File/FTPF/FTPF_taggers_output/QT/flav_AktVR30Rmax4Rmin02Tr_BTagging201903_derorigin_001.root");
  f->Add("/eos/user/m/mcentonz/File/FTPF/FTPF_taggers_output/QT/flav_AktVR30Rmax4Rmin02Tr_BTagging201903_derorigin_002.root");
  f->Add("/eos/user/m/mcentonz/File/FTPF/FTPF_taggers_output/QT/flav_AktVR30Rmax4Rmin02Tr_BTagging201903_derorigin_003.root");
  f->Add("/eos/user/m/mcentonz/File/FTPF/FTPF_taggers_output/QT/flav_AktVR30Rmax4Rmin02Tr_BTagging201903_derorigin_004.root");
  f->Add("/eos/user/m/mcentonz/File/FTPF/FTPF_taggers_output/QT/flav_AktVR30Rmax4Rmin02Tr_BTagging201903_derorigin_005.root");
  f->Add("/eos/user/m/mcentonz/File/FTPF/FTPF_taggers_output/QT/flav_AktVR30Rmax4Rmin02Tr_BTagging201903_derorigin_006.root");
  f->Add("/eos/user/m/mcentonz/File/FTPF/FTPF_taggers_output/QT/flav_AktVR30Rmax4Rmin02Tr_BTagging201903_derorigin_007.root");
  f->Add("/eos/user/m/mcentonz/File/FTPF/FTPF_taggers_output/QT/flav_AktVR30Rmax4Rmin02Tr_BTagging201903_derorigin_008.root");
  f->Add("/eos/user/m/mcentonz/File/FTPF/FTPF_taggers_output/QT/flav_AktVR30Rmax4Rmin02Tr_BTagging201903_derorigin_009.root");
  f->Add("/eos/user/m/mcentonz/File/FTPF/FTPF_taggers_output/QT/flav_AktVR30Rmax4Rmin02Tr_BTagging201903_derorigin_010.root");
  f->Add("/eos/user/m/mcentonz/File/FTPF/FTPF_taggers_output/QT/flav_AktVR30Rmax4Rmin02Tr_BTagging201903_derorigin_011.root");
  f->Add("/eos/user/m/mcentonz/File/FTPF/FTPF_taggers_output/QT/flav_AktVR30Rmax4Rmin02Tr_BTagging201903_derorigin_012.root");
  f->Add("/eos/user/m/mcentonz/File/FTPF/FTPF_taggers_output/QT/flav_AktVR30Rmax4Rmin02Tr_BTagging201903_derorigin_013.root");
  f->Add("/eos/user/m/mcentonz/File/FTPF/FTPF_taggers_output/QT/flav_AktVR30Rmax4Rmin02Tr_BTagging201903_derorigin_014.root");
  */

}
