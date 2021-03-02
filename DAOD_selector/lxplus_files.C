void lxplus_files(TChain *f, std::string jcoll)
{
  std::string fileLocation = "/eos/user/m/mcentonz/File/FTPF/FTPF_taggers_output/QT/";
  std::cout<<"laoding files from "<< fileLocation <<std::endl;
  std::cout<<"and for tree named "<< jcoll <<std::endl;

  std::string fileName;
  if      (jcoll=="bTag_AntiKt4EMPFlowJets_BTagging201903")             fileName="flav_Akt4EMPf_BTagging201903_derorigin";
  else if (jcoll=="bTag_AntiKtVR30Rmax4Rmin02TrackJets_BTagging201903") fileName="flav_AktVR30Rmax4Rmin02Tr_BTagging201903_derorigin";
  //else if (jcoll=="bTag_AntiKtVR30Rmax4Rmin02TrackJets")               fileName="flav_AktVR30Rmax4Rmin02Tr_derorigin";
  else if (jcoll=="bTag_AntiKtVR30Rmax4Rmin02TrackGhostTagJets")        fileName="flav_AktVR30Rmax4Rmin02TrackGhostTagJets_derorigin";
  else
    {
      std::cout<<"Unknown jetCollection "<<  jcoll <<" ... unable to build the chain; exit"<<std::endl;
      return;
    }
  std::cout<<"files to be read ... "<< fileName <<std::endl;
  
  f->Add((fileLocation+fileName+"_001.root").c_str());
  f->Add((fileLocation+fileName+"_002.root").c_str());
  f->Add((fileLocation+fileName+"_003.root").c_str());
  f->Add((fileLocation+fileName+"_004.root").c_str());
  f->Add((fileLocation+fileName+"_005.root").c_str());
  f->Add((fileLocation+fileName+"_006.root").c_str());
  f->Add((fileLocation+fileName+"_007.root").c_str());
  f->Add((fileLocation+fileName+"_008.root").c_str());
  f->Add((fileLocation+fileName+"_009.root").c_str());
  f->Add((fileLocation+fileName+"_010.root").c_str());
  f->Add((fileLocation+fileName+"_011.root").c_str());
  f->Add((fileLocation+fileName+"_012.root").c_str());
  f->Add((fileLocation+fileName+"_013.root").c_str());
  f->Add((fileLocation+fileName+"_014.root").c_str());
}
