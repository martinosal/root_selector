void lecce_files(TChain *f, std::string jcoll)
{


  std::string fileLocation = "/nfs/kloe/einstein3/stefania/Ftag/processedDAOD_FTAG1/mc16_13TeV.410470.PhPy8EG_A14_ttbar_hdamp258p75_nonallhad.deriv.DAOD_FTAG1.e6337_e5984_s3126_r10201_r10210_p4062/noRetag_v1";
  std::cout<<"laoding files from "<< fileLocation <<std::endl;
  std::cout<<"and for tree named "<< jcoll <<std::endl;

  std::string fileName;
  if      (jcoll=="bTag_AntiKt4EMPFlowJets_BTagging201903")             fileName="flav_Akt4EMPf_BTagging201903.root";
  else if (jcoll=="bTag_AntiKtVR30Rmax4Rmin02TrackJets_BTagging201903") fileName="flav_AktVR30Rmax4Rmin02Tr_BTagging201903.root";
  else if (jcoll=="bTag_AntiKtVR30Rmax4Rmin02TrackJets")                fileName="flav_AktVR30Rmax4Rmin02Tr.root";
  else if (jcoll=="bTag_AntiKtVR30Rmax4Rmin02TrackGhostTagJets")        fileName="flav_AktVR30Rmax4Rmin02TrackGhostTagJets.root";
  else
    {
      std::cout<<"Unknown jetCollection "<<  jcoll <<" ... unable to build the chain; exit"<<std::endl;
      return;
    }
  std::cout<<"files to be read ... "<< fileName <<std::endl;
 

  f->Add((fileLocation+ + "/DAOD_FTAG1.20480638._000001.pool.root.1/" + fileName).c_str());    
  f->Add((fileLocation+ + "/DAOD_FTAG1.20480638._000002.pool.root.1/" + fileName).c_str());  
  f->Add((fileLocation+ + "/DAOD_FTAG1.20480638._000003.pool.root.1/" + fileName).c_str());  
  f->Add((fileLocation+ + "/DAOD_FTAG1.20480638._000004.pool.root.1/" + fileName).c_str());  
  f->Add((fileLocation+ + "/DAOD_FTAG1.20480638._000005.pool.root.1/" + fileName).c_str());  
  f->Add((fileLocation+ + "/DAOD_FTAG1.20480638._000006.pool.root.1/" + fileName).c_str());  
  f->Add((fileLocation+ + "/DAOD_FTAG1.20480638._000007.pool.root.1/" + fileName).c_str());  
  f->Add((fileLocation+ + "/DAOD_FTAG1.20480638._000008.pool.root.1/" + fileName).c_str());  
  f->Add((fileLocation+ + "/DAOD_FTAG1.20480638._000009.pool.root.1/" + fileName).c_str());  
  f->Add((fileLocation+ + "/DAOD_FTAG1.20480638._000010.pool.root.1/" + fileName).c_str());  
  f->Add((fileLocation+ + "/DAOD_FTAG1.20480638._000011.pool.root.1/" + fileName).c_str());  
  f->Add((fileLocation+ + "/DAOD_FTAG1.20480638._000012.pool.root.1/" + fileName).c_str());  
  f->Add((fileLocation+ + "/DAOD_FTAG1.20480638._000013.pool.root.1/" + fileName).c_str());  
  f->Add((fileLocation+ + "/DAOD_FTAG1.20480638._000014.pool.root.1/" + fileName).c_str());  
  
  
  
}
