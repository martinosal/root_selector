void lecce_files(TChain *f)
{

/*
    f->Add(
	   "/afs/le.infn.it/project/itk/FTag/flav_Akt4EMPf_BTagging201903.root"
	   );
*/

  std::string fileLocation = "/nfs/kloe/einstein3/stefania/Ftag/processedDAOD_FTAG1/mc16_13TeV.410470.PhPy8EG_A14_ttbar_hdamp258p75_nonallhad.deriv.DAOD_FTAG1.e6337_e5984_s3126_r10201_r10210_p4062/noRetag_v1";



  /// EMPFlow 

  
  std::string fileName = "flav_Akt4EMPf_BTagging201903.root"  ;
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
  
  /*
  
  ///  flav_AktVR30Rmax4Rmin02Tr_BTagging201903
  
  std::string fileName = "flav_AktVR30Rmax4Rmin02Tr_BTagging201903.root"  ;
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
  
  
  ////   flav_AktVR30Rmax4Rmin02TrackGhostTagJets.root
  
  std::string fileName = "flav_AktVR30Rmax4Rmin02TrackGhostTagJets.root";
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

 
  /////       flav_AktVR30Rmax4Rmin02Tr.root


  std::string fileName = "flav_AktVR30Rmax4Rmin02Tr.root";
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

  */
  
}
