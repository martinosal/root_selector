# root_selector

To start:

git clone https://github.com/martinosal/root_selector.git

cd root_selector/DAOD_selector/

Input files are in .root extension, MUST be produced by the FtagPerformanceFramework tool under my branch mcentonz_lxplus_norigin in https://gitlab.cern.ch/mcentonz/FlavourTagPerformanceFramework.git. 

So if you need to start from producing these FTagPF files from DAODs, create a new folder and install it

git clone https://gitlab.cern.ch/mcentonz/FlavourTagPerformanceFramework.git

Some preprocessed FtagPF files are available under lxplus:/eos/user/m/mcentonz/File/FTPF/FTPF_taggers_output/, in particular

flav_AktVR30Rmax4Rmin02TrackGhostTagJets.root

which is a sample of 10k events produced from a ttbar_nonhallhadr DAOD and run on a VR30 ghost associated-jet collection. Let me know if there are any questions.

To run on this ntuple, if you are on lxplus just open the lxplus_files.C and add the following line

f->Add("/eos/user/m/mcentonz/File/FTPF/FTPF_taggers_output/flav_AktVR30Rmax4Rmin02TrackGhostTagJets.root");

This file should be available for external users, let me know if there is any problem on that.

Then open launch_selector.C and edit it from 

  bool laptop=true;
  bool lxplus=false;
  bool lecce=false;

to

  bool laptop=false;
  bool lxplus=true;
  bool lecce=false;

Ensure that the jetcollection variable in the same file is selecting the correct TrackGhostTagJets jet collection:

const char *jetcollection="bTag_AntiKtVR30Rmax4Rmin02TrackGhostTagJets";

Now you can finally run the selector with

root -l

.x go.C

This produces your root ntuple named debug_.root.
