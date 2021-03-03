#!/bin/bash
#
#################################### EMPflow
#
sed 's/Jcoll/bTag_AntiKt4EMPFlowJets_BTagging201903/'    lets_go.C > go1.C 
sed 's/Label/Cone/'                                      go1.C > letsgo.C 
root -l letsgo.C >& bTag_AntiKt4EMPFlowJets_BTagging201903_Cone.log 
#
#
#
sed 's/Jcoll/bTag_AntiKt4EMPFlowJets_BTagging201903/'    lets_go.C > go1.C
sed 's/Label/ConeIncl/'                                  go1.C > letsgo.C 
root -l letsgo.C >& bTag_AntiKt4EMPFlowJets_BTagging201903_ConeIncl.log 
#
#
#
sed 's/Jcoll/bTag_AntiKt4EMPFlowJets_BTagging201903/'    lets_go.C > go1.C
sed 's/Label/GhostCone/'                                 go1.C > letsgo.C 
root -l letsgo.C >& bTag_AntiKt4EMPFlowJets_BTagging201903_GhostCone.log 
#
#
#
sed 's/Jcoll/bTag_AntiKt4EMPFlowJets_BTagging201903/'    lets_go.C > go1.C
sed 's/Label/GhostConeIncl/'                             go1.C > letsgo.C 
root -l letsgo.C >& bTag_AntiKt4EMPFlowJets_BTagging201903_GhostConeIncl.log 
#
#
#
#################################### VR30
#
sed 's/Jcoll/bTag_AntiKtVR30Rmax4Rmin02TrackJets_BTagging201903/'    lets_go.C > go1.C 
sed 's/Label/Cone/'                                      go1.C > letsgo.C 
root -l letsgo.C >& bTag_AntiKtVR30Rmax4Rmin02TrackJets_BTagging201903_Cone.log 
#
#
#
sed 's/Jcoll/bTag_AntiKtVR30Rmax4Rmin02TrackJets_BTagging201903/'    lets_go.C > go1.C
sed 's/Label/ConeIncl/'                                  go1.C > letsgo.C 
root -l letsgo.C >& bTag_AntiKtVR30Rmax4Rmin02TrackJets_BTagging201903_ConeIncl.log 
#
#
#
sed 's/Jcoll/bTag_AntiKtVR30Rmax4Rmin02TrackJets_BTagging201903/'    lets_go.C > go1.C
sed 's/Label/GhostCone/'                                 go1.C > letsgo.C 
root -l letsgo.C >& bTag_AntiKtVR30Rmax4Rmin02TrackJets_BTagging201903_GhostCone.log 
#
#
#
sed 's/Jcoll/bTag_AntiKtVR30Rmax4Rmin02TrackJets_BTagging201903/'    lets_go.C > go1.C
sed 's/Label/GhostConeIncl/'                             go1.C > letsgo.C 
root -l letsgo.C >& bTag_AntiKtVR30Rmax4Rmin02TrackJets_BTagging201903_GhostConeIncl.log 
#
#
#
#################################### VR30Ghost
#
sed 's/Jcoll/bTag_AntiKtVR30Rmax4Rmin02TrackGhostTagJets/'    lets_go.C > go1.C 
sed 's/Label/Ghost/'                                      go1.C > letsgo.C 
root -l letsgo.C >& bTag_AntiKtVR30Rmax4Rmin02TrackGhostTagJets_Cone.log 
#
#
#
sed 's/Jcoll/bTag_AntiKtVR30Rmax4Rmin02TrackGhostTagJets/'    lets_go.C > go1.C
sed 's/Label/GhostIncl/'                                  go1.C > letsgo.C 
root -l letsgo.C >& bTag_AntiKtVR30Rmax4Rmin02TrackGhostTagJets_ConeIncl.log 
#
#
#
sed 's/Jcoll/bTag_AntiKtVR30Rmax4Rmin02TrackGhostTagJets/'    lets_go.C > go1.C
sed 's/Label/GhostCone/'                                 go1.C > letsgo.C 
root -l letsgo.C >& bTag_AntiKtVR30Rmax4Rmin02TrackGhostTagJets_GhostCone.log 
#
#
#
sed 's/Jcoll/bTag_AntiKtVR30Rmax4Rmin02TrackGhostTagJets/'    lets_go.C > go1.C
sed 's/Label/GhostConeIncl/'                             go1.C > letsgo.C 
root -l letsgo.C >& bTag_AntiKtVR30Rmax4Rmin02TrackGhostTagJets_GhostConeIncl.log 
#
#
#
#################################### VR30 12 GeV 
#
sed 's/Jcoll/bTag_AntiKtVR30Rmax4Rmin02TrackJets_BTagging201903/'    lets_go.C > go1.C 
sed 's/Label/Cone/'                                                  go1.C > go2.C 
sed 's/false/true/'                                                  go2.C > letsgo.C 
root -l letsgo.C >& bTag_AntiKtVR30Rmax4Rmin02TrackJets_BTagging201903_Cone_12GeV.log 
#
#
#
sed 's/Jcoll/bTag_AntiKtVR30Rmax4Rmin02TrackJets_BTagging201903/'    lets_go.C > go1.C
sed 's/Label/ConeIncl/'                                              go1.C > go2.C 
sed 's/false/true/'                                                  go2.C > letsgo.C 
root -l letsgo.C >& bTag_AntiKtVR30Rmax4Rmin02TrackJets_BTagging201903_ConeIncl_12GeV.log 
#
#
#
sed 's/Jcoll/bTag_AntiKtVR30Rmax4Rmin02TrackJets_BTagging201903/'    lets_go.C > go1.C
sed 's/Label/GhostCone/'                                             go1.C > go2.C 
sed 's/false/true/'                                                  go2.C > letsgo.C 
root -l letsgo.C >& bTag_AntiKtVR30Rmax4Rmin02TrackJets_BTagging201903_GhostCone_12GeV.log 
#
#
#
sed 's/Jcoll/bTag_AntiKtVR30Rmax4Rmin02TrackJets_BTagging201903/'    lets_go.C > go1.C
sed 's/Label/GhostConeIncl/'                                         go1.C > go2.C 
sed 's/false/true/'                                                  go2.C > letsgo.C 
root -l letsgo.C >& bTag_AntiKtVR30Rmax4Rmin02TrackJets_BTagging201903_GhostConeIncl_12GeV.log 
#
#
#
#################################### VR30Ghost 12 GeV 
#
sed 's/Jcoll/bTag_AntiKtVR30Rmax4Rmin02TrackGhostTagJets/'    lets_go.C > go1.C 
sed 's/Label/Ghost/'                                          go1.C > go2.C 	
sed 's/false/true/'                                           go2.C > letsgo.C  
root -l letsgo.C >& bTag_AntiKtVR30Rmax4Rmin02TrackGhostTagJets_Cone_12GeV.log 
#
#
#
sed 's/Jcoll/bTag_AntiKtVR30Rmax4Rmin02TrackGhostTagJets/'    lets_go.C > go1.C 
sed 's/Label/GhostIncl/'                                      go1.C > go2.C 	
sed 's/false/true/'                                           go2.C > letsgo.C  
root -l letsgo.C >& bTag_AntiKtVR30Rmax4Rmin02TrackGhostTagJets_ConeIncl_12GeV.log 
#
#
#
sed 's/Jcoll/bTag_AntiKtVR30Rmax4Rmin02TrackGhostTagJets/'    lets_go.C > go1.C 
sed 's/Label/GhostCone/'                                      go1.C > go2.C 	
sed 's/false/true/'                                           go2.C > letsgo.C  
root -l letsgo.C >& bTag_AntiKtVR30Rmax4Rmin02TrackGhostTagJets_GhostCone_12GeV.log 
#
#
#
sed 's/Jcoll/bTag_AntiKtVR30Rmax4Rmin02TrackGhostTagJets/'    lets_go.C > go1.C 
sed 's/Label/GhostConeIncl/'                                  go1.C > go2.C 	
sed 's/false/true/'                                           go2.C > letsgo.C  
root -l letsgo.C >& bTag_AntiKtVR30Rmax4Rmin02TrackGhostTagJets_GhostConeIncl_12GeV.log 
#
#
