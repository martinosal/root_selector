#!/bin/bash

declare -a variable=("jetEta" "jetPt" "nTrkAlgo_IP2D" "nTrkAlgo_IP3D" "nTrkAlgo_RNNIP" "nTrkAlgo_JF" "nTrkAlgo_SV1" "JFdR" "JFeFc" "JFmVtx" "JFnVtx" "JFsig3d" "SV1dR" "SV1eFc" "SV1mVtx" "SV1nVtx" "SV1sig3d" "SV1dR" "bHPt" "bHPtFraction" "bHjetDR")
#declare -a variable=("jetPt")
declare -a bin=("1B" "1D0B" "0B0D" "2B" "2D0B" "2B2D")

for i in ${variable[@]};
do 
    for j in ${bin[@]};
    do
	root -l -b <<EOF
.x show_variable_inLabelMatrix.C("EMPFlow","$i","$j") 
.q
EOF
    done
done

