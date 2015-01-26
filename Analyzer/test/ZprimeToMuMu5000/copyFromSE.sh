#!/bin/bash
echo Now copying Ntuple files from CASTOR to WORKDIR
CASTOR=/lustre/cms/store/user/selgamma/ZprimeToMuMu_M-5000_Tune4C_13TeV-pythia8/CMSSW720_13TeV_MC_ZprimetoMuMu_M-5000-GenReco-notrigger-pattuple50/b694b2ae17b240271d4de784b97e1a8f
WORKDIR=/cmshome/selgammal/CMSSW_7_2_0/src/MyCodeArea/Analyzer/test/ZprimeToMuMu5000
echo From: ${CASTOR}
echo To: ${WORKDIR}

rfdir ${CASTOR}

FILES=`rfdir ${CASTOR} | grep CMSSW720_MC_ZprimeMuMu5000_13TeV_tree | awk '{ print $9}'`
                             
mkdir -p ${WORKDIR}

echo $FILES

for file in $FILES 
do

echo Copying ${CASTOR}/$file
cp ${CASTOR}/$file ${WORKDIR}

done #for


echo Done.















