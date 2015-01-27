#!/bin/bash
echo Now copying Ntuple files from CASTOR to WORKDIR
CASTOR=/lustre/cms/store/user/selgamma/ZprimeToMuMu_M-5000_Tune4C_13TeV-pythia8/CMSSW720_13TeV_MC_ZprimetoMuMu_M-5000-GenReco-aod-tree12347/ec860fa037c92d70898b3ab19a8aee34
WORKDIR=/cmshome/selgammal/CMSSW_7_2_0/src/MyCodeArea/Analyzer/test/ZprimeToMuMu5000AOD
echo From: ${CASTOR}
echo To: ${WORKDIR}

rfdir ${CASTOR}

FILES=`rfdir ${CASTOR} | grep ZprimetoMuMu-MC-Mass5000-CMSSW720 | awk '{ print $9}'`
                             
mkdir -p ${WORKDIR}

echo $FILES

for file in $FILES 
do

echo Copying ${CASTOR}/$file
cp ${CASTOR}/$file ${WORKDIR}

done #for


echo Done.















