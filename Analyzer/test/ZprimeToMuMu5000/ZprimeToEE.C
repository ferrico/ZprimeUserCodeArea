//=========================================================================  
//          Analysis code for Z' boson to Mu Mu analysis                  =  
//    [1] In this code we select the high pt di-muons events              =
//    [2] It also computes Acceptance X Efficiency (Accxeff)              =                
//             Written by Sherif Elgammal                                 =
//                                                                        =
//                    07/08/2014                                          =
//=========================================================================
#define ZprimeToEE_cxx
#include "ZprimeToEE.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TFile.h"
#include <math.h>
#include "TF2.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TLorentzVector.h"
#include <time.h>
#include <iostream>
#include <vector>
#include <algorithm>    
#include "Math/VectorUtil.h"
using namespace std;
#define PI 3.14159265
bool myfunction (int i,int j) { return (i<j); }
bool picklargemass (float lhs,float rhs) { return (lhs > rhs); }
void ZprimeToEE::Loop()
{
  time_t start,end;
  double dif;
  time (&start);
  FILE * pFile;
  pFile = fopen ("myfile.txt","w");
  //values needed for AccepXeff study
  NbGen = 0;
  NbReco= 0;
  int binMass   = 100;
  float minMass = 0.0;
  float maxMass = 8000.0;
  //values needed for tag & probe study
  MassCutMin = 0.0;
  MassCutMax = 8000.0;
  TagProbeEtaCut = 0.9;
  MassResolution = 0.10;
  EtaCut = 2.1;
  RecoGenMatchingDeltaRcut = 0.15;
  RecoHLTMatchingDeltaRcut = 0.20;

  
  int binPt = 18; //30;
  minPt = 200.0;
  maxPt = 1000.0;  //1000.0;
  Etcut = 35.0;

  int binPV = 50;
  float minPV = 0.0;
  float maxPv = 50.0;
  //================================================================================== 
  //                                                                                 =
  //            Start the histograms for the mass of Z                               =  
  //                                                                                 =
  //==================================================================================
  h1_ZprimeGenmass_       = new TH1F("ZprimeGenmass","",binMass,minMass,maxMass);
  h1_ZprimeGenEta1_       = new TH1F("ZprimeGenEta1","",100,-8.0,8.0);
  h1_ZprimeGenEta2_       = new TH1F("ZprimeGenEta2","",100,-8.0,8.0);
  h1_ZprimeGenPt1_        = new TH1F("ZprimeGenPt1","",100,0.0,2000.0);
  h1_ZprimeGenPt2_        = new TH1F("ZprimeGenPt2","",100,0.0,2000.0);
  h1_ZprimeGenEn1_        = new TH1F("ZprimeGenEn1","",100,0.0,2000.0);
  h1_ZprimeGenEn2_        = new TH1F("ZprimeGenEn2","",100,0.0,2000.0);  
  h1_MassRecoGenDif_      = new TH1F("MassRecoGenDif","",100,-0.5,0.5);
  h1_ZprimeRecomass_      = new TH1F("ZprimeRecomass","",binMass,minMass,maxMass);
  h1_ZprimeRecoEta1_      = new TH1F("ZprimeRecoEta1","",100,-3.0,3.0);
  h1_ZprimeRecoEta2_      = new TH1F("ZprimeRecoEta2","",100,-3.0,3.0);
  h1_ZprimeRecoPhi1_      = new TH1F("ZprimeRecoPhi1","",100,-4.0,4.0);
  h1_ZprimeRecoPhi2_      = new TH1F("ZprimeRecoPhi2","",100,-4.0,4.0);
  h1_LeadingRecoMuPt_     = new TH1F("LeadingRecoMuPt","",50,0.0,1000.0);
  h1_SecondLeadingMuPt_   = new TH1F("SecondLeadingMuPt","",50,0.0,1000.0);
  h1_ZprimeRecoEn1_       = new TH1F("ZprimeRecoEn1","",50,0.0,2000.0);
  h1_ZprimeRecoEn2_       = new TH1F("ZprimeRecoEn2","",50,0.0,2000.0);
  h1_MassTest_            = new TH1F("MassTest","",100,0.0,2000.0);
  h1_RecoPhiMu1DifPhiMu2_ = new TH1F("RecoPhiMu1DifPhiMu2","",100,-1.2,1.2);
  h1_GenPhiMu1DifPhiMu2_  = new TH1F("GenPhiMu1DifPhiMu2","",100,-1.2,1.2);
  h1_dPToverPT_                            = new TH1F("dPToverPT","",100,0.0,0.5);
  h1_normalizedChi2_                       = new TH1F("normalizedChi2","",100,0.0,20.0);
  h1_numberOftrackerLayersWithMeasurement_ = new TH1F("numberOftrackerLayersWithMeasurement","",20,0.0,20.0);
  h1_numberOfValidPixelHits_               = new TH1F("numberOfValidPixelHits","",10,0.0,10.0);
  h1_numberOfValidMuonHits_                = new TH1F("numberOfValidMuonHits","",60,0.0,60.0);
  h1_numberOfMatchedStations_              = new TH1F("numberOfMatchedStations","",10,0.0,10.0);
  h1_trackiso_                             = new TH1F("trackiso","",50,0.0,0.3);
  h1_absdxy_                               = new TH1F("absdxy","",100,0.0,0.3);
  h1_absdz_                                = new TH1F("absdz","",100,0.0,1.0);
  h1_deltaR_                               = new TH1F("deltaR","",100,-0.01,0.01);
  h1_PtEffpterror_ = new TH1F("PtEffpterror","",binPt,minPt,maxPt);
  h1_PtEffptnumberOftrackerLayers_ = new TH1F("PtEffptnumberOftrackerLayers","",binPt,minPt,maxPt);
  h1_PtEffptnumberOfPixelHits_       = new TH1F("PtEffptnumberOfPixelHits","",binPt,minPt,maxPt);
  h1_PtEffptnumberOfMuonHits_        = new TH1F("PtEffptnumberOfMuonHits","",binPt,minPt,maxPt);
  h1_PtEffptnumberOfMatchedStations_ = new TH1F("PtEffptnumberOfMatchedStations","",binPt,minPt,maxPt);
  h1_PtEffptTrackIso_                = new TH1F("PtEffptTrackIso","",binPt,minPt,maxPt);
  h1_PtEffptabsdsy_                  = new TH1F("PtEffptabsdsy","",binPt,minPt,maxPt);
  h1_PtEffpfSumChargedHadron_ = new TH1F("PtEffpfSumChargedHadron","",binPt,minPt,maxPt);
  h1_PtEffpfSumNeutralHadron_ = new TH1F("PtEffpfSumNeutralHadron","",binPt,minPt,maxPt);
  h1_PtEffpfPhotonIso_ = new TH1F("PtEffpfPhotonIso","",binPt,minPt,maxPt);
  h1_FracpfSumChargedHadron_ = new TH1F("FracpfSumChargedHadron","",100,0.0,10.0);
  h1_FracpfSumNeutralHadron_ = new TH1F("FracpfSumNeutralHadron","",100,0.0,10.0);
  h1_FracpfPhotonIso_        = new TH1F("FracpfPhotonIso","",100,0.0,10.0);
  h1_EtaEffpterror_                  = new TH1F("EtaEffpterror","",30,0.0,3.0);
  h1_EtaEffptnumberOftrackerLayers_  = new TH1F("EtaEffptnumberOftrackerLayers","",30,0.0,3.0);
  h1_EtaEffptnumberOfPixelHits_      = new TH1F("EtaEffptnumberOfPixelHits","",30,0.0,3.0);
  h1_EtaEffptnumberOfMuonHits_       = new TH1F("EtaEffptnumberOfMuonHits","",30,0.0,3.0);
  h1_EtaEffptnumberOfMatchedStations_= new TH1F("EtaEffptnumberOfMatchedStations","",30,0.0,3.0);
  h1_EtaEffptTrackIso_               = new TH1F("EtaEffptTrackIso","",30,0.0,3.0);
  h1_EtaEffptabsdsy_                 = new TH1F("EtaEffptabsdsy","",30,0.0,3.0);
  h1_nbPVEffptTrackIso_ = new TH1F("nbPVEffptTrackIso","",binPV,minPV,maxPv);
  
  h2_RecoMass1vsMass2_ = new TH2F("RecoMass1vsMass2","",binMass,minMass,maxMass,binMass,minMass,maxMass);
  h1_test_= new TH1F("test","",50,0.0,3000.0);
  h1_nbPVID_               = new TH1F("nbPVID","",binPV,minPV,maxPv);
  h1_PtID_                 = new TH1F("PtID","",binPt,minPt,maxPt);
  h1_EtaID_                = new TH1F("EtaID","",30,0.0,3.0);
  
  h1_PtEffPFdeltaBeta_     = new TH1F("PtEffPFdeltaBeta","",binPt,minPt,maxPt);
  h1_EtaEffPFdeltaBeta_    = new TH1F("EtaEffPFdeltaBeta","",30,0.0,3.0);
  h1_nbPVEffPFdeltaBeta_   = new TH1F("nbPVEffPFdeltaBeta","",binPV,minPV,maxPv);
  h1_pfDeltaBetaID_         = new TH1F("pfDeltaBetaID","",100,0.0,10.0);

  h1_PtHeepEleEB_     = new TH1F("PtHeepEleEB","",binPt,minPt,maxPt);
  h1_EtaHeepEleEB_    = new TH1F("EtaHeepEleEB","",30,0.0,3.0);
  h1_PVHeepEleEB_     = new TH1F("PVHeepEleEB","",binPV,minPV,maxPv);
  h1_PtHeepEleEE_     = new TH1F("PtHeepEleEE","",binPt,minPt,maxPt);
  h1_EtaHeepEleEE_    = new TH1F("EtaHeepEleEE","",30,0.0,3.0);
  h1_PVHeepEleEE_     = new TH1F("PVHeepEleEE","",binPV,minPV,maxPv);
  h1_pfDeltaBetaEB_   = new TH1F("pfDeltaBetaEB","",100,0.0,2.0);
  h1_pfDeltaBetaEE_   = new TH1F("pfDeltaBetaEE","",100,0.0,2.0);

  h1_PtHeepdeltaBetaEB_     = new TH1F("PtHeepdeltaBetaEB","",binPt,minPt,maxPt);
  h1_EtaHeepdeltaBetaEB_    = new TH1F("EtaHeepdeltaBetaEB","",30,0.0,3.0);
  h1_PVHeepdeltaBetaEB_     = new TH1F("PVHeepdeltaBetaEB","",binPV,minPV,maxPv);
  h1_PtHeepdeltaBetaEE_     = new TH1F("PtHeepdeltaBetaEE","",binPt,minPt,maxPt);
  h1_EtaHeepdeltaBetaEE_    = new TH1F("EtaHeepdeltaBetaEE","",30,0.0,3.0);
  h1_PVHeepdeltaBetaEE_     = new TH1F("PVHeepdeltaBetaEE","",binPV,minPV,maxPv);

  h1_PtHeepphotonconvEB_     = new TH1F("PtHeepphotonconvEB","",binPt,minPt,maxPt);
  h1_EtaHeepphotonconvEB_    = new TH1F("EtaHeepphotonconvEB","",30,0.0,3.0);
  h1_PVHeepphotonconvEB_     = new TH1F("PVHeepphotonconvEB","",binPV,minPV,maxPv);
  h1_PtHeepphotonconvEE_     = new TH1F("PtHeepphotonconvEE","",binPt,minPt,maxPt);
  h1_EtaHeepphotonconvEE_    = new TH1F("EtaHeepphotonconvEE","",30,0.0,3.0);
  h1_PVHeepphotonconvEE_     = new TH1F("PVHeepphotonconvEE","",binPV,minPV,maxPv);

  h1_PtHeepHoeNewEB_     = new TH1F("PtHeepHoeNewEB","",binPt,minPt,maxPt);
  h1_EtaHeepHoeNewEB_    = new TH1F("EtaHeepHoeNewEB","",30,0.0,3.0);
  h1_PVHeepHoeNewEB_     = new TH1F("PVHeepHoeNewEB","",binPV,minPV,maxPv);
  h1_PtHeepHoeNewEE_     = new TH1F("PtHeepHoeNewEE","",binPt,minPt,maxPt);
  h1_EtaHeepHoeNewEE_    = new TH1F("EtaHeepHoeNewEE","",30,0.0,3.0);
  h1_PVHeepHoeNewEE_     = new TH1F("PVHeepHoeNewEE","",binPV,minPV,maxPv);


  h1_PtHeepDeltaEtaNewEB_     = new TH1F("PtHeepDeltaEtaNewEB","",binPt,minPt,maxPt);
  h1_EtaHeepDeltaEtaNewEB_    = new TH1F("EtaHeepDeltaEtaNewEB","",30,0.0,3.0);
  h1_PVHeepDeltaEtaNewEB_     = new TH1F("PVHeepDeltaEtaNewEB","",binPV,minPV,maxPv);
  h1_PtHeepDeltaEtaNewEE_     = new TH1F("PtHeepDeltaEtaNewEE","",binPt,minPt,maxPt);
  h1_EtaHeepDeltaEtaNewEE_    = new TH1F("EtaHeepDeltaEtaNewEE","",30,0.0,3.0);
  h1_PVHeepDeltaEtaNewEE_     = new TH1F("PVHeepDeltaEtaNewEE","",binPV,minPV,maxPv);

  h1_pfSumChargedHadronEB_      = new TH1F("pfSumChargedHadronEB","",100,0.0,100.0);
  h1_pfSumNeutralHadronEB_      = new TH1F("pfSumNeutralHadronEB","",100,0.0,100.0);
  h1_pfPhotonIsoEB_             = new TH1F("pfPhotonIsoEB","",100,0.0,100.0);
  h1_pfSumChargedHadronEE_      = new TH1F("pfSumChargedHadronEE","",100,0.0,100.0);
  h1_pfSumNeutralHadronEE_      = new TH1F("pfSumNeutralHadronEE","",100,0.0,100.0);
  h1_pfPhotonIsoEE_             = new TH1F("pfPhotonIsoEE","",100,0.0,100.0);

  h1_PtHeepTrackIsoNewEB_     = new TH1F("PtHeepTrackIsoNewEB","",binPt,minPt,maxPt);
  h1_EtaHeepTrackIsoNewEB_    = new TH1F("EtaHeepTrackIsoNewEB","",30,0.0,3.0);
  h1_PVHeepTrackIsoNewEB_     = new TH1F("PVHeepTrackIsoNewEB","",binPV,minPV,maxPv);
  h1_PtHeepTrackIsoNewEE_     = new TH1F("PtHeepTrackIsoNewEE","",binPt,minPt,maxPt);
  h1_EtaHeepTrackIsoNewEE_    = new TH1F("EtaHeepTrackIsoNewEE","",30,0.0,3.0);
  h1_PVHeepTrackIsoNewEE_     = new TH1F("PVHeepTrackIsoNewEE","",binPV,minPV,maxPv);

  h1_PtHeepMissHitsNewEB_    = new TH1F("PtHeepMissHitsNewEB","",binPt,minPt,maxPt);
  h1_EtaHeepMissHitsNewEB_   = new TH1F("EtaHeepMissHitsNewEB","",30,0.0,3.0);
  h1_PVHeepMissHitsNewEB_    = new TH1F("PVHeepMissHitsNewEB","",binPV,minPV,maxPv);
  h1_PtHeepMissHitsNewEE_    = new TH1F("PtHeepMissHitsNewEE","",binPt,minPt,maxPt);
  h1_EtaHeepMissHitsNewEE_   = new TH1F("EtaHeepMissHitsNewEE","",30,0.0,3.0);
  h1_PVHeepMissHitsNewEE_    = new TH1F("PVHeepMissHitsNewEE","",binPV,minPV,maxPv);

  h1_PtHeepDeltaPhiNewEB_     = new TH1F("PtHeepDeltaPhiNewEB","",binPt,minPt,maxPt);
  h1_EtaHeepDeltaPhiNewEB_    = new TH1F("EtaHeepDeltaPhiNewEB","",30,0.0,3.0);
  h1_PVHeepDeltaPhiNewEB_     = new TH1F("PVHeepDeltaPhiNewEB","",binPV,minPV,maxPv);
  h1_PtHeepDeltaPhiNewEE_     = new TH1F("PtHeepDeltaPhiNewEE","",binPt,minPt,maxPt);
  h1_EtaHeepDeltaPhiNewEE_    = new TH1F("EtaHeepDeltaPhiNewEE","",30,0.0,3.0);
  h1_PVHeepDeltaPhiNewEE_     = new TH1F("PVHeepDeltaPhiNewEE","",binPV,minPV,maxPv);

  h1_PtHeepSigmaNewEE_     = new TH1F("PtHeepSigmaNewEE","",binPt,minPt,maxPt);
  h1_EtaHeepSigmaNewEE_    = new TH1F("EtaHeepSigmaNewEE","",30,0.0,3.0);
  h1_PVHeepSigmaNewEE_     = new TH1F("PVHeepSigmaNewEE","",binPV,minPV,maxPv);


  h1_PtHeepEcalHcalIsoNewEB_     = new TH1F("PtHeepEcalHcalIsoNewEB","",binPt,minPt,maxPt);
  h1_EtaHeepEcalHcalIsoNewEB_    = new TH1F("EtaHeepEcalHcalIsoNewEB","",30,0.0,3.0);
  h1_PVHeepEcalHcalIsoNewEB_     = new TH1F("PVHeepEcalHcalIsoNewEB","",binPV,minPV,maxPv);
  h1_PtHeepEcalHcalIsoNewEE_     = new TH1F("PtHeepEcalHcalIsoNewEE","",binPt,minPt,maxPt);
  h1_EtaHeepEcalHcalIsoNewEE_    = new TH1F("EtaHeepEcalHcalIsoNewEE","",30,0.0,3.0);
  h1_PVHeepEcalHcalIsoNewEE_     = new TH1F("PVHeepEcalHcalIsoNewEE","",binPV,minPV,maxPv);


  h1_PtHeepShowerShapeNewEB_     = new TH1F("PtHeepShowerShapeNewEB","",binPt,minPt,maxPt);
  h1_EtaHeepShowerShapeNewEB_    = new TH1F("EtaHeepShowerShapeNewEB","",30,0.0,3.0);
  h1_PVHeepShowerShapeNewEB_     = new TH1F("PVHeepShowerShapeNewEB","",binPV,minPV,maxPv);
  //===================================================================================
  //
  //
  //                       Part for HEEP selection version 4.1
  //
  //
  //===================================================================================
  h1_PtHeepEleV4EB_     = new TH1F("PtHeepEleV4EB","",binPt,minPt,maxPt);
  h1_EtaHeepEleV4EB_    = new TH1F("EtaHeepEleV4EB","",30,0.0,3.0);
  h1_PVHeepEleV4EB_     = new TH1F("PVHeepEleV4EB","",binPV,minPV,maxPv);
  h1_PtHeepEleV4EE_     = new TH1F("PtHeepEleV4EE","",binPt,minPt,maxPt);
  h1_EtaHeepEleV4EE_    = new TH1F("EtaHeepEleV4EE","",30,0.0,3.0);
  h1_PVHeepEleV4EE_     = new TH1F("PVHeepEleV4EE","",binPV,minPV,maxPv);

  h1_PtHeepHoeV4EB_     = new TH1F("PtHeepHoeV4EB","",binPt,minPt,maxPt);
  h1_EtaHeepHoeV4EB_    = new TH1F("EtaHeepHoeV4EB","",30,0.0,3.0);
  h1_PVHeepHoeV4EB_     = new TH1F("PVHeepHoeV4EB","",binPV,minPV,maxPv);
  h1_PtHeepHoeV4EE_     = new TH1F("PtHeepHoeV4EE","",binPt,minPt,maxPt);               
  h1_EtaHeepHoeV4EE_    = new TH1F("EtaHeepHoeV4EE","",30,0.0,3.0);
  h1_PVHeepHoeV4EE_     = new TH1F("PVHeepHoeV4EE","",binPV,minPV,maxPv);

  h1_PtHeepDeltaEtaInV4EB_     = new TH1F("PtHeepDeltaEtaInV4EB","",binPt,minPt,maxPt);
  h1_EtaHeepDeltaEtaInV4EB_    = new TH1F("EtaHeepDeltaEtaInV4EB","",30,0.0,3.0);
  h1_PVHeepDeltaEtaInV4EB_     = new TH1F("PVHeepDeltaEtaInV4EB","",binPV,minPV,maxPv);
  h1_PtHeepDeltaEtaInV4EE_     = new TH1F("PtHeepDeltaEtaInV4EE","",binPt,minPt,maxPt);               
  h1_EtaHeepDeltaEtaInV4EE_    = new TH1F("EtaHeepDeltaEtaInV4EE","",30,0.0,3.0);
  h1_PVHeepDeltaEtaInV4EE_     = new TH1F("PVHeepDeltaEtaInV4EE","",binPV,minPV,maxPv);

  h1_DxyNewEB_          = new TH1F("DxyNewEB","",100,-0.1,0.1);
  h1_DxyNewEE_          = new TH1F("DxyNewEE","",100,-0.1,0.1);
  h1_PtHeepDxyNewEB_    = new TH1F("PtHeepDxyNewEB","",binPt,minPt,maxPt);
  h1_EtaHeepDxyNewEB_   = new TH1F("EtaHeepDxyNewEB","",30,0.0,3.0);
  h1_PVHeepDxyNewEB_    = new TH1F("PVHeepDxyNewEB","",binPV,minPV,maxPv);
  h1_PtHeepDxyNewEE_    = new TH1F("PtHeepDxyNewEE","",binPt,minPt,maxPt);
  h1_EtaHeepDxyNewEE_   = new TH1F("EtaHeepDxyNewEE","",30,0.0,3.0);
  h1_PVHeepDxyNewEE_    = new TH1F("PVHeepDxyNewEE","",binPV,minPV,maxPv);
  //==================================================================================  

  
  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    //    /*
    cout<<"=======> jentry = "<<jentry<< 
      "=======> Evt = "<<event_evtNo<< 
      "=======> Run = "<<event_runNo<< 
      "=======> Lumi = "<<event_lumi<< 
      "=======> bunch = "<<event_bunch<<endl;
    // */
    //=====================================================
    //                                                    =
    //  Calling methods to get the generated information  =
    //                                                    =
    //=====================================================
    bool firstGenMu  = SelectFirstGenMu(genET1,genPhi1,genEta1,genEn1,genID1,genStat1,flag1);
    bool secondGenMu = SelectSecondGenMu(flag1,genET1,genET2,genPhi2,genEta2,genEn2,genID2,genStat2);
    if(firstGenMu == 0 || secondGenMu == 0) continue;
    MassGen = Mass(genET1,genEta1,genPhi1,genEn1,genET2,genEta2,genPhi2,genEn2);
    PlotGenInfo(MassGen,genEta1,genEta2,genET1,genET2,genEn1,genEn2,genPhi1,genPhi2);
    NbGen ++;
    //=========================================================
    //                                                        =
    // Calling methods to get events with 2 muons passing ID  =
    //                                                        =
    //=========================================================
    float EtaTracker1,PhiTracker1;
    float EtaTracker2,PhiTracker2;
    bool firstMuFinal  = SelectFirstMuon(PtRecMu1,EnRecMu1,EtaRecMu1,PhiRecMu1,flagmu1,EtaTracker1,PhiTracker1);
    bool secondMuFinal = SelectSecondMuon(flagmu1,PtRecMu1,PtRecMu2,EnRecMu2,EtaRecMu2,PhiRecMu2,
					  EtaTracker2,PhiTracker2);
    //=========================================================
    //                  compute invariant mass                =
    //                                                        =
    //=========================================================
    if(firstMuFinal == 0 || secondMuFinal == 0) continue; 
    float mRecoMass = Mass(PtRecMu1,EtaTracker1,PhiTracker1,EnRecMu1,
			   PtRecMu2,EtaTracker2,PhiTracker2,EnRecMu2);  
    //=========================================================
    //        start doing matching between reco & gen         =
    //                                                        =
    //=========================================================
    bool GenRecoMatch1 = GenRecoMatchMu1(EtaRecMu1,PhiRecMu1,mPtGen1,mPhiGen1,mEtaGen1,mEnGen1,mGenFlag1);
    bool GenRecoMatch2 = GenRecoMatchMu2(mGenFlag1,EtaRecMu2,PhiRecMu2,mPtGen2,mPhiGen2,mEtaGen2,mEnGen2);
    //if(	GenRecoMatch1 == 0 || GenRecoMatch2 == 0 ) continue; 
    std::cout<<"I matched two muons at reco and gen in a cone 0.15"<<endl;
    //=========================================================
    //        start doing matching between reco & HLT         =
    //                                                        =
    //=========================================================
    bool RecoMuon1MatchingWithHLT = RecoHLMuonTMatching(EtaRecMu1,PhiRecMu1);
    bool RecoMuon2MatchingWithHLT = RecoHLMuonTMatching(EtaRecMu2,PhiRecMu2);
    //if(RecoMuon1MatchingWithHLT==1 || RecoMuon2MatchingWithHLT==1) 
    //{
	h1_MassRecoGenDif_->Fill((mRecoMass-MassGen)/MassGen);
	PlotleadingAndSecondleading(PtRecMu1,PtRecMu2,EnRecMu1,EnRecMu2,EtaRecMu1,EtaRecMu2,
				    PhiRecMu1,PhiRecMu2,mRecoMass);
	if( (mRecoMass-MassGen)/MassGen < MassResolution ) {NbReco++;}
	//}
  }
  //===================================================
  //                                                  =
  //          Compute Acceptance X Eff.               =
  //                                                  =
  //===================================================
  float Accxeff = float(NbReco)/NbGen;
  cout<<"Accxeff = "<<Accxeff<< endl;
  fclose (pFile);
  //================================================================================== 
  //                                                                                 =
  //                         writing Histograms to a file                            =             
  //                                                                                 =
  //==================================================================================
  TFile *output = new TFile("CMSSW720-Analyse_MC_ZprimeToEE_Mass_5000_8TeV.root","recreate");
  //TFile *output = new TFile("CMSSW539-Analyse_MC_ZprimeToMuMu_Mass_1000_8TeV_photoniso.root","recreate");
  output->cd();
  /** Write the histograms */
  h1_ZprimeGenmass_->Write();
  h1_ZprimeGenEta1_->Write();
  h1_ZprimeGenEta2_->Write();
  h1_ZprimeGenPt1_->Write();
  h1_ZprimeGenPt2_->Write();
  h1_ZprimeGenEn1_->Write();
  h1_ZprimeGenEn2_->Write();
  h1_MassRecoGenDif_->Write();
  h1_ZprimeRecomass_->Write();
  h1_LeadingRecoMuPt_->Write();
  h1_SecondLeadingMuPt_->Write();
  h1_ZprimeRecoEta1_->Write();
  h1_ZprimeRecoEta2_->Write();
  h1_ZprimeRecoPhi1_->Write();
  h1_ZprimeRecoPhi2_->Write();
  h1_ZprimeRecoEn1_->Write();
  h1_ZprimeRecoEn2_->Write();
  h1_RecoPhiMu1DifPhiMu2_->Write();
  h1_GenPhiMu1DifPhiMu2_->Write();
  
  h1_test_->Write(); 
  //------------------------------
  h1_nbPVID_->Write(); 
  h1_PtID_->Write(); 
  h1_EtaID_->Write(); 
  h1_pfDeltaBetaID_->Write();
  h1_PtEffPFdeltaBeta_->Write();
  h1_EtaEffPFdeltaBeta_->Write();
  h1_nbPVEffPFdeltaBeta_->Write();
  h1_PtHeepEleEB_->Write();
  h1_EtaHeepEleEB_->Write();
  h1_PVHeepEleEB_->Write();
  h1_PtHeepEleEE_->Write();
  h1_EtaHeepEleEE_->Write();
  h1_PVHeepEleEE_->Write();
  h1_PtHeepdeltaBetaEB_->Write();
  h1_EtaHeepdeltaBetaEB_->Write();
  h1_PVHeepdeltaBetaEB_->Write();
  h1_PtHeepdeltaBetaEE_->Write();
  h1_EtaHeepdeltaBetaEE_->Write();
  h1_PVHeepdeltaBetaEE_->Write();
  h1_pfDeltaBetaEB_->Write();
  h1_pfDeltaBetaEE_->Write();

  h1_PtHeepphotonconvEB_->Write();
  h1_EtaHeepphotonconvEB_->Write();
  h1_PVHeepphotonconvEB_->Write();
  h1_PtHeepphotonconvEE_->Write();
  h1_EtaHeepphotonconvEE_->Write();
  h1_PVHeepphotonconvEE_->Write();

  h1_PtHeepHoeNewEB_->Write();
  h1_EtaHeepHoeNewEB_->Write();
  h1_PVHeepHoeNewEB_->Write();
  h1_PtHeepHoeNewEE_->Write();
  h1_EtaHeepHoeNewEE_->Write();
  h1_PVHeepHoeNewEE_->Write();

  h1_PtHeepDeltaEtaNewEB_->Write();
  h1_EtaHeepDeltaEtaNewEB_->Write();
  h1_PVHeepDeltaEtaNewEB_->Write();
  h1_PtHeepDeltaEtaNewEE_->Write();
  h1_EtaHeepDeltaEtaNewEE_->Write();
  h1_PVHeepDeltaEtaNewEE_->Write();
  
  
  h1_PtHeepDeltaPhiNewEB_->Write();
  h1_EtaHeepDeltaPhiNewEB_->Write();
  h1_PVHeepDeltaPhiNewEB_->Write();
  h1_PtHeepDeltaPhiNewEE_->Write();
  h1_EtaHeepDeltaPhiNewEE_->Write();
  h1_PVHeepDeltaPhiNewEE_->Write();
  
  
  h1_PtHeepTrackIsoNewEB_->Write();
  h1_EtaHeepTrackIsoNewEB_->Write();
  h1_PVHeepTrackIsoNewEB_->Write();
  h1_PtHeepTrackIsoNewEE_->Write();
  h1_EtaHeepTrackIsoNewEE_->Write();
  h1_PVHeepTrackIsoNewEE_->Write();
  
  h1_PtHeepMissHitsNewEB_->Write();
  h1_EtaHeepMissHitsNewEB_->Write();
  h1_PVHeepMissHitsNewEB_->Write();
  h1_PtHeepMissHitsNewEE_->Write();
  h1_EtaHeepMissHitsNewEE_->Write();
  h1_PVHeepMissHitsNewEE_->Write();

  h1_PtHeepSigmaNewEE_->Write();
  h1_EtaHeepSigmaNewEE_->Write();
  h1_PVHeepSigmaNewEE_->Write();

  h1_PtHeepEcalHcalIsoNewEB_->Write();
  h1_EtaHeepEcalHcalIsoNewEB_->Write();
  h1_PVHeepEcalHcalIsoNewEB_->Write();
  h1_PtHeepEcalHcalIsoNewEE_->Write();
  h1_EtaHeepEcalHcalIsoNewEE_->Write();
  h1_PVHeepEcalHcalIsoNewEE_->Write();
  

  h1_PtHeepShowerShapeNewEB_->Write();
  h1_EtaHeepShowerShapeNewEB_->Write();
  h1_PVHeepShowerShapeNewEB_->Write();

  h1_pfSumChargedHadronEB_->Write();
  h1_pfSumNeutralHadronEB_->Write();
  h1_pfPhotonIsoEB_->Write();
  h1_pfSumChargedHadronEE_->Write();
  h1_pfSumNeutralHadronEE_->Write();
  h1_pfPhotonIsoEE_->Write();

  h1_PtHeepEleV4EB_->Write();
  h1_EtaHeepEleV4EB_->Write();
  h1_PVHeepEleV4EB_->Write();
  h1_PtHeepEleV4EE_->Write();
  h1_EtaHeepEleV4EE_->Write();
  h1_PVHeepEleV4EE_->Write();

  h1_PtHeepHoeV4EB_->Write();
  h1_EtaHeepHoeV4EB_->Write();
  h1_PVHeepHoeV4EB_->Write();
  h1_PtHeepHoeV4EE_->Write();
  h1_EtaHeepHoeV4EE_->Write();
  h1_PVHeepHoeV4EE_->Write();

  h1_PtHeepDeltaEtaInV4EB_->Write();
  h1_EtaHeepDeltaEtaInV4EB_->Write();
  h1_PVHeepDeltaEtaInV4EB_->Write();
  h1_PtHeepDeltaEtaInV4EE_->Write();
  h1_EtaHeepDeltaEtaInV4EE_->Write();
  h1_PVHeepDeltaEtaInV4EE_->Write();

  h1_PtHeepDxyNewEB_->Write();
  h1_EtaHeepDxyNewEB_->Write();
  h1_PVHeepDxyNewEB_->Write();
  h1_PtHeepDxyNewEE_->Write();
  h1_EtaHeepDxyNewEE_->Write();
  h1_PVHeepDxyNewEE_->Write();
  h1_DxyNewEB_->Write();
  h1_DxyNewEE_->Write();

  output->Close();
  //========================================================================
  time (&end);
  dif = difftime (end,start);
  printf ("It took you %.2lf minutes to run your program.\n", (dif/60.0) );
}
//================================================================================== 
//                                                                                 =
//                                    Start methods                                =             
//                                                                                 =
//==================================================================================
//========================== Method to select firt Gen Mu =======================
bool ZprimeToEE::SelectFirstGenMu(float &ETMu1,float &PhiSCMu1,
				  float &EtaSCMu1,float &EnMu1,
				  int &IDele1,int &Statele1,
				  unsigned &GenFlag1){
  int NbHEEPele = 0;
  int iflag = -10;
  ETMu1 = 0.0;
  for(unsigned i=0; i<iGen->size(); i++){
    if( fabs(idGen->at(i)) != 11 ) continue;
    if( statusGen->at(i) != 1 )  continue;
    if( ptGen->at(i) > ETMu1) {
      ETMu1 = ptGen->at(i);
      iflag  = i;
      NbHEEPele ++;
    }
    else continue;
  }
  if(NbHEEPele>0) {
    GenFlag1       = iflag;
    ETMu1          = ptGen->at(iflag);
    PhiSCMu1       = phiGen->at(iflag);
    EtaSCMu1       = etaGen->at(iflag);
    EnMu1          = EnergyGen->at(iflag);
    IDele1         = idGen->at(iflag);
    Statele1       = statusGen->at(iflag);
    return true;
  }         
  else return false;
}
//============================ Method to select second Gen Mu ========================
bool ZprimeToEE::SelectSecondGenMu(unsigned GenFlag1,float ETMu1,float &ETMu2,float &PhiSCMu2,
				   float &EtaSCMu2,float &EnMu2,int &IDele2,int &Statele2){
  int NbHEEPele = 0;
  int iflag = -10;
  ETMu2 = 0.0;
  for(unsigned i=0; i<iGen->size(); i++){
    if( fabs(idGen->at(i)) != 11 ) continue;
    if( statusGen->at(i) != 1 )  continue;
    if(i == GenFlag1) continue;
    if( fabs(ptGen->at(i) - ETMu1) <0.00001 ) continue;
    if( ptGen->at(i) > ETMu2) {
      ETMu2 = ptGen->at(i);
      iflag  = i;
      NbHEEPele ++;
    }
    else continue;
  }
  if(NbHEEPele>0) {
    ETMu2      = ptGen->at(iflag);
    PhiSCMu2   = phiGen->at(iflag);
    EtaSCMu2   = etaGen->at(iflag);
    EnMu2      = EnergyGen->at(iflag);
    IDele2     = idGen->at(iflag);
    Statele2   = statusGen->at(iflag); 
    return true;
  }
  else return false;
}
//============================ Method to plot Gen Mu ========================
void ZprimeToEE::PlotGenInfo(float ZprimeGenMass,float EtaGenMu1,float EtaGenMu2,float PtGenMu1,
			     float PtGenMu2,float EnGenMu1,float EnGenMu2,float PhiGenMu1,
			     float PhiGenMu2){
  h1_ZprimeGenmass_->Fill(ZprimeGenMass);
  h1_ZprimeGenEta1_->Fill(EtaGenMu1);
  h1_ZprimeGenEta2_->Fill(EtaGenMu2);
  h1_ZprimeGenPt1_->Fill(PtGenMu1);
  h1_ZprimeGenPt2_->Fill(PtGenMu2);
  h1_ZprimeGenEn1_->Fill(EnGenMu1);
  h1_ZprimeGenEn2_->Fill(EnGenMu2);
  //float GenPhiMu1DifPhiMu2 = fabs((PhiGenMu1-PhiGenMu2)*(180.0/PI));
  float CosAngle = cos((PhiGenMu1-PhiGenMu2)*(180.0/PI));
  h1_GenPhiMu1DifPhiMu2_->Fill(CosAngle);
}
//============================ Method to select first high pt muon ========================
bool ZprimeToEE::SelectFirstMuon(float &pTmuon1,float &Enmuon1,float &Etamuon1,
				 float &Phimuon1,unsigned &FlagMu1,
				 float &EtaTrak1,float &PhiTrak1){
  int NbHEEPele = 0;
  unsigned iflag = -10;
  for(unsigned i=0; i<Ele_nbElectrons->size(); i++){
    if( fabs(Ele_etaSC->at(i)) < 1.442 && 
	Ele_Et->at(i) > 35.0 && Ele_isEcalDrivenSeed->at(i) == 1 &&
	Ele_nbOfMissingHits->at(i) < 2 &&
	fabs( Ele_deltaEtaInSC->at(i) ) < std::max(0.016-1*pow(10,-4)*Ele_Et->at(i),0.004) && 
	fabs( Ele_deltaPhiInSC->at(i) ) < 0.06 &&
	Ele_hadronicOverEm->at(i) < 2.0/Ele_energySC->at(i)+0.05 &&
	(Ele_e2x5MaxOver5x5Full5x5->at(i) > 0.94 || Ele_e1x5Over5x5Full5x5->at(i) > 0.83) &&
	(Ele_EcalPlusHcald1iso->at(i)) < (2.0+0.03*(Ele_Et->at(i))+(Ele_rhoIso->at(i)*0.28)) &&
	Ele_dr03TkSumPt->at(i) < 5.0 &&
	fabs(Ele_dxy->at(i)) < 0.02 ) {
      iflag  = i;
      NbHEEPele ++;
    }
    
    else if( fabs(Ele_etaSC->at(i)) > 1.560 && fabs(Ele_etaSC->at(i)) < 2.5 &&
	     Ele_Et->at(i) > 35.0 && Ele_Et->at(i) < 50.0 && Ele_isEcalDrivenSeed->at(i) == 1 &&
	     Ele_nbOfMissingHits->at(i) < 2 &&
	     fabs( Ele_deltaEtaInSC->at(i) ) < std::max(0.015-8.5*pow(10,-5)*Ele_Et->at(i),0.006) && 
	     fabs( Ele_deltaPhiInSC->at(i) ) < 0.06 &&
	     Ele_hadronicOverEm->at(i)  < 12.5/Ele_energySC->at(i)+0.05 &&
	     Ele_sigmaIetaIeta->at(i) < 0.03 &&
	     (Ele_EcalPlusHcald1iso->at(i)) < (2.5+(Ele_rhoIso->at(i)*0.28)) &&
	     Ele_dr03TkSumPt->at(i) < 5.0 &&
	     fabs(Ele_dxy->at(i)) < 0.05 ) {
      iflag  = i;
      NbHEEPele ++;
    }
    else if( fabs(Ele_etaSC->at(i)) > 1.560 && fabs(Ele_etaSC->at(i)) < 2.5 && 
	     Ele_Et->at(i) > 50.0 && Ele_isEcalDrivenSeed->at(i) == 1 &&
	     Ele_nbOfMissingHits->at(i) < 2 &&
	     fabs( Ele_deltaEtaInSC->at(i) ) < std::max(0.015-8.5*pow(10,-5)*Ele_Et->at(i),0.006) && 
	     fabs( Ele_deltaPhiInSC->at(i) ) < 0.06 &&
	     Ele_hadronicOverEm->at(i) < 12.5/Ele_energySC->at(i)+0.05 &&
	     Ele_sigmaIetaIeta->at(i) < 0.03 &&
	     (Ele_EcalPlusHcald1iso->at(i)) < (2.5+0.03*(Ele_Et->at(i)-50.0)+(Ele_rhoIso->at(i)*0.28)) &&
	     Ele_dr03TkSumPt->at(i) < 5.0 &&
	     fabs(Ele_dxy->at(i)) < 0.05 ) {
      iflag  = i;
      NbHEEPele ++;
    }
    
    else continue;
  }
  if( NbHEEPele > 0 ){
    FlagMu1             = iflag;
    pTmuon1             = Ele_Et->at(iflag);
    Enmuon1             = Ele_energySC->at(iflag);
    Etamuon1            = Ele_etaSC->at(iflag);
    Phimuon1            = Ele_phiSC->at(iflag);
    EtaTrak1            = Ele_etaTrack->at(iflag);
    PhiTrak1            = Ele_phiTrack->at(iflag);
    //printf ("AS(%d) Et = %f Eta = %f  Phi = %f En = %f\n",iflag,Mu_pt->at(iflag),Mu_eta->at(iflag),Mu_phi->at(iflag),Mu_en->at(iflag));
    return true;  
  }
  else return false;
}
//============================ Method to select second high pt muon ========================
bool ZprimeToEE::SelectSecondMuon(unsigned FlagMu1,float pTmuon1,float &pTmuon2,float &Enmuon2,
				  float &Etamuon2,float &Phimuon2,float &EtaTrak2,float &PhiTrak2)
{
  int NbHEEPele = 0;
  unsigned iflag = -10;
  for(unsigned i=0; i<Ele_nbElectrons->size(); i++){
    if(i == FlagMu1) continue;
    if(Ele_Et->at(i) == pTmuon1) continue;
    if( fabs(Ele_etaSC->at(i)) < 1.442 && 
	Ele_Et->at(i) > 35.0 && Ele_isEcalDrivenSeed->at(i) == 1 &&
	Ele_nbOfMissingHits->at(i) < 2 &&
	fabs( Ele_deltaEtaInSC->at(i) ) < std::max(0.016-1*pow(10,-4)*Ele_Et->at(i),0.004) && 
	fabs( Ele_deltaPhiInSC->at(i) ) < 0.06 &&
	Ele_hadronicOverEm->at(i) < 2.0/Ele_energySC->at(i)+0.05 &&
	(Ele_e2x5MaxOver5x5Full5x5->at(i) > 0.94 || Ele_e1x5Over5x5Full5x5->at(i) > 0.83) &&
	(Ele_EcalPlusHcald1iso->at(i)) < (2.0+0.03*(Ele_Et->at(i))+(Ele_rhoIso->at(i)*0.28)) &&
	Ele_dr03TkSumPt->at(i) < 5.0 &&
	fabs(Ele_dxy->at(i)) < 0.02 ) {
      iflag  = i;
      NbHEEPele ++;
    }
    
    else if( fabs(Ele_etaSC->at(i)) > 1.560 && fabs(Ele_etaSC->at(i)) < 2.5 &&
	     Ele_Et->at(i) > 35.0 && Ele_Et->at(i) < 50.0 && Ele_isEcalDrivenSeed->at(i) == 1 &&
	     Ele_nbOfMissingHits->at(i) < 2 &&
	     fabs( Ele_deltaEtaInSC->at(i) ) < std::max(0.015-8.5*pow(10,-5)*Ele_Et->at(i),0.006) && 
	     fabs( Ele_deltaPhiInSC->at(i) ) < 0.06 &&
	     Ele_hadronicOverEm->at(i)  < 12.5/Ele_energySC->at(i)+0.05 &&
	     Ele_sigmaIetaIeta->at(i) < 0.03 &&
	     (Ele_EcalPlusHcald1iso->at(i)) < (2.5+(Ele_rhoIso->at(i)*0.28)) &&
	     Ele_dr03TkSumPt->at(i) < 5.0 &&
	     fabs(Ele_dxy->at(i)) < 0.05 ) {
      iflag  = i;
      NbHEEPele ++;
    }
    else if( fabs(Ele_etaSC->at(i)) > 1.560 && fabs(Ele_etaSC->at(i)) < 2.5 && 
	     Ele_Et->at(i) > 50.0 && Ele_isEcalDrivenSeed->at(i) == 1 &&
	     Ele_nbOfMissingHits->at(i) < 2 &&
	     fabs( Ele_deltaEtaInSC->at(i) ) < std::max(0.015-8.5*pow(10,-5)*Ele_Et->at(i),0.006) && 
	     fabs( Ele_deltaPhiInSC->at(i) ) < 0.06 &&
	     Ele_hadronicOverEm->at(i) < 12.5/Ele_energySC->at(i)+0.05 &&
	     Ele_sigmaIetaIeta->at(i) < 0.03 &&
	     (Ele_EcalPlusHcald1iso->at(i)) < (2.5+0.03*(Ele_Et->at(i)-50.0)+(Ele_rhoIso->at(i)*0.28)) &&
	     Ele_dr03TkSumPt->at(i) < 5.0 &&
	     fabs(Ele_dxy->at(i)) < 0.05 ) {
      iflag  = i;
      NbHEEPele ++;
    }
    
    else continue;
  }
  if( NbHEEPele > 0 ){
    pTmuon2          = Ele_Et->at(iflag);
    Enmuon2          = Ele_energySC->at(iflag);
    Etamuon2         = Ele_etaSC->at(iflag);
    Phimuon2         = Ele_phiSC->at(iflag);
    EtaTrak2         = Ele_etaTrack->at(iflag);
    PhiTrak2         = Ele_phiTrack->at(iflag);
    return true;  
  }
  else return false;
}
void ZprimeToEE::PlotleadingAndSecondleading(float ptMuon1,float ptMuon2,
					     float EnMuon1,float EnMuon2,
					     float EtaMuon1,float EtaMuon2,
					     float PhiMuon1,float PhiMuon2,
					     float momentumMassMu){
  SelectHEEPelectrons();
  SelectHEEPplusPFcut();
  SelectHEEPplusPhotonConvcut();
  SelectHEEPHoEnewcut();
  SelectHEEPDeltaEtanewcut();
  SelectHEEPDeltaPhinewcut();
  SelectHEEPTrackIsonewcut();
  SelectHEEPMissHitsnewcut();
  SelectHEEPSigmaIetaIetanewcut();
  SelectHEEPEcalHcalIsonewcut();
  SelectHEEPShowerShapenewcut();

  SelectHEEPelectronsV4();
  SelectHEEPHoEV4cut();
  SelectHEEPDeltaEtaInV4cut();
  SelectHEEPV5Dxtcut();
  h1_ZprimeRecomass_->Fill(momentumMassMu);
  float CosAngle = cos((PhiMuon1-PhiMuon2)*(180.0/PI));
  h1_RecoPhiMu1DifPhiMu2_->Fill(CosAngle);
  std::vector<float> PtVector;
  std::vector<float> EnVector;
  std::vector<float> EtaVector;
  std::vector<float> PhiVector;
  PtVector.clear();
  EnVector.clear();
  EtaVector.clear();
  PhiVector.clear();
  PtVector.push_back(ptMuon1);
  PtVector.push_back(ptMuon2);
  EnVector.push_back(EnMuon1);
  EnVector.push_back(EnMuon2);  
  EtaVector.push_back(EtaMuon1);
  EtaVector.push_back(EtaMuon2);
  PhiVector.push_back(PhiMuon1);
  PhiVector.push_back(PhiMuon2);
  std::sort (PtVector.begin(), PtVector.end(), myfunction);
  for(unsigned i=0; i<PtVector.size(); i++){
    //printf ("BS[%d] Pt = %f Eta = %f  En = %f\n",i,PtVector.at(i),EtaVector.at(i),EnVector.at(i));
    if(i==1) 
      {
	h1_LeadingRecoMuPt_->Fill(PtVector.at(i));
	h1_ZprimeRecoEn2_->Fill(EnVector.at(i));
	h1_ZprimeRecoEta2_->Fill(EtaVector.at(i));
	h1_ZprimeRecoPhi2_->Fill(PhiVector.at(i));
	//printf ("AS[%d] Pt = %f Eta = %f  En = %f\n",i,PtVector.at(i),EtaVector.at(i),EnVector.at(i));
      }
    if(i==0) 
      {
	h1_SecondLeadingMuPt_->Fill(PtVector.at(i));
	h1_ZprimeRecoEn1_->Fill(EnVector.at(i));
	h1_ZprimeRecoEta1_->Fill(EtaVector.at(i));
	h1_ZprimeRecoPhi1_->Fill(PhiVector.at(i));
      }
  }
}
//===================== Methode to calculate the mass ========================
float ZprimeToEE::Mass(float Pt1,float Eta1,float Phi1,float En1,
		       float Pt2,float Eta2,float Phi2,float En2){
  float MuMuMass = 0.0;
  TLorentzVector Mu1;
  TLorentzVector Mu2;
  Mu1.SetPtEtaPhiE(Pt1,Eta1,Phi1,En1);
  Mu2.SetPtEtaPhiE(Pt2,Eta2,Phi2,En2);
  MuMuMass = (Mu1 + Mu2).M();
  return MuMuMass;
}
//----------------------------------------------------
//                                                   -
//       Part for HLT & Reco Matching                -
//                                                   -  
//----------------------------------------------------
bool ZprimeToEE::RecoHLMuonTMatching(float RecoEta,float RecoPhi){
  int nbMatch = 0;
  for(unsigned i=0; i<MuHLTObj_nbObj->size(); i++){
    float deltaEta = RecoEta - MuHLTObj_eta->at(i);
    float deltaPhi = RecoPhi - MuHLTObj_phi->at(i);
    float deltaR   = sqrt(pow(deltaEta,2)+pow(deltaPhi,2));
    if(fabs(deltaR)>RecoHLTMatchingDeltaRcut) continue;
    nbMatch++;
  }
  if(nbMatch>0) return true;
  else return false;
}

//----------------------------------------------------
//                                                   -
//       Part for Gen & Reco Matching                -
//                                                   -  
//----------------------------------------------------
bool ZprimeToEE::GenRecoMatchMu1(float RecoEta1,float RecoPhi1,
				 float &ETMu1,float &PhiSCMu1,
				 float &EtaSCMu1,float &EnMu1,
				 unsigned &GenFlag1){
  int NbHEEPele = 0;
  unsigned iflag = -10;
  for(unsigned i=0; i<iGen->size(); i++){
    //printf ("BS[%d] genEta = %f  genPhi = %f recoEta = %f recoPhi = %f deltaR = %f\n",
    //    i,etaGen->at(i),phiGen->at(i),RecoEta1,RecoPhi1,deltaR);
    if( fabs(idGen->at(i)) != 11 ) continue;
    if( statusGen->at(i) != 1 )  continue;
    float deltaEta = RecoEta1 - etaGen->at(i);
    float deltaPhi = RecoPhi1 - phiGen->at(i);
    float deltaR   = sqrt(pow(deltaEta,2)+pow(deltaPhi,2));
    h1_deltaR_->Fill(deltaR);
    if(fabs(deltaR)>RecoGenMatchingDeltaRcut) continue;
    iflag  = i;
    NbHEEPele ++;
  }
  if(NbHEEPele > 0) {
    GenFlag1       = iflag;
    ETMu1          = ptGen->at(iflag);
    PhiSCMu1       = phiGen->at(iflag);
    EtaSCMu1       = etaGen->at(iflag);
    EnMu1          = EnergyGen->at(iflag);
    //    printf ("AS[%d] genEta = %f  genPhi = %f recoEta = %f recoPhi = %f\n",iflag,etaGen->at(iflag),phiGen->at(iflag),RecoEta1,RecoPhi1);
    return true;         
  }
  else return false;
}

bool ZprimeToEE::GenRecoMatchMu2(unsigned GenFlag1,float RecoEta2,float RecoPhi2,
				 float &ETMu2,float &PhiSCMu2,
				 float &EtaSCMu2,float &EnMu2){
  int NbHEEPele = 0;
  unsigned iflag = -10;
  for(unsigned i=0; i<iGen->size(); i++){
    //printf ("BS[%d] genEta = %f  genPhi = %f recoEta = %f recoPhi = %f deltaR = %f\n",
    //    i,etaGen->at(i),phiGen->at(i),RecoEta1,RecoPhi1,deltaR);
    if( fabs(idGen->at(i)) != 11 ) continue;
    if( statusGen->at(i) != 1 )  continue;
    if( i ==  GenFlag1)  continue;
    float deltaEta2 = RecoEta2 - etaGen->at(i);
    float deltaPhi2 = RecoPhi2 - phiGen->at(i);
    float deltaR2   = sqrt(pow(deltaEta2,2)+pow(deltaPhi2,2));
    h1_deltaR_->Fill(deltaR2);
    if(fabs(deltaR2)>RecoGenMatchingDeltaRcut) continue;
    iflag  = i;
    NbHEEPele ++;
  }
  if(NbHEEPele > 0) {
    ETMu2          = ptGen->at(iflag);
    PhiSCMu2       = phiGen->at(iflag);
    EtaSCMu2       = etaGen->at(iflag);
    EnMu2          = EnergyGen->at(iflag);
    //    printf ("AS[%d] genEta = %f  genPhi = %f recoEta = %f recoPhi = %f\n",iflag,etaGen->at(iflag),phiGen->at(iflag),RecoEta1,RecoPhi1);
    return true;         
  }
  else return false;
}

//============================ Method to select first high pt muon ========================
void ZprimeToEE::SelectHEEPelectrons(){
  for(unsigned i=0; i<Ele_nbElectrons->size(); i++){
    float charged = Ele_pfSumChargedHadronPt->at(i);
    float neutral = Ele_pfSumNeutralHadronEt->at(i);
    float pileup  = Ele_pfSumPUPt->at(i);
    float photoniso = Ele_pfSumPhotonEt->at(i);
    float deltaBeta = (charged + std::max(0.0, neutral + photoniso - 0.5 * pileup))/Ele_Et->at(i);
    if( fabs(Ele_etaSC->at(i)) < 1.442 && 
	Ele_Et->at(i) > Etcut && Ele_isEcalDrivenSeed->at(i) == 1 &&
	Ele_nbOfMissingHits->at(i) < 2 &&
	fabs( Ele_deltaEtaInSC->at(i) ) < std::max(0.016-1*pow(10,-4)*Ele_Et->at(i),0.004) && 
	fabs( Ele_deltaPhiInSC->at(i) ) < 0.06 &&
	Ele_hadronicOverEm->at(i) < 2.0/Ele_energySC->at(i)+0.05 &&
	(Ele_e2x5MaxOver5x5Full5x5->at(i) > 0.94 || Ele_e1x5Over5x5Full5x5->at(i) > 0.83) &&
	(Ele_EcalPlusHcald1iso->at(i)) < (2.0+0.03*(Ele_Et->at(i))+(Ele_rhoIso->at(i)*0.28)) &&
	Ele_dr03TkSumPt->at(i) < 5.0 &&
	fabs(Ele_dxy->at(i)) < 0.02 ) {
      if(Ele_Et->at(i)>minPt && Ele_Et->at(i)<maxPt)
	{
	  h1_PtHeepEleEB_->Fill(Ele_Et->at(i));
	  h1_EtaHeepEleEB_->Fill(fabs(Ele_etaSC->at(i)));
	  h1_PVHeepEleEB_->Fill(Ele_nbVtx->at(i));
	  h1_pfDeltaBetaEB_->Fill(deltaBeta);
	  h1_pfSumChargedHadronEB_->Fill(charged);
	  h1_pfSumNeutralHadronEB_->Fill(neutral);
	  h1_pfPhotonIsoEB_->Fill(photoniso);
	}
    }
    else if( fabs(Ele_etaSC->at(i)) > 1.560 && fabs(Ele_etaSC->at(i)) < 2.5 &&
	     Ele_Et->at(i) > Etcut && Ele_Et->at(i) < 50.0 && Ele_isEcalDrivenSeed->at(i) == 1 &&
	     Ele_nbOfMissingHits->at(i) < 2 &&
	     fabs( Ele_deltaEtaInSC->at(i) ) < std::max(0.015-8.5*pow(10,-5)*Ele_Et->at(i),0.006) && 
	     fabs( Ele_deltaPhiInSC->at(i) ) < 0.06 &&
	     Ele_hadronicOverEm->at(i)  < 12.5/Ele_energySC->at(i)+0.05 &&
	     Ele_sigmaIetaIeta->at(i) < 0.03 &&
	     (Ele_EcalPlusHcald1iso->at(i)) < (2.5+(Ele_rhoIso->at(i)*0.28)) &&
	     Ele_dr03TkSumPt->at(i) < 5.0 &&
	     fabs(Ele_dxy->at(i)) < 0.05 ) {
      if(Ele_Et->at(i)>minPt && Ele_Et->at(i)<maxPt)
	{
	  h1_PtHeepEleEE_->Fill(Ele_Et->at(i));
	  h1_EtaHeepEleEE_->Fill(fabs(Ele_etaSC->at(i)));
	  h1_PVHeepEleEE_->Fill(Ele_nbVtx->at(i));
	  h1_pfDeltaBetaEE_->Fill(deltaBeta);
	  h1_pfSumChargedHadronEE_->Fill(charged);
	  h1_pfSumNeutralHadronEE_->Fill(neutral);
	  h1_pfPhotonIsoEE_->Fill(photoniso);
	}
    }
    else if( fabs(Ele_etaSC->at(i)) > 1.560 && fabs(Ele_etaSC->at(i)) < 2.5 && 
	     Ele_Et->at(i) > 50.0 && Ele_isEcalDrivenSeed->at(i) == 1 &&
	     Ele_nbOfMissingHits->at(i) < 2 &&
	     fabs( Ele_deltaEtaInSC->at(i) ) < std::max(0.015-8.5*pow(10,-5)*Ele_Et->at(i),0.006) && 
	     fabs( Ele_deltaPhiInSC->at(i) ) < 0.06 &&
	     Ele_hadronicOverEm->at(i) < 12.5/Ele_energySC->at(i)+0.05 &&
	     Ele_sigmaIetaIeta->at(i) < 0.03 &&
	     (Ele_EcalPlusHcald1iso->at(i)) < (2.5+0.03*(Ele_Et->at(i)-50.0)+(Ele_rhoIso->at(i)*0.28)) &&
	     Ele_dr03TkSumPt->at(i) < 5.0 &&
	     fabs(Ele_dxy->at(i)) < 0.05 ) {
      if(Ele_Et->at(i)>minPt && Ele_Et->at(i)<maxPt)
	{
	  h1_PtHeepEleEE_->Fill(Ele_Et->at(i));
	  h1_EtaHeepEleEE_->Fill(fabs(Ele_etaSC->at(i)));
	  h1_PVHeepEleEE_->Fill(Ele_nbVtx->at(i));
	  h1_pfDeltaBetaEE_->Fill(deltaBeta);
	  h1_pfSumChargedHadronEE_->Fill(charged);
	  h1_pfSumNeutralHadronEE_->Fill(neutral);
	  h1_pfPhotonIsoEE_->Fill(photoniso);
	}
    }
    else continue;
  }
}
  
void ZprimeToEE::SelectHEEPplusPFcut(){
  for(unsigned i=0; i<Ele_nbElectrons->size(); i++){
    float charged = Ele_pfSumChargedHadronPt->at(i);
    float neutral = Ele_pfSumNeutralHadronEt->at(i);
    float pileup  = Ele_pfSumPUPt->at(i);
    float photoniso = Ele_pfSumPhotonEt->at(i);  //
    float deltaBeta = (charged + std::max(0.0, neutral + photoniso - 0.5 * pileup))/Ele_Et->at(i);
    if( fabs(Ele_etaSC->at(i)) < 1.442 && 
	Ele_Et->at(i) > Etcut && Ele_isEcalDrivenSeed->at(i) == 1 &&
	Ele_nbOfMissingHits->at(i) < 2 &&
	fabs( Ele_deltaEtaInSC->at(i) ) < std::max(0.016-1*pow(10,-4)*Ele_Et->at(i),0.004) && 
	fabs( Ele_deltaPhiInSC->at(i) ) < 0.06 &&
	Ele_hadronicOverEm->at(i) < 2.0/Ele_energySC->at(i)+0.05 &&
	(Ele_e2x5MaxOver5x5Full5x5->at(i) > 0.94 || Ele_e1x5Over5x5Full5x5->at(i) > 0.83) &&
	(Ele_EcalPlusHcald1iso->at(i)) < (2.0+0.03*(Ele_Et->at(i))+(Ele_rhoIso->at(i)*0.28)) &&
	Ele_dr03TkSumPt->at(i) < 5.0 &&
	fabs(Ele_dxy->at(i)) < 0.02 &&
	deltaBeta<0.10 ) {
      if(Ele_Et->at(i)>minPt && Ele_Et->at(i)<maxPt)
	{
	  h1_PtHeepdeltaBetaEB_->Fill(Ele_Et->at(i));
	  h1_EtaHeepdeltaBetaEB_->Fill(fabs(Ele_etaSC->at(i)));
	  h1_PVHeepdeltaBetaEB_->Fill(Ele_nbVtx->at(i));
	}
    }
    
    else if( fabs(Ele_etaSC->at(i)) > 1.560 && fabs(Ele_etaSC->at(i)) < 2.5 &&
	     Ele_Et->at(i) > Etcut && Ele_Et->at(i) < 50.0 && Ele_isEcalDrivenSeed->at(i) == 1 &&
	     Ele_nbOfMissingHits->at(i) < 2 &&
	     fabs( Ele_deltaEtaInSC->at(i) ) < std::max(0.015-8.5*pow(10,-5)*Ele_Et->at(i),0.006) && 
	     fabs( Ele_deltaPhiInSC->at(i) ) < 0.06 &&
	     Ele_hadronicOverEm->at(i)  < 12.5/Ele_energySC->at(i)+0.05 &&
	     Ele_sigmaIetaIeta->at(i) < 0.03 &&
	     (Ele_EcalPlusHcald1iso->at(i)) < (2.5+(Ele_rhoIso->at(i)*0.28)) &&
	     Ele_dr03TkSumPt->at(i) < 5.0 &&
	     fabs(Ele_dxy->at(i)) < 0.05  &&
	     deltaBeta<0.10 ) {
      if(Ele_Et->at(i)>minPt && Ele_Et->at(i)<maxPt)
	{
	  h1_PtHeepdeltaBetaEE_->Fill(Ele_Et->at(i));
	  h1_EtaHeepdeltaBetaEE_->Fill(fabs(Ele_etaSC->at(i)));
	  h1_PVHeepdeltaBetaEE_->Fill(Ele_nbVtx->at(i));
	}
    }
    else if( fabs(Ele_etaSC->at(i)) > 1.560 && fabs(Ele_etaSC->at(i)) < 2.5 && 
	     Ele_Et->at(i) > 50.0 && Ele_isEcalDrivenSeed->at(i) == 1 &&
	     Ele_nbOfMissingHits->at(i) < 2 &&
	     fabs( Ele_deltaEtaInSC->at(i) ) < std::max(0.015-8.5*pow(10,-5)*Ele_Et->at(i),0.006) && 
	     fabs( Ele_deltaPhiInSC->at(i) ) < 0.06 &&
	     Ele_hadronicOverEm->at(i) < 12.5/Ele_energySC->at(i)+0.05 &&
	     Ele_sigmaIetaIeta->at(i) < 0.03 &&
	     (Ele_EcalPlusHcald1iso->at(i)) < (2.5+0.03*(Ele_Et->at(i)-50.0)+(Ele_rhoIso->at(i)*0.28)) &&
	     Ele_dr03TkSumPt->at(i) < 5.0 &&
	     fabs(Ele_dxy->at(i)) < 0.05  &&
	     deltaBeta<0.10 ) {
      if(Ele_Et->at(i)>minPt && Ele_Et->at(i)<maxPt)
	{
	  h1_PtHeepdeltaBetaEE_->Fill(Ele_Et->at(i));
	  h1_EtaHeepdeltaBetaEE_->Fill(fabs(Ele_etaSC->at(i)));
	  h1_PVHeepdeltaBetaEE_->Fill(Ele_nbVtx->at(i));
	}
    }
    else continue;
  }
}
  
void ZprimeToEE::SelectHEEPplusPhotonConvcut(){
  for(unsigned i=0; i<Ele_nbElectrons->size(); i++){
    if( fabs(Ele_etaSC->at(i)) < 1.442 && 
	Ele_Et->at(i) > Etcut && Ele_isEcalDrivenSeed->at(i) == 1 &&
	Ele_nbOfMissingHits->at(i) < 2 &&
	fabs( Ele_deltaEtaInSC->at(i) ) < std::max(0.016-1*pow(10,-4)*Ele_Et->at(i),0.004) && 
	fabs( Ele_deltaPhiInSC->at(i) ) < 0.06 &&
	Ele_hadronicOverEm->at(i) < 2.0/Ele_energySC->at(i)+0.05 &&
	(Ele_e2x5MaxOver5x5Full5x5->at(i) > 0.94 || Ele_e1x5Over5x5Full5x5->at(i) > 0.83) &&
	(Ele_EcalPlusHcald1iso->at(i)) < (2.0+0.03*(Ele_Et->at(i))+(Ele_rhoIso->at(i)*0.28)) &&
	Ele_dr03TkSumPt->at(i) < 5.0 &&
	fabs(Ele_dxy->at(i)) < 0.02 &&
	Ele_isPassConversionVeto->at(i) == 1 ) {
      if(Ele_Et->at(i)>minPt && Ele_Et->at(i)<maxPt)
	{
	  h1_PtHeepphotonconvEB_->Fill(Ele_Et->at(i));
	  h1_EtaHeepphotonconvEB_->Fill(fabs(Ele_etaSC->at(i)));
	  h1_PVHeepphotonconvEB_->Fill(Ele_nbVtx->at(i));
	}
    }
    
    else if( fabs(Ele_etaSC->at(i)) > 1.560 && fabs(Ele_etaSC->at(i)) < 2.5 &&
	     Ele_Et->at(i) > Etcut && Ele_Et->at(i) < 50.0 && Ele_isEcalDrivenSeed->at(i) == 1 &&
	     Ele_nbOfMissingHits->at(i) < 2 &&
	     fabs( Ele_deltaEtaInSC->at(i) ) < std::max(0.015-8.5*pow(10,-5)*Ele_Et->at(i),0.006) && 
	     fabs( Ele_deltaPhiInSC->at(i) ) < 0.06 &&
	     Ele_hadronicOverEm->at(i)  < 12.5/Ele_energySC->at(i)+0.05 &&
	     Ele_sigmaIetaIeta->at(i) < 0.03 &&
	     (Ele_EcalPlusHcald1iso->at(i)) < (2.5+(Ele_rhoIso->at(i)*0.28)) &&
	     Ele_dr03TkSumPt->at(i) < 5.0 &&
	     fabs(Ele_dxy->at(i)) < 0.05  &&
	     Ele_isPassConversionVeto->at(i) == 1 ) {
      if(Ele_Et->at(i)>minPt && Ele_Et->at(i)<maxPt)
	{
	  h1_PtHeepphotonconvEE_->Fill(Ele_Et->at(i));
	  h1_EtaHeepphotonconvEE_->Fill(fabs(Ele_etaSC->at(i)));
	  h1_PVHeepphotonconvEE_->Fill(Ele_nbVtx->at(i));
	}
    }
    else if( fabs(Ele_etaSC->at(i)) > 1.560 && fabs(Ele_etaSC->at(i)) < 2.5 && 
	     Ele_Et->at(i) > 50.0 && Ele_isEcalDrivenSeed->at(i) == 1 &&
	     Ele_nbOfMissingHits->at(i) < 2 &&
	     fabs( Ele_deltaEtaInSC->at(i) ) < std::max(0.015-8.5*pow(10,-5)*Ele_Et->at(i),0.006) && 
	     fabs( Ele_deltaPhiInSC->at(i) ) < 0.06 &&
	     Ele_hadronicOverEm->at(i) < 12.5/Ele_energySC->at(i)+0.05 &&
	     Ele_sigmaIetaIeta->at(i) < 0.03 &&
	     (Ele_EcalPlusHcald1iso->at(i)) < (2.5+0.03*(Ele_Et->at(i)-50.0)+(Ele_rhoIso->at(i)*0.28)) &&
	     Ele_dr03TkSumPt->at(i) < 5.0 &&
	     fabs(Ele_dxy->at(i)) < 0.05  &&
	     Ele_isPassConversionVeto->at(i) == 1 ) {
      if(Ele_Et->at(i)>minPt && Ele_Et->at(i)<maxPt)
	{
	  h1_PtHeepphotonconvEE_->Fill(Ele_Et->at(i));
	  h1_EtaHeepphotonconvEE_->Fill(fabs(Ele_etaSC->at(i)));
	  h1_PVHeepphotonconvEE_->Fill(Ele_nbVtx->at(i));
	}
    }
    else continue;
  }
}

void ZprimeToEE::SelectHEEPHoEnewcut(){
  for(unsigned i=0; i<Ele_nbElectrons->size(); i++){
    if( fabs(Ele_etaSC->at(i)) < 1.442 && 
	Ele_Et->at(i) > Etcut && Ele_isEcalDrivenSeed->at(i) == 1 &&
	Ele_nbOfMissingHits->at(i) < 2 &&
	fabs( Ele_deltaEtaInSC->at(i) ) < std::max(0.016-1*pow(10,-4)*Ele_Et->at(i),0.004) && 
	fabs( Ele_deltaPhiInSC->at(i) ) < 0.06 &&
	//Ele_hadronicOverEm->at(i) < 2.0/Ele_energySC->at(i)+0.05 &&
	(Ele_e2x5MaxOver5x5Full5x5->at(i) > 0.94 || Ele_e1x5Over5x5Full5x5->at(i) > 0.83) &&
	(Ele_EcalPlusHcald1iso->at(i)) < (2.0+0.03*(Ele_Et->at(i))+(Ele_rhoIso->at(i)*0.28)) &&
	Ele_dr03TkSumPt->at(i) < 5.0 &&
	fabs(Ele_dxy->at(i)) < 0.02 ) {
      if(Ele_Et->at(i)>minPt && Ele_Et->at(i)<maxPt)
	{
	  h1_PtHeepHoeNewEB_->Fill(Ele_Et->at(i));
	  h1_EtaHeepHoeNewEB_->Fill(fabs(Ele_etaSC->at(i)));
	  h1_PVHeepHoeNewEB_->Fill(Ele_nbVtx->at(i));
	}
    }
    
    else if( fabs(Ele_etaSC->at(i)) > 1.560 && fabs(Ele_etaSC->at(i)) < 2.5 &&
	     Ele_Et->at(i) > Etcut && Ele_Et->at(i) < 50.0 && Ele_isEcalDrivenSeed->at(i) == 1 &&
	     Ele_nbOfMissingHits->at(i) < 2 &&
	     fabs( Ele_deltaEtaInSC->at(i) ) < std::max(0.015-8.5*pow(10,-5)*Ele_Et->at(i),0.006) && 
	     fabs( Ele_deltaPhiInSC->at(i) ) < 0.06 &&
	     //Ele_hadronicOverEm->at(i)  < 12.5/Ele_energySC->at(i)+0.05 &&
	     Ele_sigmaIetaIeta->at(i) < 0.03 &&
	     (Ele_EcalPlusHcald1iso->at(i)) < (2.5+(Ele_rhoIso->at(i)*0.28)) &&
	     Ele_dr03TkSumPt->at(i) < 5.0 &&
	     fabs(Ele_dxy->at(i)) < 0.05 ) {
      if(Ele_Et->at(i)>minPt && Ele_Et->at(i)<maxPt)
	{
	  h1_PtHeepHoeNewEE_->Fill(Ele_Et->at(i));
	  h1_EtaHeepHoeNewEE_->Fill(fabs(Ele_etaSC->at(i)));
	  h1_PVHeepHoeNewEE_->Fill(Ele_nbVtx->at(i));
	}
    }
    else if( fabs(Ele_etaSC->at(i)) > 1.560 && fabs(Ele_etaSC->at(i)) < 2.5 && 
	     Ele_Et->at(i) > 50.0 && Ele_isEcalDrivenSeed->at(i) == 1 &&
	     Ele_nbOfMissingHits->at(i) < 2 &&
	     fabs( Ele_deltaEtaInSC->at(i) ) < std::max(0.015-8.5*pow(10,-5)*Ele_Et->at(i),0.006) && 
	     fabs( Ele_deltaPhiInSC->at(i) ) < 0.06 &&
	     //Ele_hadronicOverEm->at(i) < 12.5/Ele_energySC->at(i)+0.05 &&
	     Ele_sigmaIetaIeta->at(i) < 0.03 &&
	     (Ele_EcalPlusHcald1iso->at(i)) < (2.5+0.03*(Ele_Et->at(i)-50.0)+(Ele_rhoIso->at(i)*0.28)) &&
	     Ele_dr03TkSumPt->at(i) < 5.0 &&
	     fabs(Ele_dxy->at(i)) < 0.05 ) {
      if(Ele_Et->at(i)>minPt && Ele_Et->at(i)<maxPt)
	{
	  h1_PtHeepHoeNewEE_->Fill(Ele_Et->at(i));
	  h1_EtaHeepHoeNewEE_->Fill(fabs(Ele_etaSC->at(i)));
	  h1_PVHeepHoeNewEE_->Fill(Ele_nbVtx->at(i));
	}
    }
    else continue;
  }
}

void ZprimeToEE::SelectHEEPDeltaEtanewcut(){
  for(unsigned i=0; i<Ele_nbElectrons->size(); i++){
    if( fabs(Ele_etaSC->at(i)) < 1.442 && 
	Ele_Et->at(i) > Etcut && Ele_isEcalDrivenSeed->at(i) == 1 &&
	Ele_nbOfMissingHits->at(i) < 2 &&
	//fabs( Ele_deltaEtaInSC->at(i) ) < std::max(0.016-1*pow(10,-4)*Ele_Et->at(i),0.004) && 
	fabs( Ele_deltaPhiInSC->at(i) ) < 0.06 &&
	Ele_hadronicOverEm->at(i) < 2.0/Ele_energySC->at(i)+0.05 &&
	(Ele_e2x5MaxOver5x5Full5x5->at(i) > 0.94 || Ele_e1x5Over5x5Full5x5->at(i) > 0.83) &&
	(Ele_EcalPlusHcald1iso->at(i)) < (2.0+0.03*(Ele_Et->at(i))+(Ele_rhoIso->at(i)*0.28)) &&
	Ele_dr03TkSumPt->at(i) < 5.0 &&
	fabs(Ele_dxy->at(i)) < 0.02 ) {
      if(Ele_Et->at(i)>minPt && Ele_Et->at(i)<maxPt)
	{
	  h1_PtHeepDeltaEtaNewEB_->Fill(Ele_Et->at(i));
	  h1_EtaHeepDeltaEtaNewEB_->Fill(fabs(Ele_etaSC->at(i)));
	  h1_PVHeepDeltaEtaNewEB_->Fill(Ele_nbVtx->at(i));
	}
    }
    
    else if( fabs(Ele_etaSC->at(i)) > 1.560 && fabs(Ele_etaSC->at(i)) < 2.5 &&
	     Ele_Et->at(i) > Etcut && Ele_Et->at(i) < 50.0 && Ele_isEcalDrivenSeed->at(i) == 1 &&
	     Ele_nbOfMissingHits->at(i) < 2 &&
	     //fabs( Ele_deltaEtaInSC->at(i) ) < std::max(0.015-8.5*pow(10,-5)*Ele_Et->at(i),0.006) && 
	     fabs( Ele_deltaPhiInSC->at(i) ) < 0.06 &&
	     Ele_hadronicOverEm->at(i)  < 12.5/Ele_energySC->at(i)+0.05 &&
	     Ele_sigmaIetaIeta->at(i) < 0.03 &&
	     (Ele_EcalPlusHcald1iso->at(i)) < (2.5+(Ele_rhoIso->at(i)*0.28)) &&
	     Ele_dr03TkSumPt->at(i) < 5.0 &&
	     fabs(Ele_dxy->at(i)) < 0.05 ) {
      if(Ele_Et->at(i)>minPt && Ele_Et->at(i)<maxPt)
	{
	  h1_PtHeepDeltaEtaNewEE_->Fill(Ele_Et->at(i));
	  h1_EtaHeepDeltaEtaNewEE_->Fill(fabs(Ele_etaSC->at(i)));
	  h1_PVHeepDeltaEtaNewEE_->Fill(Ele_nbVtx->at(i));
	}
    }
    else if( fabs(Ele_etaSC->at(i)) > 1.560 && fabs(Ele_etaSC->at(i)) < 2.5 && 
	     Ele_Et->at(i) > 50.0 && Ele_isEcalDrivenSeed->at(i) == 1 &&
	     Ele_nbOfMissingHits->at(i) < 2 &&
	     //fabs( Ele_deltaEtaInSC->at(i) ) < std::max(0.015-8.5*pow(10,-5)*Ele_Et->at(i),0.006) && 
	     fabs( Ele_deltaPhiInSC->at(i) ) < 0.06 &&
	     Ele_hadronicOverEm->at(i) < 12.5/Ele_energySC->at(i)+0.05 &&
	     Ele_sigmaIetaIeta->at(i) < 0.03 &&
	     (Ele_EcalPlusHcald1iso->at(i)) < (2.5+0.03*(Ele_Et->at(i)-50.0)+(Ele_rhoIso->at(i)*0.28)) &&
	     Ele_dr03TkSumPt->at(i) < 5.0 &&
	     fabs(Ele_dxy->at(i)) < 0.05 ) {
      if(Ele_Et->at(i)>minPt && Ele_Et->at(i)<maxPt)
	{
	  h1_PtHeepDeltaEtaNewEE_->Fill(Ele_Et->at(i));
	  h1_EtaHeepDeltaEtaNewEE_->Fill(fabs(Ele_etaSC->at(i)));
	  h1_PVHeepDeltaEtaNewEE_->Fill(Ele_nbVtx->at(i));
	}
    }
    else continue;
  }
}


void ZprimeToEE::SelectHEEPDeltaPhinewcut(){
  for(unsigned i=0; i<Ele_nbElectrons->size(); i++){
    if( fabs(Ele_etaSC->at(i)) < 1.442 && 
	Ele_Et->at(i) > Etcut && Ele_isEcalDrivenSeed->at(i) == 1 &&
	Ele_nbOfMissingHits->at(i) < 2 &&
	fabs( Ele_deltaEtaInSC->at(i) ) < std::max(0.016-1*pow(10,-4)*Ele_Et->at(i),0.004) && 
	//fabs( Ele_deltaPhiInSC->at(i) ) < 0.06 &&
	Ele_hadronicOverEm->at(i) < 2.0/Ele_energySC->at(i)+0.05 &&
	(Ele_e2x5MaxOver5x5Full5x5->at(i) > 0.94 || Ele_e1x5Over5x5Full5x5->at(i) > 0.83) &&
	(Ele_EcalPlusHcald1iso->at(i)) < (2.0+0.03*(Ele_Et->at(i))+(Ele_rhoIso->at(i)*0.28)) &&
	Ele_dr03TkSumPt->at(i) < 5.0 &&
	fabs(Ele_dxy->at(i)) < 0.02 ) {
      if(Ele_Et->at(i)>minPt && Ele_Et->at(i)<maxPt)
	{
	  h1_PtHeepDeltaPhiNewEB_->Fill(Ele_Et->at(i));
	  h1_EtaHeepDeltaPhiNewEB_->Fill(fabs(Ele_etaSC->at(i)));
	  h1_PVHeepDeltaPhiNewEB_->Fill(Ele_nbVtx->at(i));
	}
    }
    
    else if( fabs(Ele_etaSC->at(i)) > 1.560 && fabs(Ele_etaSC->at(i)) < 2.5 &&
	     Ele_Et->at(i) > Etcut && Ele_Et->at(i) < 50.0 && Ele_isEcalDrivenSeed->at(i) == 1 &&
	     Ele_nbOfMissingHits->at(i) < 2 &&
	     fabs( Ele_deltaEtaInSC->at(i) ) < std::max(0.015-8.5*pow(10,-5)*Ele_Et->at(i),0.006) && 
	     //fabs( Ele_deltaPhiInSC->at(i) ) < 0.06 &&
	     Ele_hadronicOverEm->at(i)  < 12.5/Ele_energySC->at(i)+0.05 &&
	     Ele_sigmaIetaIeta->at(i) < 0.03 &&
	     (Ele_EcalPlusHcald1iso->at(i)) < (2.5+(Ele_rhoIso->at(i)*0.28)) &&
	     Ele_dr03TkSumPt->at(i) < 5.0 &&
	     fabs(Ele_dxy->at(i)) < 0.05 ) {
      if(Ele_Et->at(i)>minPt && Ele_Et->at(i)<maxPt)
	{
	  h1_PtHeepDeltaPhiNewEE_->Fill(Ele_Et->at(i));
	  h1_EtaHeepDeltaPhiNewEE_->Fill(fabs(Ele_etaSC->at(i)));
	  h1_PVHeepDeltaPhiNewEE_->Fill(Ele_nbVtx->at(i));
	}
    }
    else if( fabs(Ele_etaSC->at(i)) > 1.560 && fabs(Ele_etaSC->at(i)) < 2.5 && 
	     Ele_Et->at(i) > 50.0 && Ele_isEcalDrivenSeed->at(i) == 1 &&
	     Ele_nbOfMissingHits->at(i) < 2 &&
	     fabs( Ele_deltaEtaInSC->at(i) ) < std::max(0.015-8.5*pow(10,-5)*Ele_Et->at(i),0.006) && 
	     //fabs( Ele_deltaPhiInSC->at(i) ) < 0.06 &&
	     Ele_hadronicOverEm->at(i) < 12.5/Ele_energySC->at(i)+0.05 &&
	     Ele_sigmaIetaIeta->at(i) < 0.03 &&
	     (Ele_EcalPlusHcald1iso->at(i)) < (2.5+0.03*(Ele_Et->at(i)-50.0)+(Ele_rhoIso->at(i)*0.28)) &&
	     Ele_dr03TkSumPt->at(i) < 5.0 &&
	     fabs(Ele_dxy->at(i)) < 0.05 ) {
      if(Ele_Et->at(i)>minPt && Ele_Et->at(i)<maxPt)
	{
	  h1_PtHeepDeltaPhiNewEE_->Fill(Ele_Et->at(i));
	  h1_EtaHeepDeltaPhiNewEE_->Fill(fabs(Ele_etaSC->at(i)));
	  h1_PVHeepDeltaPhiNewEE_->Fill(Ele_nbVtx->at(i));
	}
    }
    else continue;
  }
}

void ZprimeToEE::SelectHEEPTrackIsonewcut(){
  for(unsigned i=0; i<Ele_nbElectrons->size(); i++){
    if( fabs(Ele_etaSC->at(i)) < 1.442 && 
	Ele_Et->at(i) > Etcut && Ele_isEcalDrivenSeed->at(i) == 1 &&
	Ele_nbOfMissingHits->at(i) < 2 &&
	fabs( Ele_deltaEtaInSC->at(i) ) < std::max(0.016-1*pow(10,-4)*Ele_Et->at(i),0.004) && 
	fabs( Ele_deltaPhiInSC->at(i) ) < 0.06 &&
	Ele_hadronicOverEm->at(i) < 2.0/Ele_energySC->at(i)+0.05 &&
	(Ele_e2x5MaxOver5x5Full5x5->at(i) > 0.94 || Ele_e1x5Over5x5Full5x5->at(i) > 0.83) &&
	(Ele_EcalPlusHcald1iso->at(i)) < (2.0+0.03*(Ele_Et->at(i))+(Ele_rhoIso->at(i)*0.28)) &&
	//Ele_dr03TkSumPt->at(i) < 5.0 &&
	fabs(Ele_dxy->at(i)) < 0.02 ) {
      if(Ele_Et->at(i)>minPt && Ele_Et->at(i)<maxPt)
	{
	  h1_PtHeepTrackIsoNewEB_->Fill(Ele_Et->at(i));
	  h1_EtaHeepTrackIsoNewEB_->Fill(fabs(Ele_etaSC->at(i)));
	  h1_PVHeepTrackIsoNewEB_->Fill(Ele_nbVtx->at(i));
	}
    }
    
    else if( fabs(Ele_etaSC->at(i)) > 1.560 && fabs(Ele_etaSC->at(i)) < 2.5 &&
	     Ele_Et->at(i) > Etcut && Ele_Et->at(i) < 50.0 && Ele_isEcalDrivenSeed->at(i) == 1 &&
	     Ele_nbOfMissingHits->at(i) < 2 &&
	     fabs( Ele_deltaEtaInSC->at(i) ) < std::max(0.015-8.5*pow(10,-5)*Ele_Et->at(i),0.006) && 
	     fabs( Ele_deltaPhiInSC->at(i) ) < 0.06 &&
	     Ele_hadronicOverEm->at(i)  < 12.5/Ele_energySC->at(i)+0.05 &&
	     Ele_sigmaIetaIeta->at(i) < 0.03 &&
	     (Ele_EcalPlusHcald1iso->at(i)) < (2.5+(Ele_rhoIso->at(i)*0.28)) &&
	     //Ele_dr03TkSumPt->at(i) < 5.0 &&
	     fabs(Ele_dxy->at(i)) < 0.05 ) {
      if(Ele_Et->at(i)>minPt && Ele_Et->at(i)<maxPt)
	{
	  h1_PtHeepTrackIsoNewEE_->Fill(Ele_Et->at(i));
	  h1_EtaHeepTrackIsoNewEE_->Fill(fabs(Ele_etaSC->at(i)));
	  h1_PVHeepTrackIsoNewEE_->Fill(Ele_nbVtx->at(i));
	}
    }
    else if( fabs(Ele_etaSC->at(i)) > 1.560 && fabs(Ele_etaSC->at(i)) < 2.5 && 
	     Ele_Et->at(i) > 50.0 && Ele_isEcalDrivenSeed->at(i) == 1 &&
	     Ele_nbOfMissingHits->at(i) < 2 &&
	     fabs( Ele_deltaEtaInSC->at(i) ) < std::max(0.015-8.5*pow(10,-5)*Ele_Et->at(i),0.006) && 
	     fabs( Ele_deltaPhiInSC->at(i) ) < 0.06 &&
	     Ele_hadronicOverEm->at(i) < 12.5/Ele_energySC->at(i)+0.05 &&
	     Ele_sigmaIetaIeta->at(i) < 0.03 &&
	     (Ele_EcalPlusHcald1iso->at(i)) < (2.5+0.03*(Ele_Et->at(i)-50.0)+(Ele_rhoIso->at(i)*0.28)) &&
	     //Ele_dr03TkSumPt->at(i) < 5.0 &&
	     fabs(Ele_dxy->at(i)) < 0.05 ) {
      if(Ele_Et->at(i)>minPt && Ele_Et->at(i)<maxPt)
	{
	  h1_PtHeepTrackIsoNewEE_->Fill(Ele_Et->at(i));
	  h1_EtaHeepTrackIsoNewEE_->Fill(fabs(Ele_etaSC->at(i)));
	  h1_PVHeepTrackIsoNewEE_->Fill(Ele_nbVtx->at(i));
	}
    }
    else continue;
  }
}


void ZprimeToEE::SelectHEEPMissHitsnewcut(){
  for(unsigned i=0; i<Ele_nbElectrons->size(); i++){
    if( fabs(Ele_etaSC->at(i)) < 1.442 && 
	Ele_Et->at(i) > Etcut && Ele_isEcalDrivenSeed->at(i) == 1 &&
	//Ele_nbOfMissingHits->at(i) < 2 &&
	fabs( Ele_deltaEtaInSC->at(i) ) < std::max(0.016-1*pow(10,-4)*Ele_Et->at(i),0.004) && 
	fabs( Ele_deltaPhiInSC->at(i) ) < 0.06 &&
	Ele_hadronicOverEm->at(i) < 2.0/Ele_energySC->at(i)+0.05 &&
	(Ele_e2x5MaxOver5x5Full5x5->at(i) > 0.94 || Ele_e1x5Over5x5Full5x5->at(i) > 0.83) &&
	(Ele_EcalPlusHcald1iso->at(i)) < (2.0+0.03*(Ele_Et->at(i))+(Ele_rhoIso->at(i)*0.28)) &&
	Ele_dr03TkSumPt->at(i) < 5.0 &&
	fabs(Ele_dxy->at(i)) < 0.02 ) {
      if(Ele_Et->at(i)>minPt && Ele_Et->at(i)<maxPt)
	{
	  h1_PtHeepMissHitsNewEB_->Fill(Ele_Et->at(i));
	  h1_EtaHeepMissHitsNewEB_->Fill(fabs(Ele_etaSC->at(i)));
	  h1_PVHeepMissHitsNewEB_->Fill(Ele_nbVtx->at(i));
	}
    }
    
    else if( fabs(Ele_etaSC->at(i)) > 1.560 && fabs(Ele_etaSC->at(i)) < 2.5 &&
	     Ele_Et->at(i) > Etcut && Ele_Et->at(i) < 50.0 && Ele_isEcalDrivenSeed->at(i) == 1 &&
	     //Ele_nbOfMissingHits->at(i) < 2 &&
	     fabs( Ele_deltaEtaInSC->at(i) ) < std::max(0.015-8.5*pow(10,-5)*Ele_Et->at(i),0.006) && 
	     fabs( Ele_deltaPhiInSC->at(i) ) < 0.06 &&
	     Ele_hadronicOverEm->at(i)  < 12.5/Ele_energySC->at(i)+0.05 &&
	     Ele_sigmaIetaIeta->at(i) < 0.03 &&
	     (Ele_EcalPlusHcald1iso->at(i)) < (2.5+(Ele_rhoIso->at(i)*0.28)) &&
	     Ele_dr03TkSumPt->at(i) < 5.0 &&
	     fabs(Ele_dxy->at(i)) < 0.05 ) {
      if(Ele_Et->at(i)>minPt && Ele_Et->at(i)<maxPt)
	{
	  h1_PtHeepMissHitsNewEE_->Fill(Ele_Et->at(i));
	  h1_EtaHeepMissHitsNewEE_->Fill(fabs(Ele_etaSC->at(i)));
	  h1_PVHeepMissHitsNewEE_->Fill(Ele_nbVtx->at(i));
	}
    }
    else if( fabs(Ele_etaSC->at(i)) > 1.560 && fabs(Ele_etaSC->at(i)) < 2.5 && 
	     Ele_Et->at(i) > 50.0 && Ele_isEcalDrivenSeed->at(i) == 1 &&
	     //Ele_nbOfMissingHits->at(i) < 2 &&
	     fabs( Ele_deltaEtaInSC->at(i) ) < std::max(0.015-8.5*pow(10,-5)*Ele_Et->at(i),0.006) && 
	     fabs( Ele_deltaPhiInSC->at(i) ) < 0.06 &&
	     Ele_hadronicOverEm->at(i) < 12.5/Ele_energySC->at(i)+0.05 &&
	     Ele_sigmaIetaIeta->at(i) < 0.03 &&
	     (Ele_EcalPlusHcald1iso->at(i)) < (2.5+0.03*(Ele_Et->at(i)-50.0)+(Ele_rhoIso->at(i)*0.28)) &&
	     Ele_dr03TkSumPt->at(i) < 5.0 &&
	     fabs(Ele_dxy->at(i)) < 0.05 ) {
      if(Ele_Et->at(i)>minPt && Ele_Et->at(i)<maxPt)
	{
	  h1_PtHeepMissHitsNewEE_->Fill(Ele_Et->at(i));
	  h1_EtaHeepMissHitsNewEE_->Fill(fabs(Ele_etaSC->at(i)));
	  h1_PVHeepMissHitsNewEE_->Fill(Ele_nbVtx->at(i));
	}
    }
    else continue;
  }
}

void ZprimeToEE::SelectHEEPSigmaIetaIetanewcut(){
  for(unsigned i=0; i<Ele_nbElectrons->size(); i++){
    if( fabs(Ele_etaSC->at(i)) > 1.560 && fabs(Ele_etaSC->at(i)) < 2.5 &&
	Ele_Et->at(i) > Etcut && Ele_Et->at(i) < 50.0 && Ele_isEcalDrivenSeed->at(i) == 1 &&
	Ele_nbOfMissingHits->at(i) < 2 &&
	fabs( Ele_deltaEtaInSC->at(i) ) < std::max(0.015-8.5*pow(10,-5)*Ele_Et->at(i),0.006) && 
	fabs( Ele_deltaPhiInSC->at(i) ) < 0.06 &&
	Ele_hadronicOverEm->at(i)  < 12.5/Ele_energySC->at(i)+0.05 &&
	//Ele_sigmaIetaIeta->at(i) < 0.03 &&
	(Ele_EcalPlusHcald1iso->at(i)) < (2.5+(Ele_rhoIso->at(i)*0.28)) &&
	Ele_dr03TkSumPt->at(i) < 5.0 &&
	fabs(Ele_dxy->at(i)) < 0.05 ) {
      if(Ele_Et->at(i)>minPt && Ele_Et->at(i)<maxPt)
	{
	  h1_PtHeepSigmaNewEE_->Fill(Ele_Et->at(i));
	  h1_EtaHeepSigmaNewEE_->Fill(fabs(Ele_etaSC->at(i)));
	  h1_PVHeepSigmaNewEE_->Fill(Ele_nbVtx->at(i));
	}
    }
    else if( fabs(Ele_etaSC->at(i)) > 1.560 && fabs(Ele_etaSC->at(i)) < 2.5 && 
	     Ele_Et->at(i) > 50.0 && Ele_isEcalDrivenSeed->at(i) == 1 &&
	     Ele_nbOfMissingHits->at(i) < 2 &&
	     fabs( Ele_deltaEtaInSC->at(i) ) < std::max(0.015-8.5*pow(10,-5)*Ele_Et->at(i),0.006) && 
	     fabs( Ele_deltaPhiInSC->at(i) ) < 0.06 &&
	     Ele_hadronicOverEm->at(i) < 12.5/Ele_energySC->at(i)+0.05 &&
	     //Ele_sigmaIetaIeta->at(i) < 0.03 &&
	     (Ele_EcalPlusHcald1iso->at(i)) < (2.5+0.03*(Ele_Et->at(i)-50.0)+(Ele_rhoIso->at(i)*0.28)) &&
	     Ele_dr03TkSumPt->at(i) < 5.0 &&
	     fabs(Ele_dxy->at(i)) < 0.05 ) {
      if(Ele_Et->at(i)>minPt && Ele_Et->at(i)<maxPt)
	{
	  h1_PtHeepSigmaNewEE_->Fill(Ele_Et->at(i));
	  h1_EtaHeepSigmaNewEE_->Fill(fabs(Ele_etaSC->at(i)));
	  h1_PVHeepSigmaNewEE_->Fill(Ele_nbVtx->at(i));
	}
    }
    else continue;
  }
}



void ZprimeToEE::SelectHEEPEcalHcalIsonewcut(){
  for(unsigned i=0; i<Ele_nbElectrons->size(); i++){
    if( fabs(Ele_etaSC->at(i)) < 1.442 && 
	Ele_Et->at(i) > Etcut && Ele_isEcalDrivenSeed->at(i) == 1 &&
	Ele_nbOfMissingHits->at(i) < 2 &&
	fabs( Ele_deltaEtaInSC->at(i) ) < std::max(0.016-1*pow(10,-4)*Ele_Et->at(i),0.004) && 
	fabs( Ele_deltaPhiInSC->at(i) ) < 0.06 &&
	Ele_hadronicOverEm->at(i) < 2.0/Ele_energySC->at(i)+0.05 &&
	(Ele_e2x5MaxOver5x5Full5x5->at(i) > 0.94 || Ele_e1x5Over5x5Full5x5->at(i) > 0.83) &&
	//(Ele_EcalPlusHcald1iso->at(i)) < (2.0+0.03*(Ele_Et->at(i))+(Ele_rhoIso->at(i)*0.28)) &&
	Ele_dr03TkSumPt->at(i) < 5.0 &&
	fabs(Ele_dxy->at(i)) < 0.02 ) {
      if(Ele_Et->at(i)>minPt && Ele_Et->at(i)<maxPt)
	{
	  h1_PtHeepEcalHcalIsoNewEB_->Fill(Ele_Et->at(i));
	  h1_EtaHeepEcalHcalIsoNewEB_->Fill(fabs(Ele_etaSC->at(i)));
	  h1_PVHeepEcalHcalIsoNewEB_->Fill(Ele_nbVtx->at(i));
	}
    }
    
    else if( fabs(Ele_etaSC->at(i)) > 1.560 && fabs(Ele_etaSC->at(i)) < 2.5 &&
	     Ele_Et->at(i) > Etcut && Ele_Et->at(i) < 50.0 && Ele_isEcalDrivenSeed->at(i) == 1 &&
	     Ele_nbOfMissingHits->at(i) < 2 &&
	     fabs( Ele_deltaEtaInSC->at(i) ) < std::max(0.015-8.5*pow(10,-5)*Ele_Et->at(i),0.006) && 
	     fabs( Ele_deltaPhiInSC->at(i) ) < 0.06 &&
	     Ele_hadronicOverEm->at(i)  < 12.5/Ele_energySC->at(i)+0.05 &&
	     Ele_sigmaIetaIeta->at(i) < 0.03 &&
	     //(Ele_EcalPlusHcald1iso->at(i)) < (2.5+(Ele_rhoIso->at(i)*0.28)) &&
	     Ele_dr03TkSumPt->at(i) < 5.0 &&
	     fabs(Ele_dxy->at(i)) < 0.05 ) {
      if(Ele_Et->at(i)>minPt && Ele_Et->at(i)<maxPt)
	{
	  h1_PtHeepEcalHcalIsoNewEE_->Fill(Ele_Et->at(i));
	  h1_EtaHeepEcalHcalIsoNewEE_->Fill(fabs(Ele_etaSC->at(i)));
	  h1_PVHeepEcalHcalIsoNewEE_->Fill(Ele_nbVtx->at(i));
	}
    }
    else if( fabs(Ele_etaSC->at(i)) > 1.560 && fabs(Ele_etaSC->at(i)) < 2.5 && 
	     Ele_Et->at(i) > 50.0 && Ele_isEcalDrivenSeed->at(i) == 1 &&
	     Ele_nbOfMissingHits->at(i) < 2 &&
	     fabs( Ele_deltaEtaInSC->at(i) ) < std::max(0.015-8.5*pow(10,-5)*Ele_Et->at(i),0.006) && 
	     fabs( Ele_deltaPhiInSC->at(i) ) < 0.06 &&
	     Ele_hadronicOverEm->at(i) < 12.5/Ele_energySC->at(i)+0.05 &&
	     Ele_sigmaIetaIeta->at(i) < 0.03 &&
	     //(Ele_EcalPlusHcald1iso->at(i)) < (2.5+0.03*(Ele_Et->at(i)-50.0)+(Ele_rhoIso->at(i)*0.28)) &&
	     Ele_dr03TkSumPt->at(i) < 5.0 &&
	     fabs(Ele_dxy->at(i)) < 0.05 ) {
      if(Ele_Et->at(i)>minPt && Ele_Et->at(i)<maxPt)
	{
	  h1_PtHeepEcalHcalIsoNewEE_->Fill(Ele_Et->at(i));
	  h1_EtaHeepEcalHcalIsoNewEE_->Fill(fabs(Ele_etaSC->at(i)));
	  h1_PVHeepEcalHcalIsoNewEE_->Fill(Ele_nbVtx->at(i));
	}
    }
    else continue;
  }
}

void ZprimeToEE::SelectHEEPShowerShapenewcut(){
  for(unsigned i=0; i<Ele_nbElectrons->size(); i++){
    if( fabs(Ele_etaSC->at(i)) < 1.442 && 
	Ele_Et->at(i) > Etcut && Ele_isEcalDrivenSeed->at(i) == 1 &&
	Ele_nbOfMissingHits->at(i) < 2 &&
	fabs( Ele_deltaEtaInSC->at(i) ) < std::max(0.016-1*pow(10,-4)*Ele_Et->at(i),0.004) && 
	fabs( Ele_deltaPhiInSC->at(i) ) < 0.06 &&
	Ele_hadronicOverEm->at(i) < 2.0/Ele_energySC->at(i)+0.05 &&
	//(Ele_e2x5MaxOver5x5Full5x5->at(i) > 0.94 || Ele_e1x5Over5x5Full5x5->at(i) > 0.83) &&
	(Ele_EcalPlusHcald1iso->at(i)) < (2.0+0.03*(Ele_Et->at(i))+(Ele_rhoIso->at(i)*0.28)) &&
	Ele_dr03TkSumPt->at(i) < 5.0 &&
	fabs(Ele_dxy->at(i)) < 0.02 ) {
      if(Ele_Et->at(i)>minPt && Ele_Et->at(i)<maxPt)
	{
	  h1_PtHeepShowerShapeNewEB_->Fill(Ele_Et->at(i));
	  h1_EtaHeepShowerShapeNewEB_->Fill(fabs(Ele_etaSC->at(i)));
	  h1_PVHeepShowerShapeNewEB_->Fill(Ele_nbVtx->at(i));
	}
    }
    else continue;
  }
}





//===================================================================================
//
//
//                       Part for HEEP selection version 4.1
//
//
//============================ Method to select old heep id  ========================
void ZprimeToEE::SelectHEEPelectronsV4(){
  for(unsigned i=0; i<Ele_nbElectrons->size(); i++){
    if( fabs(Ele_etaSC->at(i)) < 1.442 && 
	Ele_Et->at(i) > Etcut && Ele_isEcalDrivenSeed->at(i) == 1 &&
	Ele_nbOfMissingHits->at(i) < 2 &&
	fabs( Ele_deltaEtaInSC->at(i) ) < 0.005 &&
	fabs( Ele_deltaPhiInSC->at(i) ) < 0.06 &&
	Ele_hadronicOverEm->at(i) < 0.05 &&
	(Ele_e2x5MaxOver5x5Full5x5->at(i) > 0.94 || Ele_e1x5Over5x5Full5x5->at(i) > 0.83) &&
	(Ele_EcalPlusHcald1iso->at(i)) < (2.0+0.03*(Ele_Et->at(i))+(Ele_rhoIso->at(i)*0.28)) &&
	Ele_dr03TkSumPt->at(i) < 5.0 &&
	fabs(Ele_dxy->at(i)) < 0.02 ) {
      if(Ele_Et->at(i)>minPt && Ele_Et->at(i)<maxPt)
	{
	  h1_PtHeepEleV4EB_->Fill(Ele_Et->at(i));
	  h1_EtaHeepEleV4EB_->Fill(fabs(Ele_etaSC->at(i)));
	  h1_PVHeepEleV4EB_->Fill(Ele_nbVtx->at(i));
	}
    }
    else if( fabs(Ele_etaSC->at(i)) > 1.560 && fabs(Ele_etaSC->at(i)) < 2.5 &&
	     Ele_Et->at(i) > Etcut && Ele_Et->at(i) < 50.0 && Ele_isEcalDrivenSeed->at(i) == 1 &&
	     Ele_nbOfMissingHits->at(i) < 2 &&
	     fabs( Ele_deltaEtaInSC->at(i) ) < 0.007 &&
	     fabs( Ele_deltaPhiInSC->at(i) ) < 0.06 &&
	     Ele_hadronicOverEm->at(i)  < 0.05 &&
	     Ele_sigmaIetaIeta->at(i) < 0.03 &&
	     (Ele_EcalPlusHcald1iso->at(i)) < (2.5+(Ele_rhoIso->at(i)*0.28)) &&
	     Ele_dr03TkSumPt->at(i) < 5.0 &&
	     fabs(Ele_dxy->at(i)) < 0.05 ) {
      if(Ele_Et->at(i)>minPt && Ele_Et->at(i)<maxPt)
	{
	  h1_PtHeepEleV4EE_->Fill(Ele_Et->at(i));
	  h1_EtaHeepEleV4EE_->Fill(fabs(Ele_etaSC->at(i)));
	  h1_PVHeepEleV4EE_->Fill(Ele_nbVtx->at(i));
	}
    }
    else if( fabs(Ele_etaSC->at(i)) > 1.560 && fabs(Ele_etaSC->at(i)) < 2.5 && 
	     Ele_Et->at(i) > 50.0 && Ele_isEcalDrivenSeed->at(i) == 1 &&
	     Ele_nbOfMissingHits->at(i) < 2 &&
	     fabs( Ele_deltaEtaInSC->at(i) ) < 0.007 &&
	     fabs( Ele_deltaPhiInSC->at(i) ) < 0.06 &&
	     Ele_hadronicOverEm->at(i) < 0.05 &&
	     Ele_sigmaIetaIeta->at(i) < 0.03 &&
	     (Ele_EcalPlusHcald1iso->at(i)) < (2.5+0.03*(Ele_Et->at(i)-50.0)+(Ele_rhoIso->at(i)*0.28)) &&
	     Ele_dr03TkSumPt->at(i) < 5.0 &&
	     fabs(Ele_dxy->at(i)) < 0.05 ) {
      if(Ele_Et->at(i)>minPt && Ele_Et->at(i)<maxPt)
	{
	  h1_PtHeepEleV4EE_->Fill(Ele_Et->at(i));
	  h1_EtaHeepEleV4EE_->Fill(fabs(Ele_etaSC->at(i)));
	  h1_PVHeepEleV4EE_->Fill(Ele_nbVtx->at(i));
	}
    }   
    else continue;
  }
}

void ZprimeToEE::SelectHEEPHoEV4cut(){
  for(unsigned i=0; i<Ele_nbElectrons->size(); i++){
    if( fabs(Ele_etaSC->at(i)) < 1.442 && 
	Ele_Et->at(i) > Etcut && Ele_isEcalDrivenSeed->at(i) == 1 &&
	Ele_nbOfMissingHits->at(i) < 2 &&
	fabs( Ele_deltaEtaInSC->at(i) ) < 0.005 &&
	fabs( Ele_deltaPhiInSC->at(i) ) < 0.06 &&
	//Ele_hadronicOverEm->at(i) < 0.05 &&
	(Ele_e2x5MaxOver5x5Full5x5->at(i) > 0.94 || Ele_e1x5Over5x5Full5x5->at(i) > 0.83) &&
	(Ele_EcalPlusHcald1iso->at(i)) < (2.0+0.03*(Ele_Et->at(i))+(Ele_rhoIso->at(i)*0.28)) &&
	Ele_dr03TkSumPt->at(i) < 5.0 &&
	fabs(Ele_dxy->at(i)) < 0.02 ) {
      if(Ele_Et->at(i)>minPt && Ele_Et->at(i)<maxPt)
	{
	  h1_PtHeepHoeV4EB_->Fill(Ele_Et->at(i));
	  h1_EtaHeepHoeV4EB_->Fill(fabs(Ele_etaSC->at(i)));
	  h1_PVHeepHoeV4EB_->Fill(Ele_nbVtx->at(i));
	}
    }
    
    else if( fabs(Ele_etaSC->at(i)) > 1.560 && fabs(Ele_etaSC->at(i)) < 2.5 &&
	     Ele_Et->at(i) > Etcut && Ele_Et->at(i) < 50.0 && Ele_isEcalDrivenSeed->at(i) == 1 &&
	     Ele_nbOfMissingHits->at(i) < 2 &&
	     fabs( Ele_deltaEtaInSC->at(i) ) < 0.007 &&
	     fabs( Ele_deltaPhiInSC->at(i) ) < 0.06 &&
	     //Ele_hadronicOverEm->at(i)  < 0.05 &&
	     Ele_sigmaIetaIeta->at(i) < 0.03 &&
	     (Ele_EcalPlusHcald1iso->at(i)) < (2.5+(Ele_rhoIso->at(i)*0.28)) &&
	     Ele_dr03TkSumPt->at(i) < 5.0 &&
	     fabs(Ele_dxy->at(i)) < 0.05 ) {
      if(Ele_Et->at(i)>minPt && Ele_Et->at(i)<maxPt)
	{
	  h1_PtHeepHoeV4EE_->Fill(Ele_Et->at(i));
	  h1_EtaHeepHoeV4EE_->Fill(fabs(Ele_etaSC->at(i)));
	  h1_PVHeepHoeV4EE_->Fill(Ele_nbVtx->at(i));
	}
    }
    else if( fabs(Ele_etaSC->at(i)) > 1.560 && fabs(Ele_etaSC->at(i)) < 2.5 && 
	     Ele_Et->at(i) > 50.0 && Ele_isEcalDrivenSeed->at(i) == 1 &&
	     Ele_nbOfMissingHits->at(i) < 2 &&
	     fabs( Ele_deltaEtaInSC->at(i) ) < 0.007 &&
	     fabs( Ele_deltaPhiInSC->at(i) ) < 0.06 &&
	     //Ele_hadronicOverEm->at(i) < 0.05 &&
	     Ele_sigmaIetaIeta->at(i) < 0.03 &&
	     (Ele_EcalPlusHcald1iso->at(i)) < (2.5+0.03*(Ele_Et->at(i)-50.0)+(Ele_rhoIso->at(i)*0.28)) &&
	     Ele_dr03TkSumPt->at(i) < 5.0 &&
	     fabs(Ele_dxy->at(i)) < 0.05 ) {
      if(Ele_Et->at(i)>minPt && Ele_Et->at(i)<maxPt)
	{
	  h1_PtHeepHoeV4EE_->Fill(Ele_Et->at(i));
	  h1_EtaHeepHoeV4EE_->Fill(fabs(Ele_etaSC->at(i)));
	  h1_PVHeepHoeV4EE_->Fill(Ele_nbVtx->at(i));
	}
    }
    else continue;
  }
}

void ZprimeToEE::SelectHEEPDeltaEtaInV4cut(){
  for(unsigned i=0; i<Ele_nbElectrons->size(); i++){
    if( fabs(Ele_etaSC->at(i)) < 1.442 && 
	Ele_Et->at(i) > Etcut && Ele_isEcalDrivenSeed->at(i) == 1 &&
	Ele_nbOfMissingHits->at(i) < 2 &&
	//fabs( Ele_deltaEtaInSC->at(i) ) < 0.005 &&
	fabs( Ele_deltaPhiInSC->at(i) ) < 0.06 &&
	Ele_hadronicOverEm->at(i) < 0.05 &&
	(Ele_e2x5MaxOver5x5Full5x5->at(i) > 0.94 || Ele_e1x5Over5x5Full5x5->at(i) > 0.83) &&
	(Ele_EcalPlusHcald1iso->at(i)) < (2.0+0.03*(Ele_Et->at(i))+(Ele_rhoIso->at(i)*0.28)) &&
	Ele_dr03TkSumPt->at(i) < 5.0 &&
	fabs(Ele_dxy->at(i)) < 0.02 ) {
      if(Ele_Et->at(i)>minPt && Ele_Et->at(i)<maxPt)
	{
	  h1_PtHeepDeltaEtaInV4EB_->Fill(Ele_Et->at(i));
	  h1_EtaHeepDeltaEtaInV4EB_->Fill(fabs(Ele_etaSC->at(i)));
	  h1_PVHeepDeltaEtaInV4EB_->Fill(Ele_nbVtx->at(i));
	}
    }
    
    else if( fabs(Ele_etaSC->at(i)) > 1.560 && fabs(Ele_etaSC->at(i)) < 2.5 &&
	     Ele_Et->at(i) > Etcut && Ele_Et->at(i) < 50.0 && Ele_isEcalDrivenSeed->at(i) == 1 &&
	     Ele_nbOfMissingHits->at(i) < 2 &&
	     //fabs( Ele_deltaEtaInSC->at(i) ) < 0.007 &&
	     fabs( Ele_deltaPhiInSC->at(i) ) < 0.06 &&
	     Ele_hadronicOverEm->at(i)  < 0.05 &&
	     Ele_sigmaIetaIeta->at(i) < 0.03 &&
	     (Ele_EcalPlusHcald1iso->at(i)) < (2.5+(Ele_rhoIso->at(i)*0.28)) &&
	     Ele_dr03TkSumPt->at(i) < 5.0 &&
	     fabs(Ele_dxy->at(i)) < 0.05 ) {
      if(Ele_Et->at(i)>minPt && Ele_Et->at(i)<maxPt)
	{
	  h1_PtHeepDeltaEtaInV4EE_->Fill(Ele_Et->at(i));
	  h1_EtaHeepDeltaEtaInV4EE_->Fill(fabs(Ele_etaSC->at(i)));
	  h1_PVHeepDeltaEtaInV4EE_->Fill(Ele_nbVtx->at(i));
	}
    }
    else if( fabs(Ele_etaSC->at(i)) > 1.560 && fabs(Ele_etaSC->at(i)) < 2.5 && 
	     Ele_Et->at(i) > 50.0 && Ele_isEcalDrivenSeed->at(i) == 1 &&
	     Ele_nbOfMissingHits->at(i) < 2 &&
	     //fabs( Ele_deltaEtaInSC->at(i) ) < 0.007 &&
	     fabs( Ele_deltaPhiInSC->at(i) ) < 0.06 &&
	     Ele_hadronicOverEm->at(i) < 0.05 &&
	     Ele_sigmaIetaIeta->at(i) < 0.03 &&
	     (Ele_EcalPlusHcald1iso->at(i)) < (2.5+0.03*(Ele_Et->at(i)-50.0)+(Ele_rhoIso->at(i)*0.28)) &&
	     Ele_dr03TkSumPt->at(i) < 5.0 &&
	     fabs(Ele_dxy->at(i)) < 0.05 ) {
      if(Ele_Et->at(i)>minPt && Ele_Et->at(i)<maxPt)
	{
	  h1_PtHeepDeltaEtaInV4EE_->Fill(Ele_Et->at(i));
	  h1_EtaHeepDeltaEtaInV4EE_->Fill(fabs(Ele_etaSC->at(i)));
	  h1_PVHeepDeltaEtaInV4EE_->Fill(Ele_nbVtx->at(i));
	}
    }
    else continue;
  }
}

void ZprimeToEE::SelectHEEPV5Dxtcut(){
  for(unsigned i=0; i<Ele_nbElectrons->size(); i++){
    if( fabs(Ele_etaSC->at(i)) < 1.442 && 
	Ele_Et->at(i) > Etcut && Ele_isEcalDrivenSeed->at(i) == 1 &&
	Ele_nbOfMissingHits->at(i) < 2 &&
	fabs( Ele_deltaEtaInSC->at(i) ) < std::max(0.016-1*pow(10,-4)*Ele_Et->at(i),0.004) && 
	fabs( Ele_deltaPhiInSC->at(i) ) < 0.06 &&
	Ele_hadronicOverEm->at(i) < 2.0/Ele_energySC->at(i)+0.05 &&
	(Ele_e2x5MaxOver5x5Full5x5->at(i) > 0.94 || Ele_e1x5Over5x5Full5x5->at(i) > 0.83) &&
	(Ele_EcalPlusHcald1iso->at(i)) < (2.0+0.03*(Ele_Et->at(i))+(Ele_rhoIso->at(i)*0.28)) &&
	//fabs(Ele_dxy->at(i)) < 0.02 &&
	Ele_dr03TkSumPt->at(i) < 5.0 ) {
      if(Ele_Et->at(i)>minPt && Ele_Et->at(i)<maxPt)
	{
	  h1_PtHeepDxyNewEB_->Fill(Ele_Et->at(i));
	  h1_EtaHeepDxyNewEB_->Fill(fabs(Ele_etaSC->at(i)));
	  h1_PVHeepDxyNewEB_->Fill(Ele_nbVtx->at(i));
	  h1_DxyNewEB_->Fill(Ele_dxy->at(i)); //*Ele_charge->at(i));
	}
    }
    
    else if( fabs(Ele_etaSC->at(i)) > 1.560 && fabs(Ele_etaSC->at(i)) < 2.5 &&
	     Ele_Et->at(i) > Etcut && Ele_Et->at(i) < 50.0 && Ele_isEcalDrivenSeed->at(i) == 1 &&
	     Ele_nbOfMissingHits->at(i) < 2 &&
	     fabs( Ele_deltaEtaInSC->at(i) ) < std::max(0.015-8.5*pow(10,-5)*Ele_Et->at(i),0.006) && 
	     fabs( Ele_deltaPhiInSC->at(i) ) < 0.06 &&
	     Ele_hadronicOverEm->at(i)  < 12.5/Ele_energySC->at(i)+0.05 &&
	     Ele_sigmaIetaIeta->at(i) < 0.03 &&
	     (Ele_EcalPlusHcald1iso->at(i)) < (2.5+(Ele_rhoIso->at(i)*0.28)) &&
	     //fabs(Ele_dxy->at(i)) < 0.05  &&
	     Ele_dr03TkSumPt->at(i) < 5.0 ) {
      if(Ele_Et->at(i)>minPt && Ele_Et->at(i)<maxPt)
	{
	  h1_PtHeepDxyNewEE_->Fill(Ele_Et->at(i));
	  h1_EtaHeepDxyNewEE_->Fill(fabs(Ele_etaSC->at(i)));
	  h1_PVHeepDxyNewEE_->Fill(Ele_nbVtx->at(i));
	  h1_DxyNewEE_->Fill(Ele_dxy->at(i)); //*Ele_charge->at(i));
	}
    }
    else if( fabs(Ele_etaSC->at(i)) > 1.560 && fabs(Ele_etaSC->at(i)) < 2.5 && 
	     Ele_Et->at(i) > 50.0 && Ele_isEcalDrivenSeed->at(i) == 1 &&
	     Ele_nbOfMissingHits->at(i) < 2 &&
	     fabs( Ele_deltaEtaInSC->at(i) ) < std::max(0.015-8.5*pow(10,-5)*Ele_Et->at(i),0.006) && 
	     fabs( Ele_deltaPhiInSC->at(i) ) < 0.06 &&
	     Ele_hadronicOverEm->at(i) < 12.5/Ele_energySC->at(i)+0.05 &&
	     Ele_sigmaIetaIeta->at(i) < 0.03 &&
	     (Ele_EcalPlusHcald1iso->at(i)) < (2.5+0.03*(Ele_Et->at(i)-50.0)+(Ele_rhoIso->at(i)*0.28)) &&
	     //fabs(Ele_dxy->at(i)) < 0.05  &&
	     Ele_dr03TkSumPt->at(i) < 5.0 ) {
      if(Ele_Et->at(i)>minPt && Ele_Et->at(i)<maxPt)
	{
	  h1_PtHeepDxyNewEE_->Fill(Ele_Et->at(i));
	  h1_EtaHeepDxyNewEE_->Fill(fabs(Ele_etaSC->at(i)));
	  h1_PVHeepDxyNewEE_->Fill(Ele_nbVtx->at(i));
	  h1_DxyNewEE_->Fill(Ele_dxy->at(i)); //*Ele_charge->at(i));
	}
    }
    else continue;
  }
}
