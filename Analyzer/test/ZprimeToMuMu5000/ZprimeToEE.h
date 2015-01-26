//=========================================================================  
//          Analysis code for Z' boson to Mu Mu analysis                  =  
//    [1] In this code we select the high pt di-muons events              =
//    [2] It also computes Acceptance X Efficiency (Accxeff)              =                
//    [3] It computes the high pt Muon Eff. using Tag & Probe             =
//             Written by Sherif Elgammal                                 =
//                                                                        =
//                    07/08/2014                                          =
//=========================================================================
//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sat Jun 28 15:22:58 2014 by ROOT version 5.30/04
// from TTree tree/tr
// found on file: CMSSW701_MC_ZprimeMuMu_13TeV_tree.root
//////////////////////////////////////////////////////////
#ifndef ZprimeToEE_h
#define ZprimeToEE_h
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
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
// Header file for the classes stored in the TTree if any.
#include <vector>

// Fixed size dimensions of array or collections stored in the TTree if any.

class ZprimeToEE {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   vector<int>     *NbEventsPassTrigger;
   vector<int>     *NbEventsPassTriggerandPVcond;
   Int_t           event_runNo;
   Int_t           event_evtNo;
   Int_t           event_lumi;
   Int_t           event_bunch;
   vector<int>     *Ele_nbElectrons;
   vector<bool>    *Ele_isEcalDrivenSeed;
   vector<bool>    *Ele_isPassConversionVeto;
   vector<int>     *Ele_charge;
   vector<int>     *Ele_nbOfMissingHits;
   vector<int>     *Ele_nbVtx;
   vector<float>   *Ele_Et;
   vector<float>   *Ele_pt;
   vector<float>   *Ele_thetaSC;
   vector<float>   *Ele_etaSC;
   vector<float>   *Ele_phiSC;
   vector<float>   *Ele_energySC;
   vector<float>   *Ele_preshowerEnergySC;
   vector<float>   *Ele_thetaTrack;
   vector<float>   *Ele_etaTrack;
   vector<float>   *Ele_phiTrack;
   vector<float>   *Ele_hadronicOverEm;
   vector<float>   *Ele_deltaEtaInSC;
   vector<float>   *Ele_deltaPhiInSC;
   vector<float>   *Ele_sigmaIetaIeta;
   vector<float>   *Ele_sigmaIetaIetaFull5x5;
   vector<float>   *Ele_e2x5MaxFull5x5;
   vector<float>   *Ele_e1x5Full5x5;
   vector<float>   *Ele_e5x5Full5x5;
   vector<float>   *Ele_e2x5MaxOver5x5Full5x5;
   vector<float>   *Ele_e1x5Over5x5Full5x5;
   vector<float>   *Ele_e2x5Right;
   vector<float>   *Ele_e2x5Left;
   vector<float>   *Ele_e2x5Top;
   vector<float>   *Ele_e2x5Bottom;
   vector<float>   *Ele_eMax;
   vector<float>   *Ele_eRight;
   vector<float>   *Ele_eLeft;
   vector<float>   *Ele_eTop;
   vector<float>   *Ele_eBottom;
   vector<float>   *Ele_dxy;
   vector<float>   *Ele_dz;
   vector<float>   *Ele_rhoIso;
   vector<float>   *Ele_fbrem;
   vector<float>   *Ele_EoverP;
   vector<float>   *Ele_Xposition;
   vector<float>   *Ele_Yposition;
   vector<float>   *Ele_EcalPlusHcald1iso;
   vector<float>   *Ele_dr03TkSumPt;
   vector<float>   *Ele_dr03EcalRecHitSumEt;
   vector<float>   *Ele_dr03HcalDepth1TowerSumEt;
   vector<float>   *Ele_dr03HcalDepth1TowerSumEtBc;
   vector<float>   *Ele_hcalDepth1OverEcal;
   vector<float>   *Ele_hcalDepth2OverEcal;
   vector<float>   *Ele_dr03HcalDepth2TowerSumEt;
   vector<float>   *Ele_hcalDepth2TowerSumEtNoVeto;
   vector<float>   *Ele_hcalDepth1TowerSumEtNoVeto;
   vector<float>   *Ele_pfSumPhotonEt;
   vector<float>   *Ele_pfSumChargedHadronPt;
   vector<float>   *Ele_pfSumNeutralHadronEt;
   vector<float>   *Ele_pfSumPUPt;
   vector<float>   *Ele_pfDeltaBeta;
   vector<float>   *Ele_pfCandidateDeltaBeta;
   vector<float>   *Ele_pfCandidateChargedHadron;
   vector<float>   *Ele_pfCandidateNeutralHadron;
   vector<int>     *iGen;
   vector<int>     *idGen;
   vector<int>     *statusGen;
   vector<float>   *ptGen;
   vector<float>   *etaGen;
   vector<float>   *phiGen;
   vector<int>     *chargeGen;
   vector<float>   *EnergyGen;
   vector<float>   *pxGen;
   vector<float>   *pyGen;
   vector<float>   *pzGen;
   vector<int>     *nbPv;
   vector<int>     *Nbdof;
   vector<float>   *PositionRho;
   vector<float>   *PositionX;
   vector<float>   *PositionY;
   vector<float>   *PositionZ;
   vector<int>     *ecal_iRechit;
   vector<int>     *ecal_RecoFlag;
   vector<bool>    *ecal_isSaturated;
   vector<float>   *ecal_EnRecHit;
   vector<int>     *MuHLTObj_nbObj;
   vector<float>   *MuHLTObj_pt;
   vector<float>   *MuHLTObj_eta;
   vector<float>   *MuHLTObj_phi;

   // List of branches
   TBranch        *b_NbEventsPassTrigger;   //!
   TBranch        *b_NbEventsPassTriggerandPVcond;   //!
   TBranch        *b_event_runNo;   //!
   TBranch        *b_event_evtNo;   //!
   TBranch        *b_event_lumi;   //!
   TBranch        *b_event_bunch;   //!
   TBranch        *b_Ele_nbElectrons;   //!
   TBranch        *b_Ele_isEcalDrivenSeed;   //!
   TBranch        *b_Ele_isPassConversionVeto;   //!
   TBranch        *b_Ele_charge;   //!
   TBranch        *b_Ele_nbOfMissingHits;   //!
   TBranch        *b_Ele_nbVtx;   //!
   TBranch        *b_Ele_Et;   //!
   TBranch        *b_Ele_pt;   //!
   TBranch        *b_Ele_thetaSC;   //!
   TBranch        *b_Ele_etaSC;   //!
   TBranch        *b_Ele_phiSC;   //!
   TBranch        *b_Ele_energySC;   //!
   TBranch        *b_Ele_preshowerEnergySC;   //!
   TBranch        *b_Ele_thetaTrack;   //!
   TBranch        *b_Ele_etaTrack;   //!
   TBranch        *b_Ele_phiTrack;   //!
   TBranch        *b_Ele_hadronicOverEm;   //!
   TBranch        *b_Ele_deltaEtaInSC;   //!
   TBranch        *b_Ele_deltaPhiInSC;   //!
   TBranch        *b_Ele_sigmaIetaIeta;   //!
   TBranch        *b_Ele_sigmaIetaIetaFull5x5;   //!
   TBranch        *b_Ele_e2x5MaxFull5x5;   //!
   TBranch        *b_Ele_e1x5Full5x5;   //!
   TBranch        *b_Ele_e5x5Full5x5;   //!
   TBranch        *b_Ele_e2x5MaxOver5x5Full5x5;   //!
   TBranch        *b_Ele_e1x5Over5x5Full5x5;   //!
   TBranch        *b_Ele_e2x5Right;   //!
   TBranch        *b_Ele_e2x5Left;   //!
   TBranch        *b_Ele_e2x5Top;   //!
   TBranch        *b_Ele_e2x5Bottom;   //!
   TBranch        *b_Ele_eMax;   //!
   TBranch        *b_Ele_eRight;   //!
   TBranch        *b_Ele_eLeft;   //!
   TBranch        *b_Ele_eTop;   //!
   TBranch        *b_Ele_eBottom;   //!
   TBranch        *b_Ele_dxy;   //!
   TBranch        *b_Ele_dz;   //!
   TBranch        *b_Ele_rhoIso;   //!
   TBranch        *b_Ele_fbrem;   //!
   TBranch        *b_Ele_EoverP;   //!
   TBranch        *b_Ele_Xposition;   //!
   TBranch        *b_Ele_Yposition;   //!
   TBranch        *b_Ele_EcalPlusHcald1iso;   //!
   TBranch        *b_Ele_dr03TkSumPt;   //!
   TBranch        *b_Ele_dr03EcalRecHitSumEt;   //!
   TBranch        *b_Ele_dr03HcalDepth1TowerSumEt;   //!
   TBranch        *b_Ele_dr03HcalDepth1TowerSumEtBc;   //!
   TBranch        *b_Ele_hcalDepth1OverEcal;   //!
   TBranch        *b_Ele_hcalDepth2OverEcal;   //!
   TBranch        *b_Ele_dr03HcalDepth2TowerSumEt;   //!
   TBranch        *b_Ele_hcalDepth2TowerSumEtNoVeto;   //!
   TBranch        *b_Ele_hcalDepth1TowerSumEtNoVeto;   //!
   TBranch        *b_Ele_pfSumPhotonEt;   //!
   TBranch        *b_Ele_pfSumChargedHadronPt;   //!
   TBranch        *b_Ele_pfSumNeutralHadronEt;   //!
   TBranch        *b_Ele_pfSumPUPt;   //!
   TBranch        *b_Ele_pfDeltaBeta;   //!
   TBranch        *b_Ele_pfCandidateDeltaBeta;   //!
   TBranch        *b_Ele_pfCandidateChargedHadron;   //!
   TBranch        *b_Ele_pfCandidateNeutralHadron;   //!
   TBranch        *b_iGen;   //!
   TBranch        *b_idGen;   //!
   TBranch        *b_statusGen;   //!
   TBranch        *b_ptGen;   //!
   TBranch        *b_etaGen;   //!
   TBranch        *b_phiGen;   //!
   TBranch        *b_chargeGen;   //!
   TBranch        *b_EnergyGen;   //!
   TBranch        *b_pxGen;   //!
   TBranch        *b_pyGen;   //!
   TBranch        *b_pzGen;   //!
   TBranch        *b_nbPv;   //!
   TBranch        *b_Nbdof;   //!
   TBranch        *b_PositionRho;   //!
   TBranch        *b_PositionX;   //!
   TBranch        *b_PositionY;   //!
   TBranch        *b_PositionZ;   //!
   TBranch        *b_ecal_iRechit;   //!
   TBranch        *b_ecal_RecoFlag;   //!
   TBranch        *b_ecal_isSaturated;   //!
   TBranch        *b_ecal_EnRecHit;   //!
   TBranch        *b_MuHLTObj_nbObj;   //!
   TBranch        *b_MuHLTObj_pt;   //!
   TBranch        *b_MuHLTObj_eta;   //!
   TBranch        *b_MuHLTObj_phi;   //!

   ZprimeToEE(TTree *tree=0);
   virtual ~ZprimeToEE();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   bool SelectFirstMuon(float &pTmuon1,float &Enmuon1,float &Etamuon1,float &Phimuon1,
			unsigned &FlagMu1,float &EtaTrak1,float &PhiTrak1);
   bool SelectSecondMuon(unsigned FlagMu1,float pTmuon1,float &pTmuon2,float &Enmuon2,
			 float &Etamuon2,float &Phimuon2,float &EtaTrak2,float &PhiTrak2);
     
			 
   float Mass(float Pt1,float Eta1,float Phi1,float En1,
	      float Pt2,float Eta2,float Phi2,float En2);
   
   bool SelectFirstGenMu(float &ETMu1,float &PhiSCMu1,
			 float &EtaSCMu1,float &EnMu1,
			 int &IDele1,int &Statele1,
			 unsigned &GenFlag1);
   bool SelectSecondGenMu(unsigned GenFlag1,float ETMu1,float &ETMu2,float &PhiSCMu2,
			  float &EtaSCMu2,float &EnMu2,int &IDele2,int &Statele2);
   
   void PlotGenInfo(float ZprimeGenMass,float EtaGenMu1,float EtaGenMu2,float PtGenMu1,
		    float PtGenMu2,float EnGenMu1,float EnGenMu2,float PhiGenMu1,
		    float PhiGenMu2);   
   void PlotleadingAndSecondleading(float ptMuon1,float ptMuon2,float EnMuon1,float EnMuon2,
				    float EtaMuon1,float EtaMuon2,float PhiMuon1,float PhiMuon2,
				    float momentumMassMu);
   bool GenRecoMatchMu1(float RecoEta1,float RecoPhi1,
			float &ETMu1,float &PhiSCMu1,
			float &EtaSCMu1,float &EnMu1,
			unsigned &GenFlag1);
   bool GenRecoMatchMu2(unsigned GenFlag1,float RecoEta2,float RecoPhi2,
			float &ETMu2,float &PhiSCMu2,
			float &EtaSCMu2,float &EnMu2);
   
   bool RecoHLMuonTMatching(float RecoEta,float RecoPhi);
   void SelectHEEPelectrons();
   void SelectHEEPplusPFcut();
   void SelectHEEPplusPhotonConvcut();
   void SelectHEEPHoEnewcut();
   void SelectHEEPDeltaEtanewcut();
   void SelectHEEPDeltaPhinewcut();
   void SelectHEEPTrackIsonewcut();
   void SelectHEEPMissHitsnewcut();
   void SelectHEEPSigmaIetaIetanewcut();
   void SelectHEEPEcalHcalIsonewcut();
   void SelectHEEPShowerShapenewcut();



   void SelectHEEPelectronsV4();
   void SelectHEEPHoEV4cut();
   void SelectHEEPDeltaEtaInV4cut();
   void SelectHEEPV5Dxtcut();


   float RecoGenMatchingDeltaRcut;
   float RecoHLTMatchingDeltaRcut;
   float vtxChi2Mu,vtxMassMu;
   float mPtGen1,mPhiGen1,mEtaGen1,mEnGen1;
   unsigned mGenFlag1;
   float mPtGen2,mPhiGen2,mEtaGen2,mEnGen2;
   int ChargeRecMu1,ChargeRecMu2;
   unsigned flagmu1;
   unsigned flag1;
   float PtRecMu1,EnRecMu1,EtaRecMu1,PhiRecMu1;
   float PtRecMu2,EnRecMu2,EtaRecMu2,PhiRecMu2;
   float genET1,genPhi1,genEta1,genEn1;
   int genID1,genStat1;
   float genET2,genPhi2,genEta2,genEn2;
   int genID2,genStat2;
   float MassGen,RecoMass;
   int NbGen,NbReco;
   int nbTP,nbTT,nbTF;
   float TagProbeEtaCut;
   float Eff;
   float MassCutMin,MassCutMax;
   float MassResolution;
   float EtaCut;
   float minPt,maxPt,Etcut;
   TH1F* h1_ZprimeGenmass_;
   TH1F* h1_ZprimeGenEta1_;
   TH1F* h1_ZprimeGenEta2_;
   TH1F* h1_ZprimeGenPt1_;
   TH1F* h1_ZprimeGenPt2_;
   TH1F* h1_ZprimeGenEn1_;
   TH1F* h1_ZprimeGenEn2_;
   TH1F* h1_MassRecoGenDif_;
   TH1F* h1_ZprimeRecomass_;
   TH1F* h1_ZprimeRecoEta1_;
   TH1F* h1_ZprimeRecoEta2_;
   TH1F* h1_ZprimeRecoPt1_;
   TH1F* h1_ZprimeRecoPt2_;
   TH1F* h1_ZprimeRecoEn1_;
   TH1F* h1_ZprimeRecoEn2_;
   TH1F* h1_RecoPhiMu1DifPhiMu2_;
   TH1F* h1_GenPhiMu1DifPhiMu2_;
   TH1F* h1_invMassTP_;
   TH1F* h1_invMassTF_;
   TH1F* h1_invMassTT_;
   TH1F* h1_PtProbeTP_;
   TH1F* h1_EtaProbeTP_;
   TH1F* h1_PtProbeTT_;
   TH1F* h1_EtaProbeTT_;
   TH1F* h1_PtProbeTF_;
   TH1F* h1_EtaProbeTF_;
   TH1F* h1_dPToverPT_;
   TH1F* h1_normalizedChi2_;
   TH1F* h1_numberOftrackerLayersWithMeasurement_;
   TH1F* h1_numberOfValidPixelHits_;
   TH1F* h1_numberOfValidMuonHits_;
   TH1F* h1_numberOfMatchedStations_;
   TH1F* h1_trackiso_;
   TH1F* h1_absdxy_;
   TH1F* h1_absdz_;
   TH1F* h1_ZprimeRecoPhi1_;
   TH1F* h1_ZprimeRecoPhi2_;
   TH1F* h1_3Dangle_;
   TH1F* h1_MassTest_;
   TH1F* h1_NbHEEPMu_;
   TH1F* h1_LeadingRecoMuPt_;
   TH1F* h1_SecondLeadingMuPt_;
   TH1F* h1_deltaR_;
   //------------------------------
   // part for N-1 eff.
   //------------------------------
   TH1F* h1_PtIDpterror_;
   TH1F* h1_PtEffpterror_;
   TH1F* h1_PtEffptnumberOftrackerLayers_;

   TH1F* h1_PtEffptnumberOfPixelHits_;
   TH1F* h1_PtEffptnumberOfMuonHits_;
   TH1F* h1_PtEffptnumberOfMatchedStations_;
   TH1F* h1_PtEffptTrackIso_;
   TH1F* h1_PtEffptabsdsy_;
   TH1F* h1_pfSumChargedHadron_;
   TH1F* h1_PtEffpfSumChargedHadron_;
   TH1F* h1_PtEffpfSumNeutralHadron_;
   TH1F* h1_pfSumNeutralHadron_;
   TH1F* h1_PtEffpfPhotonIso_;
   TH1F* h1_pfPhotonIso_; 
   TH1F* h1_FracpfSumChargedHadron_;
   TH1F* h1_FracpfSumNeutralHadron_;
   TH1F* h1_FracpfPhotonIso_;

   TH1F* h1_PtID_;
   TH1F* h1_EtaID_;
   TH1F* h1_EtaEffpterror_;
   TH1F* h1_EtaEffptnumberOftrackerLayers_;
   TH1F* h1_EtaEffptnumberOfPixelHits_;
   TH1F* h1_EtaEffptnumberOfMuonHits_;
   TH1F* h1_EtaEffptnumberOfMatchedStations_;
   TH1F* h1_EtaEffptTrackIso_;
   TH1F* h1_EtaEffptabsdsy_;
   TH1F* h1_nbPVEffptTrackIso_;
   TH1F* h1_nbPVID_;
   TH2F* h2_RecoMass1vsMass2_;
   TH1F* h1_test_;
   //-------------- new histos ---------------
   TH1F* h1_pfDeltaBetaID_;
   TH1F* h1_PtEffPFdeltaBeta_;
   TH1F* h1_EtaEffPFdeltaBeta_;
   TH1F* h1_nbPVEffPFdeltaBeta_;
   
   TH1F* h1_PtHeepEleEB_;
   TH1F* h1_EtaHeepEleEB_;
   TH1F* h1_PVHeepEleEB_;
   TH1F* h1_PtHeepEleEE_;
   TH1F* h1_EtaHeepEleEE_;
   TH1F* h1_PVHeepEleEE_;  
   TH1F* h1_pfDeltaBetaEB_;
   TH1F* h1_pfDeltaBetaEE_;
   TH1F* h1_PtHeepdeltaBetaEB_;
   TH1F* h1_EtaHeepdeltaBetaEB_;
   TH1F* h1_PVHeepdeltaBetaEB_;
   TH1F* h1_PtHeepdeltaBetaEE_;
   TH1F* h1_EtaHeepdeltaBetaEE_;
   TH1F* h1_PVHeepdeltaBetaEE_;

   TH1F* h1_PtHeepphotonconvEB_;
   TH1F* h1_EtaHeepphotonconvEB_;
   TH1F* h1_PVHeepphotonconvEB_;
   TH1F* h1_PtHeepphotonconvEE_;
   TH1F* h1_EtaHeepphotonconvEE_;
   TH1F* h1_PVHeepphotonconvEE_;

   TH1F* h1_PtHeepHoeNewEB_;
   TH1F* h1_EtaHeepHoeNewEB_;
   TH1F* h1_PVHeepHoeNewEB_;
   TH1F* h1_PtHeepHoeNewEE_;
   TH1F* h1_EtaHeepHoeNewEE_;
   TH1F* h1_PVHeepHoeNewEE_;


   TH1F* h1_PtHeepDeltaEtaNewEB_;
   TH1F* h1_EtaHeepDeltaEtaNewEB_;
   TH1F* h1_PVHeepDeltaEtaNewEB_;
   TH1F* h1_PtHeepDeltaEtaNewEE_;
   TH1F* h1_EtaHeepDeltaEtaNewEE_;
   TH1F* h1_PVHeepDeltaEtaNewEE_;

   TH1F* h1_PtHeepDeltaPhiNewEB_;
   TH1F* h1_EtaHeepDeltaPhiNewEB_;
   TH1F* h1_PVHeepDeltaPhiNewEB_;
   TH1F* h1_PtHeepDeltaPhiNewEE_;
   TH1F* h1_EtaHeepDeltaPhiNewEE_;
   TH1F* h1_PVHeepDeltaPhiNewEE_;



   TH1F* h1_PtHeepTrackIsoNewEB_;
   TH1F* h1_EtaHeepTrackIsoNewEB_;
   TH1F* h1_PVHeepTrackIsoNewEB_;
   TH1F* h1_PtHeepTrackIsoNewEE_;
   TH1F* h1_EtaHeepTrackIsoNewEE_;
   TH1F* h1_PVHeepTrackIsoNewEE_;

   TH1F* h1_PtHeepMissHitsNewEB_;
   TH1F* h1_EtaHeepMissHitsNewEB_;
   TH1F* h1_PVHeepMissHitsNewEB_;
   TH1F* h1_PtHeepMissHitsNewEE_;
   TH1F* h1_EtaHeepMissHitsNewEE_;
   TH1F* h1_PVHeepMissHitsNewEE_;

   TH1F* h1_PtHeepSigmaNewEE_;
   TH1F* h1_EtaHeepSigmaNewEE_;
   TH1F* h1_PVHeepSigmaNewEE_;

   TH1F* h1_PtHeepEcalHcalIsoNewEB_;
   TH1F* h1_EtaHeepEcalHcalIsoNewEB_;
   TH1F* h1_PVHeepEcalHcalIsoNewEB_;
   TH1F* h1_PtHeepEcalHcalIsoNewEE_;
   TH1F* h1_EtaHeepEcalHcalIsoNewEE_;
   TH1F* h1_PVHeepEcalHcalIsoNewEE_;

   TH1F* h1_PtHeepShowerShapeNewEB_;
   TH1F* h1_EtaHeepShowerShapeNewEB_;
   TH1F* h1_PVHeepShowerShapeNewEB_; 

   TH1F* h1_pfSumChargedHadronEB_;
   TH1F* h1_pfSumNeutralHadronEB_;
   TH1F* h1_pfPhotonIsoEB_;
   TH1F* h1_pfSumChargedHadronEE_;
   TH1F* h1_pfSumNeutralHadronEE_;
   TH1F* h1_pfPhotonIsoEE_; 

   TH1F* h1_PtHeepEleV4EB_;
   TH1F* h1_EtaHeepEleV4EB_;
   TH1F* h1_PVHeepEleV4EB_;
   TH1F* h1_PtHeepEleV4EE_;
   TH1F* h1_EtaHeepEleV4EE_;
   TH1F* h1_PVHeepEleV4EE_;

   TH1F* h1_PtHeepHoeV4EB_;
   TH1F* h1_EtaHeepHoeV4EB_;
   TH1F* h1_PVHeepHoeV4EB_;
   TH1F* h1_PtHeepHoeV4EE_;
   TH1F* h1_EtaHeepHoeV4EE_;
   TH1F* h1_PVHeepHoeV4EE_;

   TH1F* h1_PtHeepDeltaEtaInV4EB_;
   TH1F* h1_EtaHeepDeltaEtaInV4EB_;
   TH1F* h1_PVHeepDeltaEtaInV4EB_;
   TH1F* h1_PtHeepDeltaEtaInV4EE_;
   TH1F* h1_EtaHeepDeltaEtaInV4EE_;
   TH1F* h1_PVHeepDeltaEtaInV4EE_;

   TH1F* h1_DxyNewEB_;
   TH1F* h1_DxyNewEE_;
   TH1F* h1_PtHeepDxyNewEB_;
   TH1F* h1_EtaHeepDxyNewEB_;
   TH1F* h1_PVHeepDxyNewEB_;
   TH1F* h1_PtHeepDxyNewEE_;
   TH1F* h1_EtaHeepDxyNewEE_;
   TH1F* h1_PVHeepDxyNewEE_; 
};

#endif

#ifdef ZprimeToEE_cxx
ZprimeToEE::ZprimeToEE(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("CMSSW720_MC_Zprime_13TeV_pattuple.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("CMSSW720_MC_Zprime_13TeV_pattuple.root");
      }
      f->GetObject("tree",tree);

   }
   Init(tree);
}

ZprimeToEE::~ZprimeToEE()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t ZprimeToEE::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t ZprimeToEE::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void ZprimeToEE::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   NbEventsPassTrigger = 0;
   NbEventsPassTriggerandPVcond = 0;
   Ele_nbElectrons = 0;
   Ele_isEcalDrivenSeed = 0;
   Ele_isPassConversionVeto = 0;
   Ele_charge = 0;
   Ele_nbOfMissingHits = 0;
   Ele_nbVtx = 0;
   Ele_Et = 0;
   Ele_pt = 0;
   Ele_thetaSC = 0;
   Ele_etaSC = 0;
   Ele_phiSC = 0;
   Ele_energySC = 0;
   Ele_preshowerEnergySC = 0;
   Ele_thetaTrack = 0;
   Ele_etaTrack = 0;
   Ele_phiTrack = 0;
   Ele_hadronicOverEm = 0;
   Ele_deltaEtaInSC = 0;
   Ele_deltaPhiInSC = 0;
   Ele_sigmaIetaIeta = 0;
   Ele_sigmaIetaIetaFull5x5 = 0;
   Ele_e2x5MaxFull5x5 = 0;
   Ele_e1x5Full5x5 = 0;
   Ele_e5x5Full5x5 = 0;
   Ele_e2x5MaxOver5x5Full5x5 = 0;
   Ele_e1x5Over5x5Full5x5 = 0;
   Ele_e2x5Right = 0;
   Ele_e2x5Left = 0;
   Ele_e2x5Top = 0;
   Ele_e2x5Bottom = 0;
   Ele_eMax = 0;
   Ele_eRight = 0;
   Ele_eLeft = 0;
   Ele_eTop = 0;
   Ele_eBottom = 0;
   Ele_dxy = 0;
   Ele_dz = 0;
   Ele_rhoIso = 0;
   Ele_fbrem = 0;
   Ele_EoverP = 0;
   Ele_Xposition = 0;
   Ele_Yposition = 0;
   Ele_EcalPlusHcald1iso = 0;
   Ele_dr03TkSumPt = 0;
   Ele_dr03EcalRecHitSumEt = 0;
   Ele_dr03HcalDepth1TowerSumEt = 0;
   Ele_dr03HcalDepth1TowerSumEtBc = 0;
   Ele_hcalDepth1OverEcal = 0;
   Ele_hcalDepth2OverEcal = 0;
   Ele_dr03HcalDepth2TowerSumEt = 0;
   Ele_hcalDepth2TowerSumEtNoVeto = 0;
   Ele_hcalDepth1TowerSumEtNoVeto = 0;
   Ele_pfSumPhotonEt = 0;
   Ele_pfSumChargedHadronPt = 0;
   Ele_pfSumNeutralHadronEt = 0;
   Ele_pfSumPUPt = 0;
   Ele_pfDeltaBeta = 0;
   Ele_pfCandidateDeltaBeta = 0;
   Ele_pfCandidateChargedHadron = 0;
   Ele_pfCandidateNeutralHadron = 0;
   iGen = 0;
   idGen = 0;
   statusGen = 0;
   ptGen = 0;
   etaGen = 0;
   phiGen = 0;
   chargeGen = 0;
   EnergyGen = 0;
   pxGen = 0;
   pyGen = 0;
   pzGen = 0;
   nbPv = 0;
   Nbdof = 0;
   PositionRho = 0;
   PositionX = 0;
   PositionY = 0;
   PositionZ = 0;
   ecal_iRechit = 0;
   ecal_RecoFlag = 0;
   ecal_isSaturated = 0;
   ecal_EnRecHit = 0;
   MuHLTObj_nbObj = 0;
   MuHLTObj_pt = 0;
   MuHLTObj_eta = 0;
   MuHLTObj_phi = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("NbEventsPassTrigger", &NbEventsPassTrigger, &b_NbEventsPassTrigger);
   fChain->SetBranchAddress("NbEventsPassTriggerandPVcond", &NbEventsPassTriggerandPVcond, &b_NbEventsPassTriggerandPVcond);
   fChain->SetBranchAddress("event_runNo", &event_runNo, &b_event_runNo);
   fChain->SetBranchAddress("event_evtNo", &event_evtNo, &b_event_evtNo);
   fChain->SetBranchAddress("event_lumi", &event_lumi, &b_event_lumi);
   fChain->SetBranchAddress("event_bunch", &event_bunch, &b_event_bunch);
   fChain->SetBranchAddress("Ele_nbElectrons", &Ele_nbElectrons, &b_Ele_nbElectrons);
   fChain->SetBranchAddress("Ele_isEcalDrivenSeed", &Ele_isEcalDrivenSeed, &b_Ele_isEcalDrivenSeed);
   fChain->SetBranchAddress("Ele_isPassConversionVeto", &Ele_isPassConversionVeto, &b_Ele_isPassConversionVeto);
   fChain->SetBranchAddress("Ele_charge", &Ele_charge, &b_Ele_charge);
   fChain->SetBranchAddress("Ele_nbOfMissingHits", &Ele_nbOfMissingHits, &b_Ele_nbOfMissingHits);
   fChain->SetBranchAddress("Ele_nbVtx", &Ele_nbVtx, &b_Ele_nbVtx);
   fChain->SetBranchAddress("Ele_Et", &Ele_Et, &b_Ele_Et);
   fChain->SetBranchAddress("Ele_pt", &Ele_pt, &b_Ele_pt);
   fChain->SetBranchAddress("Ele_thetaSC", &Ele_thetaSC, &b_Ele_thetaSC);
   fChain->SetBranchAddress("Ele_etaSC", &Ele_etaSC, &b_Ele_etaSC);
   fChain->SetBranchAddress("Ele_phiSC", &Ele_phiSC, &b_Ele_phiSC);
   fChain->SetBranchAddress("Ele_energySC", &Ele_energySC, &b_Ele_energySC);
   fChain->SetBranchAddress("Ele_preshowerEnergySC", &Ele_preshowerEnergySC, &b_Ele_preshowerEnergySC);
   fChain->SetBranchAddress("Ele_thetaTrack", &Ele_thetaTrack, &b_Ele_thetaTrack);
   fChain->SetBranchAddress("Ele_etaTrack", &Ele_etaTrack, &b_Ele_etaTrack);
   fChain->SetBranchAddress("Ele_phiTrack", &Ele_phiTrack, &b_Ele_phiTrack);
   fChain->SetBranchAddress("Ele_hadronicOverEm", &Ele_hadronicOverEm, &b_Ele_hadronicOverEm);
   fChain->SetBranchAddress("Ele_deltaEtaInSC", &Ele_deltaEtaInSC, &b_Ele_deltaEtaInSC);
   fChain->SetBranchAddress("Ele_deltaPhiInSC", &Ele_deltaPhiInSC, &b_Ele_deltaPhiInSC);
   fChain->SetBranchAddress("Ele_sigmaIetaIeta", &Ele_sigmaIetaIeta, &b_Ele_sigmaIetaIeta);
   fChain->SetBranchAddress("Ele_sigmaIetaIetaFull5x5", &Ele_sigmaIetaIetaFull5x5, &b_Ele_sigmaIetaIetaFull5x5);
   fChain->SetBranchAddress("Ele_e2x5MaxFull5x5", &Ele_e2x5MaxFull5x5, &b_Ele_e2x5MaxFull5x5);
   fChain->SetBranchAddress("Ele_e1x5Full5x5", &Ele_e1x5Full5x5, &b_Ele_e1x5Full5x5);
   fChain->SetBranchAddress("Ele_e5x5Full5x5", &Ele_e5x5Full5x5, &b_Ele_e5x5Full5x5);
   fChain->SetBranchAddress("Ele_e2x5MaxOver5x5Full5x5", &Ele_e2x5MaxOver5x5Full5x5, &b_Ele_e2x5MaxOver5x5Full5x5);
   fChain->SetBranchAddress("Ele_e1x5Over5x5Full5x5", &Ele_e1x5Over5x5Full5x5, &b_Ele_e1x5Over5x5Full5x5);
   fChain->SetBranchAddress("Ele_e2x5Right", &Ele_e2x5Right, &b_Ele_e2x5Right);
   fChain->SetBranchAddress("Ele_e2x5Left", &Ele_e2x5Left, &b_Ele_e2x5Left);
   fChain->SetBranchAddress("Ele_e2x5Top", &Ele_e2x5Top, &b_Ele_e2x5Top);
   fChain->SetBranchAddress("Ele_e2x5Bottom", &Ele_e2x5Bottom, &b_Ele_e2x5Bottom);
   fChain->SetBranchAddress("Ele_eMax", &Ele_eMax, &b_Ele_eMax);
   fChain->SetBranchAddress("Ele_eRight", &Ele_eRight, &b_Ele_eRight);
   fChain->SetBranchAddress("Ele_eLeft", &Ele_eLeft, &b_Ele_eLeft);
   fChain->SetBranchAddress("Ele_eTop", &Ele_eTop, &b_Ele_eTop);
   fChain->SetBranchAddress("Ele_eBottom", &Ele_eBottom, &b_Ele_eBottom);
   fChain->SetBranchAddress("Ele_dxy", &Ele_dxy, &b_Ele_dxy);
   fChain->SetBranchAddress("Ele_dz", &Ele_dz, &b_Ele_dz);
   fChain->SetBranchAddress("Ele_rhoIso", &Ele_rhoIso, &b_Ele_rhoIso);
   fChain->SetBranchAddress("Ele_fbrem", &Ele_fbrem, &b_Ele_fbrem);
   fChain->SetBranchAddress("Ele_EoverP", &Ele_EoverP, &b_Ele_EoverP);
   fChain->SetBranchAddress("Ele_Xposition", &Ele_Xposition, &b_Ele_Xposition);
   fChain->SetBranchAddress("Ele_Yposition", &Ele_Yposition, &b_Ele_Yposition);
   fChain->SetBranchAddress("Ele_EcalPlusHcald1iso", &Ele_EcalPlusHcald1iso, &b_Ele_EcalPlusHcald1iso);
   fChain->SetBranchAddress("Ele_dr03TkSumPt", &Ele_dr03TkSumPt, &b_Ele_dr03TkSumPt);
   fChain->SetBranchAddress("Ele_dr03EcalRecHitSumEt", &Ele_dr03EcalRecHitSumEt, &b_Ele_dr03EcalRecHitSumEt);
   fChain->SetBranchAddress("Ele_dr03HcalDepth1TowerSumEt", &Ele_dr03HcalDepth1TowerSumEt, &b_Ele_dr03HcalDepth1TowerSumEt);
   fChain->SetBranchAddress("Ele_dr03HcalDepth1TowerSumEtBc", &Ele_dr03HcalDepth1TowerSumEtBc, &b_Ele_dr03HcalDepth1TowerSumEtBc);
   fChain->SetBranchAddress("Ele_hcalDepth1OverEcal", &Ele_hcalDepth1OverEcal, &b_Ele_hcalDepth1OverEcal);
   fChain->SetBranchAddress("Ele_hcalDepth2OverEcal", &Ele_hcalDepth2OverEcal, &b_Ele_hcalDepth2OverEcal);
   fChain->SetBranchAddress("Ele_dr03HcalDepth2TowerSumEt", &Ele_dr03HcalDepth2TowerSumEt, &b_Ele_dr03HcalDepth2TowerSumEt);
   fChain->SetBranchAddress("Ele_hcalDepth2TowerSumEtNoVeto", &Ele_hcalDepth2TowerSumEtNoVeto, &b_Ele_hcalDepth2TowerSumEtNoVeto);
   fChain->SetBranchAddress("Ele_hcalDepth1TowerSumEtNoVeto", &Ele_hcalDepth1TowerSumEtNoVeto, &b_Ele_hcalDepth1TowerSumEtNoVeto);
   fChain->SetBranchAddress("Ele_pfSumPhotonEt", &Ele_pfSumPhotonEt, &b_Ele_pfSumPhotonEt);
   fChain->SetBranchAddress("Ele_pfSumChargedHadronPt", &Ele_pfSumChargedHadronPt, &b_Ele_pfSumChargedHadronPt);
   fChain->SetBranchAddress("Ele_pfSumNeutralHadronEt", &Ele_pfSumNeutralHadronEt, &b_Ele_pfSumNeutralHadronEt);
   fChain->SetBranchAddress("Ele_pfSumPUPt", &Ele_pfSumPUPt, &b_Ele_pfSumPUPt);
   fChain->SetBranchAddress("Ele_pfDeltaBeta", &Ele_pfDeltaBeta, &b_Ele_pfDeltaBeta);
   fChain->SetBranchAddress("Ele_pfCandidateDeltaBeta", &Ele_pfCandidateDeltaBeta, &b_Ele_pfCandidateDeltaBeta);
   fChain->SetBranchAddress("Ele_pfCandidateChargedHadron", &Ele_pfCandidateChargedHadron, &b_Ele_pfCandidateChargedHadron);
   fChain->SetBranchAddress("Ele_pfCandidateNeutralHadron", &Ele_pfCandidateNeutralHadron, &b_Ele_pfCandidateNeutralHadron);
   fChain->SetBranchAddress("iGen", &iGen, &b_iGen);
   fChain->SetBranchAddress("idGen", &idGen, &b_idGen);
   fChain->SetBranchAddress("statusGen", &statusGen, &b_statusGen);
   fChain->SetBranchAddress("ptGen", &ptGen, &b_ptGen);
   fChain->SetBranchAddress("etaGen", &etaGen, &b_etaGen);
   fChain->SetBranchAddress("phiGen", &phiGen, &b_phiGen);
   fChain->SetBranchAddress("chargeGen", &chargeGen, &b_chargeGen);
   fChain->SetBranchAddress("EnergyGen", &EnergyGen, &b_EnergyGen);
   fChain->SetBranchAddress("pxGen", &pxGen, &b_pxGen);
   fChain->SetBranchAddress("pyGen", &pyGen, &b_pyGen);
   fChain->SetBranchAddress("pzGen", &pzGen, &b_pzGen);
   fChain->SetBranchAddress("nbPv", &nbPv, &b_nbPv);
   fChain->SetBranchAddress("Nbdof", &Nbdof, &b_Nbdof);
   fChain->SetBranchAddress("PositionRho", &PositionRho, &b_PositionRho);
   fChain->SetBranchAddress("PositionX", &PositionX, &b_PositionX);
   fChain->SetBranchAddress("PositionY", &PositionY, &b_PositionY);
   fChain->SetBranchAddress("PositionZ", &PositionZ, &b_PositionZ);
   fChain->SetBranchAddress("ecal_iRechit", &ecal_iRechit, &b_ecal_iRechit);
   fChain->SetBranchAddress("ecal_RecoFlag", &ecal_RecoFlag, &b_ecal_RecoFlag);
   fChain->SetBranchAddress("ecal_isSaturated", &ecal_isSaturated, &b_ecal_isSaturated);
   fChain->SetBranchAddress("ecal_EnRecHit", &ecal_EnRecHit, &b_ecal_EnRecHit);
   fChain->SetBranchAddress("MuHLTObj_nbObj", &MuHLTObj_nbObj, &b_MuHLTObj_nbObj);
   fChain->SetBranchAddress("MuHLTObj_pt", &MuHLTObj_pt, &b_MuHLTObj_pt);
   fChain->SetBranchAddress("MuHLTObj_eta", &MuHLTObj_eta, &b_MuHLTObj_eta);
   fChain->SetBranchAddress("MuHLTObj_phi", &MuHLTObj_phi, &b_MuHLTObj_phi);
   Notify();
}

Bool_t ZprimeToEE::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void ZprimeToEE::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t ZprimeToEE::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef ZprimeToEE_cxx
