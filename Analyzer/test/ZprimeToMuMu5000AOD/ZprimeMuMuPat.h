//=========================================================================  
//          Analysis code for Z' boson to Mu Mu analysis                  =  
//    [1] In this code we select the high pt di-muons events              =
//                        To run over MC                                  =  
//             Written by Sherif Elgammal                                 =
//                                                                        =
//                      1/1/2015                                          =
//=========================================================================
#ifndef ZprimeMuMuPat_h
#define ZprimeMuMuPat_h

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
#include <fstream>
#include <sstream>
#include "TVector3.h"
// Header file for the classes stored in the TTree if any.
#include <vector>
using namespace std;
// Fixed size dimensions of array or collections stored in the TTree if any.
class ZprimeMuMuPat {
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
   vector<int>     *Mu_nbMuon;
   vector<bool>    *Mu_isGlobalMuon;
   vector<bool>    *Mu_isPF;
   vector<bool>    *Mu_isGoodMuon;
   vector<float>   *Mu_et;
   vector<float>   *Mu_en;
   vector<float>   *Mu_ptcocktail;
   vector<float>   *Mu_etaCocktail;
   vector<float>   *Mu_phiCocktail;
   vector<float>   *Mu_thetaCocktail;
   vector<float>   *Mu_pxCocktail;
   vector<float>   *Mu_pyCocktail;
   vector<float>   *Mu_pzCocktail;
   vector<float>   *Mu_pCocktail;
   vector<float>   *Mu_dPToverPTcocktail;
   vector<int>     *Mu_chargeCocktail;
   vector<float>   *Mu_absdxy;
   vector<float>   *Mu_absdz;
   vector<float>   *Mu_normalizedChi2;
   vector<float>   *Mu_vtxMass;
   vector<float>   *Mu_vtxNormChi2;
   vector<int>     *Mu_numberOfMatchedStations;
   vector<int>     *Mu_numberOfValidPixelHits;
   vector<int>     *Mu_numberOfValidMuonHits;
   vector<int>     *Mu_numberOftrackerLayersWithMeasurement;
   vector<float>   *Mu_emIso;
   vector<float>   *Mu_hadIso;
   vector<float>   *Mu_trackiso;
   vector<float>   *Mu_pfSumChargedHadronPt;
   vector<float>   *Mu_pfSumNeutralHadronEt;
   vector<float>   *Mu_PFSumPhotonEt;
   vector<float>   *Mu_pfSumPUPt;
   vector<int>     *Mu_nbofpv;
   vector<float>   *Mu_ptInnerTrack;
   vector<float>   *Mu_ptTunePMuonBestTrack;
   vector<float>   *Mu_dPToverPTTunePMuonBestTrack;
   vector<float>   *Mu_pxTunePMuonBestTrack;
   vector<float>   *Mu_pyTunePMuonBestTrack;
   vector<float>   *Mu_pzTunePMuonBestTrack;
   vector<float>   *Mu_pTunePMuonBestTrack;
   vector<float>   *Mu_etaTunePMuonBestTrack;
   vector<float>   *Mu_phiTunePMuonBestTrack;
   vector<float>   *Mu_thetaTunePMuonBestTrack;
   vector<float>   *Mu_chargeTunePMuonBestTrack;
   vector<float>   *Mu_absdxyTunePMuonBestTrack;
   vector<float>   *Mu_ptDYTTrack;
   vector<float>   *Mu_pxDYTTrack;
   vector<float>   *Mu_pyDYTTrack;
   vector<float>   *Mu_pzDYTTrack;
   vector<float>   *Mu_pDYTTrack;
   vector<float>   *Mu_etaDYTTrack;
   vector<float>   *Mu_phiDYTTrack;
   vector<float>   *Mu_thetaDYTTrack;
   vector<float>   *Mu_chargeDYTTrack;
   vector<float>   *Mu_absdxyDYTTrack;
   vector<float>   *Mu_dPToverPTDYTTrack;
   vector<float>   *Mu_ptPickyTrack;
   vector<float>   *Mu_pxPickyTrack;
   vector<float>   *Mu_pyPickyTrack;
   vector<float>   *Mu_pzPickyTrack;
   vector<float>   *Mu_pPickyTrack;
   vector<float>   *Mu_etaPickyTrack;
   vector<float>   *Mu_phiPickyTrack;
   vector<float>   *Mu_thetaPickyTrack;
   vector<float>   *Mu_chargePickyTrack;
   vector<float>   *Mu_absdxyPickyTrack;
   vector<float>   *Mu_dPToverPTPickyTrack;
   vector<float>   *Mu_ptMuonBestTrack;
   vector<float>   *Mu_dPToverPTMuonBestTrack;
   vector<float>   *Mu_pxMuonBestTrack;
   vector<float>   *Mu_pyMuonBestTrack;
   vector<float>   *Mu_pzMuonBestTrack;
   vector<float>   *Mu_pMuonBestTrack;
   vector<float>   *Mu_etaMuonBestTrack;
   vector<float>   *Mu_phiMuonBestTrack;
   vector<float>   *Mu_thetaMuonBestTrack;
   vector<float>   *Mu_chargeMuonBestTrack;
   vector<float>   *Mu_absdxyMuonBestTrack;
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
   vector<int>     *NbOfDaughters;
   vector<float>   *McZmass;
   vector<int>     *nbPv;
   vector<int>     *Nbdof;
   vector<float>   *PositionRho;
   vector<float>   *PositionX;
   vector<float>   *PositionY;
   vector<float>   *PositionZ;
   vector<int>     *MuHLTMatch_nbMuonMatchHLT;
   vector<float>   *MuHLTMatch_pt;
   vector<float>   *MuHLTMatch_eta;
   vector<float>   *MuHLTMatch_phi;

   // List of branches
   TBranch        *b_NbEventsPassTrigger;   //!
   TBranch        *b_NbEventsPassTriggerandPVcond;   //!
   TBranch        *b_event_runNo;   //!
   TBranch        *b_event_evtNo;   //!
   TBranch        *b_event_lumi;   //!
   TBranch        *b_event_bunch;   //!
   TBranch        *b_Mu_nbMuon;   //!
   TBranch        *b_Mu_isGlobalMuon;   //!
   TBranch        *b_Mu_isPF;   //!
   TBranch        *b_Mu_isGoodMuon;   //!
   TBranch        *b_Mu_et;   //!
   TBranch        *b_Mu_en;   //!
   TBranch        *b_Mu_ptcocktail;   //!
   TBranch        *b_Mu_etaCocktail;   //!
   TBranch        *b_Mu_phiCocktail;   //!
   TBranch        *b_Mu_thetaCocktail;   //!
   TBranch        *b_Mu_pxCocktail;   //!
   TBranch        *b_Mu_pyCocktail;   //!
   TBranch        *b_Mu_pzCocktail;   //!
   TBranch        *b_Mu_pCocktail;   //!
   TBranch        *b_Mu_dPToverPTcocktail;   //!
   TBranch        *b_Mu_chargeCocktail;   //!
   TBranch        *b_Mu_absdxy;   //!
   TBranch        *b_Mu_absdz;   //!
   TBranch        *b_Mu_normalizedChi2;   //!
   TBranch        *b_Mu_vtxMass;   //!
   TBranch        *b_Mu_vtxNormChi2;   //!
   TBranch        *b_Mu_numberOfMatchedStations;   //!
   TBranch        *b_Mu_numberOfValidPixelHits;   //!
   TBranch        *b_Mu_numberOfValidMuonHits;   //!
   TBranch        *b_Mu_numberOftrackerLayersWithMeasurement;   //!
   TBranch        *b_Mu_emIso;   //!
   TBranch        *b_Mu_hadIso;   //!
   TBranch        *b_Mu_trackiso;   //!
   TBranch        *b_Mu_pfSumChargedHadronPt;   //!
   TBranch        *b_Mu_pfSumNeutralHadronEt;   //!
   TBranch        *b_Mu_PFSumPhotonEt;   //!
   TBranch        *b_Mu_pfSumPUPt;   //!
   TBranch        *b_Mu_nbofpv;   //!
   TBranch        *b_Mu_ptInnerTrack;   //!
   TBranch        *b_Mu_ptTunePMuonBestTrack;   //!
   TBranch        *b_Mu_dPToverPTTunePMuonBestTrack;   //!
   TBranch        *b_Mu_pxTunePMuonBestTrack;   //!
   TBranch        *b_Mu_pyTunePMuonBestTrack;   //!
   TBranch        *b_Mu_pzTunePMuonBestTrack;   //!
   TBranch        *b_Mu_pTunePMuonBestTrack;   //!
   TBranch        *b_Mu_etaTunePMuonBestTrack;   //!
   TBranch        *b_Mu_phiTunePMuonBestTrack;   //!
   TBranch        *b_Mu_thetaTunePMuonBestTrack;   //!
   TBranch        *b_Mu_chargeTunePMuonBestTrack;   //!
   TBranch        *b_Mu_absdxyTunePMuonBestTrack;   //!
   TBranch        *b_Mu_ptDYTTrack;   //!
   TBranch        *b_Mu_pxDYTTrack;   //!
   TBranch        *b_Mu_pyDYTTrack;   //!
   TBranch        *b_Mu_pzDYTTrack;   //!
   TBranch        *b_Mu_pDYTTrack;   //!
   TBranch        *b_Mu_etaDYTTrack;   //!
   TBranch        *b_Mu_phiDYTTrack;   //!
   TBranch        *b_Mu_thetaDYTTrack;   //!
   TBranch        *b_Mu_chargeDYTTrack;   //!
   TBranch        *b_Mu_absdxyDYTTrack;   //!
   TBranch        *b_Mu_dPToverPTDYTTrack;   //!
   TBranch        *b_Mu_ptPickyTrack;   //!
   TBranch        *b_Mu_pxPickyTrack;   //!
   TBranch        *b_Mu_pyPickyTrack;   //!
   TBranch        *b_Mu_pzPickyTrack;   //!
   TBranch        *b_Mu_pPickyTrack;   //!
   TBranch        *b_Mu_etaPickyTrack;   //!
   TBranch        *b_Mu_phiPickyTrack;   //!
   TBranch        *b_Mu_thetaPickyTrack;   //!
   TBranch        *b_Mu_chargePickyTrack;   //!
   TBranch        *b_Mu_absdxyPickyTrack;   //!
   TBranch        *b_Mu_dPToverPTPickyTrack;   //!
   TBranch        *b_Mu_ptMuonBestTrack;   //!
   TBranch        *b_Mu_dPToverPTMuonBestTrack;   //!
   TBranch        *b_Mu_pxMuonBestTrack;   //!
   TBranch        *b_Mu_pyMuonBestTrack;   //!
   TBranch        *b_Mu_pzMuonBestTrack;   //!
   TBranch        *b_Mu_pMuonBestTrack;   //!
   TBranch        *b_Mu_etaMuonBestTrack;   //!
   TBranch        *b_Mu_phiMuonBestTrack;   //!
   TBranch        *b_Mu_thetaMuonBestTrack;   //!
   TBranch        *b_Mu_chargeMuonBestTrack;   //!
   TBranch        *b_Mu_absdxyMuonBestTrack;   //!
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
   TBranch        *b_NbOfDaughters;   //!
   TBranch        *b_McZmass;   //!
   TBranch        *b_nbPv;   //!
   TBranch        *b_Nbdof;   //!
   TBranch        *b_PositionRho;   //!
   TBranch        *b_PositionX;   //!
   TBranch        *b_PositionY;   //!
   TBranch        *b_PositionZ;   //!
   TBranch        *b_MuHLTMatch_nbMuonMatchHLT;   //!
   TBranch        *b_MuHLTMatch_pt;   //!
   TBranch        *b_MuHLTMatch_eta;   //!
   TBranch        *b_MuHLTMatch_phi;   //!

   ZprimeMuMuPat(TTree *tree=0);
   virtual ~ZprimeMuMuPat();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   bool SelectFirstMuon(float &pTmuon1,float &Enmuon1,float &Etamuon1,
                        float &Phimuon1,int &ChargeMu1,unsigned &FlagMu1,
                        float &pxmuon1,float &pymuon1,float &pzmuon1,
                        float &pmuon1,float &dxymuon1,float &pTmuon1tuneP,
			float &pTmuonBestTrack1);
   bool SelectSecondMuon(int ChargeMu1,unsigned FlagMu1,float pTmuon1,float &pTmuon2,float &Enmuon2,
                         float &Etamuon2,float &Phimuon2,int &ChargeMu2,float &pxmuon2,
                         float &pymuon2,float &pzmuon2,float &pmuon2,float &dxymuon2,
			 float &pTmuon2tuneP,float &pTmuonBestTrack2);
   float Mass(float Pt1,float Eta1,float Phi1,float En1,
              float Pt2,float Eta2,float Phi2,float En2);
   void PlotRecoInfo(float CosmicMuonRejec,float vertexMassMu,float MassGenerated,
		     float PtTunePMuBestTrack,float PtTunePMu,float PtMuBestTrack,
		     float PtGenerated);
   void PickThehighestMass(float &vtxHighestMass,float &vtxHighestChi2);
   double ThreeDangle(float pxMu1,float pyMu1,float pzMu1,float pMu1,
                      float pxMu2,float pyMu2,float pzMu2,float pMu2);

   bool SelectFirstGenMu(float &ETMu1,float &PhiSCMu1,
                         float &EtaSCMu1,float &EnMu1,
                         int &IDele1,int &Statele1,
                         unsigned &GenFlag1);

   bool SelectSecondGenMu(unsigned GenFlag1,float ETMu1,float &ETMu2,float &PhiSCMu2,
                          float &EtaSCMu2,float &EnMu2,int &IDele2,int &Statele2);

   bool GenRecoMatchMu1(float RecoEta1,float RecoPhi1,
                        float &ETMu1,float &PhiSCMu1,
                        float &EtaSCMu1,float &EnMu1,
                        unsigned &GenFlag1);

    bool GenRecoMatchMu2(unsigned GenFlag1,float RecoEta2,float RecoPhi2,
                        float &ETMu2,float &PhiSCMu2,
                        float &EtaSCMu2,float &EnMu2);

    void PlotGenInfo(float ZprimeGenMass,float EtaGenMu1,float EtaGenMu2,float PtGenMu1,
		     float PtGenMu2,float EnGenMu1,float EnGenMu2);

   void MuonPassingID();
   void PlotPterror();
   void PlotNbTrackLayers();
   void PlotNBValidPixelHits();
   void PlotNbValidMuonHits();
   void PlotNbMatchedStations();
   void PlotTrackiso();
   void PlotAbsDxy();
   void plotAllHighPtMuonsID();
   void MuonPassingNewID();
   void MuonPassingTightID();
   void CosThetaCollinSoper(float Et1,float Eta1,float Phi1,float En1,
			    float Et2,float Eta2,float Phi2,float En2,
			    float ChargeEle1,float RecoMass);
   //================================================================================
   float PtDYTRecMu1,PtDYTRecMu2,PtRecTunePMu1,PtRecTunePMu2,PtRecMuBestTrack1,PtRecMuBestTrack2;
   float RecoHLTMatchingDeltaRcut,deltaRcut,minMassCut,maxMassCut;
   float vtxChi2Mu,vtxMassMu;
   float mPtGen1,mPhiGen1,mEtaGen1,mEnGen1;
   unsigned mGenFlag1;
   float mPtGen2,mPhiGen2,mEtaGen2,mEnGen2;
   int ChargeRecMu1,ChargeRecMu2;
   unsigned flagmu1;
   unsigned flag1;
   float PtRecTunePMuBestTrack1,EnRecMu1,EtaRecMu1,PhiRecMu1;
   float PtRecTunePMuBestTrack2,EnRecMu2,EtaRecMu2,PhiRecMu2;
   float pxRecMu1,pyRecMu1,pzRecMu1,pRecMu1,dxyRecMu1;
   float pxRecMu2,pyRecMu2,pzRecMu2,pRecMu2,dxyRecMu2;
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
   TH1F* h1_ZprimeRecomassBeforeTrigger_;
   TH1F* h1_ZprimeRecomass_;
   TH1F* h1_ZprimeRecomassAbove400GeV_;
   TH1F* h1_ZprimeRecomassAbove1000GeV_;
   TH1F* h1_3Dangle_;
   TH1F* h1_DxyDiff_;
   TH1F* h1_ZprimeGenmass_;
   TH1F* h1_ZprimeGenEta1_;
   TH1F* h1_ZprimeGenEta2_;
   TH1F* h1_ZprimeGenPt1_;
   TH1F* h1_ZprimeGenPt2_;
   TH1F* h1_ZprimeGenEn1_;
   TH1F* h1_ZprimeGenEn2_;
   TH1F* h1_MassRecoGenDif_;
   TH1F* h1_dPToverPT_;
   TH1F* h1_normalizedChi2_;
   TH1F* h1_numberOftrackerLayersWithMeasurement_;
   TH1F* h1_numberOfValidPixelHits_;
   TH1F* h1_numberOfValidMuonHits_;
   TH1F* h1_numberOfMatchedStations_;
   TH1F* h1_trackiso_;
   TH1F* h1_absdxy_;
   TH1F* h1_PtEffpterror_;
   TH1F* h1_PtEffptnumberOftrackerLayers_;
   TH1F* h1_PtEffptnumberOfPixelHits_;
   TH1F* h1_PtEffptnumberOfMuonHits_;
   TH1F* h1_PtEffptnumberOfMatchedStations_;
   TH1F* h1_PtEffptTrackIso_;
   TH1F* h1_PtEffptabsdsy_;
   TH1F* h1_PtEffpfSumChargedHadron_;
   TH1F* h1_PtEffpfSumNeutralHadron_;
   TH1F* h1_PtEffpfPhotonIso_;
   TH1F* h1_FracpfSumChargedHadron_;
   TH1F* h1_FracpfSumNeutralHadron_;
   TH1F* h1_FracpfPhotonIso_;
   TH1F* h1_EtaEffpterror_;
   TH1F* h1_EtaEffptnumberOftrackerLayers_;
   TH1F* h1_EtaEffptnumberOfPixelHits_;
   TH1F* h1_EtaEffptnumberOfMuonHits_;
   TH1F* h1_EtaEffptnumberOfMatchedStations_;
   TH1F* h1_EtaEffptTrackIso_;
   TH1F* h1_EtaEffptabsdsy_;
   TH1F* h1_nbPVID_;
   TH1F* h1_PtID_;
   TH1F* h1_EtaID_;
   TH1F* h1_nbPVNewID_;
   TH1F* h1_PtNewID_;
   TH1F* h1_EtaNewID_;
   TH1F* h1_nbPVTightID_;
   TH1F* h1_PtTightID_;
   TH1F* h1_EtaTightID_;
   ofstream output_txt; 
   TH1F* h1_PtResolutionMBT_;
   TH1F* h1_PtResolutionTunePMBT_;
   TH1F* h1_PtResolutiontuneP_;
   TH1F* h1_MassResultionEBEB1_;
   TH1F* h1_MassResultionEBEB2_;
   TH1F* h1_MassResultionEBEB3_;
   TH1F* h1_MassResultionEBEB4_;
   TH1F* h1_MassResultionEBEB5_;
   TH1F* h1_MassResultionEBEB6_;
   TH1F* h1_MassResultionEBEB7_;
   TH1F* h1_MassResultionEBEB8_;
   TH1F* h1_MassResultionEBEB9_;
   TH1F* h1_MassResultionEBEB10_; 
   TH1F* h1_MassGenInAccep_;
   TH1F* h1_MassRecoInAccep_;
   TH1F* h1_CosAngleCollinSoperCorrect60Mass120_;
   TH1F* h1_CosAngleCollinSoperCorrect120Mass300_;
   TH1F* h1_CosAngleCollinSoperCorrect300Mass700_;
   TH1F* h1_CosAngleCollinSoperCorrect700Mass3000_;
   TH1F* h1_CosAngleCollinSoperCorrect4900Mass5100_;
};

#endif

#ifdef ZprimeMuMuPat_cxx
ZprimeMuMuPat::ZprimeMuMuPat(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("ZprimetoMuMu-MC-Mass5000-CMSSW720-tree.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("ZprimetoMuMu-MC-Mass5000-CMSSW720-tree.root");
      }
      f->GetObject("tree",tree);

   }
   Init(tree);
}

ZprimeMuMuPat::~ZprimeMuMuPat()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t ZprimeMuMuPat::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t ZprimeMuMuPat::LoadTree(Long64_t entry)
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

void ZprimeMuMuPat::Init(TTree *tree)
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
   Mu_nbMuon = 0;
   Mu_isGlobalMuon = 0;
   Mu_isPF = 0;
   Mu_isGoodMuon = 0;
   Mu_et = 0;
   Mu_en = 0;
   Mu_ptcocktail = 0;
   Mu_etaCocktail = 0;
   Mu_phiCocktail = 0;
   Mu_thetaCocktail = 0;
   Mu_pxCocktail = 0;
   Mu_pyCocktail = 0;
   Mu_pzCocktail = 0;
   Mu_pCocktail = 0;
   Mu_dPToverPTcocktail = 0;
   Mu_chargeCocktail = 0;
   Mu_absdxy = 0;
   Mu_absdz = 0;
   Mu_normalizedChi2 = 0;
   Mu_vtxMass = 0;
   Mu_vtxNormChi2 = 0;
   Mu_numberOfMatchedStations = 0;
   Mu_numberOfValidPixelHits = 0;
   Mu_numberOfValidMuonHits = 0;
   Mu_numberOftrackerLayersWithMeasurement = 0;
   Mu_emIso = 0;
   Mu_hadIso = 0;
   Mu_trackiso = 0;
   Mu_pfSumChargedHadronPt = 0;
   Mu_pfSumNeutralHadronEt = 0;
   Mu_PFSumPhotonEt = 0;
   Mu_pfSumPUPt = 0;
   Mu_nbofpv = 0;
   Mu_ptInnerTrack = 0;
   Mu_ptTunePMuonBestTrack = 0;
   Mu_dPToverPTTunePMuonBestTrack = 0;
   Mu_pxTunePMuonBestTrack = 0;
   Mu_pyTunePMuonBestTrack = 0;
   Mu_pzTunePMuonBestTrack = 0;
   Mu_pTunePMuonBestTrack = 0;
   Mu_etaTunePMuonBestTrack = 0;
   Mu_phiTunePMuonBestTrack = 0;
   Mu_thetaTunePMuonBestTrack = 0;
   Mu_chargeTunePMuonBestTrack = 0;
   Mu_absdxyTunePMuonBestTrack = 0;
   Mu_ptDYTTrack = 0;
   Mu_pxDYTTrack = 0;
   Mu_pyDYTTrack = 0;
   Mu_pzDYTTrack = 0;
   Mu_pDYTTrack = 0;
   Mu_etaDYTTrack = 0;
   Mu_phiDYTTrack = 0;
   Mu_thetaDYTTrack = 0;
   Mu_chargeDYTTrack = 0;
   Mu_absdxyDYTTrack = 0;
   Mu_dPToverPTDYTTrack = 0;
   Mu_ptPickyTrack = 0;
   Mu_pxPickyTrack = 0;
   Mu_pyPickyTrack = 0;
   Mu_pzPickyTrack = 0;
   Mu_pPickyTrack = 0;
   Mu_etaPickyTrack = 0;
   Mu_phiPickyTrack = 0;
   Mu_thetaPickyTrack = 0;
   Mu_chargePickyTrack = 0;
   Mu_absdxyPickyTrack = 0;
   Mu_dPToverPTPickyTrack = 0;
   Mu_ptMuonBestTrack = 0;
   Mu_dPToverPTMuonBestTrack = 0;
   Mu_pxMuonBestTrack = 0;
   Mu_pyMuonBestTrack = 0;
   Mu_pzMuonBestTrack = 0;
   Mu_pMuonBestTrack = 0;
   Mu_etaMuonBestTrack = 0;
   Mu_phiMuonBestTrack = 0;
   Mu_thetaMuonBestTrack = 0;
   Mu_chargeMuonBestTrack = 0;
   Mu_absdxyMuonBestTrack = 0;
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
   NbOfDaughters = 0;
   McZmass = 0;
   nbPv = 0;
   Nbdof = 0;
   PositionRho = 0;
   PositionX = 0;
   PositionY = 0;
   PositionZ = 0;
   MuHLTMatch_nbMuonMatchHLT = 0;
   MuHLTMatch_pt = 0;
   MuHLTMatch_eta = 0;
   MuHLTMatch_phi = 0;
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
   fChain->SetBranchAddress("Mu_nbMuon", &Mu_nbMuon, &b_Mu_nbMuon);
   fChain->SetBranchAddress("Mu_isGlobalMuon", &Mu_isGlobalMuon, &b_Mu_isGlobalMuon);
   fChain->SetBranchAddress("Mu_isPF", &Mu_isPF, &b_Mu_isPF);
   fChain->SetBranchAddress("Mu_isGoodMuon", &Mu_isGoodMuon, &b_Mu_isGoodMuon);
   fChain->SetBranchAddress("Mu_et", &Mu_et, &b_Mu_et);
   fChain->SetBranchAddress("Mu_en", &Mu_en, &b_Mu_en);
   fChain->SetBranchAddress("Mu_ptcocktail", &Mu_ptcocktail, &b_Mu_ptcocktail);
   fChain->SetBranchAddress("Mu_etaCocktail", &Mu_etaCocktail, &b_Mu_etaCocktail);
   fChain->SetBranchAddress("Mu_phiCocktail", &Mu_phiCocktail, &b_Mu_phiCocktail);
   fChain->SetBranchAddress("Mu_thetaCocktail", &Mu_thetaCocktail, &b_Mu_thetaCocktail);
   fChain->SetBranchAddress("Mu_pxCocktail", &Mu_pxCocktail, &b_Mu_pxCocktail);
   fChain->SetBranchAddress("Mu_pyCocktail", &Mu_pyCocktail, &b_Mu_pyCocktail);
   fChain->SetBranchAddress("Mu_pzCocktail", &Mu_pzCocktail, &b_Mu_pzCocktail);
   fChain->SetBranchAddress("Mu_pCocktail", &Mu_pCocktail, &b_Mu_pCocktail);
   fChain->SetBranchAddress("Mu_dPToverPTcocktail", &Mu_dPToverPTcocktail, &b_Mu_dPToverPTcocktail);
   fChain->SetBranchAddress("Mu_chargeCocktail", &Mu_chargeCocktail, &b_Mu_chargeCocktail);
   fChain->SetBranchAddress("Mu_absdxy", &Mu_absdxy, &b_Mu_absdxy);
   fChain->SetBranchAddress("Mu_absdz", &Mu_absdz, &b_Mu_absdz);
   fChain->SetBranchAddress("Mu_normalizedChi2", &Mu_normalizedChi2, &b_Mu_normalizedChi2);
   fChain->SetBranchAddress("Mu_vtxMass", &Mu_vtxMass, &b_Mu_vtxMass);
   fChain->SetBranchAddress("Mu_vtxNormChi2", &Mu_vtxNormChi2, &b_Mu_vtxNormChi2);
   fChain->SetBranchAddress("Mu_numberOfMatchedStations", &Mu_numberOfMatchedStations, &b_Mu_numberOfMatchedStations);
   fChain->SetBranchAddress("Mu_numberOfValidPixelHits", &Mu_numberOfValidPixelHits, &b_Mu_numberOfValidPixelHits);
   fChain->SetBranchAddress("Mu_numberOfValidMuonHits", &Mu_numberOfValidMuonHits, &b_Mu_numberOfValidMuonHits);
   fChain->SetBranchAddress("Mu_numberOftrackerLayersWithMeasurement", &Mu_numberOftrackerLayersWithMeasurement, &b_Mu_numberOftrackerLayersWithMeasurement);
   fChain->SetBranchAddress("Mu_emIso", &Mu_emIso, &b_Mu_emIso);
   fChain->SetBranchAddress("Mu_hadIso", &Mu_hadIso, &b_Mu_hadIso);
   fChain->SetBranchAddress("Mu_trackiso", &Mu_trackiso, &b_Mu_trackiso);
   fChain->SetBranchAddress("Mu_pfSumChargedHadronPt", &Mu_pfSumChargedHadronPt, &b_Mu_pfSumChargedHadronPt);
   fChain->SetBranchAddress("Mu_pfSumNeutralHadronEt", &Mu_pfSumNeutralHadronEt, &b_Mu_pfSumNeutralHadronEt);
   fChain->SetBranchAddress("Mu_PFSumPhotonEt", &Mu_PFSumPhotonEt, &b_Mu_PFSumPhotonEt);
   fChain->SetBranchAddress("Mu_pfSumPUPt", &Mu_pfSumPUPt, &b_Mu_pfSumPUPt);
   fChain->SetBranchAddress("Mu_nbofpv", &Mu_nbofpv, &b_Mu_nbofpv);
   fChain->SetBranchAddress("Mu_ptInnerTrack", &Mu_ptInnerTrack, &b_Mu_ptInnerTrack);
   fChain->SetBranchAddress("Mu_ptTunePMuonBestTrack", &Mu_ptTunePMuonBestTrack, &b_Mu_ptTunePMuonBestTrack);
   fChain->SetBranchAddress("Mu_dPToverPTTunePMuonBestTrack", &Mu_dPToverPTTunePMuonBestTrack, &b_Mu_dPToverPTTunePMuonBestTrack);
   fChain->SetBranchAddress("Mu_pxTunePMuonBestTrack", &Mu_pxTunePMuonBestTrack, &b_Mu_pxTunePMuonBestTrack);
   fChain->SetBranchAddress("Mu_pyTunePMuonBestTrack", &Mu_pyTunePMuonBestTrack, &b_Mu_pyTunePMuonBestTrack);
   fChain->SetBranchAddress("Mu_pzTunePMuonBestTrack", &Mu_pzTunePMuonBestTrack, &b_Mu_pzTunePMuonBestTrack);
   fChain->SetBranchAddress("Mu_pTunePMuonBestTrack", &Mu_pTunePMuonBestTrack, &b_Mu_pTunePMuonBestTrack);
   fChain->SetBranchAddress("Mu_etaTunePMuonBestTrack", &Mu_etaTunePMuonBestTrack, &b_Mu_etaTunePMuonBestTrack);
   fChain->SetBranchAddress("Mu_phiTunePMuonBestTrack", &Mu_phiTunePMuonBestTrack, &b_Mu_phiTunePMuonBestTrack);
   fChain->SetBranchAddress("Mu_thetaTunePMuonBestTrack", &Mu_thetaTunePMuonBestTrack, &b_Mu_thetaTunePMuonBestTrack);
   fChain->SetBranchAddress("Mu_chargeTunePMuonBestTrack", &Mu_chargeTunePMuonBestTrack, &b_Mu_chargeTunePMuonBestTrack);
   fChain->SetBranchAddress("Mu_absdxyTunePMuonBestTrack", &Mu_absdxyTunePMuonBestTrack, &b_Mu_absdxyTunePMuonBestTrack);
   fChain->SetBranchAddress("Mu_ptDYTTrack", &Mu_ptDYTTrack, &b_Mu_ptDYTTrack);
   fChain->SetBranchAddress("Mu_pxDYTTrack", &Mu_pxDYTTrack, &b_Mu_pxDYTTrack);
   fChain->SetBranchAddress("Mu_pyDYTTrack", &Mu_pyDYTTrack, &b_Mu_pyDYTTrack);
   fChain->SetBranchAddress("Mu_pzDYTTrack", &Mu_pzDYTTrack, &b_Mu_pzDYTTrack);
   fChain->SetBranchAddress("Mu_pDYTTrack", &Mu_pDYTTrack, &b_Mu_pDYTTrack);
   fChain->SetBranchAddress("Mu_etaDYTTrack", &Mu_etaDYTTrack, &b_Mu_etaDYTTrack);
   fChain->SetBranchAddress("Mu_phiDYTTrack", &Mu_phiDYTTrack, &b_Mu_phiDYTTrack);
   fChain->SetBranchAddress("Mu_thetaDYTTrack", &Mu_thetaDYTTrack, &b_Mu_thetaDYTTrack);
   fChain->SetBranchAddress("Mu_chargeDYTTrack", &Mu_chargeDYTTrack, &b_Mu_chargeDYTTrack);
   fChain->SetBranchAddress("Mu_absdxyDYTTrack", &Mu_absdxyDYTTrack, &b_Mu_absdxyDYTTrack);
   fChain->SetBranchAddress("Mu_dPToverPTDYTTrack", &Mu_dPToverPTDYTTrack, &b_Mu_dPToverPTDYTTrack);
   fChain->SetBranchAddress("Mu_ptPickyTrack", &Mu_ptPickyTrack, &b_Mu_ptPickyTrack);
   fChain->SetBranchAddress("Mu_pxPickyTrack", &Mu_pxPickyTrack, &b_Mu_pxPickyTrack);
   fChain->SetBranchAddress("Mu_pyPickyTrack", &Mu_pyPickyTrack, &b_Mu_pyPickyTrack);
   fChain->SetBranchAddress("Mu_pzPickyTrack", &Mu_pzPickyTrack, &b_Mu_pzPickyTrack);
   fChain->SetBranchAddress("Mu_pPickyTrack", &Mu_pPickyTrack, &b_Mu_pPickyTrack);
   fChain->SetBranchAddress("Mu_etaPickyTrack", &Mu_etaPickyTrack, &b_Mu_etaPickyTrack);
   fChain->SetBranchAddress("Mu_phiPickyTrack", &Mu_phiPickyTrack, &b_Mu_phiPickyTrack);
   fChain->SetBranchAddress("Mu_thetaPickyTrack", &Mu_thetaPickyTrack, &b_Mu_thetaPickyTrack);
   fChain->SetBranchAddress("Mu_chargePickyTrack", &Mu_chargePickyTrack, &b_Mu_chargePickyTrack);
   fChain->SetBranchAddress("Mu_absdxyPickyTrack", &Mu_absdxyPickyTrack, &b_Mu_absdxyPickyTrack);
   fChain->SetBranchAddress("Mu_dPToverPTPickyTrack", &Mu_dPToverPTPickyTrack, &b_Mu_dPToverPTPickyTrack);
   fChain->SetBranchAddress("Mu_ptMuonBestTrack", &Mu_ptMuonBestTrack, &b_Mu_ptMuonBestTrack);
   fChain->SetBranchAddress("Mu_dPToverPTMuonBestTrack", &Mu_dPToverPTMuonBestTrack, &b_Mu_dPToverPTMuonBestTrack);
   fChain->SetBranchAddress("Mu_pxMuonBestTrack", &Mu_pxMuonBestTrack, &b_Mu_pxMuonBestTrack);
   fChain->SetBranchAddress("Mu_pyMuonBestTrack", &Mu_pyMuonBestTrack, &b_Mu_pyMuonBestTrack);
   fChain->SetBranchAddress("Mu_pzMuonBestTrack", &Mu_pzMuonBestTrack, &b_Mu_pzMuonBestTrack);
   fChain->SetBranchAddress("Mu_pMuonBestTrack", &Mu_pMuonBestTrack, &b_Mu_pMuonBestTrack);
   fChain->SetBranchAddress("Mu_etaMuonBestTrack", &Mu_etaMuonBestTrack, &b_Mu_etaMuonBestTrack);
   fChain->SetBranchAddress("Mu_phiMuonBestTrack", &Mu_phiMuonBestTrack, &b_Mu_phiMuonBestTrack);
   fChain->SetBranchAddress("Mu_thetaMuonBestTrack", &Mu_thetaMuonBestTrack, &b_Mu_thetaMuonBestTrack);
   fChain->SetBranchAddress("Mu_chargeMuonBestTrack", &Mu_chargeMuonBestTrack, &b_Mu_chargeMuonBestTrack);
   fChain->SetBranchAddress("Mu_absdxyMuonBestTrack", &Mu_absdxyMuonBestTrack, &b_Mu_absdxyMuonBestTrack);
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
   fChain->SetBranchAddress("NbOfDaughters", &NbOfDaughters, &b_NbOfDaughters);
   fChain->SetBranchAddress("McZmass", &McZmass, &b_McZmass);
   fChain->SetBranchAddress("nbPv", &nbPv, &b_nbPv);
   fChain->SetBranchAddress("Nbdof", &Nbdof, &b_Nbdof);
   fChain->SetBranchAddress("PositionRho", &PositionRho, &b_PositionRho);
   fChain->SetBranchAddress("PositionX", &PositionX, &b_PositionX);
   fChain->SetBranchAddress("PositionY", &PositionY, &b_PositionY);
   fChain->SetBranchAddress("PositionZ", &PositionZ, &b_PositionZ);
   fChain->SetBranchAddress("MuHLTMatch_nbMuonMatchHLT", &MuHLTMatch_nbMuonMatchHLT, &b_MuHLTMatch_nbMuonMatchHLT);
   fChain->SetBranchAddress("MuHLTMatch_pt", &MuHLTMatch_pt, &b_MuHLTMatch_pt);
   fChain->SetBranchAddress("MuHLTMatch_eta", &MuHLTMatch_eta, &b_MuHLTMatch_eta);
   fChain->SetBranchAddress("MuHLTMatch_phi", &MuHLTMatch_phi, &b_MuHLTMatch_phi);
   Notify();
}

Bool_t ZprimeMuMuPat::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void ZprimeMuMuPat::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t ZprimeMuMuPat::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef ZprimeMuMuPat_cxx
