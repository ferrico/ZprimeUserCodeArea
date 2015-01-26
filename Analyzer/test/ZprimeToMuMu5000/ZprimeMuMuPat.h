//=========================================================================  
//          Analysis code for Z' boson to Mu Mu analysis                  =  
//        In this code we select the high pt di-muons events              =
//                        To run over MC                                  =  
//               Author: Sherif Elgammal                                  =
//                                                                        =
//                     14/01/2015                                         =
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
   vector<float>   *Mu_ptBestTrack;
   vector<float>   *Mu_et;
   vector<float>   *Mu_en;
   vector<float>   *Mu_dPToverPTBestTrack;
   vector<float>   *Mu_absdxyBestTrack;
   vector<float>   *Mu_absdzBestTrack;
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
   vector<float>   *Mu_pxBestTrack;
   vector<float>   *Mu_pyBestTrack;
   vector<float>   *Mu_pzBestTrack;
   vector<float>   *Mu_pBestTrack;
   vector<float>   *Mu_etaBestTrack;
   vector<float>   *Mu_ptInnerTrack;
   vector<float>   *Mu_phiBestTrack;
   vector<float>   *Mu_thetaBestTrack;
   vector<float>   *Mu_chargeBestTrack;
   vector<float>   *Mu_pfDeltaBeta;
   vector<float>   *Mu_patDeltaBeta;
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
   TBranch        *b_Mu_nbMuon;   //!
   TBranch        *b_Mu_isGlobalMuon;   //!
   TBranch        *b_Mu_isPF;   //!
   TBranch        *b_Mu_ptBestTrack;   //!
   TBranch        *b_Mu_et;   //!
   TBranch        *b_Mu_en;   //!
   TBranch        *b_Mu_dPToverPTBestTrack;   //!
   TBranch        *b_Mu_absdxyBestTrack;   //!
   TBranch        *b_Mu_absdzBestTrack;   //!
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
   TBranch        *b_Mu_pxBestTrack;   //!
   TBranch        *b_Mu_pyBestTrack;   //!
   TBranch        *b_Mu_pzBestTrack;   //!
   TBranch        *b_Mu_pBestTrack;   //!
   TBranch        *b_Mu_etaBestTrack;   //!
   TBranch        *b_Mu_ptInnerTrack;   //!
   TBranch        *b_Mu_phiBestTrack;   //!
   TBranch        *b_Mu_thetaBestTrack;   //!
   TBranch        *b_Mu_chargeBestTrack;   //!
   TBranch        *b_Mu_pfDeltaBeta;   //!
   TBranch        *b_Mu_patDeltaBeta;   //!
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
   TBranch        *b_MuHLTObj_nbObj;   //!
   TBranch        *b_MuHLTObj_pt;   //!
   TBranch        *b_MuHLTObj_eta;   //!
   TBranch        *b_MuHLTObj_phi;   //!

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
			float &pmuon1,float &dxymuon1);
   bool SelectSecondMuon(int ChargeMu1,unsigned FlagMu1,float pTmuon1,float &pTmuon2,float &Enmuon2,
			 float &Etamuon2,float &Phimuon2,int &ChargeMu2,float &pxmuon2,
			 float &pymuon2,float &pzmuon2,float &pmuon2,float &dxymuon2);
   float Mass(float Pt1,float Eta1,float Phi1,float En1,
	      float Pt2,float Eta2,float Phi2,float En2);
   void PlotleadingAndSecondleading(float vtxMassMu);
   void PickThehighestMass(float &vtxHighestMass,float &vtxHighestChi2);
   double ThreeDangle(float pxMu1,float pyMu1,float pzMu1,float pMu1,
   		      float pxMu2,float pyMu2,float pzMu2,float pMu2);
   
   void PrintEventInformation(int runNumber, int lumiNumber, int eventNumber,
			      float vtxChi2, float vtxMass, float CosmicRejection);
   
   
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
   
   bool RecoHLMuonTMatching(float RecoEta,float RecoPhi);
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
   //================================================================================
   float RecoHLTMatchingDeltaRcut,deltaRcut,minMassCut,maxMassCut;
   float vtxChi2Mu,vtxMassMu;
   float mPtGen1,mPhiGen1,mEtaGen1,mEnGen1;
   unsigned mGenFlag1;
   float mPtGen2,mPhiGen2,mEtaGen2,mEnGen2;
   int ChargeRecMu1,ChargeRecMu2;
   unsigned flagmu1;
   unsigned flag1;
   float PtRecMu1,EnRecMu1,EtaRecMu1,PhiRecMu1;
   float PtRecMu2,EnRecMu2,EtaRecMu2,PhiRecMu2;
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
   TH1F* h1_RecoMassResolution_;
   ofstream output_txt;
};

#endif

#ifdef ZprimeMuMuPat_cxx
ZprimeMuMuPat::ZprimeMuMuPat(TTree *tree) : fChain(0) 
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
   Mu_ptBestTrack = 0;
   Mu_et = 0;
   Mu_en = 0;
   Mu_dPToverPTBestTrack = 0;
   Mu_absdxyBestTrack = 0;
   Mu_absdzBestTrack = 0;
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
   Mu_pxBestTrack = 0;
   Mu_pyBestTrack = 0;
   Mu_pzBestTrack = 0;
   Mu_pBestTrack = 0;
   Mu_etaBestTrack = 0;
   Mu_ptInnerTrack = 0;
   Mu_phiBestTrack = 0;
   Mu_thetaBestTrack = 0;
   Mu_chargeBestTrack = 0;
   Mu_pfDeltaBeta = 0;
   Mu_patDeltaBeta = 0;
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
   fChain->SetBranchAddress("Mu_nbMuon", &Mu_nbMuon, &b_Mu_nbMuon);
   fChain->SetBranchAddress("Mu_isGlobalMuon", &Mu_isGlobalMuon, &b_Mu_isGlobalMuon);
   fChain->SetBranchAddress("Mu_isPF", &Mu_isPF, &b_Mu_isPF);
   fChain->SetBranchAddress("Mu_ptBestTrack", &Mu_ptBestTrack, &b_Mu_ptBestTrack);
   fChain->SetBranchAddress("Mu_et", &Mu_et, &b_Mu_et);
   fChain->SetBranchAddress("Mu_en", &Mu_en, &b_Mu_en);
   fChain->SetBranchAddress("Mu_dPToverPTBestTrack", &Mu_dPToverPTBestTrack, &b_Mu_dPToverPTBestTrack);
   fChain->SetBranchAddress("Mu_absdxyBestTrack", &Mu_absdxyBestTrack, &b_Mu_absdxyBestTrack);
   fChain->SetBranchAddress("Mu_absdzBestTrack", &Mu_absdzBestTrack, &b_Mu_absdzBestTrack);
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
   fChain->SetBranchAddress("Mu_pxBestTrack", &Mu_pxBestTrack, &b_Mu_pxBestTrack);
   fChain->SetBranchAddress("Mu_pyBestTrack", &Mu_pyBestTrack, &b_Mu_pyBestTrack);
   fChain->SetBranchAddress("Mu_pzBestTrack", &Mu_pzBestTrack, &b_Mu_pzBestTrack);
   fChain->SetBranchAddress("Mu_pBestTrack", &Mu_pBestTrack, &b_Mu_pBestTrack);
   fChain->SetBranchAddress("Mu_etaBestTrack", &Mu_etaBestTrack, &b_Mu_etaBestTrack);
   fChain->SetBranchAddress("Mu_ptInnerTrack", &Mu_ptInnerTrack, &b_Mu_ptInnerTrack);
   fChain->SetBranchAddress("Mu_phiBestTrack", &Mu_phiBestTrack, &b_Mu_phiBestTrack);
   fChain->SetBranchAddress("Mu_thetaBestTrack", &Mu_thetaBestTrack, &b_Mu_thetaBestTrack);
   fChain->SetBranchAddress("Mu_chargeBestTrack", &Mu_chargeBestTrack, &b_Mu_chargeBestTrack);
   fChain->SetBranchAddress("Mu_pfDeltaBeta", &Mu_pfDeltaBeta, &b_Mu_pfDeltaBeta);
   fChain->SetBranchAddress("Mu_patDeltaBeta", &Mu_patDeltaBeta, &b_Mu_patDeltaBeta);
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
   fChain->SetBranchAddress("MuHLTObj_nbObj", &MuHLTObj_nbObj, &b_MuHLTObj_nbObj);
   fChain->SetBranchAddress("MuHLTObj_pt", &MuHLTObj_pt, &b_MuHLTObj_pt);
   fChain->SetBranchAddress("MuHLTObj_eta", &MuHLTObj_eta, &b_MuHLTObj_eta);
   fChain->SetBranchAddress("MuHLTObj_phi", &MuHLTObj_phi, &b_MuHLTObj_phi);
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
