//=========================================================================  
//          Analysis code for Z' boson to Mu Mu analysis                  =  
//    [1] In this code we select the high pt di-muons events              =
//                        To run over MC                                  =  
//             Written by Sherif Elgammal                                 =
//                                                                        =
//                      1/1/2015                                          =
//=========================================================================
#define ZprimeMuMuPat_cxx
#include "ZprimeMuMuPat.h"
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
#include <fstream>
#include <sstream>
using namespace std;
#define PI 3.14159265
bool myfunction (int i,int j) { return (i<j); }
bool picklargemass (float lhs,float rhs) { return (lhs > rhs); }
void ZprimeMuMuPat::Loop()
{
  time_t start,end;
  double dif;
  time (&start);
  FILE * pFile;
  pFile = fopen ("myfile.txt","w");
  //values needed for AccepXeff study
  NbGen = 0;
  NbReco= 0;
  int binMass   = 50;
  float minMass = 0.0;
  float maxMass = 8000.0;
  //values needed for tag & probe study
  MassCutMin = 0.0;
  MassCutMax = 2000.0;
  EtaCut = 2.4;
  MassResolution = 0.10;
  deltaRcut = 0.15;
  RecoHLTMatchingDeltaRcut = 0.20;
  minMassCut = 50.0;
  maxMassCut = 4000.0;
  //================================================================================== 
  //                                                                                 =
  //             Start the histograms for CollinSoper CMF                            =
  //                                                                                 =
  //==================================================================================
  int NbBins   = 10;
  float MinBin = -1.0;
  float MaxBin =  1.0;
  h1_CosAngleCollinSoperCorrect60Mass120_    = new TH1F("CosAngleCollinSoperCorrect60Mass120","",NbBins,MinBin,MaxBin);
  h1_CosAngleCollinSoperCorrect120Mass300_   = new TH1F("CosAngleCollinSoperCorrect120Mass300","",NbBins,MinBin,MaxBin);
  h1_CosAngleCollinSoperCorrect300Mass700_   = new TH1F("CosAngleCollinSoperCorrect300Mass700","",NbBins,MinBin,MaxBin);
  h1_CosAngleCollinSoperCorrect700Mass3000_  = new TH1F("CosAngleCollinSoperCorrect700Mass3000","",NbBins,MinBin,MaxBin);
  h1_CosAngleCollinSoperCorrect4900Mass5100_ = new TH1F("CosAngleCollinSoperCorrect4900Mass5100","",NbBins,MinBin,MaxBin); 
  //==================================================================================
  //================================================================================== 
  //                                                                                 =  
  //            Start a histograms for Mass resolution                               =  
  //                                                                                 =
  //==================================================================================
  h1_MassResultionEBEB1_      = new TH1F("MassResultionEBEB1","",100,-0.5,0.5);
  h1_MassResultionEBEB2_      = new TH1F("MassResultionEBEB2","",100,-0.5,0.5);
  h1_MassResultionEBEB3_      = new TH1F("MassResultionEBEB3","",100,-0.5,0.5);
  h1_MassResultionEBEB4_      = new TH1F("MassResultionEBEB4","",100,-0.5,0.5);
  h1_MassResultionEBEB5_      = new TH1F("MassResultionEBEB5","",100,-0.5,0.5);
  h1_MassResultionEBEB6_      = new TH1F("MassResultionEBEB6","",100,-0.5,0.5);
  h1_MassResultionEBEB7_      = new TH1F("MassResultionEBEB7","",100,-0.5,0.5);
  h1_MassResultionEBEB8_      = new TH1F("MassResultionEBEB8","",100,-0.5,0.5);
  h1_MassResultionEBEB9_      = new TH1F("MassResultionEBEB9","",100,-0.5,0.5);
  h1_MassResultionEBEB10_     = new TH1F("MassResultionEBEB10","",100,-0.5,0.5);
  //================================================================================== 
  //                                                                                 =
  //            Start the histograms for the mass of Z                               =  
  //                                                                                 =
  //==================================================================================
  h1_MassGenInAccep_             = new TH1F("MassGenInAccep","",58,0.0,5000.0);
  h1_MassRecoInAccep_            = new TH1F("MassRecoInAccep","",58,0.0,5000.0);
  h1_ZprimeRecomass_             = new TH1F("ZprimeRecomass","",binMass,minMass,maxMass);
  h1_ZprimeRecomassAbove1000GeV_ = new TH1F("ZprimeRecomassAbove1000GeV","",binMass,minMass,maxMass);
  h1_ZprimeGenmass_              = new TH1F("ZprimeGenmass","",binMass,minMass,maxMass);
  h1_ZprimeGenEta1_              = new TH1F("ZprimeGenEta1","",100,-8.0,8.0);
  h1_ZprimeGenEta2_              = new TH1F("ZprimeGenEta2","",100,-8.0,8.0);
  h1_ZprimeGenPt1_               = new TH1F("ZprimeGenPt1","",100,0.0,2000.0);
  h1_ZprimeGenPt2_               = new TH1F("ZprimeGenPt2","",100,0.0,2000.0);
  h1_ZprimeGenEn1_               = new TH1F("ZprimeGenEn1","",100,0.0,2000.0);
  h1_ZprimeGenEn2_               = new TH1F("ZprimeGenEn2","",100,0.0,2000.0);  
  h1_3Dangle_                    = new TH1F("3Dangle","",100,-2.0,2.0);
  h1_DxyDiff_                    = new TH1F("DxyDiff","",100,10.0,10.0);
  h1_MassRecoGenDif_             = new TH1F("MassRecoGenDif","",100,-0.5,0.5);
  h1_PtResolutionTunePMBT_       =  new TH1F("PtResolutionTunePMBT","",100,-0.5,0.5);
  h1_PtResolutiontuneP_          =  new TH1F("PtResolutiontuneP","",100,-0.5,0.5);
  h1_PtResolutionMBT_            =  new TH1F("PtResolutionMBT","",100,-0.5,0.5);
  //================================================================================== 
  //                                                                                 =
  //                 Start the histograms for N-1 dist                               =  
  //                                                                                 =
  //==================================================================================
  h1_dPToverPT_                            = new TH1F("dPToverPT","",100,0.0,0.5);
  h1_normalizedChi2_                       = new TH1F("normalizedChi2","",100,0.0,20.0);
  h1_numberOftrackerLayersWithMeasurement_ = new TH1F("numberOftrackerLayersWithMeasurement","",20,0.0,20.0);
  h1_numberOfValidPixelHits_               = new TH1F("numberOfValidPixelHits","",10,0.0,10.0);
  h1_numberOfValidMuonHits_                = new TH1F("numberOfValidMuonHits","",60,0.0,60.0);
  h1_numberOfMatchedStations_              = new TH1F("numberOfMatchedStations","",10,0.0,10.0);
  h1_trackiso_                             = new TH1F("trackiso","",50,0.0,0.3);
  h1_absdxy_                               = new TH1F("absdxy","",100,0.0,0.3);
  h1_PtEffpterror_                   = new TH1F("PtEffpterror","",20,0.0,3000.0);
  h1_PtEffptnumberOftrackerLayers_   = new TH1F("PtEffptnumberOftrackerLayers","",20,0.0,3000.0);
  h1_PtEffptnumberOfPixelHits_       = new TH1F("PtEffptnumberOfPixelHits","",20,0.0,3000.0);
  h1_PtEffptnumberOfMuonHits_        = new TH1F("PtEffptnumberOfMuonHits","",20,0.0,3000.0);
  h1_PtEffptnumberOfMatchedStations_ = new TH1F("PtEffptnumberOfMatchedStations","",20,0.0,3000.0);
  h1_PtEffptTrackIso_                = new TH1F("PtEffptTrackIso","",20,0.0,3000.0);
  h1_PtEffptabsdsy_                  = new TH1F("PtEffptabsdsy","",20,0.0,3000.0);
  h1_PtEffpfSumChargedHadron_        = new TH1F("PtEffpfSumChargedHadron","",20,0.0,3000.0);
  h1_PtEffpfSumNeutralHadron_        = new TH1F("PtEffpfSumNeutralHadron","",20,0.0,3000.0);
  h1_PtEffpfPhotonIso_               = new TH1F("PtEffpfPhotonIso","",20,0.0,3000.0);
  h1_EtaEffpterror_                  = new TH1F("EtaEffpterror","",30,0.0,3.0);
  h1_EtaEffptnumberOftrackerLayers_  = new TH1F("EtaEffptnumberOftrackerLayers","",30,0.0,3.0);
  h1_EtaEffptnumberOfPixelHits_      = new TH1F("EtaEffptnumberOfPixelHits","",30,0.0,3.0);
  h1_EtaEffptnumberOfMuonHits_       = new TH1F("EtaEffptnumberOfMuonHits","",30,0.0,3.0);
  h1_EtaEffptnumberOfMatchedStations_= new TH1F("EtaEffptnumberOfMatchedStations","",30,0.0,3.0);
  h1_EtaEffptTrackIso_               = new TH1F("EtaEffptTrackIso","",30,0.0,3.0);
  h1_EtaEffptabsdsy_                 = new TH1F("EtaEffptabsdsy","",30,0.0,3.0);
  h1_nbPVID_               = new TH1F("nbPVID","",50,0.0,50.0);
  h1_PtID_                 = new TH1F("PtID","",20,0.0,3000.0);
  h1_EtaID_                = new TH1F("EtaID","",30,0.0,3.0);
  h1_nbPVNewID_            = new TH1F("nbPVNewID","",50,0.0,50.0);
  h1_PtNewID_              = new TH1F("PtNewID","",20,0.0,3000.0);
  h1_EtaNewID_             = new TH1F("EtaNewID","",30,0.0,3.0);
  h1_nbPVTightID_          = new TH1F("nbPVTightID","",50,0.0,50.0);
  h1_PtTightID_            = new TH1F("PtTightID","",20,0.0,3000.0);
  h1_EtaTightID_           = new TH1F("EtaTightID","",30,0.0,3.0);
  //----------------------------------------------------------------------------------------------
  // Book txt file for candidate events
  Char_t txtOUT[500];
  //sprintf(txtOUT,"%s_txt.txt",datasetName.Data());
  sprintf(txtOUT,"outputD.txt");
  cout << "Opening a txt file with candidate events " << txtOUT << endl;
  //ofstream output_txt;
  output_txt.open(txtOUT);
  output_txt << "CANDIDATES Events: \n" 
  << "run: " 
  << "    lumi: " 
  << "    event: " 
  << "          dil_mass: \n \n" << endl;
  //==================================================================================  
  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    /*
      cout<<"=======> jentry = "<<jentry<< 
      "=======> Evt = "<<event_evtNo<< 
      "=======> Run = "<<event_runNo<< 
      "=======> Lumi = "<<event_lumi<< 
      "=======> bunch = "<<event_bunch<<endl;
    */
    //=====================================================
    //                                                    =
    //  Calling methods to get the generated information  =
    //                                                    =
    //=====================================================
    bool firstGenMu  = SelectFirstGenMu(genET1,genPhi1,genEta1,genEn1,genID1,genStat1,flag1);
    bool secondGenMu = SelectSecondGenMu(flag1,genET1,genET2,genPhi2,genEta2,genEn2,genID2,genStat2);
    if(firstGenMu == 0 || secondGenMu == 0) continue;
    MassGen = Mass(genET1,genEta1,genPhi1,genEn1,genET2,genEta2,genPhi2,genEn2);
    PlotGenInfo(MassGen,genEta1,genEta2,genET1,genET2,genEn1,genEn2);
    NbGen ++;
    //=========================================================
    //                                                        =
    // Calling methods to get events with 2 muons passing ID  =   
    //                                                        =
    //=========================================================
    bool firstMuFinal  = SelectFirstMuon(PtRecTunePMuBestTrack1,EnRecMu1,EtaRecMu1,PhiRecMu1,ChargeRecMu1,flagmu1,
					 pxRecMu1,pyRecMu1,pzRecMu1,pRecMu1,dxyRecMu1,PtRecTunePMu1,
					 PtRecMuBestTrack1);
    bool secondMuFinal = SelectSecondMuon(ChargeRecMu1,flagmu1,PtRecTunePMuBestTrack1,PtRecTunePMuBestTrack2,EnRecMu2,
					  EtaRecMu2,PhiRecMu2,ChargeRecMu2,pxRecMu2,pyRecMu2,pzRecMu2,pRecMu2,dxyRecMu2,
					  PtRecTunePMu2,PtRecMuBestTrack2);
    PickThehighestMass(vtxMassMu,vtxChi2Mu);
    double CosmicRejec = ThreeDangle(pxRecMu1,pyRecMu1,pzRecMu1,pRecMu1,
				     pxRecMu2,pyRecMu2,pzRecMu2,pRecMu2);
    //=========================================================
    //        call the method for N-1 plots                   =
    //                                                        =
    //=========================================================
    plotAllHighPtMuonsID();
    if(firstMuFinal == 0 || secondMuFinal == 0) continue; 
    //=========================================================
    //        start doing matching between reco & gen         =
    //                                                        =
    //=========================================================
    bool GenRecoMatch1 = GenRecoMatchMu1(EtaRecMu1,PhiRecMu1,mPtGen1,mPhiGen1,mEtaGen1,mEnGen1,mGenFlag1);
    bool GenRecoMatch2 = GenRecoMatchMu2(mGenFlag1,EtaRecMu2,PhiRecMu2,mPtGen2,mPhiGen2,mEtaGen2,mEnGen2);
    if(	GenRecoMatch1 == 0 || GenRecoMatch2 == 0 ) continue; 
    if(vtxChi2Mu<=10.0 && CosmicRejec>=-0.9998) 
      {
	if( (fabs(EtaRecMu1) < 2.1 && fabs(EtaRecMu2) < EtaCut) || (fabs(EtaRecMu2) < 2.1 && fabs(EtaRecMu1) < EtaCut) )
	  {  
	    PlotRecoInfo(CosmicRejec,vtxMassMu,MassGen,PtRecTunePMuBestTrack1,PtRecTunePMu1,PtRecMuBestTrack1,mPtGen1);
	    CosThetaCollinSoper(PtRecTunePMuBestTrack1,EtaRecMu1,PhiRecMu1,EnRecMu1,
                                PtRecTunePMuBestTrack2,EtaRecMu2,PhiRecMu2,EnRecMu2,
				ChargeRecMu1,vtxMassMu);
	    if( (vtxMassMu-MassGen)/MassGen < MassResolution ) {
	      NbReco++;
	    }
	  }
      }
  }
  //===================================================
  //                                                  =
  //          Compute Acceptance X Eff.               =
  //                                                  =
  //===================================================
  float Accxeff = float(NbReco)/NbGen;
  cout<<"Accxeff = "<<Accxeff<< endl;
  //===================================================
  fclose (pFile);
  //================================================================================== 
  //                                                                                 =
  //                         writing Histograms to a file                            =             
  //                                                                                 =
  //==================================================================================
  TFile *output = new TFile("ZprimetoMuMu-MC-Mass5000-CMSSW720.root","recreate");
  output->cd();
  /** Write the histograms */
  h1_MassGenInAccep_->Write();
  h1_MassRecoInAccep_->Write();
  h1_ZprimeRecomass_->Write();
  h1_3Dangle_->Write();
  h1_ZprimeGenmass_->Write();
  h1_ZprimeGenEta1_->Write();
  h1_ZprimeGenEta2_->Write();
  h1_ZprimeGenPt1_->Write();
  h1_ZprimeGenPt2_->Write();
  h1_ZprimeGenEn1_->Write();
  h1_ZprimeGenEn2_->Write();
  h1_dPToverPT_->Write();
  h1_normalizedChi2_->Write();
  h1_numberOftrackerLayersWithMeasurement_->Write();
  h1_numberOfValidPixelHits_->Write();
  h1_numberOfValidMuonHits_->Write();   
  h1_numberOfMatchedStations_->Write(); 
  h1_trackiso_->Write();                
  h1_absdxy_->Write();                  
  h1_PtEffpterror_->Write();                   
  h1_PtEffptnumberOftrackerLayers_->Write();   
  h1_PtEffptnumberOfPixelHits_->Write();       
  h1_PtEffptnumberOfMuonHits_->Write();        
  h1_PtEffptnumberOfMatchedStations_->Write(); 
  h1_PtEffptTrackIso_->Write();                
  h1_PtEffptabsdsy_->Write();                  
  h1_EtaEffpterror_->Write();                  
  h1_EtaEffptnumberOftrackerLayers_->Write();  
  h1_EtaEffptnumberOfPixelHits_->Write();      
  h1_EtaEffptnumberOfMuonHits_->Write();       
  h1_EtaEffptnumberOfMatchedStations_->Write();
  h1_EtaEffptTrackIso_->Write();               
  h1_EtaEffptabsdsy_->Write();
  h1_nbPVID_->Write();           
  h1_PtID_->Write();             
  h1_EtaID_->Write();            
  h1_nbPVNewID_->Write();            
  h1_PtNewID_->Write();              
  h1_EtaNewID_->Write();
  h1_nbPVTightID_->Write(); 
  h1_PtTightID_->Write();         
  h1_EtaTightID_->Write();
  h1_PtResolutionTunePMBT_->Write();
  h1_PtResolutiontuneP_->Write();
  h1_PtResolutionMBT_->Write();
  h1_MassResultionEBEB1_->Write();
  h1_MassResultionEBEB2_->Write();
  h1_MassResultionEBEB3_->Write();
  h1_MassResultionEBEB4_->Write();
  h1_MassResultionEBEB5_->Write();
  h1_MassResultionEBEB6_->Write();
  h1_MassResultionEBEB7_->Write();
  //h1_MassResultionEBEB8_->Write();
  //h1_MassResultionEBEB9_->Write();
  //h1_MassResultionEBEB10_->Write();
  h1_CosAngleCollinSoperCorrect60Mass120_->Write();
  h1_CosAngleCollinSoperCorrect120Mass300_->Write();
  h1_CosAngleCollinSoperCorrect300Mass700_->Write();
  h1_CosAngleCollinSoperCorrect700Mass3000_->Write();
  h1_CosAngleCollinSoperCorrect4900Mass5100_->Write();
  output->Close();
  output_txt.close();
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
//============================ Method to select first high pt muon ========================
bool ZprimeMuMuPat::SelectFirstMuon(float &pTmuon1,float &Enmuon1,float &Etamuon1,
				    float &Phimuon1,int &ChargeMu1,unsigned &FlagMu1,
				    float &pxmuon1,float &pymuon1,float &pzmuon1,
				    float &pmuon1,float &dxymuon1,float &pTmuon1tuneP,
				    float &pTmuonBestTrack1){
  int NbHEEPele = 0;
  unsigned iflag = -10;
  for(unsigned i=0; i<Mu_nbMuon->size(); i++){
    if( Mu_isGlobalMuon->at(i) == 1 &&
	Mu_ptTunePMuonBestTrack->at(i) > 45.0 &&
	Mu_absdxyTunePMuonBestTrack->at(i) < 0.2 &&
	(Mu_trackiso->at(i)/Mu_ptInnerTrack->at(i)) < 0.10  &&
	Mu_numberOftrackerLayersWithMeasurement->at(i) > 5 && 
	Mu_numberOfValidPixelHits->at(i) > 0 &&
	Mu_numberOfValidMuonHits->at(i) > 0 &&
	Mu_numberOfMatchedStations->at(i) > 1 &&
	Mu_dPToverPTTunePMuonBestTrack->at(i) < 0.3 ) { 
      iflag  = i;
      NbHEEPele ++;
    }
    else continue;
  }
  if( NbHEEPele > 0 ){
    FlagMu1             = iflag;
    pTmuon1             = Mu_ptTunePMuonBestTrack->at(iflag);
    Enmuon1             = Mu_en->at(iflag);
    Etamuon1            = Mu_etaTunePMuonBestTrack->at(iflag);
    Phimuon1            = Mu_phiTunePMuonBestTrack->at(iflag);
    ChargeMu1           = Mu_chargeTunePMuonBestTrack->at(iflag);
    pxmuon1             = Mu_pxTunePMuonBestTrack->at(iflag);
    pymuon1             = Mu_pyTunePMuonBestTrack->at(iflag);
    pzmuon1             = Mu_pzTunePMuonBestTrack->at(iflag);
    pmuon1              = Mu_pTunePMuonBestTrack->at(iflag);
    dxymuon1            = Mu_absdxyTunePMuonBestTrack->at(iflag);
    pTmuon1tuneP        = Mu_ptcocktail->at(iflag);
    pTmuonBestTrack1    = Mu_ptMuonBestTrack->at(iflag);
    return true;  
  }
  else return false;
}
//============================ Method to select second high pt muon ========================
bool ZprimeMuMuPat::SelectSecondMuon(int ChargeMu1,unsigned FlagMu1,float pTmuon1,float &pTmuon2,float &Enmuon2,
				     float &Etamuon2,float &Phimuon2,int &ChargeMu2,float &pxmuon2,
				     float &pymuon2,float &pzmuon2,float &pmuon2,float &dxymuon2,
				     float &pTmuon2tuneP,float &pTmuonBestTrack2)
{
  int NbHEEPele = 0;
  unsigned iflag = -10;
  for(unsigned i=0; i<Mu_nbMuon->size(); i++){
    if(i == FlagMu1) continue;
    if(Mu_ptTunePMuonBestTrack->at(i) == pTmuon1) continue;
    if(ChargeMu1*Mu_chargeTunePMuonBestTrack->at(i)>0) continue;
    if( Mu_isGlobalMuon->at(i) == 1 &&
	Mu_ptTunePMuonBestTrack->at(i) > 45.0 &&
	Mu_absdxyTunePMuonBestTrack->at(i) < 0.2 &&
	(Mu_trackiso->at(i)/Mu_ptInnerTrack->at(i)) < 0.10  &&
	Mu_numberOftrackerLayersWithMeasurement->at(i) > 5 && 
	Mu_numberOfValidPixelHits->at(i) > 0 &&
	Mu_numberOfValidMuonHits->at(i) > 0 &&
	Mu_numberOfMatchedStations->at(i) > 1 &&
	Mu_dPToverPTTunePMuonBestTrack->at(i) < 0.3 ) { 
      iflag  = i;
      NbHEEPele ++;
    }
    else continue;
  }
  if( NbHEEPele > 0 ){
    pTmuon2          = Mu_ptTunePMuonBestTrack->at(iflag);
    Enmuon2          = Mu_en->at(iflag);
    Etamuon2         = Mu_etaTunePMuonBestTrack->at(iflag);
    Phimuon2         = Mu_phiTunePMuonBestTrack->at(iflag);
    ChargeMu2        = Mu_chargeTunePMuonBestTrack->at(iflag);
    pxmuon2          = Mu_pxTunePMuonBestTrack->at(iflag);
    pymuon2          = Mu_pyTunePMuonBestTrack->at(iflag);
    pzmuon2          = Mu_pzTunePMuonBestTrack->at(iflag);
    pmuon2           = Mu_pTunePMuonBestTrack->at(iflag);
    dxymuon2         = Mu_absdxyTunePMuonBestTrack->at(iflag);
    pTmuon2tuneP     = Mu_ptcocktail->at(iflag);
    pTmuonBestTrack2 = Mu_ptMuonBestTrack->at(iflag);
    return true;  
  }
  else return false;
}
void ZprimeMuMuPat::PlotRecoInfo(float CosmicMuonRejec, float vertexMassMu,float MassGenerated,
				 float PtTunePMuBestTrack,float PtTunePMu,float PtMuBestTrack,
				 float PtGenerated){
  h1_MassRecoInAccep_->Fill(MassGenerated);
  h1_ZprimeRecomass_->Fill(vertexMassMu);
  h1_3Dangle_->Fill(CosmicMuonRejec);
  //part for Pt resolution
  h1_PtResolutionTunePMBT_->Fill((PtTunePMuBestTrack-PtGenerated)/PtGenerated);
  h1_PtResolutiontuneP_->Fill((PtTunePMu-PtGenerated)/PtGenerated);
  h1_PtResolutionMBT_->Fill((PtMuBestTrack-PtGenerated)/PtGenerated);
  //part for mass resolution
  if( vertexMassMu > 0.0 && vertexMassMu < 250.0 ){h1_MassResultionEBEB1_->Fill((vertexMassMu-MassGenerated)/MassGenerated);}
  if( vertexMassMu > 250 && vertexMassMu < 750.0 ){h1_MassResultionEBEB2_->Fill((vertexMassMu-MassGenerated)/MassGenerated);}
  if( vertexMassMu > 750 && vertexMassMu < 1250.0 ){h1_MassResultionEBEB3_->Fill((vertexMassMu-MassGenerated)/MassGenerated);}
  if( vertexMassMu > 1250 && vertexMassMu < 1750.0 ){h1_MassResultionEBEB4_->Fill((vertexMassMu-MassGenerated)/MassGenerated);}
  if( vertexMassMu > 1750 && vertexMassMu < 2250.0 ){h1_MassResultionEBEB5_->Fill((vertexMassMu-MassGenerated)/MassGenerated);}
  if( vertexMassMu > 2000 && vertexMassMu < 4000.0 ){h1_MassResultionEBEB6_->Fill((vertexMassMu-MassGenerated)/MassGenerated);}
  if( vertexMassMu > 4000 && vertexMassMu < 6000.0 ){h1_MassResultionEBEB7_->Fill((vertexMassMu-MassGenerated)/MassGenerated);} 
  //if( vertexMassMu > 3250 && vertexMassMu < 3750.0 ){h1_MassResultionEBEB8_->Fill((vertexMassMu-MassGenerated)/MassGenerated);}
  //if( vertexMassMu > 3750 && vertexMassMu < 4250.0 ){h1_MassResultionEBEB9_->Fill((vertexMassMu-MassGenerated)/MassGenerated);}
}
//===================== Methode to calculate the mass ========================
float ZprimeMuMuPat::Mass(float Pt1,float Eta1,float Phi1,float En1,
		       float Pt2,float Eta2,float Phi2,float En2){
  float MuMuMass = 0.0;
  TLorentzVector Mu1;
  TLorentzVector Mu2;
  Mu1.SetPtEtaPhiE(Pt1,Eta1,Phi1,En1);
  Mu2.SetPtEtaPhiE(Pt2,Eta2,Phi2,En2);
  MuMuMass = (Mu1 + Mu2).M();
  return MuMuMass;
}

void ZprimeMuMuPat::PickThehighestMass(float &vtxHighestMass,float &vtxHighestChi2)
{
  float Massinv = -10.0;
  unsigned iflag = -10;
  int NbMu = 0;
  for(unsigned i=0; i<Mu_vtxMass->size(); i++)
    {
      if(Mu_vtxMass->at(i)>Massinv){
	Massinv = Mu_vtxMass->at(i);
	iflag  = i;
	NbMu++;
      }
    }
  if( NbMu > 0 ){
    vtxHighestMass = Mu_vtxMass->at(iflag);
    vtxHighestChi2 = Mu_vtxNormChi2->at(iflag);
  }
}
double ZprimeMuMuPat::ThreeDangle(float pxMu1,float pyMu1,float pzMu1,float pMu1,
			       float pxMu2,float pyMu2,float pzMu2,float pMu2)
{
  TVector3 Mu1(pxMu1,pyMu1,pzMu1);
  TVector3 Mu2(pxMu2,pyMu2,pzMu2);
  double cos_angle = Mu1.Dot(Mu2) / pMu1 / pMu2;
  return cos_angle;
}

//----------------------------------------------------
//                                                   -
//       Part for Gen & Reco Matching                -
//                                                   -  
//----------------------------------------------------
//========================== Method to select firt Gen Mu =======================
bool ZprimeMuMuPat::SelectFirstGenMu(float &ETMu1,float &PhiSCMu1,
				  float &EtaSCMu1,float &EnMu1,
				  int &IDele1,int &Statele1,
				  unsigned &GenFlag1){
  int NbHEEPele = 0;
  int iflag = -10;
  ETMu1 = 0.0;
  for(unsigned i=0; i<iGen->size(); i++){
    //printf ("BS[%d] ID = %d status = %d \n",i,idGen->at(i),statusGen->at(i));
    if( fabs(idGen->at(i)) != 13 ) continue;
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
bool ZprimeMuMuPat::SelectSecondGenMu(unsigned GenFlag1,float ETMu1,float &ETMu2,float &PhiSCMu2,
				   float &EtaSCMu2,float &EnMu2,int &IDele2,int &Statele2){
  int NbHEEPele = 0;
  int iflag = -10;
  ETMu2 = 0.0;
  for(unsigned i=0; i<iGen->size(); i++){
    if( fabs(idGen->at(i)) != 13 ) continue;
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

//----------------------------------------------------
//                                                   -
//       Part for Gen & Reco Matching                -
//                                                   -  
//----------------------------------------------------
bool ZprimeMuMuPat::GenRecoMatchMu1(float RecoEta1,float RecoPhi1,
				 float &ETMu1,float &PhiSCMu1,
				 float &EtaSCMu1,float &EnMu1,
				 unsigned &GenFlag1){
  int NbHEEPele = 0;
  unsigned iflag = -10;
  for(unsigned i=0; i<iGen->size(); i++){
    //printf ("BS[%d] genEta = %f  genPhi = %f recoEta = %f recoPhi = %f deltaR = %f\n",
    //    i,etaGen->at(i),phiGen->at(i),RecoEta1,RecoPhi1,deltaR);
    if( fabs(idGen->at(i)) != 13 ) continue;
    if( statusGen->at(i) != 1 )  continue;
    float deltaEta = RecoEta1 - etaGen->at(i);
    float deltaPhi = RecoPhi1 - phiGen->at(i);
    float deltaR   = sqrt(pow(deltaEta,2)+pow(deltaPhi,2));
    if(fabs(deltaR)>deltaRcut) continue;
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

bool ZprimeMuMuPat::GenRecoMatchMu2(unsigned GenFlag1,float RecoEta2,float RecoPhi2,
				 float &ETMu2,float &PhiSCMu2,
				 float &EtaSCMu2,float &EnMu2){
  int NbHEEPele = 0;
  unsigned iflag = -10;
  for(unsigned i=0; i<iGen->size(); i++){
    //printf ("BS[%d] genEta = %f  genPhi = %f recoEta = %f recoPhi = %f deltaR = %f\n",
    //    i,etaGen->at(i),phiGen->at(i),RecoEta1,RecoPhi1,deltaR);
    if( fabs(idGen->at(i)) != 13 ) continue;
    if( statusGen->at(i) != 1 )  continue;
    if( i ==  GenFlag1)  continue;
    float deltaEta2 = RecoEta2 - etaGen->at(i);
    float deltaPhi2 = RecoPhi2 - phiGen->at(i);
    float deltaR2   = sqrt(pow(deltaEta2,2)+pow(deltaPhi2,2));
    if(fabs(deltaR2)>deltaRcut) continue;
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
//============================ Method to plot Gen Mu ========================
void ZprimeMuMuPat::PlotGenInfo(float ZprimeGenMass,float EtaGenMu1,float EtaGenMu2,float PtGenMu1,
				float PtGenMu2,float EnGenMu1,float EnGenMu2){
  h1_MassGenInAccep_->Fill(ZprimeGenMass);
  h1_ZprimeGenmass_->Fill(ZprimeGenMass);
  h1_ZprimeGenEta1_->Fill(EtaGenMu1);
  h1_ZprimeGenEta2_->Fill(EtaGenMu2);
  h1_ZprimeGenPt1_->Fill(PtGenMu1);
  h1_ZprimeGenPt2_->Fill(PtGenMu2);
  h1_ZprimeGenEn1_->Fill(EnGenMu1);
  h1_ZprimeGenEn2_->Fill(EnGenMu2);
}
//----------------------------------------------------
//                                                   -
//                           N-1 eff.                -
//                                                   -  
//----------------------------------------------------
void ZprimeMuMuPat::MuonPassingID(){
  for(unsigned i=0; i<Mu_nbMuon->size(); i++){
    if(fabs(Mu_etaTunePMuonBestTrack->at(i)) < 2.4 &&
       Mu_isGlobalMuon->at(i) == 1 &&
       (Mu_ptTunePMuonBestTrack->at(i) > 45.0 && Mu_ptTunePMuonBestTrack->at(i) < 3000.0) &&
       Mu_absdxyTunePMuonBestTrack->at(i) < 0.2 &&
       (Mu_trackiso->at(i)/Mu_ptInnerTrack->at(i)) < 0.10  &&
       Mu_numberOftrackerLayersWithMeasurement->at(i) > 5 && 
       Mu_numberOfValidPixelHits->at(i) > 0 &&
       Mu_numberOfValidMuonHits->at(i) > 0 &&
       Mu_numberOfMatchedStations->at(i) > 1 &&
       Mu_dPToverPTTunePMuonBestTrack->at(i) < 0.3 ){
      h1_PtID_->Fill(Mu_ptTunePMuonBestTrack->at(i));
      h1_EtaID_->Fill(fabs(Mu_etaTunePMuonBestTrack->at(i)));
      h1_nbPVID_->Fill(Mu_nbofpv->at(i));
    }
    else continue;
  }
}

void ZprimeMuMuPat::PlotPterror(){
  for(unsigned i=0; i<Mu_nbMuon->size(); i++){
    if(fabs(Mu_etaTunePMuonBestTrack->at(i)) < 2.4 &&
       Mu_isGlobalMuon->at(i) == 1 &&
       (Mu_ptTunePMuonBestTrack->at(i) > 45.0 && Mu_ptTunePMuonBestTrack->at(i) < 3000.0) &&
       Mu_absdxyTunePMuonBestTrack->at(i) < 0.2 &&
       (Mu_trackiso->at(i)/Mu_ptInnerTrack->at(i)) < 0.10  &&
       Mu_numberOftrackerLayersWithMeasurement->at(i) > 5 && 
       Mu_numberOfValidPixelHits->at(i) > 0 &&
       Mu_numberOfValidMuonHits->at(i) > 0 &&
       Mu_numberOfMatchedStations->at(i) > 1){
      h1_dPToverPT_->Fill(Mu_dPToverPTTunePMuonBestTrack->at(i) );
      h1_PtEffpterror_->Fill(Mu_ptTunePMuonBestTrack->at(i));
      h1_EtaEffpterror_->Fill(fabs(Mu_etaTunePMuonBestTrack->at(i)));
    }
    else continue;
  }
}

void ZprimeMuMuPat::PlotNbTrackLayers(){
  for(unsigned i=0; i<Mu_nbMuon->size(); i++){
    if(fabs(Mu_etaTunePMuonBestTrack->at(i)) < 2.4 &&
       Mu_isGlobalMuon->at(i) == 1 &&
       (Mu_ptTunePMuonBestTrack->at(i) > 45.0 && Mu_ptTunePMuonBestTrack->at(i) < 3000.0) &&
       Mu_absdxyTunePMuonBestTrack->at(i) < 0.2 &&
       (Mu_trackiso->at(i)/Mu_ptInnerTrack->at(i)) < 0.10  &&
       //Mu_numberOftrackerLayersWithMeasurement->at(i) > 5 && 
       Mu_numberOfValidPixelHits->at(i) > 0 &&
       Mu_numberOfValidMuonHits->at(i) > 0 &&
       Mu_numberOfMatchedStations->at(i) > 1 &&
       Mu_dPToverPTTunePMuonBestTrack->at(i) < 0.3 ){
      h1_numberOftrackerLayersWithMeasurement_->Fill( Mu_numberOftrackerLayersWithMeasurement->at(i) );
      h1_PtEffptnumberOftrackerLayers_->Fill(Mu_ptTunePMuonBestTrack->at(i));
      h1_EtaEffptnumberOftrackerLayers_->Fill(fabs(Mu_etaTunePMuonBestTrack->at(i)));
    }
    else continue;
  }
}


void ZprimeMuMuPat::PlotNBValidPixelHits(){
  for(unsigned i=0; i<Mu_nbMuon->size(); i++){
    if(fabs(Mu_etaTunePMuonBestTrack->at(i)) < 2.4 &&
       Mu_isGlobalMuon->at(i) == 1 &&
       (Mu_ptTunePMuonBestTrack->at(i) > 45.0 && Mu_ptTunePMuonBestTrack->at(i) < 3000.0) &&
       Mu_absdxyTunePMuonBestTrack->at(i) < 0.2 &&
       (Mu_trackiso->at(i)/Mu_ptInnerTrack->at(i)) < 0.10  &&
       Mu_numberOftrackerLayersWithMeasurement->at(i) > 5 && 
       //Mu_numberOfValidPixelHits->at(i) > 0 &&
       Mu_numberOfValidMuonHits->at(i) > 0 &&
       Mu_numberOfMatchedStations->at(i) > 1 &&
       Mu_dPToverPTTunePMuonBestTrack->at(i) < 0.3 ){
      h1_numberOfValidPixelHits_->Fill( Mu_numberOfValidPixelHits->at(i) );
      h1_PtEffptnumberOfPixelHits_->Fill(Mu_ptTunePMuonBestTrack->at(i));
      h1_EtaEffptnumberOfPixelHits_->Fill(fabs(Mu_etaTunePMuonBestTrack->at(i)));
    }
    else continue;
  }
}

void ZprimeMuMuPat::PlotNbValidMuonHits(){
  for(unsigned i=0; i<Mu_nbMuon->size(); i++){
    if(fabs(Mu_etaTunePMuonBestTrack->at(i)) < 2.4 &&
       Mu_isGlobalMuon->at(i) == 1 &&
       (Mu_ptTunePMuonBestTrack->at(i) > 45.0 && Mu_ptTunePMuonBestTrack->at(i) < 3000.0) &&
       Mu_absdxyTunePMuonBestTrack->at(i) < 0.2 &&
       (Mu_trackiso->at(i)/Mu_ptInnerTrack->at(i)) < 0.10  &&
       Mu_numberOftrackerLayersWithMeasurement->at(i) > 5 && 
       Mu_numberOfValidPixelHits->at(i) > 0 &&
       //Mu_numberOfValidMuonHits->at(i) > 0 &&
       Mu_numberOfMatchedStations->at(i) > 1 &&
       Mu_dPToverPTTunePMuonBestTrack->at(i) < 0.3 ){
      h1_numberOfValidMuonHits_->Fill( Mu_numberOfValidMuonHits->at(i) );
      h1_PtEffptnumberOfMuonHits_->Fill(Mu_ptTunePMuonBestTrack->at(i));
      h1_EtaEffptnumberOfMuonHits_->Fill(fabs(Mu_etaTunePMuonBestTrack->at(i)));
    }
    else continue;
  }
}


void ZprimeMuMuPat::PlotNbMatchedStations(){
  for(unsigned i=0; i<Mu_nbMuon->size(); i++){
    if(fabs(Mu_etaTunePMuonBestTrack->at(i)) < 2.4 &&
       Mu_isGlobalMuon->at(i) == 1 &&
       (Mu_ptTunePMuonBestTrack->at(i) > 45.0 && Mu_ptTunePMuonBestTrack->at(i) < 3000.0) &&
       Mu_absdxyTunePMuonBestTrack->at(i) < 0.2 &&
       (Mu_trackiso->at(i)/Mu_ptInnerTrack->at(i)) < 0.10  &&
       Mu_numberOftrackerLayersWithMeasurement->at(i) > 5 && 
       Mu_numberOfValidPixelHits->at(i) > 0 &&
       Mu_numberOfValidMuonHits->at(i) > 0 &&
       //Mu_numberOfMatchedStations->at(i) > 1 &&
       Mu_dPToverPTTunePMuonBestTrack->at(i) < 0.3 ){
      h1_numberOfMatchedStations_->Fill( Mu_numberOfMatchedStations->at(i) );
      h1_PtEffptnumberOfMatchedStations_->Fill(Mu_ptTunePMuonBestTrack->at(i));
      h1_EtaEffptnumberOfMatchedStations_->Fill(fabs(Mu_etaTunePMuonBestTrack->at(i)));
    }
    else continue;
  }
}


void ZprimeMuMuPat::PlotTrackiso(){
  for(unsigned i=0; i<Mu_nbMuon->size(); i++){
    if(fabs(Mu_etaTunePMuonBestTrack->at(i)) < 2.4 &&
       Mu_isGlobalMuon->at(i) == 1 &&
       (Mu_ptTunePMuonBestTrack->at(i) > 45.0 && Mu_ptTunePMuonBestTrack->at(i) < 3000.0) &&
       Mu_absdxyTunePMuonBestTrack->at(i) < 0.2 &&
       //(Mu_trackiso->at(i)/Mu_ptInnerTrack->at(i)) < 0.10  &&
       Mu_numberOftrackerLayersWithMeasurement->at(i) > 5 && 
       Mu_numberOfValidPixelHits->at(i) > 0 &&
       Mu_numberOfValidMuonHits->at(i) > 0 &&
       Mu_numberOfMatchedStations->at(i) > 1 &&
       Mu_dPToverPTTunePMuonBestTrack->at(i) < 0.3 ){
      h1_trackiso_->Fill( Mu_trackiso->at(i)/Mu_ptTunePMuonBestTrack->at(i) );
      h1_PtEffptTrackIso_->Fill(Mu_ptTunePMuonBestTrack->at(i));
      h1_EtaEffptTrackIso_->Fill(fabs(Mu_etaTunePMuonBestTrack->at(i)));
    }
    else continue;
  }
}


void ZprimeMuMuPat::PlotAbsDxy(){
  for(unsigned i=0; i<Mu_nbMuon->size(); i++){
    if(fabs(Mu_etaTunePMuonBestTrack->at(i)) < 2.4 &&
       Mu_isGlobalMuon->at(i) == 1 &&
       (Mu_ptTunePMuonBestTrack->at(i) > 45.0 && Mu_ptTunePMuonBestTrack->at(i) < 3000.0) &&
       //Mu_absdxyTunePMuonBestTrack->at(i) < 0.2 &&
       (Mu_trackiso->at(i)/Mu_ptInnerTrack->at(i)) < 0.10  &&
       Mu_numberOftrackerLayersWithMeasurement->at(i) > 5 && 
       Mu_numberOfValidPixelHits->at(i) > 0 &&
       Mu_numberOfValidMuonHits->at(i) > 0 &&
       Mu_numberOfMatchedStations->at(i) > 1 &&
       Mu_dPToverPTTunePMuonBestTrack->at(i) < 0.3 ){
      h1_absdxy_->Fill( Mu_absdxyTunePMuonBestTrack->at(i) );
      h1_PtEffptabsdsy_->Fill(Mu_ptTunePMuonBestTrack->at(i));
      h1_EtaEffptabsdsy_->Fill(fabs(Mu_etaTunePMuonBestTrack->at(i)));
    }
    else continue;
  }
}



void ZprimeMuMuPat::plotAllHighPtMuonsID(){
  MuonPassingID();
  PlotPterror();
  PlotNbTrackLayers();
  PlotNBValidPixelHits();
  PlotNbValidMuonHits();
  PlotNbMatchedStations();
  PlotTrackiso();
  PlotAbsDxy();
  MuonPassingNewID();
  MuonPassingTightID();
}

void ZprimeMuMuPat::MuonPassingNewID(){
  for(unsigned i=0; i<Mu_nbMuon->size(); i++){
    if(fabs(Mu_etaTunePMuonBestTrack->at(i)) < 2.4 &&
       Mu_isGlobalMuon->at(i) == 1 &&
       (Mu_ptTunePMuonBestTrack->at(i) > 45.0 && Mu_ptTunePMuonBestTrack->at(i) < 3000.0) &&
       Mu_absdxyTunePMuonBestTrack->at(i) < 0.02 &&
       (Mu_trackiso->at(i)/Mu_ptInnerTrack->at(i)) < 0.10  &&
       Mu_numberOftrackerLayersWithMeasurement->at(i) > 5 &&
       Mu_numberOfValidPixelHits->at(i) > 0 &&
       Mu_numberOfValidMuonHits->at(i) > 0 &&
       Mu_numberOfMatchedStations->at(i) > 1 &&
       Mu_dPToverPTTunePMuonBestTrack->at(i) < 0.3 ){
      h1_PtNewID_->Fill(Mu_ptTunePMuonBestTrack->at(i));
      h1_EtaNewID_->Fill(fabs(Mu_etaTunePMuonBestTrack->at(i)));
      h1_nbPVNewID_->Fill(Mu_nbofpv->at(i));
    }
    else continue;
  }
}



void ZprimeMuMuPat::MuonPassingTightID(){
  for(unsigned i=0; i<Mu_nbMuon->size(); i++){
    if(fabs(Mu_etaTunePMuonBestTrack->at(i)) < 2.4 &&
       Mu_isGlobalMuon->at(i) == 1 &&
       (Mu_ptTunePMuonBestTrack->at(i) > 45.0 && Mu_ptTunePMuonBestTrack->at(i) < 3000.0) &&
       Mu_absdxyTunePMuonBestTrack->at(i) < 0.01 &&
       (Mu_trackiso->at(i)/Mu_ptInnerTrack->at(i)) < 0.10  &&
       Mu_numberOftrackerLayersWithMeasurement->at(i) > 5 &&
       Mu_numberOfValidPixelHits->at(i) > 0 &&
       Mu_numberOfValidMuonHits->at(i) > 0 &&
       Mu_numberOfMatchedStations->at(i) > 1 &&
       Mu_dPToverPTTunePMuonBestTrack->at(i) < 0.3 ){
      h1_PtTightID_->Fill(Mu_ptTunePMuonBestTrack->at(i));
      h1_EtaTightID_->Fill(fabs(Mu_etaTunePMuonBestTrack->at(i)));
      h1_nbPVTightID_->Fill(Mu_nbofpv->at(i));
    }
    else continue;
  }
}


void ZprimeMuMuPat::CosThetaCollinSoper(float Et1,float Eta1,float Phi1,float En1,
					float Et2,float Eta2,float Phi2,float En2,
					float ChargeEle1,float RecoMass){
  
  TLorentzVector Ele;
  TLorentzVector Elebar;
  if(ChargeEle1<0) {
    Ele.SetPtEtaPhiE(Et1,Eta1,Phi1,En1);
    Elebar.SetPtEtaPhiE(Et2,Eta2,Phi2,En2);
  }
  if(ChargeEle1>0) {
    Ele.SetPtEtaPhiE(Et2,Eta2,Phi2,En2);
    Elebar.SetPtEtaPhiE(Et1,Eta1,Phi1,En1);
  }
  TLorentzVector Q(Ele+Elebar);
  
  //************************************************************************
  //
  // 1) cos(theta) = 2 Q^-1 (Q^2+Qt^2)^-1 (Ele^+ Elebar^- - Ele^- Elebar^+)
  //
  //
  //************************************************************************
  double Eleplus  = 1.0/sqrt(2.0) * (Ele.E() + Ele.Z());
  double Eleminus = 1.0/sqrt(2.0) * (Ele.E() - Ele.Z());
  
  double Elebarplus  = 1.0/sqrt(2.0) * (Elebar.E() + Elebar.Z());
  double Elebarminus = 1.0/sqrt(2.0) * (Elebar.E() - Elebar.Z());
  
  double costheta = 2.0 / Q.Mag() / sqrt(pow(Q.Mag(), 2) + pow(Q.Pt(), 2)) * 
    (Eleplus * Elebarminus - Eleminus * Elebarplus);
  if (Q.Pz()<0.0) costheta = -costheta;

  if( RecoMass > 60.0 && RecoMass < 120.0 ) 
    {
      h1_CosAngleCollinSoperCorrect60Mass120_->Fill(costheta);
    }

  
  if( RecoMass > 120.0 && RecoMass < 300.0 ) 
    {
      h1_CosAngleCollinSoperCorrect120Mass300_->Fill(costheta);
    }

  if( RecoMass > 300.0 && RecoMass < 700.0 ) 
    {
      h1_CosAngleCollinSoperCorrect300Mass700_->Fill(costheta);
    }

  if( RecoMass > 700.0 && RecoMass < 3000.0 ) 
    {
      h1_CosAngleCollinSoperCorrect700Mass3000_->Fill(costheta);
    }

   if( RecoMass > 4900.0 && RecoMass < 5100.0 ) 
    {
      h1_CosAngleCollinSoperCorrect4900Mass5100_->Fill(costheta);
    }

   /************************************************************************
   *
   * 2) tanphi = (Q^2 + Qt^2)^1/2 / Q (Dt dot R unit) /(Dt dot Qt unit)
   *
   ************************************************************************/
  TLorentzVector   Pbeam(0.0, 0.0,  4000., 4000.); // beam momentum in lab frame
  TLorentzVector Ptarget(0.0, 0.0, -4000., 4000.); // beam momentum in lab frame
  TLorentzVector D(Ele-Elebar);
  // unit vector on R direction
  TVector3 R = Pbeam.Vect().Cross(Q.Vect());
  TVector3 Runit = R.Unit();
  // unit vector on Qt
  TVector3 Qt = Q.Vect(); Qt.SetZ(0);
  TVector3 Qtunit = Qt.Unit();
  TVector3 Dt = D.Vect(); Dt.SetZ(0);
  double tanphi = sqrt(pow(Q.Mag(), 2) + pow(Q.Perp(), 2)) / Q.Mag() * Dt.Dot(Runit) / Dt.Dot(Qtunit);
  if (Q.Pz()<0.0) tanphi = -tanphi;
  //h1_TanPhiCollinSoperCorrect_->Fill(tanphi);
  /************************************************************************
   *
   * 3) sin2(theta) = Q^-2 Dt^2 - Q^-2 (Q^2 + Qt^2)^-1 * (Dt dot Qt)^2
   *
   ************************************************************************/
  double dt_qt = D.X()*Q.X() + D.Y()*Q.Y();
  double sin2theta = pow(D.Pt()/Q.Mag(), 2)
    - 1.0/pow(Q.Mag(), 2)/(pow(Q.Mag(), 2) + pow(Q.Pt(), 2))*pow(dt_qt, 2);
  //h1_Sin2AngleCollinSoperCorrect_->Fill(sin2theta);
}

