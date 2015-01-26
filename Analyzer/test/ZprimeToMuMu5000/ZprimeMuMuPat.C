//=========================================================================  
//          Analysis code for Z' boson to Mu Mu analysis                  =  
//        In this code we select the high pt di-muons events              =
//                        To run over MC                                  =  
//               Author: Sherif Elgammal                                  =
//                                                                        =
//                     14/01/2015                                         =
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
  int binMass   = 100;
  float minMass = 0.0;
  float maxMass = 8000.0;
  //values needed for tag & probe study
  MassCutMin = 0.0;
  MassCutMax = 2000.0;
  EtaCut = 2.4;
  deltaRcut = 0.15;
  RecoHLTMatchingDeltaRcut = 0.20;
  minMassCut = 50.0;
  maxMassCut = 4000.0;
  //================================================================================== 
  //                                                                                 =
  //            Start the histograms for the mass of Z                               =  
  //                                                                                 =
  //==================================================================================
  h1_ZprimeRecomass_             = new TH1F("ZprimeRecomass","",binMass,minMass,maxMass);
  h1_ZprimeGenmass_              = new TH1F("ZprimeGenmass","",binMass,minMass,maxMass);
  h1_ZprimeGenEta1_              = new TH1F("ZprimeGenEta1","",100,-8.0,8.0);
  h1_ZprimeGenEta2_              = new TH1F("ZprimeGenEta2","",100,-8.0,8.0);
  h1_ZprimeGenPt1_               = new TH1F("ZprimeGenPt1","",100,0.0,8000.0);
  h1_ZprimeGenPt2_               = new TH1F("ZprimeGenPt2","",100,0.0,8000.0);
  h1_ZprimeGenEn1_               = new TH1F("ZprimeGenEn1","",100,0.0,8000.0);
  h1_ZprimeGenEn2_               = new TH1F("ZprimeGenEn2","",100,0.0,8000.0);  
  h1_RecoMassResolution_         = new TH1F("RecoMassResolution","",100,-0.5,0.5);
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
  h1_EtaEffpterror_                  = new TH1F("EtaEffpterror","",60,0.0,3.0);
  h1_EtaEffptnumberOftrackerLayers_  = new TH1F("EtaEffptnumberOftrackerLayers","",60,0.0,3.0);
  h1_EtaEffptnumberOfPixelHits_      = new TH1F("EtaEffptnumberOfPixelHits","",60,0.0,3.0);
  h1_EtaEffptnumberOfMuonHits_       = new TH1F("EtaEffptnumberOfMuonHits","",60,0.0,3.0);
  h1_EtaEffptnumberOfMatchedStations_= new TH1F("EtaEffptnumberOfMatchedStations","",60,0.0,3.0);
  h1_EtaEffptTrackIso_               = new TH1F("EtaEffptTrackIso","",60,0.0,3.0);
  h1_EtaEffptabsdsy_                 = new TH1F("EtaEffptabsdsy","",60,0.0,3.0);
  h1_nbPVID_               = new TH1F("nbPVID","",50,0.0,50.0);
  h1_PtID_                 = new TH1F("PtID","",20,0.0,3000.0);
  h1_EtaID_                = new TH1F("EtaID","",60,0.0,3.0);
  h1_nbPVNewID_            = new TH1F("nbPVNewID","",50,0.0,50.0);
  h1_PtNewID_              = new TH1F("PtNewID","",20,0.0,3000.0);
  h1_EtaNewID_             = new TH1F("EtaNewID","",60,0.0,3.0);
  h1_nbPVTightID_          = new TH1F("nbPVTightID","",50,0.0,50.0);
  h1_PtTightID_            = new TH1F("PtTightID","",20,0.0,3000.0);
  h1_EtaTightID_           = new TH1F("EtaTightID","",60,0.0,3.0);
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
    //PlotGenInfo(MassGen,genEta1,genEta2,genET1,genET2,genEn1,genEn2);
    //acceptance conditions
    if(fabs(genEta1)>EtaCut || fabs(genEta2)>EtaCut) continue;
    if(MassGen<60.0) continue;
    PlotGenInfo(MassGen,genEta1,genEta2,genET1,genET2,genEn1,genEn2);
    NbGen ++;
    //=========================================================
    //                                                        =
    // Calling methods to get events with 2 muons passing ID  =
    //                                                        =
    //=========================================================
    bool firstMuFinal  = SelectFirstMuon(PtRecMu1,EnRecMu1,EtaRecMu1,PhiRecMu1,ChargeRecMu1,flagmu1,
					 pxRecMu1,pyRecMu1,pzRecMu1,pRecMu1,dxyRecMu1);
    bool secondMuFinal = SelectSecondMuon(ChargeRecMu1,flagmu1,PtRecMu1,PtRecMu2,EnRecMu2,EtaRecMu2,PhiRecMu2,
					  ChargeRecMu2,pxRecMu2,pyRecMu2,pzRecMu2,pRecMu2,dxyRecMu2);
    PickThehighestMass(vtxMassMu,vtxChi2Mu);
    double CosmicRejec = ThreeDangle(pxRecMu1,pyRecMu1,pzRecMu1,pRecMu1,
				     pxRecMu2,pyRecMu2,pzRecMu2,pRecMu2);
    //=========================================================
    //                   compute invariant mass               =
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
    //std::cout<<"I matched two muons at reco and gen in a cone 0.15"<<endl;
    if(vtxChi2Mu<=10.0 && CosmicRejec>=-0.9998) 
      {
	if( (fabs(EtaRecMu1) < 2.1 && fabs(EtaRecMu2) < EtaCut) || (fabs(EtaRecMu2) < 2.1 && fabs(EtaRecMu1) < EtaCut) )
	  {  
	    //=========================================================
	    //	      start doing matching between reco & HLT         =
	    //	                                                      =
	    //=========================================================
            bool RecoMuon1MatchingWithHLT = RecoHLMuonTMatching(EtaRecMu1,PhiRecMu1);
	    bool RecoMuon2MatchingWithHLT = RecoHLMuonTMatching(EtaRecMu2,PhiRecMu2);
	    if( (fabs(EtaRecMu1)<2.1 && RecoMuon1MatchingWithHLT==1) || (fabs(EtaRecMu2)<2.1 && RecoMuon2MatchingWithHLT==1) )
	      {
		NbReco++;
		PlotleadingAndSecondleading(vtxMassMu);
                h1_RecoMassResolution_->Fill((vtxMassMu-MassGen)/MassGen);
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
  TFile *output = new TFile("ZprimetoMuMu-MC-Mass5000-CMSSW720_miniaod.root","recreate");
  output->cd();
  /** Write the histograms */
  h1_ZprimeRecomass_->Write();
  h1_RecoMassResolution_->Write();
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
				 float &pmuon1,float &dxymuon1){
  int NbHEEPele = 0;
  unsigned iflag = -10;
  for(unsigned i=0; i<Mu_nbMuon->size(); i++){
    if( Mu_isGlobalMuon->at(i) == 1 &&
	Mu_ptBestTrack->at(i) > 45.0 &&
	Mu_absdxyBestTrack->at(i) < 0.2 &&
	(Mu_trackiso->at(i)/Mu_ptInnerTrack->at(i)) < 0.10  &&
	Mu_numberOftrackerLayersWithMeasurement->at(i) > 5 && 
	Mu_numberOfValidPixelHits->at(i) > 0 &&
	Mu_numberOfValidMuonHits->at(i) > 0 &&
	Mu_numberOfMatchedStations->at(i) > 1 &&
	Mu_dPToverPTBestTrack->at(i) < 0.3 ) { 
      iflag  = i;
      NbHEEPele ++;
    }
    else continue;
  }
  if( NbHEEPele > 0 ){
    FlagMu1             = iflag;
    pTmuon1             = Mu_ptBestTrack->at(iflag);
    Enmuon1             = Mu_en->at(iflag);
    Etamuon1            = Mu_etaBestTrack->at(iflag);
    Phimuon1            = Mu_phiBestTrack->at(iflag);
    ChargeMu1           = Mu_chargeBestTrack->at(iflag);
    pxmuon1             = Mu_pxBestTrack->at(iflag);
    pymuon1             = Mu_pyBestTrack->at(iflag);
    pzmuon1             = Mu_pzBestTrack->at(iflag);
    pmuon1              = Mu_pBestTrack->at(iflag);
    dxymuon1            = Mu_absdxyBestTrack->at(iflag);
    return true;  
  }
  else return false;
}
//============================ Method to select second high pt muon ========================
bool ZprimeMuMuPat::SelectSecondMuon(int ChargeMu1,unsigned FlagMu1,float pTmuon1,float &pTmuon2,float &Enmuon2,
				  float &Etamuon2,float &Phimuon2,int &ChargeMu2,float &pxmuon2,
				  float &pymuon2,float &pzmuon2,float &pmuon2,float &dxymuon2)
{
  int NbHEEPele = 0;
  unsigned iflag = -10;
  for(unsigned i=0; i<Mu_nbMuon->size(); i++){
    if(i == FlagMu1) continue;
    if(Mu_ptBestTrack->at(i) == pTmuon1) continue;
    if(ChargeMu1*Mu_chargeBestTrack->at(i)>0) continue;
    if( Mu_isGlobalMuon->at(i) == 1 &&
	Mu_ptBestTrack->at(i) > 45.0 &&
	Mu_absdxyBestTrack->at(i) < 0.2 &&
	(Mu_trackiso->at(i)/Mu_ptInnerTrack->at(i)) < 0.10  &&
	Mu_numberOftrackerLayersWithMeasurement->at(i) > 5 && 
	Mu_numberOfValidPixelHits->at(i) > 0 &&
	Mu_numberOfValidMuonHits->at(i) > 0 &&
	Mu_numberOfMatchedStations->at(i) > 1 &&
	Mu_dPToverPTBestTrack->at(i) < 0.3 ) { 
      iflag  = i;
      NbHEEPele ++;
    }
    else continue;
  }
  if( NbHEEPele > 0 ){
    pTmuon2          = Mu_ptBestTrack->at(iflag);
    Enmuon2          = Mu_en->at(iflag);
    Etamuon2         = Mu_etaBestTrack->at(iflag);
    Phimuon2         = Mu_phiBestTrack->at(iflag);
    ChargeMu2        = Mu_chargeBestTrack->at(iflag);
    pxmuon2          = Mu_pxBestTrack->at(iflag);
    pymuon2          = Mu_pyBestTrack->at(iflag);
    pzmuon2          = Mu_pzBestTrack->at(iflag);
    pmuon2           = Mu_pBestTrack->at(iflag);
    dxymuon2         = Mu_absdxyBestTrack->at(iflag);
    return true;  
  }
  else return false;
}
void ZprimeMuMuPat::PlotleadingAndSecondleading(float vertexMassMu){
  h1_ZprimeRecomass_->Fill(vertexMassMu);
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
void ZprimeMuMuPat::PrintEventInformation(int runNumber, int lumiNumber, int eventNumber,
				       float vtxChi2, float vtxMass, float CosmicRejection)
{
  if(event_runNo == runNumber && event_lumi == lumiNumber && event_evtNo == eventNumber)
    {
      output_txt << event_runNo
		 << "        " << event_lumi
		 << "        " << event_evtNo
		 << "        " << vtxChi2
		 << "        " << vtxMass << endl;
      for(unsigned i=0; i<Mu_nbMuon->size(); i++){
	//if( fabs(Mu_etaBestTrack->at(i)) < EtaCut ){cout<<"[1] eta="<<Mu_etaBestTrack->at(i)<<endl;}
	cout<<"[0] phi="<<Mu_phiBestTrack->at(i)<<endl;
	cout<<"[1] eta="<<Mu_etaBestTrack->at(i)<<endl;
	if(Mu_isGlobalMuon->at(i) == 1) {cout<<"[2] isGlobal="<<Mu_isGlobalMuon->at(i)<<endl;}
	//if(Mu_ptc->at(i) > 45.0) {cout<<"[3] ptcocktail="<<Mu_ptcocktail->at(i)<<endl;}
	if(Mu_absdxyBestTrack->at(i) < 0.2) {cout<<"[4] absdxy="<<Mu_absdxyBestTrack->at(i)<<endl;}
	if(Mu_trackiso->at(i)/Mu_ptInnerTrack->at(i) < 0.10) {cout<<"[5] trackiso="<<Mu_trackiso->at(i)/Mu_ptInnerTrack->at(i)<<endl;}
	if(Mu_numberOftrackerLayersWithMeasurement->at(i) > 5) {cout<<"[6] nbTrackerLayer="<<Mu_numberOftrackerLayersWithMeasurement->at(i)<<endl;}
	if(Mu_numberOfValidPixelHits->at(i) > 0) {cout<<"[7] nbPixelHits="<<Mu_numberOfValidPixelHits->at(i)<<endl;}
	if(Mu_numberOfValidMuonHits->at(i) > 0) {cout<<"[8] nbMuonHits="<<Mu_numberOfValidMuonHits->at(i)<<endl;}
	if(Mu_numberOfMatchedStations->at(i) > 1) {cout<<"[9] nbStation="<<Mu_numberOfMatchedStations->at(i)<<endl;}
	//if(Mu_dPToverPTcocktail->at(i) < 0.3) {cout<<"[10] DeltaPterror="<<Mu_dPToverPTcocktail->at(i)<<endl;}
	cout<<"[11] Charge="<<Mu_chargeBestTrack->at(i)<<endl;
      }
      cout<<"[000] vtxMassMu="<<vtxMass<<endl;
      cout<<"[000] vtxChi2Mu="<<vtxChi2<<endl;
      cout<<"[000] CosAngle="<<CosmicRejection<<endl;
    }
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
//----------------------------------------------------
//                                                   -
//       Part for HLT & Reco Matching                -
//                                                   -  
//----------------------------------------------------
bool ZprimeMuMuPat::RecoHLMuonTMatching(float RecoEta,float RecoPhi){
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
//============================ Method to plot Gen Mu ========================
void ZprimeMuMuPat::PlotGenInfo(float ZprimeGenMass,float EtaGenMu1,float EtaGenMu2,float PtGenMu1,
				float PtGenMu2,float EnGenMu1,float EnGenMu2){
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
    if(fabs(Mu_etaBestTrack->at(i)) < 2.4 &&
       Mu_isGlobalMuon->at(i) == 1 &&
       (Mu_ptBestTrack->at(i) > 45.0 && Mu_ptBestTrack->at(i) < 3000.0) &&
       Mu_absdxyBestTrack->at(i) < 0.2 &&
       (Mu_trackiso->at(i)/Mu_ptInnerTrack->at(i)) < 0.10  &&
       Mu_numberOftrackerLayersWithMeasurement->at(i) > 5 && 
       Mu_numberOfValidPixelHits->at(i) > 0 &&
       Mu_numberOfValidMuonHits->at(i) > 0 &&
       Mu_numberOfMatchedStations->at(i) > 1 &&
       Mu_dPToverPTBestTrack->at(i) < 0.3 ){
      h1_PtID_->Fill(Mu_ptBestTrack->at(i));
      h1_EtaID_->Fill(fabs(Mu_etaBestTrack->at(i)));
      h1_nbPVID_->Fill(Mu_nbofpv->at(i));
    }
    else continue;
  }
}

void ZprimeMuMuPat::PlotPterror(){
  for(unsigned i=0; i<Mu_nbMuon->size(); i++){
    if(fabs(Mu_etaBestTrack->at(i)) < 2.4 &&
       Mu_isGlobalMuon->at(i) == 1 &&
       (Mu_ptBestTrack->at(i) > 45.0 && Mu_ptBestTrack->at(i) < 3000.0) &&
       Mu_absdxyBestTrack->at(i) < 0.2 &&
       (Mu_trackiso->at(i)/Mu_ptInnerTrack->at(i)) < 0.10  &&
       Mu_numberOftrackerLayersWithMeasurement->at(i) > 5 && 
       Mu_numberOfValidPixelHits->at(i) > 0 &&
       Mu_numberOfValidMuonHits->at(i) > 0 &&
       Mu_numberOfMatchedStations->at(i) > 1){
      h1_dPToverPT_->Fill(Mu_dPToverPTBestTrack->at(i) );
      h1_PtEffpterror_->Fill(Mu_ptBestTrack->at(i));
      h1_EtaEffpterror_->Fill(fabs(Mu_etaBestTrack->at(i)));
    }
    else continue;
  }
}

void ZprimeMuMuPat::PlotNbTrackLayers(){
  for(unsigned i=0; i<Mu_nbMuon->size(); i++){
    if(fabs(Mu_etaBestTrack->at(i)) < 2.4 &&
       Mu_isGlobalMuon->at(i) == 1 &&
       (Mu_ptBestTrack->at(i) > 45.0 && Mu_ptBestTrack->at(i) < 3000.0) &&
       Mu_absdxyBestTrack->at(i) < 0.2 &&
       (Mu_trackiso->at(i)/Mu_ptInnerTrack->at(i)) < 0.10  &&
       //Mu_numberOftrackerLayersWithMeasurement->at(i) > 5 && 
       Mu_numberOfValidPixelHits->at(i) > 0 &&
       Mu_numberOfValidMuonHits->at(i) > 0 &&
       Mu_numberOfMatchedStations->at(i) > 1 &&
       Mu_dPToverPTBestTrack->at(i) < 0.3 ){
      h1_numberOftrackerLayersWithMeasurement_->Fill( Mu_numberOftrackerLayersWithMeasurement->at(i) );
      h1_PtEffptnumberOftrackerLayers_->Fill(Mu_ptBestTrack->at(i));
      h1_EtaEffptnumberOftrackerLayers_->Fill(fabs(Mu_etaBestTrack->at(i)));
    }
    else continue;
  }
}


void ZprimeMuMuPat::PlotNBValidPixelHits(){
  for(unsigned i=0; i<Mu_nbMuon->size(); i++){
    if(fabs(Mu_etaBestTrack->at(i)) < 2.4 &&
       Mu_isGlobalMuon->at(i) == 1 &&
       (Mu_ptBestTrack->at(i) > 45.0 && Mu_ptBestTrack->at(i) < 3000.0) &&
       Mu_absdxyBestTrack->at(i) < 0.2 &&
       (Mu_trackiso->at(i)/Mu_ptInnerTrack->at(i)) < 0.10  &&
       Mu_numberOftrackerLayersWithMeasurement->at(i) > 5 && 
       //Mu_numberOfValidPixelHits->at(i) > 0 &&
       Mu_numberOfValidMuonHits->at(i) > 0 &&
       Mu_numberOfMatchedStations->at(i) > 1 &&
       Mu_dPToverPTBestTrack->at(i) < 0.3 ){
      h1_numberOfValidPixelHits_->Fill( Mu_numberOfValidPixelHits->at(i) );
      h1_PtEffptnumberOfPixelHits_->Fill(Mu_ptBestTrack->at(i));
      h1_EtaEffptnumberOfPixelHits_->Fill(fabs(Mu_etaBestTrack->at(i)));
    }
    else continue;
  }
}

void ZprimeMuMuPat::PlotNbValidMuonHits(){
  for(unsigned i=0; i<Mu_nbMuon->size(); i++){
    if(fabs(Mu_etaBestTrack->at(i)) < 2.4 &&
       Mu_isGlobalMuon->at(i) == 1 &&
       (Mu_ptBestTrack->at(i) > 45.0 && Mu_ptBestTrack->at(i) < 3000.0) &&
       Mu_absdxyBestTrack->at(i) < 0.2 &&
       (Mu_trackiso->at(i)/Mu_ptInnerTrack->at(i)) < 0.10  &&
       Mu_numberOftrackerLayersWithMeasurement->at(i) > 5 && 
       Mu_numberOfValidPixelHits->at(i) > 0 &&
       //Mu_numberOfValidMuonHits->at(i) > 0 &&
       Mu_numberOfMatchedStations->at(i) > 1 &&
       Mu_dPToverPTBestTrack->at(i) < 0.3 ){
      h1_numberOfValidMuonHits_->Fill( Mu_numberOfValidMuonHits->at(i) );
      h1_PtEffptnumberOfMuonHits_->Fill(Mu_ptBestTrack->at(i));
      h1_EtaEffptnumberOfMuonHits_->Fill(fabs(Mu_etaBestTrack->at(i)));
    }
    else continue;
  }
}


void ZprimeMuMuPat::PlotNbMatchedStations(){
  for(unsigned i=0; i<Mu_nbMuon->size(); i++){
    if(fabs(Mu_etaBestTrack->at(i)) < 2.4 &&
       Mu_isGlobalMuon->at(i) == 1 &&
       (Mu_ptBestTrack->at(i) > 45.0 && Mu_ptBestTrack->at(i) < 3000.0) &&
       Mu_absdxyBestTrack->at(i) < 0.2 &&
       (Mu_trackiso->at(i)/Mu_ptInnerTrack->at(i)) < 0.10  &&
       Mu_numberOftrackerLayersWithMeasurement->at(i) > 5 && 
       Mu_numberOfValidPixelHits->at(i) > 0 &&
       Mu_numberOfValidMuonHits->at(i) > 0 &&
       //Mu_numberOfMatchedStations->at(i) > 1 &&
       Mu_dPToverPTBestTrack->at(i) < 0.3 ){
      h1_numberOfMatchedStations_->Fill( Mu_numberOfMatchedStations->at(i) );
      h1_PtEffptnumberOfMatchedStations_->Fill(Mu_ptBestTrack->at(i));
      h1_EtaEffptnumberOfMatchedStations_->Fill(fabs(Mu_etaBestTrack->at(i)));
    }
    else continue;
  }
}


void ZprimeMuMuPat::PlotTrackiso(){
  for(unsigned i=0; i<Mu_nbMuon->size(); i++){
    if(fabs(Mu_etaBestTrack->at(i)) < 2.4 &&
       Mu_isGlobalMuon->at(i) == 1 &&
       (Mu_ptBestTrack->at(i) > 45.0 && Mu_ptBestTrack->at(i) < 3000.0) &&
       Mu_absdxyBestTrack->at(i) < 0.2 &&
       //(Mu_trackiso->at(i)/Mu_ptInnerTrack->at(i)) < 0.10  &&
       Mu_numberOftrackerLayersWithMeasurement->at(i) > 5 && 
       Mu_numberOfValidPixelHits->at(i) > 0 &&
       Mu_numberOfValidMuonHits->at(i) > 0 &&
       Mu_numberOfMatchedStations->at(i) > 1 &&
       Mu_dPToverPTBestTrack->at(i) < 0.3 ){
      h1_trackiso_->Fill( Mu_trackiso->at(i)/Mu_ptInnerTrack->at(i) );
      h1_PtEffptTrackIso_->Fill(Mu_ptBestTrack->at(i));
      h1_EtaEffptTrackIso_->Fill(fabs(Mu_etaBestTrack->at(i)));
    }
    else continue;
  }
}


void ZprimeMuMuPat::PlotAbsDxy(){
  for(unsigned i=0; i<Mu_nbMuon->size(); i++){
    if(fabs(Mu_etaBestTrack->at(i)) < 2.4 &&
       Mu_isGlobalMuon->at(i) == 1 &&
       (Mu_ptBestTrack->at(i) > 45.0 && Mu_ptBestTrack->at(i) < 3000.0) &&
       //Mu_absdxyBestTrack->at(i) < 0.2 &&
       (Mu_trackiso->at(i)/Mu_ptInnerTrack->at(i)) < 0.10  &&
       Mu_numberOftrackerLayersWithMeasurement->at(i) > 5 && 
       Mu_numberOfValidPixelHits->at(i) > 0 &&
       Mu_numberOfValidMuonHits->at(i) > 0 &&
       Mu_numberOfMatchedStations->at(i) > 1 &&
       Mu_dPToverPTBestTrack->at(i) < 0.3 ){
      h1_absdxy_->Fill( Mu_absdxyBestTrack->at(i) );
      h1_PtEffptabsdsy_->Fill(Mu_ptBestTrack->at(i));
      h1_EtaEffptabsdsy_->Fill(fabs(Mu_etaBestTrack->at(i)));
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
    if(fabs(Mu_etaBestTrack->at(i)) < 2.4 &&
       Mu_isGlobalMuon->at(i) == 1 &&
       (Mu_ptBestTrack->at(i) > 45.0 && Mu_ptBestTrack->at(i) < 3000.0) &&
       Mu_absdxyBestTrack->at(i) < 0.02 &&
       (Mu_trackiso->at(i)/Mu_ptInnerTrack->at(i)) < 0.10  &&
       Mu_numberOftrackerLayersWithMeasurement->at(i) > 5 &&
       Mu_numberOfValidPixelHits->at(i) > 0 &&
       Mu_numberOfValidMuonHits->at(i) > 0 &&
       Mu_numberOfMatchedStations->at(i) > 1 &&
       Mu_dPToverPTBestTrack->at(i) < 0.3 ){
      h1_PtNewID_->Fill(Mu_ptBestTrack->at(i));
      h1_EtaNewID_->Fill(fabs(Mu_etaBestTrack->at(i)));
      h1_nbPVNewID_->Fill(Mu_nbofpv->at(i));
    }
    else continue;
  }
}



void ZprimeMuMuPat::MuonPassingTightID(){
  for(unsigned i=0; i<Mu_nbMuon->size(); i++){
    if(fabs(Mu_etaBestTrack->at(i)) < 2.4 &&
       Mu_isGlobalMuon->at(i) == 1 &&
       (Mu_ptBestTrack->at(i) > 45.0 && Mu_ptBestTrack->at(i) < 3000.0) &&
       Mu_absdxyBestTrack->at(i) < 0.01 &&
       (Mu_trackiso->at(i)/Mu_ptInnerTrack->at(i)) < 0.10  &&
       Mu_numberOftrackerLayersWithMeasurement->at(i) > 5 &&
       Mu_numberOfValidPixelHits->at(i) > 0 &&
       Mu_numberOfValidMuonHits->at(i) > 0 &&
       Mu_numberOfMatchedStations->at(i) > 1 &&
       Mu_dPToverPTBestTrack->at(i) < 0.3 ){
      h1_PtTightID_->Fill(Mu_ptBestTrack->at(i));
      h1_EtaTightID_->Fill(fabs(Mu_etaBestTrack->at(i)));
      h1_nbPVTightID_->Fill(Mu_nbofpv->at(i));
    }
    else continue;
  }
}


