//=========================================================================  
//             Make the pat tuple for Ele & Mu analysis                   =
//                (works only on miniaod formate)                         =   
//                         CMSSW_7_2_0                                    =
//                                                                        =
//                    Author: Sherif Elgammal                             =
//                                                                        =
//                         14/01/2015                                     =
//=========================================================================
#include "MyCodeArea/Analyzer/interface/MakeMuonElectronPatTuple.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/GeometryVector/interface/Phi.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/Common/interface/View.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/MuonReco/interface/MuonCocktails.h"
#include "DataFormats/Common/interface/RefCore.h"
#include "RecoVertex/VertexTools/interface/InvariantMassFromVertex.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/Measurement1D.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/Common/interface/RefToBase.h"
#include "HLTrigger/HLTfilters/interface/HLTHighLevel.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHit.h"
using namespace edm;
using namespace std;
using namespace reco;

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidateFwd.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
//====================== end a part for photons ==========================
#include <vector>
#include <set>
#include <stdio.h>
#include "TFile.h"
#include <math.h>
#include "TF2.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TLorentzVector.h"
#include "Math/VectorUtil.h"
//========================================================================
MakeMuonElectronPatTuple::MakeMuonElectronPatTuple( const edm::ParameterSet& ps )
//========================================================================
{
  BosonID_                       = ps.getParameter<int>("GenBosonID");
  prunedGenToken_                = ps.getParameter<edm::InputTag>("pruned");
  packedGenToken_                = ps.getParameter<edm::InputTag>("packed");  
  vertexSrc                      = ps.getParameter<edm::InputTag>("vertexCollection");
  maxAbsZ_                       = ps.getParameter<double>("maxAbsZ");
  maxd0_                         = ps.getParameter<double>("maxd0");
  minNdof_                       = ps.getParameter<int>("minndof");
  NbGoodPv_                      = ps.getParameter<int>("NbGoodPv");
  muonToken_                     = ps.getParameter<edm::InputTag>("muons");
  triggerBits_                   = ps.getParameter<edm::InputTag>("bits");
  triggerObjects_                = ps.getParameter<edm::InputTag>("objects");
  triggerPrescales_              = ps.getParameter<edm::InputTag>("prescales");
  pfToken_                       = ps.getParameter<edm::InputTag>("pfCands");
  jetToken_                      = ps.getParameter<edm::InputTag>("jets");
  //------------- start a part for electrons -------------------------------
  electronToken_  = ps.getParameter<edm::InputTag>("electrons");
  tok_EB_         = ps.getParameter<edm::InputTag>("EBrecHitCollectionLabel");
  tok_EE_         = ps.getParameter<edm::InputTag>("EErecHitCollectionLabel");
  rhoIsoInputTag_ = ps.getParameter<edm::InputTag>("rhoIsoInputTag");
  
  // reducedBarrelRecHitToken_ = ps.getParameter<edm::InputTag>("reducedBarrelRecHitCollection");
   //reducedEndcapRecHitToken_ = ps.getParameter<edm::InputTag>("reducedEndcapRecHitCollection");

  //------------- end a part for electrons -------------------------------

  //===============================================================================================
  outputFile_ = ps.getParameter<std::string>("outputFile");
  rootFile_   = TFile::Open(outputFile_.c_str(),"RECREATE"); // open output file to store histograms
}

 
//========================================================================
MakeMuonElectronPatTuple::~MakeMuonElectronPatTuple()
//========================================================================    
{
  delete rootFile_;
}

//========================================================================
void MakeMuonElectronPatTuple::beginJob() {
  //======================================================================== 
  // go to *OUR* rootfile and book histograms
  rootFile_->cd();
  //========================================================
  mytree  = new TTree("tree","tr");
  //=============================================================
  //
  //           Create Branchs for Nb of event,run,lumi
  //
  //=============================================================
  TotalNbEvents = 0;
  NbEventsPass  = 0;
  Run       = 0;
  Event     = 0;
  lumi      = 0;
  bunch     = 0;
  mytree->Branch("NbEventsPassTrigger",&NbEventsPassTrigger);
  mytree->Branch("NbEventsPassTriggerandPVcond",&NbEventsPassTriggerandPVcond);
  mytree->Branch("event_runNo",  &Run,   "event_runNo/I");
  mytree->Branch("event_evtNo",  &Event, "event_evtNo/I");
  mytree->Branch("event_lumi",   &lumi,  "event_lumi/I");
  mytree->Branch("event_bunch",  &bunch, "event_bunch/I");
  //==================================================
  //
  //      Create Branchs for Electrons variables
  //
  //==================================================
  mytree->Branch("Ele_nbElectrons",&Ele_nbElectrons);
  mytree->Branch("Ele_isEcalDrivenSeed",&Ele_isEcalDrivenSeed);
  mytree->Branch("Ele_isPassConversionVeto",&Ele_isPassConversionVeto);
  mytree->Branch("Ele_charge",&Ele_charge);
  mytree->Branch("Ele_nbOfMissingHits",&Ele_nbOfMissingHits);
  mytree->Branch("Ele_nbVtx",&Ele_nbVtx);
  mytree->Branch("Ele_Et",&Ele_Et);
  mytree->Branch("Ele_EtFromCaloEn",&Ele_EtFromCaloEn);
  mytree->Branch("Ele_pt",&Ele_pt);
  mytree->Branch("Ele_thetaSC",&Ele_thetaSC);
  mytree->Branch("Ele_etaSC",&Ele_etaSC);
  mytree->Branch("Ele_phiSC",&Ele_phiSC);
  mytree->Branch("Ele_energySC",&Ele_energySC);
  mytree->Branch("Ele_preshowerEnergySC",&Ele_preshowerEnergySC);
  mytree->Branch("Ele_thetaTrack",&Ele_thetaTrack);
  mytree->Branch("Ele_etaTrack",&Ele_etaTrack);
  mytree->Branch("Ele_phiTrack",&Ele_phiTrack);
  mytree->Branch("Ele_hadronicOverEm",&Ele_hadronicOverEm);
  mytree->Branch("Ele_deltaEtaInSC",&Ele_deltaEtaInSC);
  mytree->Branch("Ele_deltaPhiInSC",&Ele_deltaPhiInSC);
  mytree->Branch("Ele_sigmaIetaIeta",&Ele_sigmaIetaIeta);
  mytree->Branch("Ele_e2x5Max",&Ele_e2x5Max);
  mytree->Branch("Ele_e1x5",&Ele_e1x5);
  mytree->Branch("Ele_e5x5",&Ele_e5x5);
  mytree->Branch("Ele_e2x5MaxOver5x5",&Ele_e2x5MaxOver5x5);
  mytree->Branch("Ele_e1x5Over5x5",&Ele_e1x5Over5x5);
  mytree->Branch("Ele_sigmaIetaIetaFull5x5",&Ele_sigmaIetaIetaFull5x5);
  mytree->Branch("Ele_e2x5MaxFull5x5",&Ele_e2x5MaxFull5x5);
  mytree->Branch("Ele_e1x5Full5x5",&Ele_e1x5Full5x5);
  mytree->Branch("Ele_e5x5Full5x5",&Ele_e5x5Full5x5);
  mytree->Branch("Ele_e2x5MaxOver5x5Full5x5",&Ele_e2x5MaxOver5x5Full5x5);
  mytree->Branch("Ele_e1x5Over5x5Full5x5",&Ele_e1x5Over5x5Full5x5);
  mytree->Branch("Ele_e2x5Right",&Ele_e2x5Right);
  mytree->Branch("Ele_e2x5Left",&Ele_e2x5Left);
  mytree->Branch("Ele_e2x5Top",&Ele_e2x5Top);
  mytree->Branch("Ele_e2x5Bottom",&Ele_e2x5Bottom);
  mytree->Branch("Ele_eMax",&Ele_eMax);
  mytree->Branch("Ele_eRight",&Ele_eRight);
  mytree->Branch("Ele_eLeft",&Ele_eLeft);
  mytree->Branch("Ele_eTop",&Ele_eTop);
  mytree->Branch("Ele_eBottom",&Ele_eBottom);
  mytree->Branch("Ele_dxy",&Ele_dxy);
  mytree->Branch("Ele_dz",&Ele_dz);
  mytree->Branch("Ele_rhoIso",&Ele_rhoIso);
  mytree->Branch("Ele_fbrem",&Ele_fbrem);
  mytree->Branch("Ele_EoverP",&Ele_EoverP);
  mytree->Branch("Ele_Xposition",&Ele_Xposition);
  mytree->Branch("Ele_Yposition",&Ele_Yposition);
  mytree->Branch("Ele_EcalPlusHcald1iso",&Ele_EcalPlusHcald1iso);
  mytree->Branch("Ele_dr03TkSumPt",&Ele_dr03TkSumPt);
  mytree->Branch("Ele_dr03EcalRecHitSumEt",&Ele_dr03EcalRecHitSumEt);
  mytree->Branch("Ele_dr03HcalDepth1TowerSumEt",&Ele_dr03HcalDepth1TowerSumEt);
  mytree->Branch("Ele_dr03HcalDepth1TowerSumEtBc",&Ele_dr03HcalDepth1TowerSumEtBc);
  mytree->Branch("Ele_hcalDepth1OverEcal",&Ele_hcalDepth1OverEcal);
  mytree->Branch("Ele_hcalDepth2OverEcal",&Ele_hcalDepth2OverEcal);
  mytree->Branch("Ele_dr03HcalDepth2TowerSumEt",&Ele_dr03HcalDepth2TowerSumEt);
  mytree->Branch("Ele_hcalDepth2TowerSumEtNoVeto",&Ele_hcalDepth2TowerSumEtNoVeto);
  mytree->Branch("Ele_hcalDepth1TowerSumEtNoVeto",&Ele_hcalDepth1TowerSumEtNoVeto);
  mytree->Branch("Ele_pfSumPhotonEt",&Ele_pfSumPhotonEt);
  mytree->Branch("Ele_pfSumChargedHadronPt",&Ele_pfSumChargedHadronPt);
  mytree->Branch("Ele_pfSumNeutralHadronEt",&Ele_pfSumNeutralHadronEt);
  mytree->Branch("Ele_pfSumPUPt",&Ele_pfSumPUPt);
  mytree->Branch("Ele_pfDeltaBeta",&Ele_pfDeltaBeta);
  //=============================================================
  //
  //           Create Branchs for Muons variables
  //
  //=============================================================
  mytree->Branch("Mu_nbMuon",&Mu_nbMuon);
  mytree->Branch("Mu_isGlobalMuon",&Mu_isGlobalMuon);
  mytree->Branch("Mu_isPF",&Mu_isPF);
  mytree->Branch("Mu_isTrackerMuon",&Mu_isTrackerMuon);
  mytree->Branch("Mu_et",&Mu_et);
  mytree->Branch("Mu_en",&Mu_en);
  mytree->Branch("Mu_pt",&Mu_pt);
  mytree->Branch("Mu_eta",&Mu_eta);
  mytree->Branch("Mu_phi",&Mu_phi);
  mytree->Branch("Mu_charge",&Mu_charge);
  mytree->Branch("Mu_ptBestTrack",&Mu_ptBestTrack);
  mytree->Branch("Mu_pxBestTrack",&Mu_pxBestTrack);
  mytree->Branch("Mu_pyBestTrack",&Mu_pyBestTrack);
  mytree->Branch("Mu_pzBestTrack",&Mu_pzBestTrack);
  mytree->Branch("Mu_pBestTrack",&Mu_pBestTrack);
  mytree->Branch("Mu_etaBestTrack",&Mu_etaBestTrack);
  mytree->Branch("Mu_phiBestTrack",&Mu_phiBestTrack);
  mytree->Branch("Mu_thetaBestTrack",&Mu_thetaBestTrack);
  mytree->Branch("Mu_chargeBestTrack",&Mu_chargeBestTrack);
  mytree->Branch("Mu_dPToverPTBestTrack",&Mu_dPToverPTBestTrack);
  mytree->Branch("Mu_absdxyBestTrack",&Mu_absdxyBestTrack);
  mytree->Branch("Mu_absdzBestTrack",&Mu_absdzBestTrack);
  mytree->Branch("Mu_ptInnerTrack",&Mu_ptInnerTrack);
  mytree->Branch("Mu_normalizedChi2",&Mu_normalizedChi2);
  mytree->Branch("Mu_vtxMass",&Mu_vtxMass);
  mytree->Branch("Mu_vtxNormChi2",&Mu_vtxNormChi2);
  mytree->Branch("Mu_numberOfMatchedStations",&Mu_numberOfMatchedStations);
  mytree->Branch("Mu_numberOfValidPixelHits",&Mu_numberOfValidPixelHits);
  mytree->Branch("Mu_numberOfValidMuonHits",&Mu_numberOfValidMuonHits);
  mytree->Branch("Mu_numberOftrackerLayersWithMeasurement",&Mu_numberOftrackerLayersWithMeasurement);
  mytree->Branch("Mu_emIso",&Mu_emIso);
  mytree->Branch("Mu_hadIso",&Mu_hadIso);
  mytree->Branch("Mu_trackiso",&Mu_trackiso);
  mytree->Branch("Mu_pfSumChargedHadronPt",&Mu_pfSumChargedHadronPt);
  mytree->Branch("Mu_pfSumNeutralHadronEt",&Mu_pfSumNeutralHadronEt);
  mytree->Branch("Mu_PFSumPhotonEt",&Mu_PFSumPhotonEt);
  mytree->Branch("Mu_pfSumPUPt",&Mu_pfSumPUPt);
  mytree->Branch("Mu_nbofpv",&Mu_nbofpv);
  mytree->Branch("Mu_pfDeltaBeta",&Mu_pfDeltaBeta);
  mytree->Branch("Mu_patDeltaBeta",&Mu_patDeltaBeta);
  mytree->Branch("Mu_ptTunePMuonBestTrack",&Mu_ptTunePMuonBestTrack);
  //=============================================================
  //
  //           Create Branchs for gen particles variables
  //
  //=============================================================
  mytree->Branch("iGen",&iGen);
  mytree->Branch("idGen",&idGen);
  mytree->Branch("statusGen",&statusGen);
  mytree->Branch("ptGen",&ptGen);
  mytree->Branch("etaGen",&etaGen);
  mytree->Branch("phiGen",&phiGen);
  mytree->Branch("chargeGen",&chargeGen);
  mytree->Branch("EnergyGen",&EnergyGen);
  mytree->Branch("pxGen",&pxGen);
  mytree->Branch("pyGen",&pyGen);
  mytree->Branch("pzGen",&pzGen);
  //=============================================================
  //
  //           Create Branchs for Pimary Vertice variables
  //
  //=============================================================
  mytree->Branch("nbPv",&nbPv);
  mytree->Branch("Nbdof",&Nbdof);
  mytree->Branch("PositionRho",&PositionRho);
  mytree->Branch("PositionX",&PositionX);
  mytree->Branch("PositionY",&PositionY);
  mytree->Branch("PositionZ",&PositionZ);
  //=============================================================
  //
  //           Create Branchs for Muons match HLT variables
  //
  //=============================================================
  mytree->Branch("MuHLTObj_nbObj",&MuHLTObj_nbObj);
  mytree->Branch("MuHLTObj_pt",&MuHLTObj_pt);
  mytree->Branch("MuHLTObj_eta",&MuHLTObj_eta);
  mytree->Branch("MuHLTObj_phi",&MuHLTObj_phi);
  //=============================================================
  //  
  //    Create Branchs for ecal rec hits variables
  //      
  //=============================================================
  mytree->Branch("ecal_iRechit",&ecal_iRechit);
  mytree->Branch("ecal_RecoFlag",&ecal_RecoFlag);
  mytree->Branch("ecal_isSaturated",&ecal_isSaturated);
  //mytree->Branch("ecal_EtaRecHit",&ecal_EtaRecHit);
  //mytree->Branch("ecal_PhiRecHit",&ecal_PhiRecHit);
  mytree->Branch("ecal_EnRecHit",&ecal_EnRecHit);
}
  

//========================================================================================
void MakeMuonElectronPatTuple::analyze( const edm::Event& evt, const edm::EventSetup& es ) {
//======================================================================================
  using namespace edm; // needed for all fwk related classe
  //==============================================
  //=        Begin of the main program           =
  //============================================== 
  TotalNbEvents ++;
  NbEventsPassTrigger.clear();
  NbEventsPassTrigger.push_back(TotalNbEvents);
  //===================== Handle For Primary Vertics ===========
  edm::Handle<reco::VertexCollection> pvHandle; 
  evt.getByLabel(vertexSrc,pvHandle);
  const reco::VertexCollection &vertices = *pvHandle.product();
  bool GoodPv  = PrimaryVertex(vertices);
  if( GoodPv == 0 ) return;    //the aim of that command is to select only events
                               // with one good reconstructed primary vertics
  NbEventsPass ++;
  NbEventsPassTriggerandPVcond.clear();
  NbEventsPassTriggerandPVcond.push_back(NbEventsPass);
  Run   = evt.id().run();
  Event = evt.id().event();
  lumi  = evt.luminosityBlock();
  bunch = evt.bunchCrossing();
  IntialValues();
  PrimaryVertexTree(vertices);
  PatMuonTree(evt,es);
  MuonParticleFlowTree(evt,es);
  accessGenInfo(evt,es);
  PatElectronTree(evt,es);
  EBecalRecHitsTree(evt,es);
  EEecalRecHitsTree(evt,es);
  //hlt matching
  TriggerMatchingTree(evt,es);
  //==============================================
  //=        End of the main program             =
  //============================================== 
  mytree->Fill();
}

//========================================================================
void MakeMuonElectronPatTuple::endJob() {
//========================================================================
  // go to *OUR* root file and store histograms
  rootFile_->cd();

  mytree->Write();

  rootFile_->Close();
}

//define this as a plug-in
DEFINE_FWK_MODULE(MakeMuonElectronPatTuple);
//===================================================================================
//=============================================================
//
//                Method for Pat Electron Tree
//
//=============================================================
void MakeMuonElectronPatTuple::PatElectronTree(const edm::Event& iEvent,const edm::EventSetup& es)
{  
  Ele_nbElectrons.clear();
  Ele_isEcalDrivenSeed.clear();
  Ele_isPassConversionVeto.clear();
  Ele_charge.clear();
  Ele_nbOfMissingHits.clear();
  Ele_nbVtx.clear();
  Ele_Et.clear();
  Ele_EtFromCaloEn.clear();
  Ele_pt.clear();
  Ele_thetaSC.clear();
  Ele_etaSC.clear();
  Ele_phiSC.clear();
  Ele_energySC.clear();
  Ele_preshowerEnergySC.clear();
  Ele_thetaTrack.clear();
  Ele_etaTrack.clear();
  Ele_phiTrack.clear();
  Ele_hadronicOverEm.clear();
  Ele_deltaEtaInSC.clear();
  Ele_deltaPhiInSC.clear();
  Ele_sigmaIetaIeta.clear();
  Ele_e2x5Max.clear();
  Ele_e1x5.clear();
  Ele_e5x5.clear();
  Ele_e2x5MaxOver5x5.clear();
  Ele_e1x5Over5x5.clear();
  Ele_sigmaIetaIetaFull5x5.clear();
  Ele_e2x5MaxFull5x5.clear();
  Ele_e1x5Full5x5.clear();
  Ele_e5x5Full5x5.clear();
  Ele_e2x5MaxOver5x5Full5x5.clear();
  Ele_e1x5Over5x5Full5x5.clear();
  Ele_e2x5Right.clear();
  Ele_e2x5Left.clear();
  Ele_e2x5Top.clear();
  Ele_e2x5Bottom.clear();
  Ele_eMax.clear();
  Ele_eRight.clear();
  Ele_eLeft.clear();
  Ele_eTop.clear();
  Ele_eBottom.clear();
  Ele_dxy.clear();
  Ele_dz.clear();
  Ele_rhoIso.clear();
  Ele_fbrem.clear();
  Ele_EoverP.clear();
  Ele_Xposition.clear();
  Ele_Yposition.clear();
  //------------- detector isolation -------------------------
  Ele_EcalPlusHcald1iso.clear();
  Ele_dr03TkSumPt.clear();
  Ele_dr03EcalRecHitSumEt.clear();
  Ele_dr03HcalDepth1TowerSumEt.clear();
  Ele_dr03HcalDepth1TowerSumEtBc.clear();
  Ele_hcalDepth1OverEcal.clear();
  Ele_hcalDepth2OverEcal.clear();
  Ele_dr03HcalDepth2TowerSumEt.clear();
  Ele_hcalDepth2TowerSumEtNoVeto.clear();
  Ele_hcalDepth1TowerSumEtNoVeto.clear();                       
  //------------- PF isolation from pat::ele -------------------------
  Ele_pfSumPhotonEt.clear();
  Ele_pfSumChargedHadronPt.clear();
  Ele_pfSumNeutralHadronEt.clear();
  Ele_pfSumPUPt.clear();
  Ele_pfDeltaBeta.clear();
  // rho for isolation
  edm::Handle<double> rhoIso_h;
  iEvent.getByLabel(rhoIsoInputTag_, rhoIso_h);
  double rhoIso = *(rhoIso_h.product());
  // primary vertex candidate collection
  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByLabel(vertexSrc,vertices);
  const reco::Vertex &PV = vertices->front();
  edm::Handle<pat::ElectronCollection> electrons;
  iEvent.getByLabel(electronToken_, electrons);
  //------------ Rec Hits in EB ----------------
  Handle<EcalRecHitCollection> ecalEB;
  iEvent.getByLabel(tok_EB_,ecalEB);
  //------------ Rec Hits in EE ----------------
  Handle<EcalRecHitCollection> ecalEE;
  iEvent.getByLabel(tok_EE_,ecalEE);
  //---------------- Topology ------------------
  //edm::ESHandle<CaloTopology> pTopology;
  //es.get<CaloTopologyRecord>().get(pTopology);
  //const CaloTopology *topology = pTopology.product();
  //prepare electron cluster shapes extraction
  //EcalClusterLazyTools lazyTools( iEvent, es, reducedBarrelRecHitToken_, reducedEndcapRecHitToken_ );  
  //std::auto_ptr<EcalClusterLazyTools> lazyTools_;
  //lazyTools_ .reset(new EcalClusterLazyTools( iEvent , es , tok_EB_ , tok_EE_ ));  
  for (const pat::Electron &el : *electrons) {
    if (el.pt() < 5) continue;
    NbElectrons++;
    Ele_nbElectrons.push_back(NbElectrons); 
    Ele_Et.push_back(el.superCluster()->energy() * sin(el.p4().theta()));
    Ele_EtFromCaloEn.push_back(el.caloEnergy() * sin(el.p4().theta()));
    Ele_pt.push_back(el.pt()); 
    Ele_thetaSC.push_back(el.caloPosition().theta()); //theta SC
    Ele_etaSC.push_back(el.superCluster()->eta());    //eta SC
    Ele_phiSC.push_back(el.superCluster()->phi());    //phi SC
    Ele_energySC.push_back(el.superCluster()->energy()); //energy SC
    Ele_preshowerEnergySC.push_back(el.superCluster()->preshowerEnergy()); 
    Ele_thetaTrack.push_back(el.p4().theta()); //theta track
    Ele_etaTrack.push_back(el.p4().eta());     //eta track
    Ele_phiTrack.push_back(el.p4().phi());     //phi track
    Ele_hadronicOverEm.push_back(el.hadronicOverEm());
    Ele_deltaEtaInSC.push_back(el.deltaEtaSuperClusterTrackAtVtx());
    Ele_deltaPhiInSC.push_back(el.deltaPhiSuperClusterTrackAtVtx());
    Ele_sigmaIetaIeta.push_back(el.sigmaIetaIeta());
    Ele_sigmaIetaIetaFull5x5.push_back(el.full5x5_sigmaIetaIeta()); 
    Ele_e2x5Max.push_back(el.e2x5Max());
    Ele_e1x5.push_back(el.e1x5());
    Ele_e5x5.push_back(el.e5x5());
    Ele_e2x5MaxOver5x5.push_back(el.e2x5Max()/el.e5x5());
    Ele_e1x5Over5x5.push_back(el.e1x5()/el.e5x5());
    Ele_e2x5MaxFull5x5.push_back(el.full5x5_e2x5Max());
    Ele_e1x5Full5x5.push_back(el.full5x5_e1x5());
    Ele_e5x5Full5x5.push_back(el.full5x5_e5x5());
    Ele_e2x5MaxOver5x5Full5x5.push_back(el.full5x5_e2x5Max()/el.full5x5_e5x5());
    Ele_e1x5Over5x5Full5x5.push_back(el.full5x5_e1x5()/el.full5x5_e5x5());
   /* 
    Ele_e2x5Right.push_back(lazyTools_->e2x5Right(*(el.superCluster()->seed())));
    Ele_e2x5Left.push_back(lazyTools_->e2x5Left(*(el.superCluster()->seed())));
    Ele_e2x5Top.push_back(lazyTools_->e2x5Top(*(el.superCluster()->seed())));
    Ele_e2x5Bottom.push_back(lazyTools_->e2x5Bottom(*(el.superCluster()->seed())));
    Ele_eMax.push_back(lazyTools_->eMax(*(el.superCluster()->seed())));
    Ele_eRight.push_back(lazyTools_->eRight(*(el.superCluster()->seed())));
    Ele_eLeft.push_back(lazyTools_->eLeft(*(el.superCluster()->seed())));
    Ele_eTop.push_back(lazyTools_->eTop(*(el.superCluster()->seed())));
    Ele_eBottom.push_back(lazyTools_->eBottom(*(el.superCluster()->seed())));
    */
    Ele_nbVtx.push_back(vertices->size());
    if(vertices->size()>0){
      Ele_dxy.push_back(el.gsfTrack()->dxy(PV.position()));   
      Ele_dz.push_back(el.gsfTrack()->dz(PV.position())); 
    }
    else{
      Ele_dxy.push_back(el.gsfTrack()->dxy());   
      Ele_dz.push_back(el.gsfTrack()->dz());
    }
    Ele_isEcalDrivenSeed.push_back(el.ecalDrivenSeed());
    Ele_isPassConversionVeto.push_back(el.passConversionVeto());
    Ele_charge.push_back(el.gsfTrack()->charge());
    Ele_rhoIso.push_back(rhoIso);
    Ele_nbOfMissingHits.push_back(el.gsfTrack()->numberOfLostHits()); 
    Ele_fbrem.push_back(el.fbrem());
    Ele_EoverP.push_back(el.eSeedClusterOverP());
    Ele_Xposition.push_back(el.caloPosition().x());   
    Ele_Yposition.push_back(el.caloPosition().y()); 
    //------------- detector isolation -------------------------
    Ele_hcalDepth1OverEcal.push_back(el.hcalDepth1OverEcal());
    Ele_hcalDepth2OverEcal.push_back(el.hcalDepth2OverEcal());
    Ele_dr03HcalDepth2TowerSumEt.push_back(el.dr03HcalDepth2TowerSumEt());
    Ele_hcalDepth2TowerSumEtNoVeto.push_back(el.isolationVariables03().hcalDepth2TowerSumEt);// hcaldepht2 iso deposit with 
                                                                                             // electron footprint removed
  
    Ele_hcalDepth1TowerSumEtNoVeto.push_back(el.isolationVariables03().hcalDepth1TowerSumEt);// hcaldepht1 iso deposit with 
                                                                                             // electron footprint removed
    
    Ele_EcalPlusHcald1iso.push_back(el.dr03EcalRecHitSumEt() + el.dr03HcalDepth1TowerSumEt());
    Ele_dr03TkSumPt.push_back(el.dr03TkSumPt());
    Ele_dr03EcalRecHitSumEt.push_back(el.dr03EcalRecHitSumEt());
    Ele_dr03HcalDepth1TowerSumEt.push_back(el.dr03HcalDepth1TowerSumEt());
    Ele_dr03HcalDepth1TowerSumEtBc.push_back(el.dr03HcalDepth1TowerSumEtBc());
    //------------- PF isolation from pat::ele -------------------------
    Ele_pfSumPhotonEt.push_back(el.pfIsolationVariables().sumPhotonEt);
    Ele_pfSumChargedHadronPt.push_back(el.pfIsolationVariables().sumChargedHadronPt); 
    Ele_pfSumNeutralHadronEt.push_back(el.pfIsolationVariables().sumNeutralHadronEt);
     Ele_pfSumPUPt.push_back(el.pfIsolationVariables().sumPUPt);  
    // do deltaBeta
    double charged   = el.pfIsolationVariables().sumPhotonEt;
    double neutral   = el.pfIsolationVariables().sumNeutralHadronEt;
    double pileup    = el.pfIsolationVariables().sumPUPt;
    double deltaBeta = charged + std::max(0.0, neutral-0.5*pileup);
    Ele_pfDeltaBeta.push_back(deltaBeta);
  }
}
//=============================================================
//
//                Method for Pat Muon Tree
//
//=============================================================
void MakeMuonElectronPatTuple::PatMuonTree(const edm::Event& evt,const edm::EventSetup& es)
{  
  Mu_nbMuon.clear();
  Mu_isGlobalMuon.clear();
  Mu_isPF.clear();
  Mu_isTrackerMuon.clear();
  Mu_isHighPtMuon.clear();
  Mu_ptBestTrack.clear();
  Mu_dPToverPTBestTrack.clear();
  Mu_en.clear();
  Mu_pt.clear();
  Mu_eta.clear();
  Mu_phi.clear();
  Mu_et.clear();
  Mu_charge.clear();
  Mu_numberOfMatchedStations.clear();
  Mu_absdxyBestTrack.clear();
  Mu_absdzBestTrack.clear();
  Mu_absdxy.clear();
  Mu_absdz.clear();
  Mu_trackiso.clear();
  Mu_pfSumChargedHadronPt.clear();
  Mu_pfSumNeutralHadronEt.clear();
  Mu_PFSumPhotonEt.clear();
  Mu_pfSumPUPt.clear();
  Mu_normalizedChi2.clear();
  Mu_numberOfValidPixelHits.clear();
  Mu_numberOftrackerLayersWithMeasurement.clear();
  Mu_numberOfValidMuonHits.clear();
  Mu_calEnergy.clear();
  Mu_emIso.clear();
  Mu_hadIso.clear();
  Mu_vtxMass.clear();
  Mu_vtxNormChi2.clear();
  Mu_nbofpv.clear();
  Mu_pxBestTrack.clear();
  Mu_pyBestTrack.clear();
  Mu_pzBestTrack.clear();
  Mu_pBestTrack.clear();
  Mu_etaBestTrack.clear();
  Mu_ptInnerTrack.clear();
  Mu_phiBestTrack.clear();
  Mu_thetaBestTrack.clear();
  Mu_chargeBestTrack.clear();
  Mu_ptTunePMuonBestTrack.clear();
  float vtxNormChi2,DiMass;
  // Get TransientTracks (for use in e.g. the vertex fit) for each of
  // the muon tracks, using e.g. the cocktail momentum.
  edm::ESHandle<TransientTrackBuilder> ttkb;
  es.get<TransientTrackRecord>().get("TransientTrackBuilder", ttkb);
  std::vector<reco::TransientTrack> ttv;
  // primary vertex candidate collection
  edm::Handle<reco::VertexCollection> vertices;
  evt.getByLabel(vertexSrc,vertices);
  const reco::Vertex &PV = vertices->front();
  // pat candidate collection
  edm::Handle<pat::MuonCollection> muons;
  evt.getByLabel(muonToken_, muons);
  for (const pat::Muon &mu : *muons) {
    if( mu.pt() < 5.0 ) continue;
    //============= Parameters related to matched Gen info =====================
    if( !mu.innerTrack().isNonnull() ) continue;
    if( !mu.globalTrack().isNonnull() ) continue;
    //if( mu.isTrackerMuon()==false ) continue;
    NbMuons++;
    Mu_nbMuon.push_back(NbMuons);
    Mu_isHighPtMuon.push_back(mu.isHighPtMuon(PV));
    Mu_isGlobalMuon.push_back(mu.isGlobalMuon());
    Mu_isPF.push_back(mu.isPFMuon());
    Mu_isTrackerMuon.push_back(mu.isTrackerMuon());
    Mu_nbofpv.push_back(vertices->size());
    //============== Parameters related to Kinematics =====================
    const reco::TrackRef& tunePTrack = mu.muonBestTrack();
    const reco::TrackRef& tunePTrack2 = mu.tunePMuonBestTrack(); 
    Mu_ptBestTrack.push_back(tunePTrack1->pt());
    Mu_ptTunePMuonBestTrack.push_back(tunePTrack2->pt());
    Mu_dPToverPTBestTrack.push_back(tunePTrack->ptError()/mu.pt());
    Mu_pxBestTrack.push_back(tunePTrack->px()); //px component of the track 
    Mu_pyBestTrack.push_back(tunePTrack->py()); //py component of the track 
    Mu_pzBestTrack.push_back(tunePTrack->pz()); //pz component of the track 
    Mu_pBestTrack.push_back(tunePTrack->p());   //magnitude of momentum vector 
    Mu_etaBestTrack.push_back(tunePTrack->eta());
    Mu_phiBestTrack.push_back(tunePTrack->phi());
    Mu_thetaBestTrack.push_back(tunePTrack->theta());
    Mu_chargeBestTrack.push_back(tunePTrack->charge());
    Mu_absdxyBestTrack.push_back(fabs(tunePTrack->dxy(PV.position())));
    Mu_absdzBestTrack.push_back(fabs(tunePTrack->dz(PV.position())));
    Mu_en.push_back(mu.energy());
    Mu_et.push_back(mu.et());
    Mu_pt.push_back(mu.pt());
    Mu_eta.push_back(mu.eta());
    Mu_phi.push_back(mu.phi());
    Mu_charge.push_back(mu.charge());
    Mu_ptInnerTrack.push_back(mu.innerTrack()->pt());    
    //====================== Parameters related to track quality =====================
    Mu_normalizedChi2.push_back(mu.globalTrack()->normalizedChi2());
    Mu_numberOfValidPixelHits.push_back(mu.globalTrack()->hitPattern().numberOfValidPixelHits());
    Mu_numberOfValidMuonHits.push_back(mu.globalTrack()->hitPattern().numberOfValidMuonHits());
    Mu_numberOftrackerLayersWithMeasurement.push_back(mu.globalTrack()->hitPattern().trackerLayersWithMeasurement());
    Mu_numberOfMatchedStations.push_back(mu.numberOfMatchedStations());
    //============= Parameters related to detector isolation =====================
    Mu_emIso.push_back(mu.isolationR03().emEt);
    Mu_hadIso.push_back(mu.isolationR03().hadEt);
    Mu_trackiso.push_back(mu.isolationR03().sumPt);
    //============= Parameters related to PF isolation =====================
    Mu_pfSumChargedHadronPt.push_back(mu.pfIsolationR03().sumChargedHadronPt);
    Mu_pfSumNeutralHadronEt.push_back(mu.pfIsolationR03().sumNeutralHadronEt);
    Mu_PFSumPhotonEt.push_back(mu.pfIsolationR03().sumPhotonEt);
    Mu_pfSumPUPt.push_back(mu.pfIsolationR03().sumPUPt);
    // do deltaBeta
    double charged   = mu.pfIsolationR03().sumChargedHadronPt;
    double neutral   = mu.pfIsolationR03().sumNeutralHadronEt;
    double pileup    = mu.pfIsolationR03().sumPUPt;
    double deltaBeta = (charged + std::max(0.0, neutral-0.5*pileup))/mu.pt();
    Mu_patDeltaBeta.push_back(deltaBeta);
    //------------------------------------------------------
    //   part to compute mass from vertex info             -   
    //   Given to me by Piotr                              -
    //------------------------------------------------------
    if( NbMuons > 0 && mu.muonBestTrack()->pt()>45.0 && 
        (mu.isolationR03().sumPt/mu.innerTrack()->pt())<0.10 &&
        (mu.muonBestTrack()->ptError()/mu.muonBestTrack()->pt())<0.3 && 
        fabs(mu.muonBestTrack()->dxy(PV.position()))<0.2 &&
        mu.globalTrack()->hitPattern().trackerLayersWithMeasurement()>5 &&
        mu.globalTrack()->hitPattern().numberOfValidPixelHits()>0 &&
        mu.globalTrack()->hitPattern().numberOfValidMuonHits()>0 &&
        mu.numberOfMatchedStations()>1 ){
      ttv.push_back(ttkb->build(tunePTrack));
      //if(ttv.size()==2){
      if(ttv.size()>1){  
        KalmanVertexFitter kvf(true); // false means no smoothing which means no track re-fit
        CachingVertex<5> vtx = kvf.vertex(ttv);
        vtxNormChi2 = vtx.totalChiSquared()/vtx.degreesOfFreedom();
        InvariantMassFromVertex imfv;
        static const double muon_mass = 0.1056583;
        Measurement1D mass = imfv.invariantMass(vtx, muon_mass);
        DiMass = mass.value();
        Mu_vtxMass.push_back(DiMass);
        Mu_vtxNormChi2.push_back(vtxNormChi2);
        nbTk++;
      }
    }
  }
}
//=============================================================
//
//     Method for accessing gen info
//
//=============================================================
bool MakeMuonElectronPatTuple::isAncestor(const reco::Candidate* ancestor, const reco::Candidate * particle)
{
  //particle is already the ancestor
  if(ancestor == particle ) return true;
  
  //otherwise loop on mothers, if any and return true if the ancestor is found
  for(size_t i=0;i< particle->numberOfMothers();i++)
    {
      if(isAncestor(ancestor,particle->mother(i))) return true;
    }
  //if we did not return yet, then particle and ancestor are not relatives
  return false;
}

void MakeMuonElectronPatTuple::accessGenInfo(const edm::Event& evt,const edm::EventSetup& es)
{ 
  iGen.clear();
  idGen.clear();
  statusGen.clear();
  ptGen.clear();
  etaGen.clear();
  phiGen.clear();
  chargeGen.clear();
  EnergyGen.clear();
  pxGen.clear();
  pyGen.clear();
  pzGen.clear();
  // Pruned particles are the one containing "important" stuff
  Handle<edm::View<reco::GenParticle> > pruned;
  evt.getByLabel(prunedGenToken_,pruned);
  // Packed particles are all the status 1, so usable to remake jets
  // The navigation from status 1 to pruned is possible (the other direction should be made by hand)
  Handle<edm::View<pat::PackedGenParticle> > packed;
  evt.getByLabel(packedGenToken_,packed);
  //let's try to find all status1 originating directly from a B meson decay 
  for(size_t i=0; i<pruned->size();i++){
    if(abs((*pruned)[i].pdgId()) == BosonID_){
      const Candidate * Zprime = &(*pruned)[i];
      //std::cout << "PdgID: " << Zprime->pdgId() << " pt " << Zprime->pt() << " eta: " << Zprime->eta() << " phi: " << Zprime->phi() << std::endl;
      //std::cout << "  found daugthers: " << std::endl;
      for(size_t j=0; j<packed->size();j++){
	//get the pointer to the first survied ancestor of a given packed GenParticle in the prunedCollection 
	const Candidate * motherInPrunedCollection = (*packed)[j].mother(0) ;
	if(motherInPrunedCollection != nullptr && isAncestor( Zprime , motherInPrunedCollection)){
	  NbGenMuons++;
	  iGen.push_back(NbGenMuons);
	  idGen.push_back((*packed)[j].pdgId());
	  statusGen.push_back((*packed)[j].status());
	  ptGen.push_back((*packed)[j].pt());
	  etaGen.push_back((*packed)[j].eta());
	  phiGen.push_back((*packed)[j].phi());
	  chargeGen.push_back((*packed)[j].charge());
	  EnergyGen.push_back((*packed)[j].energy());
	  pxGen.push_back((*packed)[j].px());
	  pyGen.push_back((*packed)[j].py());
	  pzGen.push_back((*packed)[j].pz());
	}
      }
    }
  }
}
//=============================================================
//
//         Method for finding good Primary Vertex
//
//=============================================================
bool MakeMuonElectronPatTuple::PrimaryVertex(const reco::VertexCollection &vertices)
{
  int nbGoodPv = 0;
  bool result = false; 
  for(reco::VertexCollection::const_iterator it=vertices.begin() ; it!=vertices.end() ; ++it)
    {
      if(it->ndof() > minNdof_ && 
	 ( (maxAbsZ_ <= 0.0) || fabs(it->position().z()) <= maxAbsZ_ )  && 
	 ( (maxd0_ <= 0.0) || fabs(it->position().rho()) <= maxd0_ ) ) nbGoodPv++;
    }
  if( nbGoodPv>=NbGoodPv_ ) result = true;
  return result;
}
//=============================================================
//
//            Method for Pimary Vertices Tree
//
//=============================================================
void MakeMuonElectronPatTuple::PrimaryVertexTree(const reco::VertexCollection &vertices)
{
  nbPv.clear();
  Nbdof.clear();
  PositionX.clear();
  PositionY.clear();
  PositionZ.clear();
  PositionRho.clear();
  for(reco::VertexCollection::const_iterator it=vertices.begin() ; it!=vertices.end() ; ++it){
    value3_ ++;
    nbPv.push_back(value3_);
    Nbdof.push_back((*it).ndof());
    PositionX.push_back((*it).position().x());
    PositionY.push_back((*it).position().y());
    PositionZ.push_back((*it).position().z());
    PositionRho.push_back((*it).position().rho());
  }
}
//=============================================================
//
//            Method for Trigger Matching Tree
//
//=============================================================
void MakeMuonElectronPatTuple::TriggerMatchingTree(const edm::Event& iEvent,const edm::EventSetup& es)
{ 
  MuHLTObj_nbObj.clear();
  MuHLTObj_pt.clear();
  MuHLTObj_eta.clear();
  MuHLTObj_phi.clear();
  edm::Handle<edm::TriggerResults> triggerBits;
  edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
  edm::Handle<pat::PackedTriggerPrescales> triggerPrescales;
  iEvent.getByLabel(triggerBits_, triggerBits);
  iEvent.getByLabel(triggerObjects_, triggerObjects);
  iEvent.getByLabel(triggerPrescales_, triggerPrescales);
  const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);
  //std::cout << "\n === TRIGGER PATHS === " << std::endl;
  for (unsigned int i = 0, n = triggerBits->size(); i < n; ++i) {
    std::cout << "Trigger " << names.triggerName(i) << 
    ", prescale " << triggerPrescales->getPrescaleForIndex(i) <<
    ": " << (triggerBits->accept(i) ? "PASS" : "fail (or not run)") 
    	      << std::endl;
  }
  //std::cout << "\n === TRIGGER OBJECTS === " << std::endl;
  for (pat::TriggerObjectStandAlone obj : *triggerObjects) { // note: not "const &" since we want to call unpackPathNames
    obj.unpackPathNames(names);
    std::cout << "\tTrigger object:  pt " << obj.pt() << ", eta " << obj.eta() << ", phi " << obj.phi() << std::endl;
    // Print trigger object collection and type
    std::cout << "\t   Collection: " << obj.collection() << std::endl;
    std::cout << "\t   Type IDs:   ";
    nMuHLT++;
    MuHLTObj_nbObj.push_back(nMuHLT);
    MuHLTObj_pt.push_back(obj.pt()); 
    MuHLTObj_eta.push_back(obj.eta());
    MuHLTObj_phi.push_back(obj.phi());
    //for (unsigned h = 0; h < obj.filterIds().size(); ++h) std::cout << " " << obj.filterIds()[h] ;
    //std::cout << std::endl;
    // Print associated trigger filters
    std::cout << "\t   Filters:    ";
    for (unsigned h = 0; h < obj.filterLabels().size(); ++h) std::cout << " " << obj.filterLabels()[h];
    std::cout << std::endl;
    std::vector<std::string> pathNamesAll  = obj.pathNames(false);
    std::vector<std::string> pathNamesLast = obj.pathNames(true);
    // Print all trigger paths, for each one record also if the object is associated to a 'l3' filter (always true for the
    // definition used in the PAT trigger producer) and if it's associated to the last filter of a successfull path (which
    // means that this object did cause this trigger to succeed; however, it doesn't work on some multi-object triggers)
    //std::cout << "\t   Paths (" << pathNamesAll.size()<<"/"<<pathNamesLast.size()<<"):    ";
    
    for (unsigned h = 0, n = pathNamesAll.size(); h < n; ++h) {
      bool isBoth = obj.hasPathName( pathNamesAll[h], true, true ); 
      bool isL3   = obj.hasPathName( pathNamesAll[h], false, true ); 
      bool isLF   = obj.hasPathName( pathNamesAll[h], true, false ); 
      bool isNone = obj.hasPathName( pathNamesAll[h], false, false ); 
      //std::cout << "   " << pathNamesAll[h];
      if (isBoth) std::cout << "(L,3)";
      if (isL3 && !isBoth) std::cout << "(*,3)";
      if (isLF && !isBoth) std::cout << "(L,*)";
      if (isNone && !isBoth && !isL3 && !isLF) std::cout << "(*,*)";
    }
    std::cout << std::endl;
   
  }
  //std::cout << std::endl;
}

//=============================================================
//
//            Method for particle Flow Tree
//
//=============================================================
void MakeMuonElectronPatTuple::MuonParticleFlowTree(const edm::Event& iEvent,const edm::EventSetup& es)
{ 
  Mu_pfDeltaBeta.clear();
  edm::Handle<pat::MuonCollection> muons;
  iEvent.getByLabel(muonToken_, muons);
  //edm::Handle<pat::ElectronCollection> electrons;
  //iEvent.getByLabel(electronToken_, electrons);
  edm::Handle<pat::PackedCandidateCollection> pfs;
  iEvent.getByLabel(pfToken_, pfs);
  edm::Handle<pat::JetCollection> jets;
  iEvent.getByLabel(jetToken_, jets);
  std::vector<const reco::Candidate *> leptons;
  for (const pat::Muon &mu : *muons) leptons.push_back(&mu);
  //for (const pat::Electron &el : *electrons) leptons.push_back(&el);
  for (const reco::Candidate *lep : leptons) {
    if (lep->pt() < 40.0) continue;
    // initialize sums
    double charged = 0, neutral = 0, pileup  = 0;
    // now get a list of the PF candidates used to build this lepton, so to exclude them
    std::vector<reco::CandidatePtr> footprint;
    for (unsigned int i = 0, n = lep->numberOfSourceCandidatePtrs(); i < n; ++i) {
      footprint.push_back(lep->sourceCandidatePtr(i));
    }
    // now loop on pf candidates
    for (unsigned int i = 0, n = pfs->size(); i < n; ++i) {
      const pat::PackedCandidate &pf = (*pfs)[i];
      if (deltaR(pf,*lep) < 0.2) {
	// pfcandidate-based footprint removal
	if (std::find(footprint.begin(), footprint.end(), reco::CandidatePtr(pfs,i)) != footprint.end()) {
	  continue;
	}
	if (pf.charge() == 0) {
	  if (pf.pt() > 0.5) neutral += pf.pt();
	} else if (pf.fromPV() >= 2) {
	  charged += pf.pt();
	} else {
	  if (pf.pt() > 0.5) pileup += pf.pt();
	}
      }
    }
    // do deltaBeta
    double iso = charged + std::max(0.0, neutral-0.5*pileup);
    //printf("%-8s of pt %6.1f, eta %+4.2f: relIso = %5.2f\n",
    //	   abs(lep->pdgId())==13 ? "muon" : "electron",
    //	   lep->pt(), lep->eta(), iso/lep->pt());
    Mu_pfDeltaBeta.push_back(iso/lep->pt());
  }

  // Let's compute the fraction of charged pt from particles with dz < 0.1 cm
  for (const pat::Jet &j :  *jets) {
    if (j.pt() < 40 || fabs(j.eta()) > 2.4) continue;
    double in = 0, out = 0; 
    for (unsigned int id = 0, nd = j.numberOfDaughters(); id < nd; ++id) {
      const pat::PackedCandidate &dau = dynamic_cast<const pat::PackedCandidate &>(*j.daughter(id));
      if (dau.charge() == 0) continue;
      (fabs(dau.dz())<0.1 ? in : out) += dau.pt();
    }
    //double sum = in + out;
    //printf("Jet with pt %6.1f, eta %+4.2f, beta(0.1) = %+5.3f, pileup mva disc %+.2f\n",
    //	   j.pt(),j.eta(), sum ? in/sum : 0, j.userFloat("pileupJetId:fullDiscriminant"));
  }
}
//=============================================================
//
//     Method for initializing values for the variables
//
//=============================================================
void MakeMuonElectronPatTuple::IntialValues()
{
  
  NbMuons     = 0;
  value3_     = 0;
  NbGenMuons  = 0;
  nMuHLT      = 0;
  NbElectrons = 0;
}
//=============================================================
////
////                Method for ECAL Rec Hits Tree
////
////=============================================================
void MakeMuonElectronPatTuple::EBecalRecHitsTree(const edm::Event& iEvent,const edm::EventSetup& es)
{ 
    ecal_iRechit.clear();
    ecal_EnRecHit.clear();
    ecal_RecoFlag.clear();
    ecal_isSaturated.clear();
    //------------ Rec Hits in EB ----------------
    Handle<EcalRecHitCollection> ecalEB;
    iEvent.getByLabel(tok_EB_,ecalEB);
    const EcalRecHitCollection *ebhits = ecalEB.product();
    for(EcalRecHitCollection::const_iterator it = ebhits->begin(); it != ebhits->end(); ++it) {
     if( (*it).energy() < 10.0 ) continue;
    nbRecHits++;
    ecal_iRechit.push_back(nbRecHits);
    ecal_EnRecHit.push_back( (*it).energy() );
    ecal_RecoFlag.push_back( (*it).recoFlag() );
    if( (*it).recoFlag() == EcalRecHit::kSaturated ) ecal_isSaturated.push_back(true);
    else ecal_isSaturated.push_back(false);
    }
}
void MakeMuonElectronPatTuple::EEecalRecHitsTree(const edm::Event& iEvent,const edm::EventSetup& es)
{ 
  ecal_iRechit.clear();
  ecal_EnRecHit.clear();
  ecal_RecoFlag.clear();
  ecal_isSaturated.clear();
  //------------ Rec Hits in EE ----------------
  Handle<EcalRecHitCollection> ecalEE;
  iEvent.getByLabel(tok_EE_,ecalEE);
  const EcalRecHitCollection *eehits = ecalEE.product();
  for(EcalRecHitCollection::const_iterator it = eehits->begin(); it != eehits->end(); ++it) {
    if( (*it).energy() < 10.0 ) continue;
    nbRecHits++;
    ecal_iRechit.push_back(nbRecHits);
    ecal_EnRecHit.push_back( (*it).energy() );
    ecal_RecoFlag.push_back( (*it).recoFlag() );
    if( (*it).recoFlag() == EcalRecHit::kSaturated ) ecal_isSaturated.push_back(true);
    else ecal_isSaturated.push_back(false);
  }
}







    









