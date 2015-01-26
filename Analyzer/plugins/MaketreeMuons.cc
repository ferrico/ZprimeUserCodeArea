//=========================================================================  
//      Make the root tree for Z boson to Mu Mu analysis                  =  
//                                                                        =
//                         CMSSW_7_2_0                                    =
//                                                                        =
//                 Written by Sherif Elgammal                             =
//                                                                        =
//                         20/12/2014                                     =
//=========================================================================
#include "MyCodeArea/Analyzer/interface/MaketreeMuons.h"
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
using namespace edm;
using namespace std;
using namespace reco;

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

//========================================================================
MaketreeMuons::MaketreeMuons( const edm::ParameterSet& ps )
//========================================================================
{
  vertexSrc                      = ps.getParameter<edm::InputTag>("vertexCollection");
  maxAbsZ_                       = ps.getParameter<double>("maxAbsZ");
  maxd0_                         = ps.getParameter<double>("maxd0");
  minNdof_                       = ps.getParameter<int>("minndof");
  NbGoodPv_                      = ps.getParameter<int>("NbGoodPv");
  globalMuons_                   = ps.getParameter<edm::InputTag>("globalMuons");
  genParticlesColl_              = ps.getParameter<edm::InputTag>("genparticleCollection");
  token_globalMuons              = ps.getParameter<edm::InputTag>( "globalMuonTracks" );
  //inputTag_                      = ps.getParameter<edm::InputTag> ("TriggerResultsTag");
  triggerEvent                   = ps.getParameter<edm::InputTag>("triggerEvent");
  triggerFilter                  = ps.getParameter<std::string>("triggerFilter");
  triggerFilter_asym             = ps.getParameter<std::vector<std::string> >("triggerFilterAsym");


  // get names from module parameters, then derive slot numbers
  n_                  = 0;
  firstevent_         = true;  
  //produces<vector<std::string> >();

  // PG and FRC 06-07-11
  debug	=	ps.getUntrackedParameter<bool> ("debug", false);
  

  //pfMuons_                       = ps.getParameter<edm::InputTag>("pfMuons");
  //pfMuonToken_                   = ps.getParameter<edm::InputTag>("PFjetsCollection");
  //tokentracks_                   = ps.getParameter<edm::InputTag>("TrackCollectionTag");
  //theRecoLabel_                  = ps.getParameter<edm::InputTag>("inputTagMuonReco");
  //===============================================================================================
  outputFile_ = ps.getParameter<std::string>("outputFile");
  rootFile_   = TFile::Open(outputFile_.c_str(),"RECREATE"); // open output file to store histograms
}

 
//========================================================================
MaketreeMuons::~MaketreeMuons()
//========================================================================    
{
  delete rootFile_;
}

//========================================================================
void MaketreeMuons::beginJob() {
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
  //=============================================================
  //
  //           Create Branchs for Muons variables
  //
  //=============================================================
  mytree->Branch("Mu_nbMuon",&Mu_nbMuon);
  mytree->Branch("Mu_isGlobalMuon",&Mu_isGlobalMuon);
  mytree->Branch("Mu_isPF",&Mu_isPF);
  mytree->Branch("Mu_isGoodMuon",&Mu_isGoodMuon);
  mytree->Branch("Mu_et",&Mu_et);
  mytree->Branch("Mu_en",&Mu_en);
  // part for old TuneP track
  mytree->Branch("Mu_ptcocktail",&Mu_ptcocktail);
  mytree->Branch("Mu_etaCocktail",&Mu_etaCocktail);
  mytree->Branch("Mu_phiCocktail",&Mu_phiCocktail);
  mytree->Branch("Mu_thetaCocktail",&Mu_thetaCocktail);  
  mytree->Branch("Mu_pxCocktail",&Mu_pxCocktail);
  mytree->Branch("Mu_pyCocktail",&Mu_pyCocktail);
  mytree->Branch("Mu_pzCocktail",&Mu_pzCocktail);
  mytree->Branch("Mu_pCocktail",&Mu_pCocktail);
  mytree->Branch("Mu_dPToverPTcocktail",&Mu_dPToverPTcocktail);
  mytree->Branch("Mu_chargeCocktail",&Mu_chargeCocktail);
  mytree->Branch("Mu_absdxy",&Mu_absdxy);
  mytree->Branch("Mu_absdz",&Mu_absdz);
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
  mytree->Branch("Mu_ptInnerTrack",&Mu_ptInnerTrack);  
  // part for old TuneP Muon Best Track
  mytree->Branch("Mu_ptTunePMuonBestTrack",&Mu_ptTunePMuonBestTrack);
  mytree->Branch("Mu_dPToverPTTunePMuonBestTrack",&Mu_dPToverPTTunePMuonBestTrack);
  mytree->Branch("Mu_pxTunePMuonBestTrack",&Mu_pxTunePMuonBestTrack);
  mytree->Branch("Mu_pyTunePMuonBestTrack",&Mu_pyTunePMuonBestTrack);
  mytree->Branch("Mu_pzTunePMuonBestTrack",&Mu_pzTunePMuonBestTrack);
  mytree->Branch("Mu_pTunePMuonBestTrack",&Mu_pTunePMuonBestTrack);
  mytree->Branch("Mu_etaTunePMuonBestTrack",&Mu_etaTunePMuonBestTrack);
  mytree->Branch("Mu_phiTunePMuonBestTrack",&Mu_phiTunePMuonBestTrack);
  mytree->Branch("Mu_thetaTunePMuonBestTrack",&Mu_thetaTunePMuonBestTrack);
  mytree->Branch("Mu_chargeTunePMuonBestTrack",&Mu_chargeTunePMuonBestTrack);
  mytree->Branch("Mu_absdxyTunePMuonBestTrack",&Mu_absdxyTunePMuonBestTrack);
  // part for DYT track
  mytree->Branch("Mu_ptDYTTrack",&Mu_ptDYTTrack);
  mytree->Branch("Mu_pxDYTTrack",&Mu_pxDYTTrack);
  mytree->Branch("Mu_pyDYTTrack",&Mu_pyDYTTrack);
  mytree->Branch("Mu_pzDYTTrack",&Mu_pzDYTTrack);
  mytree->Branch("Mu_pDYTTrack",&Mu_pDYTTrack);
  mytree->Branch("Mu_etaDYTTrack",&Mu_etaDYTTrack);
  mytree->Branch("Mu_phiDYTTrack",&Mu_phiDYTTrack);
  mytree->Branch("Mu_thetaDYTTrack",&Mu_thetaDYTTrack);
  mytree->Branch("Mu_chargeDYTTrack",&Mu_chargeDYTTrack);
  mytree->Branch("Mu_absdxyDYTTrack",&Mu_absdxyDYTTrack);
  mytree->Branch("Mu_dPToverPTDYTTrack",&Mu_dPToverPTDYTTrack);
  // part for Picky Track
  mytree->Branch("Mu_ptPickyTrack",&Mu_ptPickyTrack);
  mytree->Branch("Mu_pxPickyTrack",&Mu_pxPickyTrack);
  mytree->Branch("Mu_pyPickyTrack",&Mu_pyPickyTrack);
  mytree->Branch("Mu_pzPickyTrack",&Mu_pzPickyTrack);
  mytree->Branch("Mu_pPickyTrack",&Mu_pPickyTrack);
  mytree->Branch("Mu_etaPickyTrack",&Mu_etaPickyTrack);
  mytree->Branch("Mu_phiPickyTrack",&Mu_phiPickyTrack);
  mytree->Branch("Mu_thetaPickyTrack",&Mu_thetaPickyTrack);
  mytree->Branch("Mu_chargePickyTrack",&Mu_chargePickyTrack);
  mytree->Branch("Mu_absdxyPickyTrack",&Mu_absdxyPickyTrack);
  mytree->Branch("Mu_dPToverPTPickyTrack",&Mu_dPToverPTPickyTrack);
  // part for Muon Best Track
  mytree->Branch("Mu_ptMuonBestTrack",&Mu_ptMuonBestTrack);
  mytree->Branch("Mu_dPToverPTMuonBestTrack",&Mu_dPToverPTMuonBestTrack);
  mytree->Branch("Mu_pxMuonBestTrack",&Mu_pxMuonBestTrack);
  mytree->Branch("Mu_pyMuonBestTrack",&Mu_pyMuonBestTrack);
  mytree->Branch("Mu_pzMuonBestTrack",&Mu_pzMuonBestTrack);
  mytree->Branch("Mu_pMuonBestTrack",&Mu_pMuonBestTrack);
  mytree->Branch("Mu_etaMuonBestTrack",&Mu_etaMuonBestTrack);
  mytree->Branch("Mu_phiMuonBestTrack",&Mu_phiMuonBestTrack);
  mytree->Branch("Mu_thetaMuonBestTrack",&Mu_thetaMuonBestTrack);
  mytree->Branch("Mu_chargeMuonBestTrack",&Mu_chargeMuonBestTrack);
  mytree->Branch("Mu_absdxyMuonBestTrack",&Mu_absdxyMuonBestTrack);
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
  mytree->Branch("NbOfDaughters",&NbOfDaughters);
  mytree->Branch("McZmass",&McZmass);
  //mytree->Branch("McZpt",&McZpt);
  //mytree->Branch("McZpx",&McZpx);
  //mytree->Branch("McZpy",&McZpy);
  //mytree->Branch("McZpz",&McZpz);
  //mytree->Branch("McZen",&McZen);
  //mytree->Branch("massGen",&massGen);
  //mytree->Branch("vxGen",&vxGen);
  //mytree->Branch("vyGen",&vyGen);
  //mytree->Branch("vzGen",&vzGen);
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
  mytree->Branch("MuHLTMatch_nbMuonMatchHLT",&MuHLTMatch_nbMuonMatchHLT);
  mytree->Branch("MuHLTMatch_pt",&MuHLTMatch_pt);
  mytree->Branch("MuHLTMatch_eta",&MuHLTMatch_eta);
  mytree->Branch("MuHLTMatch_phi",&MuHLTMatch_phi);
}
  

//========================================================================================
void MaketreeMuons::analyze( const edm::Event& evt, const edm::EventSetup& es ) {
//======================================================================================
  using namespace edm; // needed for all fwk related classe
  //==============================================
  //=        Begin of the main program           =
  //============================================== 
  TotalNbEvents ++;
  NbEventsPassTrigger.clear();
  NbEventsPassTrigger.push_back(TotalNbEvents);
  //===================== Handle For Muons =====================
  edm::Handle<reco::MuonCollection> muons;
  evt.getByLabel(globalMuons_, muons);
  const reco::MuonCollection* realMuons = muons.product();
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
  eventVertexPosition_ = GlobalPoint(0., 0., 0.);
  const reco::Vertex& thePrimaryEventVertex = (*pvHandle->begin());
  eventVertexPosition_ = GlobalPoint(thePrimaryEventVertex.x(), thePrimaryEventVertex.y(), thePrimaryEventVertex.z());
  Run   = evt.id().run();
  Event = evt.id().event();
  lumi  = evt.luminosityBlock();
  bunch = evt.bunchCrossing();
  IntialValues();
  GenParticleTree(evt);
  PrimaryVertexTree(vertices);
  MuonTree(evt,es,realMuons,thePrimaryEventVertex,vertices);  
  //hlt matching
  TriggerMatchingTree(evt,realMuons);
  //==============================================
  //=        End of the main program             =
  //============================================== 
  mytree->Fill();
}

//========================================================================
void MaketreeMuons::endJob() {
//========================================================================
  // go to *OUR* root file and store histograms
  rootFile_->cd();

  mytree->Write();

  rootFile_->Close();
}

//define this as a plug-in
DEFINE_FWK_MODULE(MaketreeMuons);
//=============================================================
//
//         Method for finding good Primary Vertex
//
//=============================================================
bool MaketreeMuons::PrimaryVertex(const reco::VertexCollection &vertices)
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
//                Method for Reco Muon Tree
//
//=============================================================
void MaketreeMuons::MuonTree(const edm::Event& evt,const edm::EventSetup& es,
			     const reco::MuonCollection* muons,const reco::Vertex &vertex,
			     const reco::VertexCollection &verticess)
{  
  Mu_nbMuon.clear();
  Mu_isGlobalMuon.clear();
  Mu_isPF.clear();
  Mu_isGoodMuon.clear();
  Mu_etaCocktail.clear();
  Mu_ptTunePMuonBestTrack.clear();
  Mu_dPToverPTTunePMuonBestTrack.clear();
  Mu_et.clear();
  Mu_phiCocktail.clear();	
  Mu_en.clear();
  Mu_numberOfMatchedStations.clear();
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
  Mu_chargeCocktail.clear();
  Mu_emIso.clear();
  Mu_hadIso.clear();
  Mu_thetaCocktail.clear();
  Mu_vtxMass.clear();
  Mu_vtxNormChi2.clear();
  Mu_ptcocktail.clear();
  Mu_dPToverPTcocktail.clear();
  Mu_nbofpv.clear();
  Mu_pxCocktail.clear();
  Mu_pyCocktail.clear();
  Mu_pzCocktail.clear();
  Mu_pCocktail.clear();
  Mu_pxTunePMuonBestTrack.clear();
  Mu_pyTunePMuonBestTrack.clear();
  Mu_pzTunePMuonBestTrack.clear();
  Mu_pTunePMuonBestTrack.clear();
  Mu_etaTunePMuonBestTrack.clear();
  Mu_ptInnerTrack.clear();
  Mu_phiTunePMuonBestTrack.clear();
  Mu_thetaTunePMuonBestTrack.clear();
  Mu_chargeTunePMuonBestTrack.clear();
  Mu_absdxyTunePMuonBestTrack.clear();
  Mu_ptDYTTrack.clear();
  Mu_pxDYTTrack.clear();
  Mu_pyDYTTrack.clear();
  Mu_pzDYTTrack.clear();
  Mu_pDYTTrack.clear(); 
  Mu_etaDYTTrack.clear();
  Mu_phiDYTTrack.clear();
  Mu_thetaDYTTrack.clear();
  Mu_chargeDYTTrack.clear();    
  Mu_absdxyDYTTrack.clear();
  Mu_dPToverPTDYTTrack.clear();
  Mu_ptPickyTrack.clear();
  Mu_pxPickyTrack.clear();
  Mu_pyPickyTrack.clear();
  Mu_pzPickyTrack.clear();
  Mu_pPickyTrack.clear(); 
  Mu_etaPickyTrack.clear();
  Mu_phiPickyTrack.clear();
  Mu_thetaPickyTrack.clear();
  Mu_chargePickyTrack.clear();    
  Mu_absdxyPickyTrack.clear();
  Mu_dPToverPTPickyTrack.clear();
  float vtxNormChi2,DiMass;
  Mu_ptMuonBestTrack.clear();
  Mu_dPToverPTMuonBestTrack.clear();
  Mu_pxMuonBestTrack.clear();
  Mu_pyMuonBestTrack.clear();
  Mu_pzMuonBestTrack.clear();
  Mu_pMuonBestTrack.clear(); 
  Mu_etaMuonBestTrack.clear();
  Mu_phiMuonBestTrack.clear();
  Mu_thetaMuonBestTrack.clear();
  Mu_chargeMuonBestTrack.clear();
  Mu_absdxyMuonBestTrack.clear();
  // Get TransientTracks (for use in e.g. the vertex fit) for each of
  // the muon tracks, using e.g. the cocktail momentum.
  edm::ESHandle<TransientTrackBuilder> ttkb;
  es.get<TransientTrackRecord>().get("TransientTrackBuilder", ttkb);
  std::vector<reco::TransientTrack> ttv;
  for (reco::MuonCollection::const_iterator muonIt = muons->begin(); muonIt!=muons->end(); ++muonIt){
    //Varaibles needed for high pt Muons are taken from the twiki below
    //https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonId#HighPT_Muon
    reco::Muon recoMu = *muonIt;
    if (recoMu.innerTrack().isNull()) continue;
    if (recoMu.globalTrack().isNull()) continue;
    if(recoMu.isTrackerMuon()==false) continue;
    NbMuons++;
    Mu_nbMuon.push_back(NbMuons);
    Mu_isGlobalMuon.push_back(recoMu.isGlobalMuon());
    Mu_isPF.push_back(recoMu.isPFMuon());
    Mu_nbofpv.push_back(verticess.size());
    //============== Parameters related to Kinematics =====================
    reco::TrackRef cktTrack = (muon::tevOptimized(recoMu, 200, 17., 40., 0.25)).first; 
    Mu_en.push_back(recoMu.energy());
    Mu_et.push_back(recoMu.et());
    Mu_ptInnerTrack.push_back(recoMu.innerTrack()->pt());
    Mu_ptcocktail.push_back(cktTrack->pt());
    Mu_dPToverPTcocktail.push_back(cktTrack->ptError()/cktTrack->pt());
    Mu_absdxy.push_back(fabs(cktTrack->dxy(vertex.position())));
    Mu_absdz.push_back(fabs(cktTrack->dz(vertex.position())));
    Mu_etaCocktail.push_back(cktTrack->eta());
    Mu_phiCocktail.push_back(cktTrack->phi());
    Mu_thetaCocktail.push_back(cktTrack->theta());
    Mu_chargeCocktail.push_back(cktTrack->charge());
    Mu_pxCocktail.push_back(cktTrack->px()); //px component of the track 
    Mu_pyCocktail.push_back(cktTrack->py()); //py component of the track 
    Mu_pzCocktail.push_back(cktTrack->pz()); //pz component of the track 
    Mu_pCocktail.push_back(cktTrack->p());   //magnitude of momentum vector
    // part for TuneP Muon Best track
    Mu_ptTunePMuonBestTrack.push_back(recoMu.tunePMuonBestTrack()->pt());
    Mu_pxTunePMuonBestTrack.push_back(recoMu.tunePMuonBestTrack()->px()); //px component of the track 
    Mu_pyTunePMuonBestTrack.push_back(recoMu.tunePMuonBestTrack()->py()); //py component of the track 
    Mu_pzTunePMuonBestTrack.push_back(recoMu.tunePMuonBestTrack()->pz()); //pz component of the track 
    Mu_pTunePMuonBestTrack.push_back(recoMu.tunePMuonBestTrack()->p());   //magnitude of momentum vector 
    Mu_etaTunePMuonBestTrack.push_back(recoMu.tunePMuonBestTrack()->eta());
    Mu_phiTunePMuonBestTrack.push_back(recoMu.tunePMuonBestTrack()->phi());
    Mu_thetaTunePMuonBestTrack.push_back(recoMu.tunePMuonBestTrack()->theta());
    Mu_chargeTunePMuonBestTrack.push_back(recoMu.tunePMuonBestTrack()->charge());
    Mu_dPToverPTTunePMuonBestTrack.push_back(recoMu.tunePMuonBestTrack()->ptError()/recoMu.tunePMuonBestTrack()->pt());
    Mu_absdxyTunePMuonBestTrack.push_back(fabs(recoMu.tunePMuonBestTrack()->dxy(vertex.position())));
    // part for Muon Best track
    Mu_ptMuonBestTrack.push_back(recoMu.muonBestTrack()->pt());
    Mu_pxMuonBestTrack.push_back(recoMu.muonBestTrack()->px()); //px component of the track 
    Mu_pyMuonBestTrack.push_back(recoMu.muonBestTrack()->py()); //py component of the track 
    Mu_pzMuonBestTrack.push_back(recoMu.muonBestTrack()->pz()); //pz component of the track 
    Mu_pMuonBestTrack.push_back(recoMu.muonBestTrack()->p());   //magnitude of momentum vector 
    Mu_etaMuonBestTrack.push_back(recoMu.muonBestTrack()->eta());
    Mu_phiMuonBestTrack.push_back(recoMu.muonBestTrack()->phi());
    Mu_thetaMuonBestTrack.push_back(recoMu.muonBestTrack()->theta());
    Mu_chargeMuonBestTrack.push_back(recoMu.muonBestTrack()->charge());
    Mu_dPToverPTMuonBestTrack.push_back(recoMu.muonBestTrack()->ptError()/recoMu.muonBestTrack()->pt());
    Mu_absdxyMuonBestTrack.push_back(fabs(recoMu.muonBestTrack()->dxy(vertex.position())));
    // part for DYT track
    /*
    Mu_ptDYTTrack.push_back(recoMu.dytTrack()->pt());
    Mu_pxDYTTrack.push_back(recoMu.dytTrack()->px());
    Mu_pyDYTTrack.push_back(recoMu.dytTrack()->py());
    Mu_pzDYTTrack.push_back(recoMu.dytTrack()->pz());
    Mu_pDYTTrack.push_back(recoMu.dytTrack()-> p());
    Mu_etaDYTTrack.push_back(recoMu.dytTrack()->eta());
    Mu_phiDYTTrack.push_back(recoMu.dytTrack()->phi());
    Mu_thetaDYTTrack.push_back(recoMu.dytTrack()->theta());
    Mu_chargeDYTTrack.push_back(recoMu.dytTrack()->charge());    
    Mu_absdxyDYTTrack.push_back(fabs(recoMu.dytTrack()->dxy(vertex.position())));
    Mu_dPToverPTDYTTrack.push_back(recoMu.dytTrack()->ptError()/recoMu.dytTrack()->pt());
    */
    // part for Picky track                          
    /*
    Mu_ptDYTTrack.push_back(recoMu.pickyTrack()->pt());
    Mu_pxDYTTrack.push_back(recoMu.pickyTrack()->px());
    Mu_pyDYTTrack.push_back(recoMu.pickyTrack()->py());
    Mu_pzDYTTrack.push_back(recoMu.pickyTrack()->pz());
    Mu_pDYTTrack.push_back(recoMu.pickyTrack()-> p());
    Mu_etaDYTTrack.push_back(recoMu.pickyTrack()->eta());
    Mu_phiDYTTrack.push_back(recoMu.pickyTrack()->phi());
    Mu_thetaDYTTrack.push_back(recoMu.pickyTrack()->theta());
    Mu_chargeDYTTrack.push_back(recoMu.pickyTrack()->charge());    
    Mu_absdxyDYTTrack.push_back(fabs(recoMu.pickyTrack()->dxy(vertex.position())));
    Mu_dPToverPTDYTTrack.push_back(recoMu.pickyTrack()->ptError()/recoMu.pickyTrack()->pt());
    */
    //====================== Parameters related to track quality =====================
    Mu_normalizedChi2.push_back(recoMu.globalTrack()->normalizedChi2());
    Mu_numberOfValidPixelHits.push_back(recoMu.globalTrack()->hitPattern().numberOfValidPixelHits());
    Mu_numberOfValidMuonHits.push_back(recoMu.globalTrack()->hitPattern().numberOfValidMuonHits());
    Mu_numberOftrackerLayersWithMeasurement.push_back(recoMu.globalTrack()->hitPattern().trackerLayersWithMeasurement());
    Mu_numberOfMatchedStations.push_back(recoMu.numberOfMatchedStations());
    //============= Parameters related to detector isolation =====================
    Mu_emIso.push_back(recoMu.isolationR03().emEt);
    Mu_hadIso.push_back(recoMu.isolationR03().hadEt);
    Mu_trackiso.push_back(recoMu.isolationR03().sumPt);
    //============= Parameters related to PF isolation =====================
    Mu_pfSumChargedHadronPt.push_back(recoMu.pfIsolationR03().sumChargedHadronPt);
    Mu_pfSumNeutralHadronEt.push_back(recoMu.pfIsolationR03().sumNeutralHadronEt);
    Mu_PFSumPhotonEt.push_back(recoMu.pfIsolationR03().sumPhotonEt);
    Mu_pfSumPUPt.push_back(recoMu.pfIsolationR03().sumPUPt);
    //------------------------------------------------------
    //   part to compute mass from vertex info             -   
    //   Given to me by Piotr                              -
    //------------------------------------------------------
    //printf ("[muonsNb %d] pt = %f\n",NbMuons,cktTrack->pt());
    if( NbMuons > 0 && cktTrack->pt()>45.0 && (recoMu.isolationR03().sumPt/recoMu.innerTrack()->pt()<0.10) 
    	&& (cktTrack->ptError()/cktTrack->pt()<0.3) && 
        fabs(cktTrack->dxy(vertex.position()))<0.2 &&
        recoMu.globalTrack()->hitPattern().trackerLayersWithMeasurement()>5 && 
        recoMu.globalTrack()->hitPattern().numberOfValidPixelHits()>0 &&
        recoMu.globalTrack()->hitPattern().numberOfValidMuonHits()>0 &&
        recoMu.numberOfMatchedStations()>1 ){
      ttv.push_back(ttkb->build(cktTrack));
      //cout<<"size(outside) = "<<ttv.size()<<endl;
      if(ttv.size()>1){
      //if(ttv.size()==2){
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
	//cout<<"size(inside) = "<<ttv.size()<<endl;
	//printf ("normChi2 = %f mass = %f\n",vtxNormChi2,DiMass);
      }
    }
  }
}
//=============================================================
//
//            Method for Genrated Particles Tree
//
//=============================================================
void MaketreeMuons::GenParticleTree(const edm::Event& evt){
  edm::Handle<GenParticleCollection> genParticles;
  evt.getByLabel(genParticlesColl_, genParticles);

   if (!(genParticles.isValid())) return;  

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
  NbOfDaughters.clear();
  McZmass.clear();
  McZpt.clear();
  //massGen.clear();
  //vxGen.clear();
  //vyGen.clear();
  //vzGen.clear();
  //McZpx.clear();
  //McZpy.clear();
  //McZpz.clear();
  //McZen.clear();
  for(size_t i = 0; i < genParticles->size(); ++ i) {
    const GenParticle & p = (*genParticles)[i];
    if( fabs(p.pdgId()) != 13 ) continue; //consider only mu
    if( p.status() > 3 ) continue;        //consider only status 1,3
    //if( p.numberOfDaughters() < 2 ) continue;
    value2_++;
    iGen.push_back(value2_);
    math::XYZTLorentzVector Zboson;
    Zboson = p.p4();
    idGen.push_back(p.pdgId());           //int id = 
    statusGen.push_back(p.status());      //int st = 
    ptGen.push_back(p.pt());
    etaGen.push_back(p.eta());
    phiGen.push_back(p.phi());
    chargeGen.push_back(p.charge());
    EnergyGen.push_back(p.energy());
    pxGen.push_back(p.px());
    pyGen.push_back(p.py());
    pzGen.push_back(p.pz());
    NbOfDaughters.push_back(p.numberOfDaughters());
    McZmass.push_back(Zboson.M()); 
    //ROOT::Math::Boost CoMBoost(Zboson.BoostToCM());
    //const Candidate * mom = p.mother();
    //massGen.push_back(p.mass());
    //vxGen.push_back(p.vx());
    //vyGen.push_back(p.vy());
    //vzGen.push_back(p.vz());
    //McZpt.push_back(Zboson.Pt()); 
    //McZpx.push_back(Zboson.Px()); 
    //McZpy.push_back(Zboson.Py()); 
    //McZpz.push_back(Zboson.Pz());
    //McZen.push_back(Zboson.E());
  }
}
//=============================================================
//
//            Method for Pimary Vertices Tree
//
//=============================================================
void MaketreeMuons::PrimaryVertexTree(const reco::VertexCollection &vertices)
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
bool MaketreeMuons::IsMuMatchedToHLTMu ( const reco::Muon &mu, std::vector<reco::Particle> HLTMu , double DR, double DPtRel) {
  size_t dim =  HLTMu.size();
  size_t nPass=0;
  if (dim==0) return false;
  for (size_t k =0; k< dim; k++ ) {
    if (  (deltaR(HLTMu[k], mu) < DR)   && (fabs(HLTMu[k].pt() - mu.innerTrack()->pt())/ HLTMu[k].pt()<DPtRel)){ 
      nPass++ ;
    }
  }
  return (nPass>0);
}

void MaketreeMuons::TriggerMatchingTree(const edm::Event& iEvent,const reco::MuonCollection* muons){
  double maxDeltaR_ = 0.2;
  double maxDPtRel_ = 1.0;
  //cout << "Start Trigger matching for muon" << endl;
  MuHLTMatch_pt.clear();
  MuHLTMatch_eta.clear();
  MuHLTMatch_phi.clear();
  MuHLTMatch_nbMuonMatchHLT.clear();
  std::vector<reco::Particle>  HLTMuMatched;
  edm::Handle<trigger::TriggerEvent> handleTriggerEvent;
  iEvent.getByLabel(triggerEvent, handleTriggerEvent );
  const trigger::TriggerObjectCollection & toc(handleTriggerEvent->getObjects());
  size_t nMuHLT =0;
  for ( size_t ia = 0; ia < handleTriggerEvent->sizeFilters(); ++ ia) {
    const trigger::Keys & k = handleTriggerEvent->filterKeys(ia);
    for (trigger::Keys::const_iterator ki = k.begin(); ki !=k.end(); ++ki ) {
      HLTMuMatched.push_back(toc[*ki].particle());
    }
  }
  
  for (reco::MuonCollection::const_iterator muonIt = muons->begin(); muonIt!=muons->end(); ++muonIt){
    reco::Muon mu = *muonIt;
    if (mu.innerTrack().isNull()) continue;
    if (mu.globalTrack().isNull()) continue;
    reco::TrackRef cktTrack = (muon::tevOptimized(mu, 200, 17., 40., 0.25)).first;
    if (cktTrack->pt()<40.0) continue;
    if (fabs(cktTrack->eta())>2.1) continue;
    if (IsMuMatchedToHLTMu(mu,HLTMuMatched,maxDeltaR_,maxDPtRel_)==true){
      nMuHLT++;
      if(nMuHLT > 1) continue;
      //cout << "Muon HLT Matched"  << endl;
      MuHLTMatch_nbMuonMatchHLT.push_back(nMuHLT);
      MuHLTMatch_pt.push_back(cktTrack->pt());
      MuHLTMatch_eta.push_back(cktTrack->eta());
      MuHLTMatch_phi.push_back(cktTrack->phi());
    }
  }
}
//=============================================================
//
//     Method for initializing values for the variables
//
//=============================================================
void MaketreeMuons::IntialValues()
{
  
  NbMuons   = 0;
  value2_   = 0;
  value3_   = 0;
  nbTk      = 0;
}
    
