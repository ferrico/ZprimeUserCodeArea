import FWCore.ParameterSet.Config as cms
from RecoMuon.TrackingTools.MuonServiceProxy_cff import *
from PhysicsTools.PatAlgos.triggerLayer1.triggerProducer_cfi import *
process = cms.Process("tree")
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.options = cms.untracked.PSet(
    Rethrow = cms.untracked.vstring('ProductNotFound'),
    wantSummary = cms.untracked.bool(True) 
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10000),
    SkipEvent = cms.untracked.vstring('ProductNotFound')
)

# reduce verbosity
#process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(10000)  

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(

'/store/mc/Phys14DR/DYToMuMu_M-50_Tune4C_13TeV-pythia8/MINIAODSIM/PU40bx25_tsg_castor_PHYS14_25_V1-v2/00000/622CAFBA-BD9A-E411-BE11-002481E14FFC.root'



#'/store/mc/Phys14DR/ZprimeToMuMu_M-5000_Tune4C_13TeV-pythia8/MINIAODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/00000/788396C0-9D6F-E411-97DF-002590494E34.root',
#'/store/mc/Phys14DR/ZprimeToMuMu_M-5000_Tune4C_13TeV-pythia8/MINIAODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/00000/90575C8B-D56F-E411-A39B-0025904B1420.root',
#'/store/mc/Phys14DR/ZprimeToMuMu_M-5000_Tune4C_13TeV-pythia8/MINIAODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/00000/FCEF82CA-9D6F-E411-A2FC-002481E0D5CE.root',
#'/store/mc/Phys14DR/ZprimeToMuMu_M-5000_Tune4C_13TeV-pythia8/MINIAODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/10000/44B8EEE7-9A6F-E411-8858-002590DB0640.root',
#'/store/mc/Phys14DR/ZprimeToMuMu_M-5000_Tune4C_13TeV-pythia8/MINIAODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/10000/F03F87D0-9A6F-E411-8BC4-00266CFFA5E0.root'
#'/store/mc/Spring14miniaod/TT_Tune4C_13TeV-pythia8-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/063013AD-9907-E411-8135-0026189438BD.root'

 )
)

##-------- Electron events of interest --------
process.HLTEle =cms.EDFilter("HLTHighLevel",
     TriggerResultsTag = cms.InputTag("TriggerResults","","HLT"),
     HLTPaths = cms.vstring("HLT_Mu40_v*"),
     eventSetupPathsKey = cms.string(''),
     andOr = cms.bool(True), #----- True = OR, False = AND between the HLTPaths
     throw = cms.bool(False) # throw exception on unknown path names
 )

#
# rho value for isolation
#
#from RecoJets.JetProducers.kt4PFJets_cfi import *
#process.kt6PFJetsForIsolation = kt4PFJets.clone( rParam = 0.6, doRhoFastjet = True )
#process.kt6PFJetsForIsolation.Rho_EtaMax = cms.double(2.5)

# Global tag
#from Configuration.AlCa.autoCond import autoCond
#process.GlobalTag.globaltag = cms.string( autoCond[ 'com10' ] )
process.GlobalTag.globaltag = cms.string('PHYS14_25_V1::All')
process.makeMuonElectronPatTuple = cms.EDAnalyzer("MakeMuonElectronPatTuple",
    outputFile              = cms.string('CMSSW720_MC_Zprime_13TeV_pattupleNew.root'),
    electrons               = cms.InputTag("slimmedElectrons"),
    EBrecHitCollectionLabel = cms.InputTag("reducedEgamma","reducedEBRecHits"),
    EErecHitCollectionLabel = cms.InputTag("reducedEgamma","reducedEERecHits"),
    muons                   = cms.InputTag('slimmedMuons'),
    GenBosonID              = cms.int32(32),
    packed                  = cms.InputTag("packedGenParticles"),
    pruned                  = cms.InputTag("prunedGenParticles"),
    vertexCollection        = cms.InputTag('offlineSlimmedPrimaryVertices'),
    maxAbsZ                 = cms.double(24),	
    maxd0                   = cms.double(2),
    minndof                 = cms.int32(4),
    NbGoodPv                = cms.int32(1),
    jets                    = cms.InputTag("slimmedJets"),
    rhoIsoInputTag          = cms.InputTag("fixedGridRhoFastjetCentralCalo"),
    pfCands                 = cms.InputTag("packedPFCandidates"),
    # Trigger matching
    bits                    = cms.InputTag("TriggerResults","","HLT"),
    prescales               = cms.InputTag("patTrigger"),
    objects                 = cms.InputTag("selectedPatTrigger"),
    # selections of trigger objects
    matchedCuts = cms.string( 'type( "TriggerMuon" ) && path( "HLT_Mu40_v*" )' )

    #TrackCollectionTag     = cms.InputTag("generalTracks"),
    #triggerEvent           = cms.InputTag("hltTriggerSummaryAOD","","HLT"),
    #triggerFilter          = cms.string('HLT_Mu40_eta2p1_*'),
    #triggerFilterAsym      =  cms.vstring('hltDiMuonL3PreFiltered8','hltDiMuonL3p5PreFiltered8'),
    
)



#process.p = cms.Path(process.HLTEle*process.makeMuonElectronPatTuple)

process.p = cms.Path(process.makeMuonElectronPatTuple)

 
