import FWCore.ParameterSet.Config as cms
from RecoMuon.TrackingTools.MuonServiceProxy_cff import *
from PhysicsTools.PatAlgos.triggerLayer1.triggerProducer_cfi import *


process = cms.Process("tree")
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
process.load("Configuration.StandardSequences.MagneticField_cff")
#process.load("Configuration.StandardSequences.Geometry_cff") #old one, to use for old releases
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")

process.options = cms.untracked.PSet(
    Rethrow = cms.untracked.vstring('ProductNotFound'),
    wantSummary = cms.untracked.bool(True) 
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

# reduce verbosity
#process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(10000)  

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(

'/store/mc/Phys14DR/ZprimeToMuMu_M-5000_Tune4C_13TeV-pythia8/AODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/00000/02BE6C0D-E46E-E411-89C3-003048F0E1B0.root',
#'/store/mc/Phys14DR/ZprimeToMuMu_M-5000_Tune4C_13TeV-pythia8/AODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/00000/02F656DA-B76F-E411-A0E8-002590AC4C74.root',
#'/store/mc/Phys14DR/ZprimeToMuMu_M-5000_Tune4C_13TeV-pythia8/AODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/00000/26FDF897-C96E-E411-A5F8-002590AC4CB8.root',
#'/store/mc/Phys14DR/ZprimeToMuMu_M-5000_Tune4C_13TeV-pythia8/AODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/00000/3A34D96D-D76E-E411-83DC-003048D43960.root',
#'/store/mc/Phys14DR/ZprimeToMuMu_M-5000_Tune4C_13TeV-pythia8/AODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/00000/3A66B45B-D76E-E411-BBDF-002590AC5012.root',
#'/store/mc/Phys14DR/ZprimeToMuMu_M-5000_Tune4C_13TeV-pythia8/AODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/00000/4A9D3E72-DB6E-E411-8AE9-002590DB9286.root',
#'/store/mc/Phys14DR/ZprimeToMuMu_M-5000_Tune4C_13TeV-pythia8/AODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/00000/64A3B071-C66F-E411-9876-002481E0D9BC.root',
#'/store/mc/Phys14DR/ZprimeToMuMu_M-5000_Tune4C_13TeV-pythia8/AODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/00000/68000025-E16E-E411-B132-0025907DC9BA.root',
#'/store/mc/Phys14DR/ZprimeToMuMu_M-5000_Tune4C_13TeV-pythia8/AODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/00000/768CE23C-7E6F-E411-8B35-002590DB92BE.root',
#'/store/mc/Phys14DR/ZprimeToMuMu_M-5000_Tune4C_13TeV-pythia8/AODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/00000/A88E2F73-C66F-E411-BB27-0025901D493E.root',
#'/store/mc/Phys14DR/ZprimeToMuMu_M-5000_Tune4C_13TeV-pythia8/AODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/00000/C2677DF5-E06E-E411-89C3-0025907DC9D0.root',
#'/store/mc/Phys14DR/ZprimeToMuMu_M-5000_Tune4C_13TeV-pythia8/AODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/10000/2C88A05F-826F-E411-AB46-002481E0D9BC.root',
#'/store/mc/Phys14DR/ZprimeToMuMu_M-5000_Tune4C_13TeV-pythia8/AODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/10000/3435D464-FE6E-E411-A209-002590AC4C80.root',
#'/store/mc/Phys14DR/ZprimeToMuMu_M-5000_Tune4C_13TeV-pythia8/AODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/10000/48A8D072-8F6F-E411-85CC-0025904B1366.root',
#'/store/mc/Phys14DR/ZprimeToMuMu_M-5000_Tune4C_13TeV-pythia8/AODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/10000/BCA7188A-8F6F-E411-B365-002590DB9178.root'

 )
)

##-------- Electron events of interest --------
process.HLTEle =cms.EDFilter("HLTHighLevel",
     TriggerResultsTag = cms.InputTag("TriggerResults","","HLT"),
     HLTPaths = cms.vstring("HLT_Mu40_eta2p1_*"),
     eventSetupPathsKey = cms.string(''),
     andOr = cms.bool(True), #----- True = OR, False = AND between the HLTPaths
     throw = cms.bool(False) # throw exception on unknown path names
 )


# Global tag
#from Configuration.AlCa.autoCond import autoCond
#process.GlobalTag.globaltag = cms.string( autoCond[ 'com10' ] )
process.GlobalTag.globaltag = cms.string('PHYS14_25_V1::All')
process.maketreeMuon = cms.EDAnalyzer("MaketreeMuons",
    outputFile             = cms.string('ZprimetoMuMu-MC-Mass5000-CMSSW720-tree.root'),
    globalMuons            = cms.InputTag('muons'),
    globalMuonTracks       = cms.InputTag('globalMuons'),
    vertexCollection       = cms.InputTag('offlinePrimaryVertices'),
    genparticleCollection  = cms.InputTag("genParticles"),                                   
    TrackCollectionTag     = cms.InputTag("generalTracks"),
    # Trigger matching                                           
    triggerEvent          = cms.InputTag("hltTriggerSummaryAOD","","HLT"),
    triggerFilter         = cms.string('HLT_Mu40_eta2p1_*'),
    triggerFilterAsym     =  cms.vstring('hltDiMuonL3PreFiltered8','hltDiMuonL3p5PreFiltered8'),
    maxAbsZ  = cms.double(24),	
    maxd0    = cms.double(2),
    minndof  = cms.int32(4),
    NbGoodPv = cms.int32(1)
)



#process.p = cms.Path(process.HLTEle*process.maketreeMuon)

process.p = cms.Path(process.maketreeMuon)

 
