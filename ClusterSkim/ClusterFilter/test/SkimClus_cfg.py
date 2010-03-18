# Import configurations
import FWCore.ParameterSet.Config as cms

# set up process
process = cms.Process("GMSBSkimWgt")

# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'


process.load('Configuration/StandardSequences/Services_cff')
process.load('FWCore/MessageService/MessageLogger_cfi')
process.load('Configuration/StandardSequences/GeometryExtended_cff')
process.load('Configuration/StandardSequences/MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration/StandardSequences/RawToDigi_Data_cff')
process.load('Configuration/StandardSequences/Reconstruction_cff')
process.load('Configuration/StandardSequences/EndOfProcess_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.load('RecoEcal/EgammaClusterProducers/islandClusteringSequence_cff')
process.load('Configuration/EventContent/EventContent_cff')
process.GlobalTag.globaltag = cms.string('GR09_R_34X_V4::All')


# this defines the input files
#from PhysicsTools.StarterKit.RecoInput_cfi import *
process.source = cms.Source("PoolSource",fileNames = cms.untracked.vstring(
'/store/data/BeamCommissioning09/MinimumBias/RECO/Feb9ReReco_v2/0030/06D02670-8F18-DF11-93F9-002618943856.root','/store/data/BeamCommissioning09/MinimumBias/RECO/Feb9ReReco_v2/0027/FEF0D80F-6816-DF11-8CCB-00261894398A.root','/store/data/BeamCommissioning09/MinimumBias/RECO/Feb9ReReco_v2/0027/F26FEF19-8E16-DF11-8AD9-003048678B38.root'

))
# this inputs the input files from the previous function
#process.source = RecoInput()

# set the number of events
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )



process.clusSkim = cms.EDFilter("ClusterFilter",
                                SCTag = cms.InputTag("correctedHybridSuperClusters"),
                                nHitSel = cms.int32(40),
                                nClusSel = cms.int32(3),
                                ESel = cms.double(100)
                                )
                                

process.p = cms.Path(process.clusSkim)


# setup event content: drop everything before 
process.patEventContent = cms.PSet(
    outputCommands = cms.untracked.vstring('keep *')
)

# define output event selection to be that which satisfies 'p'
process.patEventSelection = cms.PSet(
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('p')
    )
)

# talk to output module
process.out = cms.OutputModule("PoolOutputModule",
    process.patEventSelection,
    process.patEventContent,
    verbose = cms.untracked.bool(False),
    fileName = cms.untracked.string('ClusterSkim.root')
)

# define output path
process.outpath = cms.EndPath(process.out)
