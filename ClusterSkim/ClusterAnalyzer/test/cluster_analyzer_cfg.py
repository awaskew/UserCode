import FWCore.ParameterSet.Config as cms

process = cms.Process("PhoTuple")

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
process.load('Configuration/EventContent/EventContent_cff')
process.GlobalTag.globaltag = cms.string('GR09_R_34X_V4::All')

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
#these are my original CRAB skimmed files.  I don't know why you'd want to go back to them
#but here they are for reference.    
#'dcache:/pnfs/cms/WAX/11/store/user/askew/minbias/Dec14/SkimPhotons_notrig_1.root',
#'dcache:/pnfs/cms/WAX/11/store/user/askew/minbias/Dec14/SkimPhotons_notrig_2.root',
#'dcache:/pnfs/cms/WAX/11/store/user/askew/minbias/Dec14/SkimPhotons_notrig_3.root',
#'dcache:/pnfs/cms/WAX/11/store/user/askew/minbias/Dec14/SkimPhotons_notrig_4.root',
#'dcache:/pnfs/cms/WAX/11/store/user/askew/minbias/Dec14/SkimPhotons_notrig_5.root',
#'dcache:/pnfs/cms/WAX/11/store/user/askew/minbias/Dec14/SkimPhotons_notrig_6.root',
#'dcache:/pnfs/cms/WAX/11/store/user/askew/minbias/Dec14/SkimPhotons_notrig_7.root',
#'dcache:/pnfs/cms/WAX/11/store/user/askew/minbias/Dec14/SkimPhotons_notrig_8.root',
#'dcache:/pnfs/cms/WAX/11/store/user/askew/minbias/Dec14/SkimPhotons_notrig_9.root',
#'dcache:/pnfs/cms/WAX/11/store/user/askew/minbias/Dec14/SkimPhotons_notrig_10.root',
#These are the same files as the other exercises use:
#'dcache:/pnfs/cms/WAX/11/store/user/askew/minbias/EJTermFile/Photest000.root',
#'dcache:/pnfs/cms/WAX/11/store/user/askew/minbias/EJTermFile/Photest001.root',
#'dcache:/pnfs/cms/WAX/11/store/user/askew/minbias/EJTermFile/Photest002.root',
#'dcache:/pnfs/cms/WAX/11/store/user/askew/minbias/EJTermFile/Photest003.root',
#'dcache:/pnfs/cms/WAX/11/store/user/askew/minbias/EJTermFile/Photest004.root',
#'dcache:/pnfs/cms/WAX/11/store/user/askew/minbias/EJTermFile/Photest005.root',
#'dcache:/pnfs/cms/WAX/11/store/user/askew/minbias/EJTermFile/Photest006.root',
#'dcache:/pnfs/cms/WAX/11/store/user/askew/minbias/EJTermFile/Photest007.root',
#'dcache:/pnfs/cms/WAX/11/store/user/askew/minbias/EJTermFile/Photest008.root',
#'dcache:/pnfs/cms/WAX/11/store/user/askew/minbias/EJTermFile/Photest009.root'
#'file:/uscms_data/d2/askew/Bester/CMSSW_3_3_4/src/PhotestLowPTMC.root',
#'file:/uscms_data/d2/askew/Bester/CMSSW_3_3_4/src/PhotestLowPTMC001.root',
#'file:/uscms_data/d2/askew/Bester/CMSSW_3_3_4/src/PhotestLowPTMC002.root',
#'file:/uscms_data/d2/askew/Bester/CMSSW_3_3_4/src/PhotestWithRAW.root',
#'file:/uscms_data/d2/askew/Bester/CMSSW_3_3_4/src/PhotestWithRAW001.root',
#'file:/uscms_data/d2/askew/Bester/CMSSW_3_3_4/src/PhotestWithRAW002.root',
#'file:/uscms_data/d2/askew/Bester/CMSSW_3_3_4/src/PhotestWithRAW003.root',
#'file:/uscms_data/d2/askew/Bester/CMSSW_3_3_4/src/PhotestWithRAW004.root',
#'file:/uscms_data/d2/askew/Bester/CMSSW_3_3_4/src/PhotestWithRAW005.root',
#'file:/uscms_data/d2/askew/Bester/CMSSW_3_3_4/src/PhotestWithRAW006.root',
#'file:/uscms_data/d2/askew/Bester/CMSSW_3_3_4/src/PhotestWithRAW007.root',
#'file:/uscms_data/d2/askew/Bester/CMSSW_3_3_4/src/PhotestWithRAW008.root',
#'file:/uscms_data/d2/askew/Bester/CMSSW_3_3_4/src/PhotestWithRAW009.root'
'file:/uscms_data/d2/askew/ClusSkimComp/CMSSW_3_4_2_patch1/src/EJTermPhotonAnalysis/PhotonLongExercise/test/ClusterSkim.root'
    )
)
process.options = cms.untracked.PSet(
  fileMode = cms.untracked.string('NOMERGE')
)
#if you wanted to drop stuff on input, here's how you could do it.
#process.source.inputCommands = cms.untracked.vstring("keep *", "drop L1GlobalTriggerObjectMapRecord_hltL1GtObjectMap__HLT","drop edmErrorSummaryEntrys_logErrorHarvester__RECO","drop L1GctEmCands_gctDigis_nonIsoEm_RECO","drop L1Gct*_gctDigis_*_RECO","drop L1MuGMTReadoutCollection_gtDigis__RECO","drop l1extraL1EmParticles_l1extraParticles_*_*","drop *_l1extraParticles_*_*","drop recoPhotons_photons__RECO")
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

#this is defined, but turned off by default.
process.Out = cms.OutputModule("PoolOutputModule",
    outputCommands = cms.untracked.vstring('keep *', 
        'keep edmHepMCProduct_*_*_*', 
        'keep recoBasicClusters_*_*_*', 
        'keep recoSuperClusters_*_*_*', 
        'keep *_PhotonIDProd_*_*', 
        'keep *_PhotonIDProd_*_*', 
        'keep recoPhotons_*_*_*'),
    fileName = cms.untracked.string('Photest.root')
)

process.photonIDAna = cms.EDAnalyzer("ClusterAnalyzer",
    outputFile  = cms.string('PhoIDHists2.root'),
    hasRaw = cms.bool(False)
)

#Filter on trigger if you wish.  Also not in the path, but can be changed.
process.hltHighLevel = cms.EDFilter("HLTHighLevel",
    TriggerResultsTag = cms.InputTag("TriggerResults","","HLT"),
    #HLTPaths = cms.vstring("HLT_L1SingleEG2","HLT_L1SingleEG5_NoBPTX"), # provide list of HLT paths (or patterns) you want
    HLTPaths = cms.vstring("HLT_L1SingleEG2"), # provide list of HLT paths (or patterns) you want
    eventSetupPathsKey = cms.string(''),
    andOr = cms.bool(True),             # how to deal with multiple triggers: True (OR) accept if ANY is true, False (AND) accept if ALL are true
    throw = cms.bool(True)    # throw exception on unknown path names
)

#process.p = cms.Path(process.hltHighLevel * process.photonIDAna)
process.p = cms.Path(process.photonIDAna)
#process.e = cms.EndPath(process.Out)
