import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.FakeConditions_cff import *
process = cms.Process("PhotonIDProc")
process.load("FWCore.MessageLogger.MessageLogger_cfi")

process.load("RecoEgamma.EgammaSelectDigis.select_digis_cff")

process.load("Geometry.CaloEventSetup.CaloGeometry_cfi")

process.load("Geometry.CMSCommonData.cmsIdealGeometryXML_cfi")

process.load("Geometry.CaloEventSetup.CaloTopology_cfi")

process.load("Geometry.EcalMapping.EcalMappingRecord_cfi")
process.load("Geometry.EcalMapping.EcalMapping_cfi")

process.load("SimCalorimetry.EcalTrigPrimProducers.ecalTriggerPrimitiveDigis_cff")

import EventFilter.EcalRawToDigiDev.EcalUnpackerData_cfi
process.ecalDigis = EventFilter.EcalRawToDigiDev.EcalUnpackerData_cfi.ecalEBunpacker.clone()
import EventFilter.ESRawToDigi.esRawToDigi_cfi
process.ecalPreshowerDigis = EventFilter.ESRawToDigi.esRawToDigi_cfi.esRawToDigi.clone()
process.load("EventFilter.EcalRawToDigiDev.EcalUnpackerMapping_cfi")

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:/uscmst1b_scratch/lpc1/lpceg/askew/DigiSel/CMSSW_2_1_0_pre5/src/RecoEgamma/EgammaSelectDigis/test/D4E57D72-ED33-DD11-AA68-000423D6CA42.root')
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1000)
)
process.Out = cms.OutputModule("PoolOutputModule",
    outputCommands = cms.untracked.vstring('drop *', 
        'keep *_selectDigi_*_*'),
    fileName = cms.untracked.string('Photest.root')
)

process.digiAna = cms.EDAnalyzer("EcalDigiSelectAnalyzer",
    outputFile = cms.string('digiHists.root')
)

process.p = cms.Path(process.ecalDigis*process.seldigis*process.digiAna)
process.e = cms.EndPath(process.Out)
process.PoolSource.fileNames = ['/store/relval/CMSSW_2_1_10/RelValTTbar/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/IDEAL_V9_v1/0001/1E04FC31-F99A-DD11-94EE-0018F3D096DE.root',
        '/store/relval/CMSSW_2_1_10/RelValTTbar/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/IDEAL_V9_v1/0001/30CF1A2C-F99A-DD11-9E3B-001A928116FC.root',
        '/store/relval/CMSSW_2_1_10/RelValTTbar/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/IDEAL_V9_v1/0001/4826BC56-319B-DD11-9280-0017312B5F3F.root',
        '/store/relval/CMSSW_2_1_10/RelValTTbar/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/IDEAL_V9_v1/0001/4C720416-F89A-DD11-B745-003048769FDF.root',
        '/store/relval/CMSSW_2_1_10/RelValTTbar/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/IDEAL_V9_v1/0001/4E40F71D-009B-DD11-9350-001731AF684D.root',
        '/store/relval/CMSSW_2_1_10/RelValTTbar/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/IDEAL_V9_v1/0001/54767785-FA9A-DD11-977E-001A92810AE2.root',
        '/store/relval/CMSSW_2_1_10/RelValTTbar/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/IDEAL_V9_v1/0001/569F692D-F99A-DD11-B0B1-0018F3D095FA.root']
process.ecalDigis.DoRegional = False
process.ecalDigis.InputLabel = 'rawDataCollector'
process.ecalPreshowerDigis.sourceTag = 'rawDataCollector'


