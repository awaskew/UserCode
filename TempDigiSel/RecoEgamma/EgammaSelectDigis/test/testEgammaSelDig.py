import FWCore.ParameterSet.Config as cms

process = cms.Process("PhotonIDProc")
process.load("FWCore.MessageLogger.MessageLogger_cfi")

process.load("RecoEgamma.EgammaSelectDigis.select_digis_cff")

process.load("Geometry.CaloEventSetup.CaloGeometry_cfi")

process.load("Geometry.CMSCommonData.cmsIdealGeometryXML_cfi")

process.load("Geometry.CaloEventSetup.CaloTopology_cfi")

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
    input = cms.untracked.int32(-1)
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
process.PoolSource.fileNames = ['/store/relval/2008/6/6/RelVal-RelValTTbar-1212531852-IDEAL_V1-2nd-02/0000/081018D5-EC33-DD11-A623-000423D6CA42.root', '/store/relval/2008/6/6/RelVal-RelValTTbar-1212531852-IDEAL_V1-2nd-02/0000/5233E17D-EF33-DD11-8B13-000423D98DB4.root', '/store/relval/2008/6/6/RelVal-RelValTTbar-1212531852-IDEAL_V1-2nd-02/0000/6C3F68D8-EC33-DD11-A624-000423D98804.root', '/store/relval/2008/6/6/RelVal-RelValTTbar-1212531852-IDEAL_V1-2nd-02/0000/6E4050D4-EC33-DD11-85C8-000423D6B48C.root', '/store/relval/2008/6/6/RelVal-RelValTTbar-1212531852-IDEAL_V1-2nd-02/0000/70C3BD68-ED33-DD11-BBCD-001617E30D40.root', 
    '/store/relval/2008/6/6/RelVal-RelValTTbar-1212531852-IDEAL_V1-2nd-02/0000/783C8A02-ED33-DD11-A1DD-00161757BF42.root', '/store/relval/2008/6/6/RelVal-RelValTTbar-1212531852-IDEAL_V1-2nd-02/0000/8006D44B-EC33-DD11-B534-001617DBD230.root', '/store/relval/2008/6/6/RelVal-RelValTTbar-1212531852-IDEAL_V1-2nd-02/0000/84C08904-ED33-DD11-93B0-000423D9870C.root', '/store/relval/2008/6/6/RelVal-RelValTTbar-1212531852-IDEAL_V1-2nd-02/0000/941F9165-ED33-DD11-973E-001617C3B64C.root', '/store/relval/2008/6/6/RelVal-RelValTTbar-1212531852-IDEAL_V1-2nd-02/0000/96AEE6C5-F133-DD11-B794-000423D6CAF2.root', 
    '/store/relval/2008/6/6/RelVal-RelValTTbar-1212531852-IDEAL_V1-2nd-02/0000/AAD8494E-EB33-DD11-B34F-00161757BF42.root', '/store/relval/2008/6/6/RelVal-RelValTTbar-1212531852-IDEAL_V1-2nd-02/0000/C22848D7-EC33-DD11-8BF0-000423D98DB4.root', '/store/relval/2008/6/6/RelVal-RelValTTbar-1212531852-IDEAL_V1-2nd-02/0000/D63078D6-EC33-DD11-A9D7-000423D9939C.root', '/store/relval/2008/6/6/RelVal-RelValTTbar-1212531852-IDEAL_V1-2nd-02/0000/D881F5B8-EE33-DD11-9679-000423D9939C.root', '/store/relval/2008/6/6/RelVal-RelValTTbar-1212531852-IDEAL_V1-2nd-02/0000/DAF2BD43-EC33-DD11-B85A-000423DD2F34.root', 
    '/store/relval/2008/6/6/RelVal-RelValTTbar-1212531852-IDEAL_V1-2nd-02/0000/EA41DDD9-EB33-DD11-9520-001617DBD5AC.root', '/store/relval/2008/6/6/RelVal-RelValTTbar-1212531852-IDEAL_V1-2nd-02/0000/ECCE3FE3-EB33-DD11-80A9-000423D986A8.root', '/store/relval/2008/6/6/RelVal-RelValTTbar-1212531852-IDEAL_V1-2nd-02/0000/F20BA363-ED33-DD11-B6D8-001617DF785A.root', '/store/relval/2008/6/6/RelVal-RelValTTbar-1212531852-IDEAL_V1-2nd-02/0000/F2B9DA5F-ED33-DD11-B82D-001617DBD5AC.root', '/store/relval/2008/6/6/RelVal-RelValTTbar-1212531852-IDEAL_V1-2nd-02/0000/FC8A8764-ED33-DD11-90DB-001617DBD5B2.root']
process.ecalDigis.DoRegional = False
process.ecalDigis.InputLabel = 'rawDataCollector'
process.ecalPreshowerDigis.Label = 'rawDataCollector'


