import FWCore.ParameterSet.Config as cms

process = cms.Process("ClustHits")
process.load('Configuration/StandardSequences/Reconstruction_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')

process.GlobalTag.globaltag = 'IDEAL_30X::All'
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.load("Configuration.StandardSequences.Services_cff")
process.load("Configuration.StandardSequences.Geometry_cff")

process.load("RecoEgamma.MaterialConversionModules.ClusterAndHits_cfi")

process.load("RecoTracker.GeometryESProducer.TrackerRecoGeometryESProducer_cfi")

process.load("RecoLocalTracker.SiStripRecHitConverter.SiStripRecHitMatcher_cfi")

process.load("TrackingTools.KalmanUpdators.KFUpdatorESProducer_cfi")

process.load("TrackingTools.TrackFitters.KFTrajectorySmootherESProducer_cfi")

process.load("TrackingTools.TrackFitters.KFTrajectoryFitterESProducer_cfi")

process.load("TrackingTools.TrackFitters.KFFittingSmootherESProducer_cfi")

process.load("RecoTracker.TransientTrackingRecHit.TransientTrackingRecHitBuilder_cfi")

process.load("Configuration.StandardSequences.MagneticField_cff")

process.load("Geometry.CMSCommonData.cmsIdealGeometryXML_cfi")

process.load("Geometry.TrackerGeometryBuilder.trackerGeometry_cfi")

process.load("Geometry.TrackerNumberingBuilder.trackerNumberingGeometry_cfi")

process.load("RecoLocalTracker.SiStripRecHitConverter.StripCPEfromTrackAngle_cfi")

process.load("SimTracker.SiStripDigitizer.SiStripDigi_cfi")

process.load("SimTracker.SiPixelDigitizer.PixelDigi_cfi")

process.load("RecoLocalTracker.SiPixelClusterizer.SiPixelClusterizer_cfi")

process.load("RecoLocalTracker.SiStripClusterizer.SiStripClusterizer_cfi")

process.load("RecoLocalTracker.SiPixelRecHits.SiPixelRecHits_cfi")

process.load("RecoLocalTracker.SiStripRecHitConverter.StripCPE_cfi")

process.load("RecoLocalTracker.SiStripRecHitConverter.SiStripRecHitConverter_cfi")
process.load("RecoLocalTracker.SiPixelRecHits.SiPixelRecHits_cfi")
process.load("RecoLocalTracker.SiPixelRecHits.PixelCPEParmError_cfi")

process.load("RecoTracker.MeasurementDet.MeasurementTrackerESProducer_cfi")

process.ClusterAndHitsAna = cms.EDAnalyzer("ClusterAndHitsAnalyzer",
    outputFile = cms.string('blarg.root'),
    pixelRecHits = cms.InputTag("siPixelRecHits"),
    matchedStripRecHits = cms.InputTag("siStripMatchedRecHits","matchedRecHit"),
    stereoStripRecHits = cms.InputTag("siStripMatchedRecHits","stereoRecHit"),
    rphiStripRecHits = cms.InputTag("siStripMatchedRecHits","rphiRecHit"),
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:/uscmst1b_scratch/lpc1/lpceg/askew/TrackDev/CMSSW_1_8_0_pre9/src/RecoEgamma/MaterialConversionModules/test/00B473D5-A2DB-DC11-99C4-001617C3B6C6.root')
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
)
process.Out = cms.OutputModule("PoolOutputModule",
    outputCommands = cms.untracked.vstring('drop *', 
        'keep edmHepMCProduct_*_*_*', 
        'keep recoBasicClusters_*_*_*', 
        'keep recoSuperClusters_*_*_*', 
        'keep *_Pi0ConversionModule_*_*', 
        'keep *_ClusterAndHitsProd_*_*'),
    fileName = cms.untracked.string('Pi0_test.root')
)

process.p = cms.Path(process.siStripMatchedRecHits*process.siPixelRecHits*process.ClusterAndHitsProd*process.ClusterAndHitsAna)

process.PoolSource.fileNames = ['file:1841D7BE-43E8-DD11-B0B6-001D09F25401.root']


