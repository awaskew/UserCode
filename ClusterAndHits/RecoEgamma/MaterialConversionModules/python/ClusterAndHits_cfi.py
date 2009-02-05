import FWCore.ParameterSet.Config as cms

ClusterAndHitsProd = cms.EDProducer("ClusterAndHitsProducer",
    clusterMatchedRecHitsColl = cms.string('savedMatchwithEcalSuperClus'),
    clusterandhitsCollection = cms.string('ClusterAndHits'),
    clusterRPhiRecHitsColl = cms.string('savedRPhiwithEcalSuperClus'),
    pixelRecHits = cms.InputTag("siPixelRecHits"),
    siClusterColl = cms.string('savedSiStripClusters'),
    hybridsuperclusterCollection = cms.string(''),
    siPixClusterColl = cms.string('savedSiPixelClusters'),
    matchedStripRecHits = cms.InputTag("siStripMatchedRecHits","matchedRecHit"),
    debugLevel = cms.string('INFO'),
    clusterStereoRecHitsColl = cms.string('savedSterwithEcalSuperClus'),
    rphiStripRecHits = cms.InputTag("siStripMatchedRecHits","rphiRecHit"),
    clusterPixelRecHitsColl = cms.string('savedPixwithEcalSuperClus'),
    hybridsuperclusterProducer = cms.string('correctedHybridSuperClusters'),
    stereoStripRecHits = cms.InputTag("siStripMatchedRecHits","stereoRecHit"),
    cluster_pt_thresh = cms.double(15.0)
)




