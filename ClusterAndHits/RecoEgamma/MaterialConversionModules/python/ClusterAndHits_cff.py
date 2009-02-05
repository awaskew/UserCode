import FWCore.ParameterSet.Config as cms

from RecoEgamma.MaterialConversionModules.ClusterAndHits_cfi import *
clusterAndHitsSequence = cms.Sequence(ClusterAndHitsProd*ClusterAndHitsAna)

