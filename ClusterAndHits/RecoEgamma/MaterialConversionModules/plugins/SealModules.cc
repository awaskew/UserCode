#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/PluginManager/interface/ModuleDef.h"
DEFINE_SEAL_MODULE();
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "RecoEgamma/MaterialConversionModules/plugins/ClusterAndHitsAnalyzer.h"
#include "RecoEgamma/MaterialConversionModules/plugins/ClusterAndHitsProducer.h"
DEFINE_ANOTHER_FWK_MODULE(ClusterAndHitsProducer);
DEFINE_ANOTHER_FWK_MODULE(ClusterAndHitsAnalyzer);

