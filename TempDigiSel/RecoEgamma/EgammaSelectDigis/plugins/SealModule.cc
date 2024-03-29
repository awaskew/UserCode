#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/PluginManager/interface/ModuleDef.h"
DEFINE_SEAL_MODULE();

#include "PhysicsTools/UtilAlgos/interface/ObjectSelector.h"
#include "PhysicsTools/UtilAlgos/interface/EventSetupInitTrait.h"
#include "RecoEgamma/EgammaSelectDigis/plugins/EcalDigiSelector.h" 
#include "RecoEgamma/EgammaSelectDigis/plugins/EcalDigiSelectAnalyzer.h"

DEFINE_ANOTHER_FWK_MODULE(EcalDigiSelector);
DEFINE_ANOTHER_FWK_MODULE(EcalDigiSelectAnalyzer);
