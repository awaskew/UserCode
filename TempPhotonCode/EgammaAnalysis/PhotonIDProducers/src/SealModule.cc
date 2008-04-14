#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "PhysicsTools/UtilAlgos/interface/ObjectSelector.h"
#include "PhysicsTools/UtilAlgos/interface/EventSetupInitTrait.h"

#include "EgammaAnalysis/PhotonIDProducers/interface/PhotonIDProducer.h"
//#include "EgammaAnalysis/PhotonIDProducers/interface/PhotonIDSelector.h"
//#include "RecoEgamma/PhotonIDProducers/plugins/PhotonIDSelectorCutBased.h"

//typedef PhotonIDSelector<PhotonIDSelectorCutBased>   PhoIdCutBasedSel;

//typedef ObjectSelector<EleIdCutBasedSel> EleIdCutBased ;
//typedef ObjectSelector<
//          EleIdCutBasedSel, 
//          edm::RefVector<reco::PixelMatchGsfElectronCollection> 
//         > EleIdCutBasedRef ;
//typedef ObjectSelector<
//          EleIdNeuralNetSel, 
//          edm::RefVector<reco::PixelMatchGsfElectronCollection> 
//         > EleIdNeuralNetRef ;
//typedef ObjectSelector<
//          EleIdLikelihoodSel, 
//          edm::RefVector<reco::PixelMatchGsfElectronCollection> 
//         > EleIdLikelihoodRef ;


DEFINE_SEAL_MODULE();

DEFINE_ANOTHER_FWK_MODULE(PhotonIDProducer);

//DEFINE_ANOTHER_FWK_MODULE(EleIdCutBased);
//DEFINE_ANOTHER_FWK_MODULE(EleIdCutBasedRef);
//DEFINE_ANOTHER_FWK_MODULE(EleIdNeuralNetRef);
//DEFINE_ANOTHER_FWK_MODULE(EleIdLikelihoodRef);
