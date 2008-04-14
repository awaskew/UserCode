#include "EgammaAnalysis/PhotonIDAlgos/interface/CutBasedPhotonIDAlgo.h"
#include "AnalysisDataFormats/Egamma/interface/PhotonID.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaReco/interface/BasicCluster.h"

void CutBasedPhotonIDAlgo::setup(const edm::ParameterSet& conf) {
  
  // Get all the parameters
  baseSetup(conf);
  
  //get cuts here
  photonBasicClusterIsolationCut_ = conf.getParameter<double>("PhotonBCIso");
  photonHollowConeTrkIsolationCut_ = conf.getParameter<double>("PhotonHollowTrk");
  photonSolidConeTrkIsolationCut_ = conf.getParameter<double>("PhotonSolidTrk");
  photonCutElectronDupes_ = conf.getParameter<bool>("DiscardAlsoElectrons");
  photonCutNonFiducialClus_ = conf.getParameter<bool>("RequireFiducial");
  
  
}

reco::PhotonID CutBasedPhotonIDAlgo::calculate(const reco::Photon*, const edm::Event&){

  //need to do the following things here:
  //1.)  Call base class methods to calculate photonID variables like fiducial and
  //     isolations.
  //2.)  Decide whether this particular photon passes the cuts that are set forth in the ps.
  //3.)  Create a new PhotonID object, complete with decision and return it.

  reco::PhotonID temp;
  return temp;

}
