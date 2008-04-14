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
  photonBasicClusterConeOuterRadius_ = conf.getParameter<double>("BasicClusterConeOuterRadius");
  photonBasicClusterConeInnerRadius_ = conf.getParameter<double>("BasicClusterConeInnerRadius");
  isolationbasicclusterThreshold_ = conf.getParameter<double>("isolationbasicclusterThreshold");

  trackConeOuterRadius_ = conf.getParameter<double>("TrackConeOuterRadius");
  trackConeInnerRadius_ = conf.getParameter<double>("TrackConeInnerRadius");
  isolationtrackThreshold_ = conf.getParameter<double>("isolationtrackThreshold");
  
}

reco::PhotonID CutBasedPhotonIDAlgo::calculate(const reco::Photon* pho, const edm::Event& e){

  //need to do the following things here:
  //1.)  Call base class methods to calculate photonID variables like fiducial and
  //     isolations.
  //2.)  Decide whether this particular photon passes the cuts that are set forth in the ps.
  //3.)  Create a new PhotonID object, complete with decision and return it.


  //Get fiducial information
  bool isEBPho   = false;
  bool isEEPho   = false;
  bool isEBGap   = false;
  bool isEEGap   = false;
  bool isEBEEGap = false;
  classify(pho, isEBPho, isEEPho, isEBGap, isEEGap, isEBEEGap);

  //Calculate hollow cone track isolation
  int ntrk=0;
  double trkiso=0;
  calculateTrackIso(pho, e, trkiso, ntrk, isolationtrackThreshold_,    
		    trackConeOuterRadius_, trackConeInnerRadius_);

  //Calculate solid cone track isolation
  int sntrk=0;
  double strkiso=0;
  calculateTrackIso(pho, e, strkiso, sntrk, isolationtrackThreshold_,    
		    trackConeOuterRadius_, 0.);
  
  double bc_iso = calculateBasicClusterIso(pho, e, 
					   photonBasicClusterConeOuterRadius_,
					   photonBasicClusterConeInnerRadius_,
					   isolationbasicclusterThreshold_);

  bool isElec = isAlsoElectron(pho, e);
  
  reco::PhotonID temp(false, bc_iso, strkiso,
		      trkiso, sntrk, ntrk,
		      isEBPho, isEEPho, isEBGap, isEEGap, isEBEEGap,
		      isElec);

  return temp;

}
