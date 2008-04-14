#ifndef CutBasedPhotonID_H
#define CutBasedPhotonID_H

#include "EgammaAnalysis/PhotonIDAlgos/interface/PhotonIDAlgo.h"
#include "AnalysisDataFormats/Egamma/interface/PhotonID.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

class CutBasedPhotonIDAlgo : PhotonIDAlgo {

public:

  CutBasedPhotonIDAlgo(){};

  virtual ~CutBasedPhotonIDAlgo(){};

  void setup(const edm::ParameterSet& conf);
  reco::PhotonID calculate(const reco::Photon*, const edm::Event&);
  void decide(const reco::Photon* photon,
	      const edm::Event& e, reco::PhotonID phID);
 private:
  double photonBasicClusterIsolationCut_;
  double photonHollowConeTrkIsolationCut_;
  double photonSolidConeTrkIsolationCut_;
  bool photonCutElectronDupes_;
  bool photonCutNonFiducialClus_;

  //Isolation parameters
  double photonBasicClusterConeOuterRadius_;
  double photonBasicClusterConeInnerRadius_;
  double isolationbasicclusterThreshold_;
  double trackConeOuterRadius_;
  double trackConeInnerRadius_;
  double isolationtrackThreshold_;
};

#endif // PTDRElectronID_H
