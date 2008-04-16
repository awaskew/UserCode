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
  void decide(reco::PhotonID &phID);
 private:
  
  //Which cuts to do?
  bool dophotonBCIsolationCut_;
  bool dophotonHCTrkIsolationCut_;
  bool dophotonSCTrkIsolationCut_;
  bool dophotonHCNTrkCut_;
  bool dophotonSCNTrkCut_;
  bool dorequireNotElectron_;
  bool dorequireFiducial_;

  //Actual cut values
  double photonBasicClusterIsolationCut_;
  double photonHollowConeTrkIsolationCut_;
  double photonSolidConeTrkIsolationCut_;
  int photonSolidConeNTrkCut_;
  int photonHollowConeNTrkCut_;

  //Isolation parameters
  double photonBasicClusterConeOuterRadius_;
  double photonBasicClusterConeInnerRadius_;
  double isolationbasicclusterThreshold_;
  double trackConeOuterRadius_;
  double trackConeInnerRadius_;
  double isolationtrackThreshold_;
};

#endif // PTDRElectronID_H
