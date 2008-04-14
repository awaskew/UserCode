#ifndef PhotonIDAlgo_H
#define PhotonIDAlgo_H

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "AnalysisDataFormats/Egamma/interface/PhotonID.h"
#include <string>
/*
Base class for doing necessary calculations for Photon ID.  Contains
isolation calculation hacked out of the ElectronTkIsolation class.
Contains basic fiducial cuts for gap and crack cuts
*/

class PhotonIDAlgo {

public:

  PhotonIDAlgo(){};

  virtual ~PhotonIDAlgo(){};

  void baseSetup(const edm::ParameterSet& conf);
  void classify(const reco::Photon* photon, 
		bool &isEBPho,
		bool &isEEPho,
		bool &isEBGap,
		bool &isEEGap,
		bool &isEBEEGap);
  void calculateTrackIso(const reco::Photon* photon,
			 const edm::Event &e,
			 double &trkCone,
			 int &ntrkCone,
			 double pTThresh=0,
			 double RCone=.4,
			 double RinnerCone=.1);
  double calculateBasicClusterIso(const reco::Photon* photon,
				  const edm::Event& iEvent,
				  double RCone=0.4,
				  double RConeInner=0,
				  double etMin=0);

  
 private:

  double photonBasicClusterConeRadius_;
  double trackConeOuterRadius_;
  double trackConeInnerRadius_;
  std::string barrelSuperClusterProducer_;
  std::string endcapSuperClusterProducer_;      
  std::string barrelhybridsuperclusterCollection_;
  std::string endcapsuperclusterCollection_;
  std::string barrelbasicclusterCollection_;
  std::string barrelbasicclusterProducer_;
  std::string endcapbasicclusterCollection_;
  std::string endcapbasicclusterProducer_;
  edm::InputTag trackInputTag_;

  };

#endif // PhotonIDAlgo_H
