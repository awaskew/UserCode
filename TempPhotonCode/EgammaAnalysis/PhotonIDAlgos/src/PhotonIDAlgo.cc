#include "EgammaAnalysis/PhotonIDAlgos/interface/PhotonIDAlgo.h"
#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include <string>
#include <TMath.h>

void PhotonIDAlgo::baseSetup(const edm::ParameterSet& conf) {
  //Need to get track collection, super cluster and basic cluster collection here.
  //Also cone Radii for tracks and bc
  photonBasicClusterConeRadius_ = conf.getParameter<double>("BasicClusterConeRadius");
  trackConeOuterRadius_ = conf.getParameter<double>("TrackConeOuterRadius");
  trackConeInnerRadius_ = conf.getParameter<double>("TrackConeInnerRadius");
  trackInputTag_ = conf.getParameter<edm::InputTag>("trackProducer");
  barrelSuperClusterProducer_ = conf.getParameter<std::string>("barrelSuperClusterProducer");
  endcapSuperClusterProducer_ = conf.getParameter<std::string>("endcapSuperClustersProducer");      
  barrelhybridsuperclusterCollection_ = conf.getParameter<std::string>("barrelhybridsuperclusterCollection");
  endcapsuperclusterCollection_ = conf.getParameter<std::string>("endcapsuperclusterCollection");
  barrelbasicclusterCollection_ = conf.getParameter<std::string>("barrelbasiccluterCollection");
  barrelbasicclusterProducer_ = conf.getParameter<std::string>("barrelbasicclusterProducer");
  endcapbasicclusterCollection_ = conf.getParameter<std::string>("endcapbasicclusterCollection");
  endcapbasicclusterProducer_ = conf.getParameter<std::string>("endcapbasicclusterProducer");
  
}



void PhotonIDAlgo::classify(const reco::Photon* photon, 
			    bool &isEBPho,
			    bool &isEEPho,
			    bool &isEBGap,
			    bool &isEEGap,
			    bool &isEBEEGap){

  
  double eta = photon->p4().Eta();
  double phi = photon->p4().Phi();
  double feta = fabs(eta);
  if(feta>1.479) 
    isEEPho = true;
  else 
    isEBPho = true;
  if (fabs(feta-1.479)<.1) isEBEEGap=true; 
  
  
  //fiducial cuts, currently only for EB, since I don't know
  //EE yet.
  float phigap = fabs(phi-int(phi*9/3.1416)*3.1416/9.);
  if(phigap > 1.65 && phigap <1.85) isEBGap=true;
  if(fabs(eta)<.05) isEBGap=true;
  if(fabs(eta)>.4 && fabs(eta)<.5) isEBGap=true;
  if(fabs(eta)>.75 && fabs(eta)<.85) isEBGap=true;
  if(fabs(eta)>1.1 && fabs(eta)<1.2) isEBGap=true;
  if(fabs(eta)>1.43) isEBGap=true;
  
}

void PhotonIDAlgo::calculateTrackIso(const reco::Photon* photon,
				     const edm::Event& e,
				     double &trkCone,
				     int &ntrkCone,
				     double pTThresh,
				     double RCone,
				     double RinnerCone){
  //Track isolation calculator goes here.
  //Not my code, I stole it from:
  //RecoEgamma/EgammaIsolationAlgos/src/ElectronTkIsolation.
  //and hacked it to my own purposes.  Therefore, consider mistakes mine (AA).
  
  int counter  =0;
  double ptSum =0.;
  
  //Photon Eta and Phi.  Hope these are correct.
  double peta = photon->p4().Eta();
  double pphi = photon->p4().Phi();
  
  //get the tracks
  edm::Handle<reco::TrackCollection> tracks;
  e.getByLabel(trackInputTag_,tracks);
  const reco::TrackCollection* trackCollection = tracks.product();
  
  for ( reco::TrackCollection::const_iterator itrTr  = (*trackCollection).begin(); 
	itrTr != (*trackCollection).end(); 
	++itrTr){
    math::XYZVector tmpTrackMomentumAtVtx = (*itrTr).momentum(); 
    double this_pt  = (*itrTr).pt();
    if ( this_pt < pTThresh ) 
      continue;  
    
    //This is vertex checking, I'll need to substitute PV somehow...
    //	if (fabs( (*itrTr).dz() - (*tmpTrack).dz() ) > lip_ )
    //  continue ;
    double trEta = (*itrTr).eta();
    double trPhi = (*itrTr).phi();
    double deta2 = (trEta-peta)*(trEta-peta);
    double dphi = fabs(trPhi-pphi);
    if (dphi > TMath::Pi()) dphi = TMath::Pi()*2 - dphi;
    double dphi2 = dphi*dphi;
    double dr = sqrt(deta2 + dphi2);
    if ( dr < RCone && dr > RinnerCone ){
      ++counter;
      ptSum += this_pt;
    }//In cone? 
  }//end loop over tracks                 
  ntrkCone = counter;
  trkCone = ptSum;
}

				     

double PhotonIDAlgo::calculateBasicClusterIso(const reco::Photon* photon,
					    const edm::Event& iEvent,
					    double RCone,
					    double RConeInner,
					    double etMin){
  //This is not my code, I stole this almost entirely from 
  //RecoEgamma/EgammaIsolationAlgos/src/EgammaEcalIsolation.cc
  //Any mistakes in adaptation for here are mine (AA).

  //Need to change this to make appropriate to EB-EE
  edm::Handle<reco::SuperClusterCollection> superClusterH;
  edm::Handle<reco::BasicClusterCollection> basicClusterH;

  double peta = photon->p4().Eta();
  if (fabs(peta) > 1.479){
    iEvent.getByLabel(endcapSuperClusterProducer_,endcapsuperclusterCollection_,superClusterH);
    iEvent.getByLabel(endcapbasicclusterProducer_,endcapbasicclusterCollection_,basicClusterH);
  }
  else{
    iEvent.getByLabel(barrelSuperClusterProducer_,barrelhybridsuperclusterCollection_,superClusterH);
    iEvent.getByLabel(barrelbasicclusterProducer_,barrelbasicclusterCollection_,basicClusterH);
  }
  const reco::SuperClusterCollection* superClusterCollection_ = superClusterH.product();
  const reco::BasicClusterCollection* basicClusterCollection_ = basicClusterH.product();


  double ecalIsol=0.;
  reco::SuperClusterRef sc = photon->superCluster();
  math::XYZVector position(sc.get()->position().x(),
			   sc.get()->position().y(),
			   sc.get()->position().z());
  
  // match the photon hybrid supercluster with those with Algo==0 (island)
  double delta1=1000.;
  const reco::SuperCluster *matchedsupercluster=0;
  bool MATCHEDSC = false;
  
  for(reco::SuperClusterCollection::const_iterator scItr = superClusterCollection_->begin(); scItr != superClusterCollection_->end(); ++scItr){
     
    const reco::SuperCluster *supercluster = &(*scItr);
    
    math::XYZVector currentPosition(supercluster->position().x(),
				    supercluster->position().y(),
				    supercluster->position().z());
  
     
    if(supercluster->seed()->algo() == 0){    
      double trEta = currentPosition.eta();
      double trPhi = currentPosition.phi();
      double peta = position.eta();
      double pphi = position.phi();
      double deta2 = (trEta-peta)*(trEta-peta);
      double dphi = fabs(trPhi-pphi);
      if (dphi > TMath::Pi()) dphi = TMath::Pi()*2 - dphi;
      double dphi2 = dphi*dphi;
      double dr = sqrt(deta2 + dphi2);
      if (dr < delta1) {
	delta1=dr;
	matchedsupercluster = supercluster;
	MATCHEDSC = true;
      }
    }
  }
 
  const reco::BasicCluster *cluster= 0;
  
  //loop over basic clusters
  for(reco::BasicClusterCollection::const_iterator cItr = basicClusterCollection_->begin(); cItr != basicClusterCollection_->end(); ++cItr){
    
    cluster = &(*cItr);
    double ebc_bcchi2 = cluster->chi2();
    int   ebc_bcalgo = cluster->algo();
    double ebc_bce    = cluster->energy();
    double ebc_bceta  = cluster->eta();
    double ebc_bcet   = ebc_bce*sin(2*atan(exp(ebc_bceta)));
    double newDelta = 0.;
 
 
    if (ebc_bcet > etMin && ebc_bcalgo == 0) {
      if (ebc_bcchi2 < 30.) {
	
	if(MATCHEDSC){
	  bool inSuperCluster = false;
	  
	  reco::basicCluster_iterator theEclust = matchedsupercluster->clustersBegin();
	  // loop over the basic clusters of the matched supercluster
	  for(;theEclust != matchedsupercluster->clustersEnd();
	      theEclust++) {
	    if ((**theEclust) ==  (*cluster) ) inSuperCluster = true;
	  }
	  if (!inSuperCluster) {
	    
	    math::XYZVector basicClusterPosition(cluster->position().x(),
						 cluster->position().y(),
						 cluster->position().z());
	    double trEta = basicClusterPosition.eta();
	    double trPhi = basicClusterPosition.phi();
	    double peta = position.eta();
	    double pphi = position.phi();
	    double deta2 = (trEta-peta)*(trEta-peta);
	    double dphi = fabs(trPhi-pphi);
	    if (dphi > TMath::Pi()) dphi = TMath::Pi()*2 - dphi;
	    double dphi2 = dphi*dphi;
	    double dr = sqrt(deta2 + dphi2);
	    if(dr < RCone
	       && newDelta > RConeInner) {
	      ecalIsol+=ebc_bcet;
	    }
	  }
	}
      } // matches ebc_bcchi2
    } // matches ebc_bcet && ebc_bcalgo
    
  }
  
  //  std::cout << "Will return ecalIsol = " << ecalIsol << std::endl; 
  return ecalIsol;
  

}
