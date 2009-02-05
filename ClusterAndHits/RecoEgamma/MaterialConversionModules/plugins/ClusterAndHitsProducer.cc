#include <vector>
#include "RecoEgamma/MaterialConversionModules/plugins/ClusterAndHitsProducer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "Geometry/CommonDetUnit/interface/TrackingGeometry.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "Geometry/CommonDetUnit/interface/GeomDetUnit.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h" 
#include "TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHitBuilder.h"
#include "TrackingTools/Records/interface/TransientRecHitRecord.h" 
#include "RecoTracker/Record/interface/CkfComponentsRecord.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit2D.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"
#include "DataFormats/SiStripDetId/interface/StripSubdetector.h"
#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"
#include "DataFormats/SiStripDetId/interface/TIBDetId.h"
#include "DataFormats/SiStripDetId/interface/TOBDetId.h"
#include "DataFormats/SiStripDetId/interface/TIDDetId.h"
#include "DataFormats/SiStripDetId/interface/TECDetId.h"
#include "DataFormats/SiPixelDetId/interface/PXBDetId.h"
#include "DataFormats/SiPixelDetId/interface/PXFDetId.h" 
#include "DataFormats/TrackerRecHit2D/interface/SiStripMatchedRecHit2DCollection.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit2DCollection.h" 
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHitCollection.h"
#include "DataFormats/SiStripCluster/interface/SiStripClusterCollection.h"
#include "DataFormats/SiStripCluster/interface/SiStripCluster.h"
#include "DataFormats/SiPixelCluster/interface/SiPixelCluster.h"
#include "TrackingTools/RoadSearchHitAccess/interface/DetHitAccess.h"
#include "RecoTracker/MeasurementDet/interface/MeasurementTracker.h"
#include <TMath.h>
#include <iostream>
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "DataFormats/Common/interface/DetSetNew.h"
#include "DataFormats/Common/interface/DetSetVectorNew.h"

using namespace reco;
ClusterAndHitsProducer::ClusterAndHitsProducer(const edm::ParameterSet& ps)
{
  
  hybridsuperclusterCollection_ = ps.getParameter<std::string>("hybridsuperclusterCollection");
  hybridsuperclusterProducer_ = ps.getParameter<std::string>("hybridsuperclusterProducer");

  cluster_pt_thresh_ = ps.getParameter<double>("cluster_pt_thresh");
 
  matchedStripRecHitsInputTag_ = ps.getParameter<edm::InputTag>("matchedStripRecHits");
  rphiStripRecHitsInputTag_    = ps.getParameter<edm::InputTag>("rphiStripRecHits");
  stereoStripRecHitsInputTag_  = ps.getParameter<edm::InputTag>("stereoStripRecHits");
  pixelRecHitsInputTag_  = ps.getParameter<edm::InputTag>("pixelRecHits");  
  
  clusterMatchedRecHitsColl_ = ps.getParameter<std::string>("clusterMatchedRecHitsColl");
  clusterRPhiRecHitsColl_ = ps.getParameter<std::string>("clusterRPhiRecHitsColl");
  clusterStereoRecHitsColl_ = ps.getParameter<std::string>("clusterStereoRecHitsColl");
  siClusterColl_ = ps.getParameter<std::string>("siClusterColl");
  siPixClusterColl_ = ps.getParameter<std::string>("siPixClusterColl");

  produces<SiStripMatchedRecHit2DCollection>(clusterMatchedRecHitsColl_);
  produces<SiStripRecHit2DCollection>(clusterRPhiRecHitsColl_);
  produces<SiStripRecHit2DCollection>(clusterStereoRecHitsColl_);
  produces<SiPixelRecHitCollection>(clusterPixelRecHitsColl_);
  produces<edmNew::DetSetVector<SiStripCluster> >(siClusterColl_);
  produces<edmNew::DetSetVector<SiPixelCluster> >(siPixClusterColl_);

  nEvt_ = 0;
}


ClusterAndHitsProducer::~ClusterAndHitsProducer()
{
}

void ClusterAndHitsProducer::produce(edm::Event& evt, const edm::EventSetup& es)
{
  //This is where the products from our event will go.
  //clear the decks.
  rphiMappa_.clear();
  rphidetIdSet_.clear();
  rphiMappaClus_.clear();
  pixMappa_.clear();
  pixdetIdSet_.clear();
  pixMappaClus_.clear();
  sterMappa_.clear();
  sterdetIdSet_.clear();
  

  //Get HybridSuperClusters to start.
  edm::Handle<reco::SuperClusterCollection> pBarrelHybridSuperClusters;
  
  evt.getByLabel(hybridsuperclusterProducer_, hybridsuperclusterCollection_, pBarrelHybridSuperClusters);
  if (!pBarrelHybridSuperClusters.isValid()){
    std::cout << "Error! can't get collection with label " << hybridsuperclusterCollection_.c_str()<<std::endl; ;
  }
  reco::SuperClusterCollection BarrelHybridSuperClusters = *pBarrelHybridSuperClusters;
  //Got the HybridSuperClusters

  //Loop over superclusters, and apply threshold
  for (int loop1=0;loop1<int(BarrelHybridSuperClusters.size());loop1++){
    SuperCluster clus1 = BarrelHybridSuperClusters[loop1];
    float eta1 = clus1.eta();
    float phi1 = clus1.phi();
    if (phi1<0) phi1 +=TMath::Pi()*2.;
    float energy1 = clus1.energy();
    float theta1 = 2*atan(exp(-1.*eta1));
    float cluspt1 = energy1 * sin(theta1);
    if (cluspt1 > cluster_pt_thresh_){
      
      edm::Handle<SiStripRecHit2DCollection> prphiRecHits;
      evt.getByLabel(rphiStripRecHitsInputTag_ ,prphiRecHits);
      edm::Handle<SiStripRecHit2DCollection> pstereoRecHits;
      evt.getByLabel(stereoStripRecHitsInputTag_ ,pstereoRecHits);
 
      const SiStripRecHit2DCollection *rphiRecHits = prphiRecHits.product();
      const SiStripRecHit2DCollection *stereoRecHits = pstereoRecHits.product();
      
      const SiPixelRecHitCollection *pixelRecHitCollection = 0;
      edm::Handle<SiPixelRecHitCollection> pixelRecHits;
      evt.getByLabel(pixelRecHitsInputTag_, pixelRecHits);
      pixelRecHitCollection = pixelRecHits.product();

      //Need tracker geometry to tell me where stuff is.
      edm::ESHandle<TrackerGeometry> tracker;
      es.get<TrackerDigiGeometryRecord>().get(tracker);
      const TrackerGeometry& geometry = *tracker;
      
      //Okay, what we'd actually like to do is sift through the different
      //DetId and select those which are geometrically close to the
      //passage of electron/photon.

      for (SiStripRecHit2DCollection::const_iterator detvec = rphiRecHits->begin();
	   detvec!=rphiRecHits->end();detvec++){
	
	unsigned int detid = detvec->detId();
	DetId detIdObject( detid ); 
	const GeomDetUnit *det = tracker->idToDetUnit(detIdObject);
	
	// calculate global position of center of detector
	GlobalPoint center = det->surface().toGlobal(LocalPoint(0,0,0)); 
	double phi = center.phi();
	if (phi<0) phi+=TMath::Pi()*2.;
	double zed = center.z();
	double dphi = 0;
	double dphi1 = fabs(phi-phi1);
	if (dphi1 > TMath::Pi()) dphi = 2.*TMath::Pi() - dphi1;
	else dphi = dphi1;
	
	bool testDetZ = TestZPosition(clus1.x(), clus1.y(), clus1.z(),
				      0.,
				      center.x(), center.y(), zed,
				      0., 20.);

	if (dphi < 0.3 && testDetZ){

	  for (edmNew::DetSet<SiStripRecHit2D>::const_iterator rphits = detvec->begin();
	       rphits!=detvec->end();++rphits){
	    
	    //Now only select likely hits. 
	    GlobalPoint position = geometry.idToDet( 
						    rphits->geographicalId()
						    )->surface().toGlobal(
									  rphits->localPosition());
	    float position_err_zz = 0;
	    if ( (unsigned int)detIdObject.subdetId() == StripSubdetector::TEC)
	      position_err_zz = 1.;
	    if ( (unsigned int)detIdObject.subdetId() == StripSubdetector::TOB){
	      TOBDetId tobid(detIdObject.rawId()); 
	      if ( !tobid.glued() ) {
		position_err_zz=12.;
	      }
	    }
	    if ( (unsigned int)detIdObject.subdetId() == StripSubdetector::TIB){
	      TIBDetId tobid(detIdObject.rawId()); 
	      if ( !tobid.glued() ) {
		position_err_zz=12.;
	      }
	    }
	    
	    bool TestZ = TestZPosition(clus1.x(), clus1.y(), clus1.z(),
				       0.,
				       position.x(), position.y(), position.z(),
				       0., position_err_zz);
	    double HitR = position.perp();
	    double PhiDiffMax = 0.01 + 0.0025 * HitR;
	    double HitPhi = position.phi();
	    if (HitPhi<0) HitPhi+=TMath::Pi()*2.;
	    double HitPhiDiff = fabs(phi1 - HitPhi);
	    if (HitPhiDiff > TMath::Pi()) HitPhiDiff = TMath::Pi()*2. - HitPhiDiff;
	    if (TestZ && HitPhiDiff < PhiDiffMax){
	      std::map<DetId, std::vector<SiStripRecHit2D> >::iterator it = rphiMappa_.find(detIdObject);
	      std::map<DetId, std::vector<SiStripCluster> >::iterator sit = rphiMappaClus_.find(detIdObject);
	      if (it == rphiMappa_.end() && sit == rphiMappaClus_.end()){
		std::vector <SiStripRecHit2D> vecca;
		std::vector <SiStripCluster> Cvecca;
		vecca.push_back(*rphits);
		SiStripRecHit2D::ClusterRef reffer = rphits->cluster();
		const SiStripCluster *clussy = reffer.get();
		Cvecca.push_back(*clussy);
		rphiMappa_.insert(make_pair(detIdObject, vecca));
		rphiMappaClus_.insert(make_pair(detIdObject, Cvecca));
		rphidetIdSet_.push_back(detIdObject);
	      }
	      else if (it!=rphiMappa_.end() && sit!=rphiMappaClus_.end()){
		//This is where I would look to see if I already added the hit.
		it->second.push_back(*rphits);
		SiStripRecHit2D::ClusterRef reffer = rphits->cluster();
		const SiStripCluster *clussy = reffer.get();
		sit->second.push_back(*clussy);
	      }
	    }   
	  }
	}
      }   

      for (SiPixelRecHitCollection::const_iterator detvec = pixelRecHitCollection->begin();
	   detvec!=pixelRecHitCollection->end();detvec++){
	unsigned int detid = detvec->detId();
	DetId detIdObject( detid );
	const GeomDetUnit *det = tracker->idToDetUnit(detIdObject);
	
	// calculate global position of center
	GlobalPoint center = det->surface().toGlobal(LocalPoint(0,0,0)); 
	double phi = center.phi();
	if (phi<0) phi+=TMath::Pi()*2.;
	double zed = center.z();
	double dphi = 0;
	double dphi1 = fabs(phi-phi1);
	if (dphi1 > TMath::Pi()) dphi = 2.*TMath::Pi() - dphi1;
	else dphi = dphi1;
	
	bool testDetZ = TestZPosition(clus1.x(), clus1.y(), clus1.z(),
				      0.,
				      center.x(), center.y(), zed,
				      0., 20.);
	
	if (dphi < .3 && testDetZ){
	  for (edmNew::DetSet<SiPixelRecHit>::const_iterator pixhits = detvec->begin();
	       pixhits!=detvec->end();++pixhits){
	    
	    //Now only select likely hits. 
	    GlobalPoint position = geometry.idToDet( 
						    pixhits->geographicalId()
						    )->surface().toGlobal(
									  pixhits->localPosition());

	    float position_err_zz = 0;
	    
	    
	    bool TestZ = TestZPosition(clus1.x(), clus1.y(), clus1.z(),
				       0.,
				       position.x(), position.y(), position.z(),
				       0., position_err_zz);
	    double HitR = position.perp();
	    double PhiDiffMax = 0.01 + 0.0025 * HitR;
	    double HitPhi = position.phi();
	    if (HitPhi<0) HitPhi+=TMath::Pi()*2.;
	    double HitPhiDiff = fabs(phi1 - HitPhi);
	    if (HitPhiDiff > TMath::Pi()) HitPhiDiff = TMath::Pi()*2. - HitPhiDiff;
	    if (TestZ && HitPhiDiff < PhiDiffMax){
	      std::map<DetId, std::vector<SiPixelRecHit> >::iterator it = pixMappa_.find(detIdObject);
	      std::map<DetId, std::vector<SiPixelCluster> >::iterator sit = pixMappaClus_.find(detIdObject);
	      if (it == pixMappa_.end() && sit == pixMappaClus_.end()){
		std::vector <SiPixelRecHit> vecca;
		std::vector <SiPixelCluster> Cvecca;
		vecca.push_back(*pixhits);
		SiPixelRecHit::ClusterRef reffer = pixhits->cluster();
		const SiPixelCluster *clussy = reffer.get();

		Cvecca.push_back(*clussy);
		pixMappa_.insert(make_pair(detIdObject, vecca));
		pixMappaClus_.insert(make_pair(detIdObject, Cvecca));
		pixdetIdSet_.push_back(detIdObject);
	      }
	      else if (it!=pixMappa_.end() && sit!=pixMappaClus_.end()){
		it->second.push_back(*pixhits);
		SiPixelRecHit::ClusterRef reffer = pixhits->cluster();
		const SiPixelCluster *clussy = reffer.get();
		sit->second.push_back(*clussy);
		
	      }
	    }
	  }
	}
      }  
      for (SiStripRecHit2DCollection::const_iterator detvec = stereoRecHits->begin();
	   detvec!=stereoRecHits->end();detvec++){
	unsigned int detid = detvec->detId();
	DetId detIdObject( detid ); 
	const GeomDetUnit *det = tracker->idToDetUnit(detIdObject);
	
	// calculate global position of center
	GlobalPoint center = det->surface().toGlobal(LocalPoint(0,0,0)); 
	double phi = center.phi();
	if (phi<0) phi+=TMath::Pi()*2.;
	double zed = center.z();
	double dphi = 0;
	double dphi1 = fabs(phi-phi1);
	if (dphi1 > TMath::Pi()) dphi = 2.*TMath::Pi() - dphi1;
	else dphi = dphi1;
	
	bool testDetZ = TestZPosition(clus1.x(), clus1.y(), clus1.z(),
				      0.,
				      center.x(), center.y(), zed,
				      0., 20.);
	if (dphi < 0.3 && testDetZ){
	  for (edmNew::DetSet<SiStripRecHit2D>::const_iterator sthits = detvec->begin();
	       sthits!=detvec->end();++sthits){
	    
	    
	    //Now only select likely hits. 
	    GlobalPoint position = geometry.idToDet( 
						    sthits->geographicalId()
						    )->surface().toGlobal(
									  sthits->localPosition());
	    float position_err_zz = 0;
	    
	    
	    bool TestZ = TestZPosition(clus1.x(), clus1.y(), clus1.z(),
				       0.,
				       position.x(), position.y(), position.z(),
				       0., position_err_zz);
	    double HitR = position.perp();
	    double PhiDiffMax = 0.01 + 0.0025 * HitR;
	    double HitPhi = position.phi();
	    if (HitPhi<0) HitPhi+=TMath::Pi()*2.;
	    double HitPhiDiff = fabs(phi1 - HitPhi);
	    if (HitPhiDiff > TMath::Pi()) HitPhiDiff = TMath::Pi()*2. - HitPhiDiff;
	    if (TestZ && HitPhiDiff < PhiDiffMax){

	      std::map<DetId, std::vector<SiStripRecHit2D> >::iterator it = sterMappa_.find(detIdObject);
	      std::map<DetId, std::vector<SiStripCluster> >::iterator sit = rphiMappaClus_.find(detIdObject);
	      if (it == sterMappa_.end() && sit == rphiMappaClus_.end()){
		std::vector <SiStripRecHit2D> vecca;
		std::vector <SiStripCluster> Cvecca;
		vecca.push_back(*sthits);
		SiStripRecHit2D::ClusterRef reffer = sthits->cluster();
		const SiStripCluster *clussy = reffer.get();
		Cvecca.push_back(*clussy);
		sterMappa_.insert(make_pair(detIdObject, vecca));
		rphiMappaClus_.insert(make_pair(detIdObject, Cvecca));
		sterdetIdSet_.push_back(detIdObject);
	      }
	      else if (it!=sterMappa_.end() && sit!=rphiMappaClus_.end()){
		it->second.push_back(*sthits);
		SiStripRecHit2D::ClusterRef reffer = sthits->cluster();
		const SiStripCluster *clussy = reffer.get();
		sit->second.push_back(*clussy);
	      }
	    }//test2	   
	  }//hit loop
	}//test 1
      }//stereo hits
    }//PT threshold
  }//SuperClusters loop
  

  std::auto_ptr<SiStripRecHit2DCollection > rpCol(new SiStripRecHit2DCollection);
  //Done with supercluster loop, go through maps and place into the RangeMap objects for storage.
  if (rphidetIdSet_.size() > 0){
    for (int ik=0;ik<int(rphidetIdSet_.size());++ik){
      std::map<DetId, std::vector<SiStripRecHit2D> >::iterator mpa = rphiMappa_.find(rphidetIdSet_[ik]);
      if (mpa != rphiMappa_.end()){
	edmNew::DetSetVector<SiStripRecHit2D>::FastFiller ffiller(*rpCol, mpa->first);
	std::vector<SiStripRecHit2D> vecerr = mpa->second;
	for (int qu=0;qu<int(vecerr.size());++qu){
	  ffiller.push_back(vecerr[qu]);
	}
      }
    }
  }

  std::auto_ptr<SiPixelRecHitCollection > pxCol(new SiPixelRecHitCollection);
                                                      
  if (pixdetIdSet_.size() > 0){
    for (int ik=0;ik<int(pixdetIdSet_.size());++ik){
      std::map<DetId, std::vector<SiPixelRecHit> >::iterator mpa = pixMappa_.find(pixdetIdSet_[ik]);
      if (mpa != pixMappa_.end()){
	edmNew::DetSetVector<SiPixelRecHit>::FastFiller ffiller(*pxCol, mpa->first);
	std::vector<SiPixelRecHit> vecerr = mpa->second;
	for (int qu=0;qu<int(vecerr.size());++qu){
	  ffiller.push_back(vecerr[qu]);
	}
      }
    }
  }

  std::auto_ptr<SiStripRecHit2DCollection > stCol(new SiStripRecHit2DCollection);
  if (sterdetIdSet_.size() > 0){
    for (int ik=0;ik<int(sterdetIdSet_.size());++ik){
      std::map<DetId, std::vector<SiStripRecHit2D> >::iterator mpa = sterMappa_.find(sterdetIdSet_[ik]);
      if (mpa != sterMappa_.end()){
	edmNew::DetSetVector<SiStripRecHit2D>::FastFiller ffiller(*stCol, mpa->first);
	std::vector<SiStripRecHit2D> vecerr = mpa->second;
	for (int qu=0;qu<int(vecerr.size());++qu){
	  ffiller.push_back(vecerr[qu]);
	}
      }
    }
  }                                                   

  //  edmNew::DetSetVector<SiStripCluster> 
  std::auto_ptr<edmNew::DetSetVector<SiStripCluster> > theDetSetVector(new edmNew::DetSetVector<SiStripCluster>);
  //  if (rphidetIdSet_.size() > 0){
  //  for (int ik=0;ik<int(rphidetIdSet_.size());++ik){
  for (std::map<DetId, std::vector<SiStripCluster> >::const_iterator mpa = rphiMappaClus_.begin();
       mpa !=rphiMappaClus_.end();
       ++mpa){
    //std::map<DetId, std::vector<SiStripCluster> >::iterator mpa = rphiMappaClus_.find(rphidetIdSet_[ik]);
    //if (mpa != rphiMappaClus_.end()){

    DetId blerg = mpa->first;
    std::vector<SiStripCluster> vecerr = mpa->second;
    edmNew::DetSetVector<SiStripCluster>::FastFiller ffiller(*theDetSetVector, blerg);

    for (int qu=0;qu<int(vecerr.size());++qu){
      ffiller.push_back(vecerr[qu]);
    }
     
  }
  std::auto_ptr<edmNew::DetSetVector<SiPixelCluster> > thePDetSetVector(new edmNew::DetSetVector<SiPixelCluster>);

  if (pixdetIdSet_.size() > 0){
    for (int ik=0;ik<int(pixdetIdSet_.size());++ik){
      std::map<DetId, std::vector<SiPixelCluster> >::iterator mpa = pixMappaClus_.find(pixdetIdSet_[ik]);
      if (mpa != pixMappaClus_.end()){
	edmNew::DetSetVector<SiPixelCluster>::FastFiller ffiller(*thePDetSetVector, mpa->first);
	std::vector<SiPixelCluster> vecerr = mpa->second;
	for (int qu=0;qu<int(vecerr.size());++qu){
	  ffiller.push_back(vecerr[qu]);
	}
      }
    }
  }


  evt.put(theDetSetVector, siClusterColl_);
  evt.put(thePDetSetVector, siPixClusterColl_);
  evt.put(rpCol,clusterRPhiRecHitsColl_);
  evt.put(pxCol,clusterPixelRecHitsColl_);
  evt.put(stCol,clusterStereoRecHitsColl_);
  
}

bool ClusterAndHitsProducer::TestZPosition(double _SeedX, double _SeedY, double _SeedZ,
				   double _VtxZ,
				   double _TestX, double _TestY, double _TestZ,
				   double _TestXXerr, double _TestYYerr
				   )
{
  //Testing the Z coordinate position of hit.
  double _ZVTXCONSTERR = 10.;
  double _CALZERR=2.;
  //Define the vertex error on Z
  double ZVtxMin = _VtxZ - _ZVTXCONSTERR;
  double ZVtxMax = _VtxZ + _ZVTXCONSTERR;
  
  double SeedZMin = _SeedZ - _CALZERR;
  double SeedZMax = _SeedZ + _CALZERR;
  
  double RCAL = sqrt(_SeedX*_SeedX + _SeedY*_SeedY);
  double slopemin = (SeedZMin - ZVtxMin) / RCAL;
  double slopemax = (SeedZMax - ZVtxMax) / RCAL;
  
  double interceptmin = ZVtxMin;
  double interceptmax = ZVtxMax;
  
  double TESTR = sqrt(_TestX*_TestX + _TestY*_TestY);
  
  double REGIONZMax = slopemax * TESTR + interceptmax + 2.*_TestYYerr;
  double REGIONZMin = slopemin * TESTR + interceptmin - 2.*_TestYYerr;
  
  if (  REGIONZMax > _TestZ && REGIONZMin < _TestZ){
    return true; 
  }
  return false;
}

