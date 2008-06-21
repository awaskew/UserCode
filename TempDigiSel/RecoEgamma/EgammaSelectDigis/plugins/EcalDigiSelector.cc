#include <vector>
#include "RecoEgamma/EgammaSelectDigis/plugins/EcalDigiSelector.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
#include "DataFormats/EcalDigi/interface/EcalDataFrame.h"
#include "DataFormats/EcalDigi/interface/EEDataFrame.h"
#include "DataFormats/EcalDigi/interface/EBDataFrame.h"
#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHit.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include <TMath.h>
#include <iostream>


using namespace reco;
EcalDigiSelector::EcalDigiSelector(const edm::ParameterSet& ps)
{
  selectedEcalEBDigiCollection_ = ps.getParameter<std::string>("selectedEcalEBDigiCollection");
  selectedEcalEEDigiCollection_ = ps.getParameter<std::string>("selectedEcalEEDigiCollection");

  barrelSuperClusterCollection_ = ps.getParameter<std::string>("barrelSuperClusterCollection");
  barrelSuperClusterProducer_ = ps.getParameter<std::string>("barrelSuperClusterProducer");

  endcapSuperClusterCollection_ = ps.getParameter<std::string>("endcapSuperClusterCollection");
  endcapSuperClusterProducer_ = ps.getParameter<std::string>("endcapSuperClusterProducer");

  EcalEBDigiTag_ = ps.getParameter<edm::InputTag>("EcalEBDigiTag");
  EcalEEDigiTag_ = ps.getParameter<edm::InputTag>("EcalEEDigiTag");

  EcalEBRecHitTag_ = ps.getParameter<edm::InputTag>("EcalEBRecHitTag");
  EcalEERecHitTag_ = ps.getParameter<edm::InputTag>("EcalEERecHitTag");
  
  cluster_pt_thresh_ = ps.getParameter<double>("cluster_pt_thresh");
  nclus_sel_ = ps.getParameter<int>("nclus_sel");
  produces<EBDigiCollection>(selectedEcalEBDigiCollection_);
  produces<EEDigiCollection>(selectedEcalEEDigiCollection_);

}


EcalDigiSelector::~EcalDigiSelector()
{
}

void EcalDigiSelector::produce(edm::Event& evt, const edm::EventSetup& es)
{
  
  //Get BarrelSuperClusters to start.
  edm::Handle<reco::SuperClusterCollection> pBarrelSuperClusters;
  
  evt.getByLabel(barrelSuperClusterProducer_, barrelSuperClusterCollection_, pBarrelSuperClusters);
  if (!pBarrelSuperClusters.isValid()){
    std::cout << "Error! can't get collection with label " << barrelSuperClusterCollection_.c_str()<<std::endl; ;
  }
  reco::SuperClusterCollection BarrelSuperClusters = *pBarrelSuperClusters;
  //Got BarrelSuperClusters

  //Get BarrelSuperClusters to start.
  edm::Handle<reco::SuperClusterCollection> pEndcapSuperClusters;
  
  evt.getByLabel(endcapSuperClusterProducer_, endcapSuperClusterCollection_, pEndcapSuperClusters);
  if (!pEndcapSuperClusters.isValid()){
    std::cout << "Error! can't get collection with label " << endcapSuperClusterCollection_.c_str()<<std::endl; ;
  }
  reco::SuperClusterCollection EndcapSuperClusters = *pEndcapSuperClusters;
  //Got EndcapSuperClusters

  reco::SuperClusterCollection saveBarrelSuperClusters;
  reco::SuperClusterCollection saveEndcapSuperClusters;

  //Loop over barrel superclusters, and apply threshold
  for (int loop=0;loop<int(BarrelSuperClusters.size());loop++){
    SuperCluster clus1 = BarrelSuperClusters[loop];
    float eta1 = clus1.eta();
    float energy1 = clus1.energy();
    float theta1 = 2*atan(exp(-1.*eta1));
    float cluspt1 = energy1 * sin(theta1);
    if (cluspt1 > cluster_pt_thresh_){
      saveBarrelSuperClusters.push_back(clus1);
    }
  }

  //Loop over endcap superclusters, and apply threshold
  for (int loop=0;loop<int(EndcapSuperClusters.size());loop++){
    SuperCluster clus1 = EndcapSuperClusters[loop];
    float eta1 = clus1.eta();
    float energy1 = clus1.energy();
    float theta1 = 2*atan(exp(-1.*eta1));
    float cluspt1 = energy1 * sin(theta1);
    if (cluspt1 > cluster_pt_thresh_){
      saveEndcapSuperClusters.push_back(clus1);
    }
  }
  
  std::auto_ptr<EBDigiCollection> SEBDigiCol(new EBDigiCollection);
  std::auto_ptr<EEDigiCollection> SEEDigiCol(new EEDigiCollection);
  int TotClus = saveBarrelSuperClusters.size() + saveEndcapSuperClusters.size();
  std::cout << "Barrel Clusters: " << saveBarrelSuperClusters.size();
  std::cout << " Endcap Clusters: " << saveEndcapSuperClusters.size();
  std::cout << " Total: " << TotClus << std::endl;
  if (TotClus >= nclus_sel_){
    EcalClusterLazyTools tooly(evt, es, EcalEBRecHitTag_, EcalEERecHitTag_);
    if (saveBarrelSuperClusters.size() > 0){
      
      //Get barrel rec hit collection
//       edm::Handle<EcalRecHitCollection> ecalhitsCollH;
//       evt.getByLabel(EcalEBRecHitTag_, ecalhitsCollH);
//       const EcalRecHitCollection* rechitsCollection = ecalhitsCollH.product();
      
      //get barrel digi collection
      edm::Handle<EBDigiCollection> pdigis;
      const EBDigiCollection* digis=0;
      evt.getByLabel(EcalEBDigiTag_,pdigis);
      digis = pdigis.product(); // get a ptr to the product
      std::vector<DetId> saveTheseDetIds;
      //pick out the detids for the 3x3 in each of the selected superclusters
      for (int loop = 0;loop < int(saveBarrelSuperClusters.size());loop++){
	SuperCluster clus1 = saveBarrelSuperClusters[loop];
	BasicClusterRef bcref = clus1.seed();
	const BasicCluster *bc = bcref.get();
	//Get the maximum detid
	std::pair<DetId, float> EDetty = tooly.getMaximum(*bc);
	//get the 3x3 array centered on maximum detid.
	std::vector<DetId> detvec = tooly.matrixDetId(EDetty.first, -1, 1, -1, 1);
	//Loop over the 3x3
	for (int ik = 0;ik<int(detvec.size());++ik)
	  saveTheseDetIds.push_back(detvec[ik]);
      }
      for (int detloop=0; detloop < int(saveTheseDetIds.size());++detloop){
	EBDetId detL = EBDetId(saveTheseDetIds[detloop]);
	for (EBDigiCollection::const_iterator blah = digis->begin();
	     blah!=digis->end();blah++){
	  if (detL == blah->id()){
	    EBDataFrame myDigi = (*blah);
	    SEBDigiCol->push_back(detL);
	    EBDataFrame df( SEBDigiCol->back());
	    for (int iq =0;iq<myDigi.size();++iq){
	      df.setSample(iq, myDigi.sample(iq).raw());
	      
	    }
	    break;
	  }
	}
      }//loop over detids
      std::cout << "size of new digi container contents: " << SEBDigiCol->size() << std::endl;
    }//If barrel superclusters need saving.
    
    
    if (saveEndcapSuperClusters.size() > 0){
      //Get endcap rec hit collection
      //get endcap digi collection
      edm::Handle<EEDigiCollection> pdigis;
      const EEDigiCollection* digis=0;
      evt.getByLabel(EcalEEDigiTag_,pdigis);
      digis = pdigis.product(); // get a ptr to the product
      std::vector<DetId> saveTheseDetIds;
      //pick out the digis for the 3x3 in each of the selected superclusters
      for (int loop = 0;loop < int(saveEndcapSuperClusters.size());loop++){
	SuperCluster clus1 = saveEndcapSuperClusters[loop];
	BasicClusterRef bcref = clus1.seed();
	const BasicCluster *bc = bcref.get();
	//Get the maximum detid
	std::pair<DetId, float> EDetty = tooly.getMaximum(*bc);
	//get the 3x3 array centered on maximum detid.
	std::vector<DetId> detvec = tooly.matrixDetId(EDetty.first, -1, 1, -1, 1);
	//Loop over the 3x3
	for (int ik = 0;ik<int(detvec.size());++ik)
	  saveTheseDetIds.push_back(detvec[ik]);
      }
      for (int detloop=0; detloop < int(saveTheseDetIds.size());++detloop){
	EEDetId detL = EEDetId(saveTheseDetIds[detloop]);
	for (EEDigiCollection::const_iterator blah = digis->begin();
	     blah!=digis->end();blah++){
	  if (detL == blah->id()){
	    EEDataFrame myDigi = (*blah);
	    SEEDigiCol->push_back(detL);
	    EEDataFrame df( SEEDigiCol->back());
	    for (int iq =0;iq<myDigi.size();++iq){
	      df.setSample(iq, myDigi.sample(iq).raw());	      
	    }
	    break;
	  }
	}
      }//loop over detids
      std::cout << "Current new digi container contents: " << SEEDigiCol->size() << std::endl;
    }//If endcap superclusters need saving.
    
  }//If we're actually saving stuff
  
  //Okay, either my collections have been filled with the requisite Digis, or they haven't.
  
  std::cout << "Saving: " << SEBDigiCol->size() << " barrel digis and " << SEEDigiCol->size() << " endcap digis." << std::endl;

  //Empty collection, or full, still put in event.
  evt.put(SEBDigiCol, selectedEcalEBDigiCollection_);
  evt.put(SEEDigiCol, selectedEcalEEDigiCollection_);
  
}



