#include <vector>
#include "RecoEgamma/EgammaSelectDigis/interface/EcalDigiSelector.h"
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
  reco::SuperClusterCollection EndcapSuperClusters = *pBarrelSuperClusters;
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
 //EBDigiCollection SEBDigiCol;
 //EEDigiCollection SEEDigiCol;
 std::auto_ptr<EEDigiCollection> SEEDigiCol(new EEDigiCollection);
  int TotClus = saveBarrelSuperClusters.size() + saveEndcapSuperClusters.size();
  if (TotClus >= nclus_sel_){
    if (saveBarrelSuperClusters.size() > 0){
      //Get barrel rec hit collection
      edm::Handle<EcalRecHitCollection> ecalhitsCollH;
      evt.getByLabel(EcalEBRecHitTag_, ecalhitsCollH);
      const EcalRecHitCollection* rechitsCollection = ecalhitsCollH.product();

      //get barrel digi collection
      edm::Handle<EBDigiCollection> pdigis;
      const EBDigiCollection* digis=0;
      evt.getByLabel(EcalEBDigiTag_,pdigis);
      digis = pdigis.product(); // get a ptr to the product
      
      //Ready the topology for the cluster tools
      edm::ESHandle<CaloTopology> pTopology;
      es.get<CaloTopologyRecord>().get(pTopology);
      const CaloTopology *topology = pTopology.product();
      
      //pick out the digis for the 3x3 in each of the selected superclusters
      for (int loop = 0;loop < int(saveBarrelSuperClusters.size());loop++){
	SuperCluster clus1 = saveBarrelSuperClusters[loop];
	BasicClusterRef bcref = clus1.seed();
	const BasicCluster *bc = bcref.get();
	//Get the maximum detid
	std::pair<DetId, float> EDetty =  EcalClusterTools::getMaximum(*bc, rechitsCollection);
	//get the 3x3 array centered on maximum detid.
	std::vector<DetId> detvec = EcalClusterTools::matrixDetId(topology, EDetty.first, -1, 1, -1, 1);
	//Loop over the 3x3
	for (int detloop=0;detloop < int(detvec.size());detloop++){
	  EBDetId detL = detvec[detloop];
	  EBDigiCollection::const_iterator myDg = digis->find(detL);
	  EBDataFrame myDigi = (*myDg);
	  SEBDigiCol->push_back(detL);
	  EBDataFrame df( SEBDigiCol->back());
	  for (int iq=0;iq<10;++iq)
	    df.setSample(iq, myDigi.sample(iq).adc());
	    //SEBDigiCol.push_back(detL.rawId(), *myDg);
	}//loop over detids
	
      }//Loop over barrel superclusters
    }//If barrel superclusters need saving.


    if (saveEndcapSuperClusters.size() > 0){
      //Get endcap rec hit collection
      edm::Handle<EcalRecHitCollection> ecalhitsCollH;
      evt.getByLabel(EcalEERecHitTag_, ecalhitsCollH);
      const EcalRecHitCollection* rechitsCollection = ecalhitsCollH.product();

      //get endcap digi collection
      edm::Handle<EEDigiCollection> pdigis;
      const EEDigiCollection* digis=0;
      evt.getByLabel(EcalEEDigiTag_,pdigis);
      digis = pdigis.product(); // get a ptr to the product
      
      //Ready the topology for the cluster tools
      edm::ESHandle<CaloTopology> pTopology;
      es.get<CaloTopologyRecord>().get(pTopology);
      const CaloTopology *topology = pTopology.product();
      
      //pick out the digis for the 3x3 in each of the selected superclusters
      for (int loop = 0;loop < int(saveEndcapSuperClusters.size());loop++){
	SuperCluster clus1 = saveEndcapSuperClusters[loop];
	BasicClusterRef bcref = clus1.seed();
	const BasicCluster *bc = bcref.get();
	//Get the maximum detid
	std::pair<DetId, float> EDetty =  EcalClusterTools::getMaximum(*bc, rechitsCollection);
	//get the 3x3 array centered on maximum detid.
	std::vector<DetId> detvec = EcalClusterTools::matrixDetId(topology, EDetty.first, -1, 1, -1, 1);
	//Loop over the 3x3
	for (int detloop=0;detloop < int(detvec.size());detloop++){
	  EEDetId detL = detvec[detloop];
	  EEDigiCollection::const_iterator myDg = digis->find(detL);
	  EEDataFrame myDigi = (*myDg);
	  SEEDigiCol->push_back(detL);
	  EEDataFrame df( SEEDigiCol->back());
	  for (int iq=0;iq<10;++iq)
	    df.setSample(iq, myDigi.sample(iq).adc());
	    //SEBDigiCol.push_back(detL.rawId(), *myDg);
	}//loop over detids
	
      }//Loop over endcap superclusters
    }//If endcap superclusters need saving.

  }//If we're actually saving stuff

  //Okay, either my collections have been filled with the requisite Digis, or they haven't.

  //  std::auto_ptr<EBDigiCollection> epEBDigiCol(&SEBDigiCol);
  //epEBDigiCol->assign(SEBDigiCol.begin(), SEBDigiCol.end());
  evt.put(SEBDigiCol, selectedEcalEBDigiCollection_);
  
  //std::auto_ptr<EEDigiCollection> epEEDigiCol(&SEEDigiCol);
  //epEEDigiCol->assign(SEEDigiCol.begin(), SEEDigiCol.end());
  evt.put(SEEDigiCol, selectedEcalEEDigiCollection_);
  
}



