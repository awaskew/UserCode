#include "ClusterSkim/ClusterFilter/interface/ClusterFilter.h"
#include <vector>
#include <TMath.h>

ClusterFilter::ClusterFilter(const edm::ParameterSet& iConfig){

    //now do what ever initialization is needed
  _SCTag = iConfig.getParameter<edm::InputTag>("SCTag");
  _nHitSel = iConfig.getParameter<int>("nHitSel");
  _nClusSel = iConfig.getParameter<int>("nClusSel");
  _ESel = iConfig.getParameter<double>("ESel");
}
 

ClusterFilter::~ClusterFilter()
{
  
  // do anything here that needs to be done at desctruction time
    // (e.g. close files, deallocate resources etc.)
 
}
 
 
 //
 // member functions
 //
 
bool ClusterFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
 
  
  edm::Handle<reco::SuperClusterCollection> SCHandle;
  iEvent.getByLabel(_SCTag,SCHandle);
  const reco::SuperClusterCollection hclus = *SCHandle;
  
  for (reco::SuperClusterCollection::const_iterator *clus = hclus.begin();
       clus!=hclus.end();++clus){
    if (clus->clustersSize()>=_nClusSel) return true;
    const std::vector< std::pair<DetId, float> > hitsel = clus->hitsAndFractions();
    if (hitsel.size()>=_nHitSel) return true;
    if (clus->energy()>=_ESel) return true;

  }//Loop over reco::SuperCluster collection

  return false;
}

// ------------ method called once each job just before starting event loop  ------------
void ClusterFilter::beginJob(const edm::EventSetup&)
{
}

// ------------ method called once each job just after ending the event loop  ------------
void ClusterFilter::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(ClusterFilter);
