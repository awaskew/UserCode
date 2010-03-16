#include <memory>
 
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include <string>
 
//
// class declaration
//
 
using namespace edm;
 
class ClusterFilter : public edm::EDFilter {
 public:
  explicit ClusterFilter(const edm::ParameterSet&);
  ~ClusterFilter();
 
 private:
  virtual void beginJob(const edm::EventSetup&) ;
  virtual bool filter(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
 
  edm::InputTag _SCTag;
  Int_t _nHitSel;
  Int_t _nClusSel;
  Float_t _ESel;

  
  
};
