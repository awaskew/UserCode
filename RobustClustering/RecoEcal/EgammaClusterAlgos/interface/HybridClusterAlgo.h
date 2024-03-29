#ifndef RecoEcal_EgammaClusterAlgos_HybridClusterAlgo_h
#define RecoEcal_EgammaClusterAlgos_HybridClusterAlgo_h

#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/EgammaReco/interface/BasicClusterFwd.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "RecoCaloTools/Navigation/interface/EcalBarrelNavigator.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHit.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloTopology/interface/EcalBarrelHardcodedTopology.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalEtaPhiRegion.h"
#include "RecoEcal/EgammaCoreTools/interface/PositionCalc.h"
#include "RecoEcal/EgammaCoreTools/interface/BremRecoveryPhiRoadAlgo.h"
#include "RecoEcal/EgammaCoreTools/interface/SuperClusterShapeAlgo.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/EcalSeverityLevelAlgo.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include <vector>
#include <set>

//Less than operator for sorting EcalRecHits according to energy.
struct less_mag : public std::binary_function<EcalRecHit, EcalRecHit, bool> {
  bool operator()(EcalRecHit x, EcalRecHit y) { return x.energy() > y.energy() ; }
};

class HybridClusterAlgo
{
 private:
  //Quick typedef for position calculation.
  typedef math::XYZPoint Point;

  //Thresholds for seeds.
  double eb_st;

  //Number of steps in phi that the Hybrid algorithm will take
  //when clustering.  Remember, uses phiSteps_ in positive direction
  //and then phiSteps_ in negative direction.
  int phiSteps_;

  // et in 25
  double et25(EcalBarrelNavigator &navigator,
                const EcalRecHitCollection *hits,
                const CaloSubdetectorGeometry *geometry);


  BremRecoveryPhiRoadAlgo *phiRoadAlgo_;

  //Threshold for basic cluster.
  double eThres_;
  double eThresA_;
  double eThresB_;

  //Threshold for becoming a sub-peak in the supercluster.
  double Eseed;

  //Threshold for adding the additional two 'wing' cells to domino. 
  double Ewing;

  // do dynamic phi road
  bool dynamicPhiRoad_;

  // do dynamic ethres
  bool dynamicEThres_;

  //Map of DetId, RecHit relationship.  EcalRecHit knows what DetId it is,
  //but DetId doesn't  know what EcalRecHit it is. 
  //  std::map<DetId, EcalRecHit>  rechits_m;

  // colection of all rechits
  const EcalRecHitCollection *recHits_;

  // topology
  EcalBarrelHardcodedTopology *topo_;

  //  SuperClusterShapeAlgo* SCShape_;

  //Set of DetIds that have already been used.
  std::set<DetId> useddetids;

  // The vector of seeds:
  std::vector<EcalRecHit> seeds;

  //The vector of seed clusters:
  std::vector<reco::BasicCluster> seedClus_;

  //Map of basicclusters and what supercluster they will belong to.
  std::map<int, std::vector<reco::BasicCluster> > clustered_;

  //Control the verbosity.
  int debugLevel_;

  //algo to calulate position of clusters
  PositionCalc posCalculator_;

  // channels not to be used for seeding 
  std::vector<int> v_chstatus_; 

  // severity levels to discriminate against
  std::vector<int> v_severitylevel_;
  float severityRecHitThreshold_;
  EcalSeverityLevelAlgo::SpikeId spId_;
  float severitySpikeThreshold_;

  bool excludeFromCluster_;
  std::set<DetId> excludedCrys_;

  //New options for more robust clustering
  //and seeding
  bool robustSeeding_;
  bool robustClustering_;
  double robustTimeCons_;

 public:
   enum DebugLevel { pDEBUG = 0, pINFO = 1, pERROR = 2 }; 
  
  //The default constructor
   HybridClusterAlgo(){ };
  
  //The real constructor
  HybridClusterAlgo(double eb_str, 
		    int step,
		    double ethres,
		    double eseed,
                    double ewing,
		    std::vector<int> v_chstatus,
                    const PositionCalc& posCalculator,
                    DebugLevel debugLevel = pINFO,
		    bool dynamicEThres = false,
                    double eThresA = 0,
                    double eThresB = 0.1,
		    std::vector<int> severityToExclude=std::vector<int>(999),
		    double severityRecHitThreshold=0.08,
		    int severitySpikeId=1,
		    double severitySpikeThreshold=0,
		    bool excludeFromCluster=false,
		    bool robustSeeding = false,
		    bool robustClustering = false,
        	    double robustTimeCons = 999.
		    );
//                    const edm::ParameterSet &bremRecoveryPset,

  // destructor
  ~HybridClusterAlgo() 
  {
     if (dynamicPhiRoad_) delete phiRoadAlgo_;
     delete topo_; 
    //     delete SCShape_;
  } 

  void setDynamicPhiRoad(const edm::ParameterSet &bremRecoveryPset)
  {
     dynamicPhiRoad_ = true;
     phiRoadAlgo_ = new BremRecoveryPhiRoadAlgo(bremRecoveryPset);
  }

  //Hand over the map, the geometry, and I'll hand you back clusters.
  void makeClusters(const EcalRecHitCollection*,
		    const CaloSubdetectorGeometry * geometry,
		    reco::BasicClusterCollection &basicClusters,
		    bool regional = false,
		    const std::vector<EcalEtaPhiRegion>& regions = std::vector<EcalEtaPhiRegion>(),
		    const EcalChannelStatus *chStatus = new EcalChannelStatus());

  //Make superclusters from the references to the BasicClusters in the event.
  reco::SuperClusterCollection makeSuperClusters(const reco::CaloClusterPtrVector&);

  //The routine doing the real work.
  void mainSearch(const EcalRecHitCollection* hits, const CaloSubdetectorGeometry * geometry);
  
  //Make dominos for the hybrid method.
  double makeDomino(EcalBarrelNavigator &navigator, std::vector <EcalRecHit> &cells, float SeedCrysTime, float SeedCrysTErr);

};

#endif
