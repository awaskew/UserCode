// -*- C++ -*-
//
// Package:    ClusterAnalyzer
// Class:      ClusterAnalyzer
// 
/**

 Description: Generate one big tree with all the photons, Ecal RecHits
HE rec hits, tracks, etc, for use in the EJTerm exercise.  This is more as
a reference, people can cut and paste what they need from this if they
need the functionality.  Alternatively, we COULD actually just make these
trees.

 Implementation:
     Probably shouldn't be imitated.  This isn't the most wonderful code
     in the world
*/
//
// Original Author:  A. Askew
//
//         Created:  Fri Dec 18 11:03:51 CDT 2009

///////////////////////////////////////////////////////////////////////
//                    header file for this analyzer                  //
///////////////////////////////////////////////////////////////////////
#include "../interface/ClusterAnalyzer.h"

///////////////////////////////////////////////////////////////////////
//                        CMSSW includes                             //
///////////////////////////////////////////////////////////////////////
#include "DataFormats/Provenance/interface/EventID.h"
#include "DataFormats/Provenance/interface/LuminosityBlockID.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"
#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/CaloTopology/interface/EcalBarrelTopology.h"
#include "Geometry/CaloTopology/interface/EcalEndcapTopology.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHit.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"

#include "HLTrigger/HLTfilters/interface/HLTHighLevel.h"
#include "DataFormats/Common/interface/TriggerResults.h"
///////////////////////////////////////////////////////////////////////
//                      Root include files                           //
///////////////////////////////////////////////////////////////////////
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TMath.h"
#include "TTree.h"

using namespace std;

///////////////////////////////////////////////////////////////////////
//                           Constructor                             //
///////////////////////////////////////////////////////////////////////

ClusterAnalyzer::ClusterAnalyzer(const edm::ParameterSet& ps)
{
  // Read Parameters from configuration file

  // output filename
  outputFile_   = ps.getParameter<std::string>("outputFile");
  // open output file to store histograms
  rootFile_ = TFile::Open(outputFile_.c_str(),"RECREATE");
  _hasRaw = ps.getParameter<bool>("hasRaw");
}

///////////////////////////////////////////////////////////////////////
//                            Destructor                             //
///////////////////////////////////////////////////////////////////////
ClusterAnalyzer::~ClusterAnalyzer()
{

// do anything here that needs to be done at desctruction time
// (e.g. close files, deallocate resources etc.)

  delete rootFile_;

}

///////////////////////////////////////////////////////////////////////
//    method called once each job just before starting event loop    //
///////////////////////////////////////////////////////////////////////
void 
ClusterAnalyzer::beginJob(edm::EventSetup const&)
{

  // go to *OUR* rootfile
  rootFile_->cd();
  ClusterTree = new TTree("ClusterTree","Photon tree");
  ClusterTree->Branch("Run",&prun,"Run/I");
  ClusterTree->Branch("Event", &pevt,"Event/I");
  ClusterTree->Branch("LumiBlk", &LumiBlk, "LumiBlk/I");

  ClusterTree->Branch("BX",&bx,"BX/I");

  ClusterTree->Branch("nSC", &nSC, "nSC/I");
  ClusterTree->Branch("SCeta", SCeta, "SCeta[nSC]/F");
  ClusterTree->Branch("SCphi", SCphi, "SCphi[nSC]/F");
  ClusterTree->Branch("SCX", SCX, "SCX[nSC]/F");
  ClusterTree->Branch("SCY", SCY, "SCY[nSC]/F");
  ClusterTree->Branch("SCZ", SCZ, "SCZ[nSC]/F");
  ClusterTree->Branch("SCet", SCet, "SCet[nSC]/F");
  ClusterTree->Branch("SCe", SCe,"SCe[nSC]/F");
  ClusterTree->Branch("SCrawE", SCrawE,"SCrawE[nSC]/F");

  ClusterTree->Branch("nAllCells", &nAllCells, "nAllCells/I");
  ClusterTree->Branch("AllCellsEta", AllCellsEta, "AllCellsEta[nAllCells]/F");
  ClusterTree->Branch("AllCellsPhi", AllCellsPhi, "AllCellsPhi[nAllCells]/F");
  ClusterTree->Branch("AllCellsIEta", AllCellsIEta, "AllCellsIEta[nAllCells]/I");
  ClusterTree->Branch("AllCellsIPhi", AllCellsIPhi, "AllCellsIPhi[nAllCells]/I");  
  ClusterTree->Branch("AllCellsE", AllCellsE, "AllCellsE[nAllCells]/F");
  ClusterTree->Branch("AllCellsEt", AllCellsEt, "AllCellsEt[nAllCells]/F");
  ClusterTree->Branch("AllClustered", AllClustered, "AllClustered[nAllCells]/I");
  ClusterTree->Branch("AllFlag", AllFlag,"AllFlag[nAllCells]/I");
  ClusterTree->Branch("AllTime", AllTime,"AllTime[nAllCells]/F");
 
}

///////////////////////////////////////////////////////////////////////
//                method called to for each event                    //
///////////////////////////////////////////////////////////////////////
void
ClusterAnalyzer::analyze(const edm::Event& evt, const edm::EventSetup& es)
{
  
  using namespace std;
  using namespace edm;


  bx = evt.bunchCrossing();


  Handle<EcalRecHitCollection> ecalhitsCollH;
  //evt.getByLabel("reducedEcalRecHitsEB","", ecalhitsCollH);
  evt.getByLabel("ecalRecHit","EcalRecHitsEB", ecalhitsCollH);
  const EcalRecHitCollection* rechitsCollection_ = ecalhitsCollH.product();

  edm::Handle<reco::SuperClusterCollection> SCHandle;
  evt.getByLabel(edm::InputTag("correctedHybridSuperClusters"),SCHandle);
  const reco::SuperClusterCollection hclus = *SCHandle;
  
  //Clear the Decks
  prun=0;
  pevt=0;
  //Always save the run/event/luminosityid.  Seriuosly.
  //ALWAYS!  The three together uniquely identify an event.
  //I don't know that it has happened yet, but event number is
  //allowed to cycle during a run.  You have been warned!
  edm::EventID eid = evt.id();
  prun=eid.run();
  pevt=eid.event();
  edm::LuminosityBlockID lid = (evt.getLuminosityBlock()).id();
  LumiBlk = lid.luminosityBlock();



  nSC=0;
  //Hm.  What are these?  You'll see in a minute:  some cleverness
  //is required to make sure we know what photon owns what crystal.     
  std::map<DetId, int> crysclus;
 for (reco::SuperClusterCollection::const_iterator clus = hclus.begin();
       clus!=hclus.end() && nSC<100;++clus){
   const std::vector< std::pair<DetId, float> > hitsel = clus->hitsAndFractions();
   for (int ik=0;ik<int(hitsel.size());++ik){ 
     crysclus.insert(make_pair(hitsel[ik].first, nSC));
   }
   SCeta[nSC]=clus->position().eta();
   SCphi[nSC]=clus->position().phi();
   SCX[nSC]=clus->position().x();
   SCY[nSC]=clus->position().y();
   SCZ[nSC]=clus->position().z();
   SCet[nSC]=clus->energy()*sin(clus->position().theta());
   SCe[nSC]=clus->energy();
   SCrawE[nSC]=clus->rawEnergy();

   nSC++;
 }
 


  //Clear ECAL cells.
  nAllCells=0;
  for (int i=0;i<30000;++i){
    AllCellsEta[i]=0;
    AllCellsPhi[i]=0;
    AllCellsIEta[i]=0;
    AllCellsIPhi[i]=0;
    AllCellsE[i]=0;
    AllCellsEt[i]=0;
    AllClustered[i]=0;
  }
  //Cells loop
  edm::ESHandle<CaloGeometry> geoHandle;
  es.get<CaloGeometryRecord>().get(geoHandle);
  const CaloGeometry* caloGeom = geoHandle.product();
  EcalRecHitCollection::const_iterator it;
  for (it = rechitsCollection_->begin();it!=rechitsCollection_->end()&&nAllCells<30000;++it){
    DetId blarg = it->detid();
    //Here's where the map is used to put what photon this crystal goes with. 
    std::map<DetId, int>::const_iterator selcrys = crysclus.find(blarg);
    if (selcrys==crysclus.end()){
      AllClustered[nAllCells]=-1;
    }
    else{
      AllClustered[nAllCells]=selcrys->second;
    }
    const CaloCellGeometry *this_cell = caloGeom->getGeometry(it->id());
    GlobalPoint position = this_cell->getPosition();
    float ET = it->energy() * sin(position.theta()); 
    AllCellsEta[nAllCells]=position.eta();
    AllCellsPhi[nAllCells]=position.phi();
    EBDetId dit = it->detid();
    AllCellsIEta[nAllCells]=dit.ieta();
    AllCellsIPhi[nAllCells]=dit.iphi();
    AllCellsE[nAllCells]=it->energy();
    AllCellsEt[nAllCells]=ET;
    AllFlag[nAllCells]=it->recoFlag();
    AllTime[nAllCells]=it->time();
    nAllCells++;
  }

  //We're not interested unless there's a photon present.
  if (nSC>0) ClusterTree->Fill();
  
}

///////////////////////////////////////////////////////////////////////
//    method called once each job just after ending the event loop   //
///////////////////////////////////////////////////////////////////////
void 
ClusterAnalyzer::endJob()
{

  // go to *OUR* root file and store histograms
  rootFile_->cd();

  ClusterTree->Write();
  // Write the root file (really writes the TTree)
  rootFile_->Write();
  rootFile_->Close();


}

//define this as a plug-in
DEFINE_FWK_MODULE(ClusterAnalyzer);
