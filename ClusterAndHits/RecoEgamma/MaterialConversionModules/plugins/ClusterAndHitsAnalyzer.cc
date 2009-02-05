// -*- C++ -*-

///////////////////////////////////////////////////////////////////////
//                    header file for this analyzer                  //
///////////////////////////////////////////////////////////////////////
#include "RecoEgamma/MaterialConversionModules/plugins/ClusterAndHitsAnalyzer.h"

///////////////////////////////////////////////////////////////////////
//                        CMSSW includes                             //
///////////////////////////////////////////////////////////////////////
//#include "DataFormats/EgammaCandidates/interface/PhotonIDFwd.h"
//#include "DataFormats/EgammaCandidates/interface/PhotonID.h"
//#include "DataFormats/EgammaCandidates/interface/PhotonIDAssociation.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"
#include "Geometry/CaloTopology/interface/EcalBarrelTopology.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHit.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"

#include "DataFormats/Common/interface/DetSetNew.h"
#include "DataFormats/Common/interface/DetSetVectorNew.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripMatchedRecHit2DCollection.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit2DCollection.h" 
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHitCollection.h"
#include "DataFormats/SiStripCluster/interface/SiStripClusterCollection.h"
#include "DataFormats/SiStripCluster/interface/SiStripCluster.h"
#include "DataFormats/SiPixelCluster/interface/SiPixelCluster.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit2D.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"
#include "DataFormats/SiStripDetId/interface/StripSubdetector.h"
#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"
#include "Geometry/CommonDetUnit/interface/TrackingGeometry.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "Geometry/CommonDetUnit/interface/GeomDetUnit.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
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
ClusterAndHitsAnalyzer::ClusterAndHitsAnalyzer(const edm::ParameterSet& ps)
{
  // Read Parameters from configuration file

  // output filename
  outputFile_   = ps.getParameter<std::string>("outputFile");

  // open output file to store histograms
  rootFile_ = TFile::Open(outputFile_.c_str(),"RECREATE");
}

///////////////////////////////////////////////////////////////////////
//                            Destructor                             //
///////////////////////////////////////////////////////////////////////
ClusterAndHitsAnalyzer::~ClusterAndHitsAnalyzer()
{

// do anything here that needs to be done at desctruction time
// (e.g. close files, deallocate resources etc.)

  delete rootFile_;

}

///////////////////////////////////////////////////////////////////////
//    method called once each job just before starting event loop    //
///////////////////////////////////////////////////////////////////////
void 
ClusterAndHitsAnalyzer::beginJob(edm::EventSetup const&)
{

  // go to *OUR* rootfile
  rootFile_->cd();

  PhotonIDTree = new TTree("PhotonIDTree","PhotonID Quantities");
  PhotonIDTree->Branch("prun", &prun, "prun/I");
  PhotonIDTree->Branch("pevt", &pevt, "pevt/I");

  PhotonIDTree->Branch("nMCpho", &nMCpho, "nMCpho/I");
  PhotonIDTree->Branch("MCphoEt", MCphoEt, "MCphoEt[nMCpho]/F");
  PhotonIDTree->Branch("MCphoEta", MCphoEta, "MCphoEta[nMCpho]/F");
  PhotonIDTree->Branch("MCphoPhi", MCphoPhi, "MCphoPhi[nMCpho]/F");

  PhotonIDTree->Branch("photon_physeta", &photon_physeta,"photon_physeta/F");
  PhotonIDTree->Branch("photon_physphi", &photon_physphi,"photon_physphi/F");
  PhotonIDTree->Branch("photonSCeta", &photonSCeta,"photonSCeta/F");
  PhotonIDTree->Branch("photonSCphi", &photonSCphi,"photonSCphi/F");
  PhotonIDTree->Branch("photonSCX", &photonSCX, "photonSCX/F");
  PhotonIDTree->Branch("photonSCY", &photonSCY, "photonSCY/F");
  PhotonIDTree->Branch("photonSCZ", &photonSCZ, "photonSCZ/F");
  PhotonIDTree->Branch("photonet", &photonet,"photonet/F");
  PhotonIDTree->Branch("photone", &photone,"photone/F");
  PhotonIDTree->Branch("pho_seedeta", &pho_seedeta, "pho_seedeta/F");
  PhotonIDTree->Branch("pho_seedphi", &pho_seedphi, "pho_seedphi/F");
  PhotonIDTree->Branch("pho_seedet", &pho_seedet, "pho_seedet/F");
  PhotonIDTree->Branch("pho_seeddR", &pho_seeddR, "pho_seeddR/F");
  PhotonIDTree->Branch("photon_rawE", &photon_rawE,"photon_rawET/F");
  PhotonIDTree->Branch("photon_ecaliso", &photon_ecaliso,"photon_ecaliso/F");
  PhotonIDTree->Branch("photon_hcaliso", &photon_hcaliso,"photon_hcaliso/F");
  PhotonIDTree->Branch("photon_trkisoptsol",&photon_trkisoptsol,"photon_trkisoptsol/F");
  PhotonIDTree->Branch("photon_trkisopthol",&photon_trkisopthol,"photon_trkisopthol/F");
  PhotonIDTree->Branch("photon_ntrksol", &photon_ntrksol,"photon_ntrksol/I");
  PhotonIDTree->Branch("photon_ntrkhol", &photon_ntrkhol,"photon_ntrkhol/I");

  PhotonIDTree->Branch("photon_ecaliso03", &photon_ecaliso03,"photon_ecaliso03/F");
  PhotonIDTree->Branch("photon_hcaliso03", &photon_hcaliso03,"photon_hcaliso03/F");
  PhotonIDTree->Branch("photon_trkisoptsol03",&photon_trkisoptsol03,"photon_trkisoptsol03/F");
  PhotonIDTree->Branch("photon_trkisopthol03",&photon_trkisopthol03,"photon_trkisopthol03/F");
  PhotonIDTree->Branch("photon_ntrksol03", &photon_ntrksol03,"photon_ntrksol03/I");
  PhotonIDTree->Branch("photon_ntrkhol03", &photon_ntrkhol03,"photon_ntrkhol03/I");

  PhotonIDTree->Branch("photon_etawid", &photon_etawid,"photon_etawid/F");
  PhotonIDTree->Branch("photon_phiwid", &photon_phiwid,"photon_phiwid/F");
  PhotonIDTree->Branch("photon_r9", &photon_r9, "photon_r9/F");
  PhotonIDTree->Branch("photon_hadoverem", &photon_hadoverem,"photon_hadoverem/F");
  PhotonIDTree->Branch("photon_ebgap", &photon_ebgap,"photon_ebgap/I");
  PhotonIDTree->Branch("photon_isEB", &photon_isEB,"photon_isEB/I");
  PhotonIDTree->Branch("photon_ebeegap", &photon_ebeegap,"photon_ebeegap/I");

  PhotonIDTree->Branch("photon_isloosepho", &photon_isloosepho,"photon_isloosepho/I");
  PhotonIDTree->Branch("photon_istightpho", &photon_istightpho,"photon_istightpho/I");
  PhotonIDTree->Branch("photon_isConv", &photon_isConv, "photon_isConv/I");
  PhotonIDTree->Branch("photon_hasPixSeed", &photon_hasPixSeed,"photon_hasPixSeed/I");
  PhotonIDTree->Branch("photon_VtxX", &photon_VtxX, "photon_VtxX/F");
  PhotonIDTree->Branch("photon_VtxY", &photon_VtxY, "photon_VtxY/F");
  PhotonIDTree->Branch("photon_VtxZ", &photon_VtxZ, "photon_VtxZ/F");
  PhotonIDTree->Branch("photon_drminjet", &photon_drminjet,"photon_drminjet/F");
  PhotonIDTree->Branch("photon_drMCpho", &photon_drMCpho,"photon_drMCpho/F");

  PhotonIDTree->Branch("nATrackerStripHits", &nATrackerStripHits,"nATrackerStripHits/I");
  PhotonIDTree->Branch("ATrackerStripPosX", ATrackerStripPosX, "ATrackerStripPosX[nATrackerStripHits]/F");
  PhotonIDTree->Branch("ATrackerStripPosY", ATrackerStripPosY, "ATrackerStripPosY[nATrackerStripHits]/F");
  PhotonIDTree->Branch("ATrackerStripPosZ", ATrackerStripPosZ, "ATrackerStripPosZ[nATrackerStripHits]/F");
  PhotonIDTree->Branch("ATrackerStripPosXErr", ATrackerStripPosXErr, "ATrackerStripPosXErr[nATrackerStripHits]/F");
  PhotonIDTree->Branch("ATrackerStripPosYErr", ATrackerStripPosYErr, "ATrackerStripPosYErr[nATrackerStripHits]/F");
  PhotonIDTree->Branch("ATrackerStripPosZErr", ATrackerStripPosZErr, "ATrackerStripPosZErr[nATrackerStripHits]/F");
  PhotonIDTree->Branch("ATrackerStripDetXPos", ATrackerStripDetXPos, "ATrackerStripDetXPos[nATrackerStripHits]/F");
  PhotonIDTree->Branch("ATrackerStripDetYPos", ATrackerStripDetYPos, "ATrackerStripDetYPos[nATrackerStripHits]/F");
  PhotonIDTree->Branch("ATrackerStripDetZPos", ATrackerStripDetZPos, "ATrackerStripDetZPos[nATrackerStripHits]/F");

  PhotonIDTree->Branch("nCTrackerStripHits", &nCTrackerStripHits,"nCTrackerStripHits/I");
  PhotonIDTree->Branch("CTrackerStripPosX", CTrackerStripPosX, "CTrackerStripPosX[nCTrackerStripHits]/F");
  PhotonIDTree->Branch("CTrackerStripPosY", CTrackerStripPosY, "CTrackerStripPosY[nCTrackerStripHits]/F");
  PhotonIDTree->Branch("CTrackerStripPosZ", CTrackerStripPosZ, "CTrackerStripPosZ[nCTrackerStripHits]/F");
  PhotonIDTree->Branch("CTrackerStripPosXErr", CTrackerStripPosXErr, "CTrackerStripPosXErr[nCTrackerStripHits]/F");
  PhotonIDTree->Branch("CTrackerStripPosYErr", CTrackerStripPosYErr, "CTrackerStripPosYErr[nCTrackerStripHits]/F");
  PhotonIDTree->Branch("CTrackerStripPosZErr", CTrackerStripPosZErr, "CTrackerStripPosZErr[nCTrackerStripHits]/F");
  PhotonIDTree->Branch("CTrackerStripDetXPos", CTrackerStripDetXPos, "CTrackerStripDetXPos[nCTrackerStripHits]/F");
  PhotonIDTree->Branch("CTrackerStripDetYPos", CTrackerStripDetYPos, "CTrackerStripDetYPos[nCTrackerStripHits]/F");
  PhotonIDTree->Branch("CTrackerStripDetZPos", CTrackerStripDetZPos, "CTrackerStripDetZPos[nCTrackerStripHits]/F");
  
  PhotonIDTree->Branch("nATrackerPixelHits", &nATrackerPixelHits,"nATrackerPixelHits/I");
  PhotonIDTree->Branch("ATrackerPixelPosX", ATrackerPixelPosX, "ATrackerPixelPosX[nATrackerPixelHits]/F");
  PhotonIDTree->Branch("ATrackerPixelPosY", ATrackerPixelPosY, "ATrackerPixelPosY[nATrackerPixelHits]/F");
  PhotonIDTree->Branch("ATrackerPixelPosZ", ATrackerPixelPosZ, "ATrackerPixelPosZ[nATrackerPixelHits]/F");
  PhotonIDTree->Branch("ATrackerPixelPosXErr", ATrackerPixelPosXErr, "ATrackerPixelPosXErr[nATrackerPixelHits]/F");
  PhotonIDTree->Branch("ATrackerPixelPosYErr", ATrackerPixelPosYErr, "ATrackerPixelPosYErr[nATrackerPixelHits]/F");
  PhotonIDTree->Branch("ATrackerPixelPosZErr", ATrackerPixelPosZErr, "ATrackerPixelPosZErr[nATrackerPixelHits]/F");
  PhotonIDTree->Branch("ATrackerPixelDetXPos", ATrackerPixelDetXPos, "ATrackerPixelDetXPos[nATrackerPixelHits]/F");
  PhotonIDTree->Branch("ATrackerPixelDetYPos", ATrackerPixelDetYPos, "ATrackerPixelDetYPos[nATrackerPixelHits]/F");
  PhotonIDTree->Branch("ATrackerPixelDetZPos", ATrackerPixelDetZPos, "ATrackerPixelDetZPos[nATrackerPixelHits]/F");

  PhotonIDTree->Branch("nCTrackerPixelHits", &nCTrackerPixelHits,"nCTrackerPixelHits/I");
  PhotonIDTree->Branch("CTrackerPixelPosX", CTrackerPixelPosX, "CTrackerPixelPosX[nCTrackerPixelHits]/F");
  PhotonIDTree->Branch("CTrackerPixelPosY", CTrackerPixelPosY, "CTrackerPixelPosY[nCTrackerPixelHits]/F");
  PhotonIDTree->Branch("CTrackerPixelPosZ", CTrackerPixelPosZ, "CTrackerPixelPosZ[nCTrackerPixelHits]/F");
  PhotonIDTree->Branch("CTrackerPixelPosXErr", CTrackerPixelPosXErr, "CTrackerPixelPosXErr[nCTrackerPixelHits]/F");
  PhotonIDTree->Branch("CTrackerPixelPosYErr", CTrackerPixelPosYErr, "CTrackerPixelPosYErr[nCTrackerPixelHits]/F");
  PhotonIDTree->Branch("CTrackerPixelPosZErr", CTrackerPixelPosZErr, "CTrackerPixelPosZErr[nCTrackerPixelHits]/F");
  PhotonIDTree->Branch("CTrackerPixelDetXPos", CTrackerPixelDetXPos, "CTrackerPixelDetXPos[nCTrackerPixelHits]/F");
  PhotonIDTree->Branch("CTrackerPixelDetYPos", CTrackerPixelDetYPos, "CTrackerPixelDetYPos[nCTrackerPixelHits]/F");
  PhotonIDTree->Branch("CTrackerPixelDetZPos", CTrackerPixelDetZPos, "CTrackerPixelDetZPos[nCTrackerPixelHits]/F");

}

///////////////////////////////////////////////////////////////////////
//                method called to for each event                    //
///////////////////////////////////////////////////////////////////////
void
ClusterAndHitsAnalyzer::analyze(const edm::Event& evt, const edm::EventSetup& es)
{
  
  using namespace std;
  using namespace edm;
  
  // grab photons
  Handle<reco::PhotonCollection> photonColl;
  evt.getByLabel("photons", "", photonColl);

  Handle<edm::ValueMap<Bool_t> > loosePhotonQual;
  evt.getByLabel("PhotonIDProd", "PhotonCutBasedIDLoose", loosePhotonQual);
  Handle<edm::ValueMap<Bool_t> > tightPhotonQual;
  evt.getByLabel("PhotonIDProd", "PhotonCutBasedIDTight", tightPhotonQual);

  edm::Handle<SiStripRecHit2DCollection> prphiRecHits;
  evt.getByLabel(edm::InputTag("siStripMatchedRecHits:rphiRecHit"),prphiRecHits);
  edm::Handle<SiStripRecHit2DCollection> pstereoRecHits;
  evt.getByLabel(edm::InputTag("siStripMatchedRecHits:stereoRecHit"),pstereoRecHits);
  
  const SiStripRecHit2DCollection *ALLrphiRecHits = prphiRecHits.product();
  const SiStripRecHit2DCollection *ALLstereoRecHits = pstereoRecHits.product();
  
  const SiPixelRecHitCollection *ALLpixelRecHitCollection = 0;
  edm::Handle<SiPixelRecHitCollection> pixelRecHits;
  evt.getByLabel(edm::InputTag("siPixelRecHits"), pixelRecHits);
  ALLpixelRecHitCollection = pixelRecHits.product();


  ////////////
  edm::Handle<SiStripRecHit2DCollection> cprphiRecHits;
  evt.getByLabel("ClusterAndHitsProd","savedRPhiwithEcalSuperClus",cprphiRecHits);
  edm::Handle<SiStripRecHit2DCollection> cpstereoRecHits;
  evt.getByLabel("ClusterAndHitsProd","savedSterwithEcalSuperClus",cpstereoRecHits);
  
  const SiStripRecHit2DCollection *ChosenrphiRecHits = cprphiRecHits.product();
  const SiStripRecHit2DCollection *ChosenstereoRecHits = cpstereoRecHits.product();
  
  const SiPixelRecHitCollection *ChosenpixelRecHitCollection = 0;
  edm::Handle<SiPixelRecHitCollection> cpixelRecHits;
  evt.getByLabel("ClusterAndHitsProd","savedPixwithEcalSuperClus", cpixelRecHits);
  ChosenpixelRecHitCollection = pixelRecHits.product();

  
  

  //Need tracker geometry to tell me where stuff is.
  edm::ESHandle<TrackerGeometry> tracker;
  es.get<TrackerDigiGeometryRecord>().get(tracker);
  const TrackerGeometry& geometry = *tracker;


  // grab PhotonId objects  
//   Handle<reco::PhotonIDAssociationCollection> photonIDMapColl;
//   evt.getByLabel("PhotonIDProd", "PhotonAssociatedID", photonIDMapColl);
  
  // create reference to the object types we are interested in
  const reco::PhotonCollection *photons = photonColl.product();  
  const edm::ValueMap<Bool_t> *phoMap = loosePhotonQual.product();
  const edm::ValueMap<Bool_t> *phoMapT = tightPhotonQual.product();

  prun=0;
  pevt=0;

  nMCpho=0;
   for (int i=0;i<100;++i){
    MCphoEt[i]=0;
    MCphoEta[i]=999;
    MCphoPhi[i]=999;
  }

  Handle<reco::CaloJetCollection> jetColl;
  evt.getByLabel("sisCone5CaloJets","",jetColl);
  const reco::CaloJetCollection *jets = jetColl.product();


  //This is all about getting the MC photon information.
  edm::Handle<reco::GenParticleCollection> particles;
  evt.getByLabel ("genParticles", particles);
  //Obtain generator level photons.
  std::vector<reco::Particle> mcphotons;
  float EtCut = 5; // Don't forget this is here.
  int con =0;
  int phcon=0;
  for(int iq=0; iq < int(particles->size()); ++iq){
    const reco::GenParticle & p = (*particles)[iq];
    if (p.pdgId()==22){  //if you're a photon, come take a ride!
      const reco::Candidate *mom = p.mother();
      if (mom!=NULL){//Who's your momma?
	int motherid = mom->pdgId();
	if (abs(p.pdgId())==22 && p.et()>=EtCut){
	  con++;//count all of them.
	}
	if(abs(p.pdgId())==22 
	   && p.et()>= EtCut 
	   && (fabs(motherid) < 10 || fabs(motherid)==22 || fabs(motherid)==21)//quarks or placeholder...
	   ){
      // 	  //	  const reco::Candidate *mom = p.mother();
	  reco::Particle canner(p.charge(),p.p4(), p.vertex(),p.pdgId(),1);
	  mcphotons.push_back(canner);  
	  phcon++;

	  // 	}  
	  //}
	}
      }
    }
  }//End genParticle Loop

  for (int q=0;q<int(mcphotons.size());++q){
    MCphoEt[nMCpho]  = mcphotons[q].et();
    MCphoEta[nMCpho] = mcphotons[q].eta();
    MCphoPhi[nMCpho] = mcphotons[q].phi();
    nMCpho++;
  }

  nATrackerStripHits=0;
  nCTrackerStripHits=0;
  nATrackerPixelHits=0;
  nCTrackerPixelHits=0;
  for (int ik=0;ik<10000;++ik){
    ATrackerStripPosX[ik]=0;
    ATrackerStripPosY[ik]=0;
    ATrackerStripPosZ[ik]=0;
    ATrackerStripPosXErr[ik]=0;
    ATrackerStripPosYErr[ik]=0;
    ATrackerStripPosZErr[ik]=0;
    ATrackerStripDetXPos[ik]=0;
    ATrackerStripDetYPos[ik]=0;
    ATrackerStripDetZPos[ik]=0;

    CTrackerStripPosX[ik]=0;
    CTrackerStripPosY[ik]=0;
    CTrackerStripPosZ[ik]=0;
    CTrackerStripPosXErr[ik]=0;
    CTrackerStripPosYErr[ik]=0;
    CTrackerStripPosZErr[ik]=0;
    CTrackerStripDetXPos[ik]=0;
    CTrackerStripDetYPos[ik]=0;
    CTrackerStripDetZPos[ik]=0;
      
    ATrackerPixelPosX[ik]=0;
    ATrackerPixelPosY[ik]=0;
    ATrackerPixelPosZ[ik]=0;
    ATrackerPixelPosXErr[ik]=0;
    ATrackerPixelPosYErr[ik]=0;
    ATrackerPixelPosZErr[ik]=0;
    ATrackerPixelDetXPos[ik]=0;
    ATrackerPixelDetYPos[ik]=0;
    ATrackerPixelDetZPos[ik]=0;
    
    CTrackerPixelPosX[ik]=0;
    CTrackerPixelPosY[ik]=0;
    CTrackerPixelPosZ[ik]=0;
    CTrackerPixelPosXErr[ik]=0;
    CTrackerPixelPosYErr[ik]=0;
    CTrackerPixelPosZErr[ik]=0;
    CTrackerPixelDetXPos[ik]=0;
    CTrackerPixelDetYPos[ik]=0;
    CTrackerPixelDetZPos[ik]=0;
  }

  for (SiStripRecHit2DCollection::const_iterator detvec = ALLrphiRecHits->begin();
       detvec!=ALLrphiRecHits->end() && nATrackerStripHits<10000;detvec++){
    unsigned int detid = detvec->detId();
    DetId detIdObject( detid ); 
    const GeomDetUnit *det = tracker->idToDetUnit(detIdObject);
    GlobalPoint detcenter = det->surface().toGlobal(LocalPoint(0,0,0)); 


    for (edmNew::DetSet<SiStripRecHit2D>::const_iterator rphits = detvec->begin();
	 rphits!=detvec->end() && nATrackerStripHits<10000;++rphits){

      ATrackerStripDetXPos[nATrackerStripHits]=detcenter.x();
      ATrackerStripDetYPos[nATrackerStripHits]=detcenter.y();
      ATrackerStripDetZPos[nATrackerStripHits]=detcenter.z();
      //Now only select likely hits. 
      GlobalPoint position = geometry.idToDet( 
					      rphits->geographicalId()
					      )->surface().toGlobal(
								    rphits->localPosition());
      
      ATrackerStripPosX[nATrackerStripHits]=position.x();
      ATrackerStripPosY[nATrackerStripHits]=position.y();
      ATrackerStripPosZ[nATrackerStripHits]=position.z();
      ATrackerStripPosXErr[nATrackerStripHits]=0;
      ATrackerStripPosYErr[nATrackerStripHits]=0;
      ATrackerStripPosZErr[nATrackerStripHits]=0;
      nATrackerStripHits++;
    }
  }


  for (SiStripRecHit2DCollection::const_iterator detvec = ChosenrphiRecHits->begin();
       detvec!=ChosenrphiRecHits->end() && nCTrackerStripHits<10000;detvec++){
    unsigned int detid = detvec->detId();
    DetId detIdObject( detid ); 
    const GeomDetUnit *det = tracker->idToDetUnit(detIdObject);
    GlobalPoint detcenter = det->surface().toGlobal(LocalPoint(0,0,0)); 


    for (edmNew::DetSet<SiStripRecHit2D>::const_iterator rphits = detvec->begin();
	 rphits!=detvec->end() && nCTrackerStripHits<10000;++rphits){

      CTrackerStripDetXPos[nCTrackerStripHits]=detcenter.x();
      CTrackerStripDetYPos[nCTrackerStripHits]=detcenter.y();
      CTrackerStripDetZPos[nCTrackerStripHits]=detcenter.z();
      //Now only select likely hits. 
      GlobalPoint position = geometry.idToDet( 
					      rphits->geographicalId()
					      )->surface().toGlobal(
								    rphits->localPosition());
      
      CTrackerStripPosX[nCTrackerStripHits]=position.x();
      CTrackerStripPosY[nCTrackerStripHits]=position.y();
      CTrackerStripPosZ[nCTrackerStripHits]=position.z();
      CTrackerStripPosXErr[nCTrackerStripHits]=0;
      CTrackerStripPosYErr[nCTrackerStripHits]=0;
      CTrackerStripPosZErr[nCTrackerStripHits]=0;
      nCTrackerStripHits++;
    }
  }

  for (SiStripRecHit2DCollection::const_iterator detvec = ALLstereoRecHits->begin();
       detvec!=ALLstereoRecHits->end() && nATrackerStripHits<10000;detvec++){
    unsigned int detid = detvec->detId();
    DetId detIdObject( detid ); 
    const GeomDetUnit *det = tracker->idToDetUnit(detIdObject);
    GlobalPoint detcenter = det->surface().toGlobal(LocalPoint(0,0,0)); 


    for (edmNew::DetSet<SiStripRecHit2D>::const_iterator rphits = detvec->begin();
	 rphits!=detvec->end() && nATrackerStripHits<10000;++rphits){

      ATrackerStripDetXPos[nATrackerStripHits]=detcenter.x();
      ATrackerStripDetYPos[nATrackerStripHits]=detcenter.y();
      ATrackerStripDetZPos[nATrackerStripHits]=detcenter.z();
      //Now only select likely hits. 
      GlobalPoint position = geometry.idToDet( 
					      rphits->geographicalId()
					      )->surface().toGlobal(
								    rphits->localPosition());
      
      ATrackerStripPosX[nATrackerStripHits]=position.x();
      ATrackerStripPosY[nATrackerStripHits]=position.y();
      ATrackerStripPosZ[nATrackerStripHits]=position.z();
      ATrackerStripPosXErr[nATrackerStripHits]=0;
      ATrackerStripPosYErr[nATrackerStripHits]=0;
      ATrackerStripPosZErr[nATrackerStripHits]=0;
      nATrackerStripHits++;
    }
  }


  for (SiStripRecHit2DCollection::const_iterator detvec = ChosenstereoRecHits->begin();
       detvec!=ChosenstereoRecHits->end() && nCTrackerStripHits<10000;detvec++){
    unsigned int detid = detvec->detId();
    DetId detIdObject( detid ); 
    const GeomDetUnit *det = tracker->idToDetUnit(detIdObject);
    GlobalPoint detcenter = det->surface().toGlobal(LocalPoint(0,0,0)); 


    for (edmNew::DetSet<SiStripRecHit2D>::const_iterator rphits = detvec->begin();
	 rphits!=detvec->end() && nCTrackerStripHits<10000;++rphits){

      CTrackerStripDetXPos[nCTrackerStripHits]=detcenter.x();
      CTrackerStripDetYPos[nCTrackerStripHits]=detcenter.y();
      CTrackerStripDetZPos[nCTrackerStripHits]=detcenter.z();
      //Now only select likely hits. 
      GlobalPoint position = geometry.idToDet( 
					      rphits->geographicalId()
					      )->surface().toGlobal(
								    rphits->localPosition());
      
      CTrackerStripPosX[nCTrackerStripHits]=position.x();
      CTrackerStripPosY[nCTrackerStripHits]=position.y();
      CTrackerStripPosZ[nCTrackerStripHits]=position.z();
      CTrackerStripPosXErr[nCTrackerStripHits]=0;
      CTrackerStripPosYErr[nCTrackerStripHits]=0;
      CTrackerStripPosZErr[nCTrackerStripHits]=0;
      nCTrackerStripHits++;
    }
  }


  for (SiPixelRecHitCollection::const_iterator detvec = ALLpixelRecHitCollection->begin();
	   detvec!=ALLpixelRecHitCollection->end() && nATrackerPixelHits<10000;detvec++){
    unsigned int detid = detvec->detId();
    DetId detIdObject( detid );
    const GeomDetUnit *det = tracker->idToDetUnit(detIdObject);
    GlobalPoint detcenter = det->surface().toGlobal(LocalPoint(0,0,0)); 

    for (edmNew::DetSet<SiPixelRecHit>::const_iterator pixhits = detvec->begin();
	 pixhits!=detvec->end() && nATrackerPixelHits<10000;++pixhits){
	    
      ATrackerPixelDetXPos[nATrackerPixelHits]=detcenter.x();
      ATrackerPixelDetYPos[nATrackerPixelHits]=detcenter.y();
      ATrackerPixelDetZPos[nATrackerPixelHits]=detcenter.z();
      //Now only select likely hits. 
      GlobalPoint position = geometry.idToDet( 
					      pixhits->geographicalId()
					      )->surface().toGlobal(
								    pixhits->localPosition());
      
      ATrackerPixelPosX[nATrackerPixelHits]=position.x();
      ATrackerPixelPosY[nATrackerPixelHits]=position.y();
      ATrackerPixelPosZ[nATrackerPixelHits]=position.z();
      ATrackerPixelPosXErr[nATrackerPixelHits]=0;
      ATrackerPixelPosYErr[nATrackerPixelHits]=0;
      ATrackerPixelPosZErr[nATrackerPixelHits]=0;
      nATrackerPixelHits++;

    }
  } 

  for (SiPixelRecHitCollection::const_iterator detvec = ChosenpixelRecHitCollection->begin();
	   detvec!=ChosenpixelRecHitCollection->end() && nCTrackerPixelHits<10000;detvec++){
    unsigned int detid = detvec->detId();
    DetId detIdObject( detid );
    const GeomDetUnit *det = tracker->idToDetUnit(detIdObject);
    GlobalPoint detcenter = det->surface().toGlobal(LocalPoint(0,0,0)); 

    for (edmNew::DetSet<SiPixelRecHit>::const_iterator pixhits = detvec->begin();
	 pixhits!=detvec->end() && nCTrackerPixelHits<10000;++pixhits){
	    
      CTrackerPixelDetXPos[nCTrackerPixelHits]=detcenter.x();
      CTrackerPixelDetYPos[nCTrackerPixelHits]=detcenter.y();
      CTrackerPixelDetZPos[nCTrackerPixelHits]=detcenter.z();
      //Now only select likely hits. 
      GlobalPoint position = geometry.idToDet( 
					      pixhits->geographicalId()
					      )->surface().toGlobal(
								    pixhits->localPosition());
      
      CTrackerPixelPosX[nCTrackerPixelHits]=position.x();
      CTrackerPixelPosY[nCTrackerPixelHits]=position.y();
      CTrackerPixelPosZ[nCTrackerPixelHits]=position.z();
      CTrackerPixelPosXErr[nCTrackerPixelHits]=0;
      CTrackerPixelPosYErr[nCTrackerPixelHits]=0;
      CTrackerPixelPosZErr[nCTrackerPixelHits]=0;
      nCTrackerPixelHits++;

    }
  } 
  

  int idxpho=0;
  reco::PhotonCollection::const_iterator pho;
  for (pho = (*photons).begin(); pho!= (*photons).end(); pho++){   
    
    edm::Ref<reco::PhotonCollection> photonref(photonColl, idxpho);
    photonSCeta=999;
    photonSCphi=999;
    photonSCX=999;
    photonSCY=999;
    photonSCZ=999;
    photonet=999;
    photon_physeta=999;
    photon_physphi=999;
    photone=999;
    pho_seedeta=999;
    pho_seedphi=999;
    pho_seedet=999;
    pho_seeddR=999;
    photon_rawE=999;
    photon_ecaliso=999;
    photon_hcaliso=999;
    photon_trkisoptsol=999;
    photon_trkisopthol=999;
    photon_ntrksol=999;
    photon_ntrkhol=999;
    
    photon_ecaliso03=999;
    photon_hcaliso03=999;
    photon_trkisoptsol03=999;
    photon_trkisopthol03=999;
    photon_ntrksol03=999;
    photon_ntrkhol03=999;
    
    photon_etawid=999;
    photon_r9=999;
    photon_hadoverem=999;
    photon_ebgap=999;
    photon_isEB=999;
    photon_ebeegap=999;
    photon_isloosepho=999;
    photon_istightpho=999;
    photon_isConv=999;
    photon_hasPixSeed=999;
    photon_VtxX=999;
    photon_VtxY=999;
    photon_VtxZ=999;
    photon_phiwid=999;
    photon_drminjet=999;
    photon_drMCpho=999;

    idxpho++;
    //Done resetting everything.
    Float_t drmcmin = 999;
    for (int q=0;q<nMCpho;++q){
      
      float deta = fabs(pho->eta() - MCphoEta[q]);
      float dphi = fabs(pho->phi() - MCphoPhi[q]);
      if (dphi > TMath::Pi()) dphi = TMath::Pi()*2. - dphi;
      float dr = sqrt(deta*deta + dphi*dphi);
      if (dr < drmcmin){
	drmcmin = dr;
      }
    }
    photon_drMCpho=drmcmin;
    photonSCeta=pho->superCluster()->position().eta();
    photonSCphi=pho->superCluster()->position().phi();
    photonSCX=pho->superCluster()->position().x();
    photonSCY=pho->superCluster()->position().y();
    photonSCZ=pho->superCluster()->position().z();
    photonet=pho->et();
    photon_physeta=pho->eta();
    photon_physphi=pho->phi();
    photone=pho->energy();
    pho_seedeta=pho->superCluster()->seed()->position().eta();
    pho_seedphi=pho->superCluster()->seed()->position().phi();
    pho_seedet=(pho->superCluster()->seed()->energy())/(cosh(pho->superCluster()->seed()->position().eta()));;
    
    photon_rawE=pho->superCluster()->rawEnergy();
    photon_ecaliso=pho->ecalRecHitSumConeDR04();
    photon_hcaliso=pho->hcalTowerSumConeDR04();
    photon_trkisoptsol=pho->isolationTrkSolidConeDR04();
    photon_trkisopthol=pho->isolationTrkHollowConeDR04();
    photon_ntrksol=pho->nTrkSolidConeDR04();
    photon_ntrkhol=pho->nTrkHollowConeDR04();
    
    photon_ecaliso03=pho->ecalRecHitSumConeDR03();
    photon_hcaliso03=pho->hcalRecHitSumConeDR03();
    photon_trkisoptsol03=pho->isolationTrkSolidConeDR03();
    photon_trkisopthol03=pho->isolationTrkHollowConeDR03();
    photon_ntrksol03=pho->nTrkSolidConeDR03();
    photon_ntrkhol03=pho->nTrkHollowConeDR03();
    
    photon_etawid=pho->covIetaIeta();
    photon_r9=pho->r9();
    photon_hadoverem=pho->hadronicOverEm();
    photon_ebgap=pho->isEBGap();
    photon_isEB=pho->isEB();
    photon_ebeegap=pho->isEBEEGap();
    
    photon_isloosepho= Int_t( (*phoMap)[photonref]);
    photon_istightpho= Int_t( (*phoMapT)[photonref]);
    photon_isConv=pho->hasConversionTracks();
    photon_hasPixSeed=pho->hasPixelSeed();
    photon_VtxX=pho->vertex().x();
    photon_VtxY=pho->vertex().y();
    photon_VtxZ=pho->vertex().z();
    photon_phiwid=pho->superCluster()->phiWidth();
    //Where's my jets at?
    double dRJmin=999;
    for (int i=0;i<int(jets->size());++i){
      float tet = (*jets)[i].et();
      float teta = (*jets)[i].eta();
      float tphi = (*jets)[i].phi();
      double deltaPhiJT = pho->superCluster()->position().phi()-tphi;
      double etaCurrent = pho->superCluster()->position().eta();
      if(deltaPhiJT > TMath::Pi()) deltaPhiJT -= 2.*TMath::Pi();
      if(deltaPhiJT < -1.*TMath::Pi()) deltaPhiJT += 2.*TMath::Pi();
      double dRt = std::sqrt(std::pow(etaCurrent-teta,2)+std::pow(deltaPhiJT,2));
      if (tet > 10. && (*jets)[i].emEnergyFraction()<.8
	  && dRt < dRJmin){
	dRJmin = dRt;
      } 
    }
    photon_drminjet=dRJmin;
    PhotonIDTree->Fill();
      
  } // End Loop over photons
  

  
}

///////////////////////////////////////////////////////////////////////
//    method called once each job just after ending the event loop   //
///////////////////////////////////////////////////////////////////////
void 
ClusterAndHitsAnalyzer::endJob()
{

  // go to *OUR* root file and store histograms
  rootFile_->cd();
  PhotonIDTree->Write();
  // Write the root file (really writes the TTree)
  rootFile_->Write();
  rootFile_->Close();

}

//define this as a plug-in
// DEFINE_FWK_MODULE(ClusterAndHitsAnalyzer);
