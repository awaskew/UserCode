#ifndef RecoEgamma_MaterialConversionModules_ClusterAndHitsAnalyzer_H
#define RecoEgamma_MaterialConversionModules_ClusterAndHitsAnalyzer_H


//

//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include <string>
#include "TH1.h"
#include "TTree.h"


class TFile;

//
// class declaration
//
class ClusterAndHitsAnalyzer : public edm::EDAnalyzer {
   public:
      explicit ClusterAndHitsAnalyzer( const edm::ParameterSet& );
      ~ClusterAndHitsAnalyzer();


      virtual void analyze( const edm::Event&, const edm::EventSetup& );
      virtual void beginJob(edm::EventSetup const&);
      virtual void endJob();
 private:

      std::string outputFile_;   // output file

      // root file to store histograms
      TFile*  rootFile_;

      // data members for histograms to be filled


      TTree *PhotonIDTree;
      Int_t prun;
      Int_t pevt;

      Int_t nMCpho;
      Float_t MCphoEt[100];
      Float_t MCphoEta[100];
      Float_t MCphoPhi[100];

      Float_t photonSCeta;
      Float_t photonSCphi;
      Float_t photonSCX;
      Float_t photonSCY;
      Float_t photonSCZ;
      Float_t photonet;
      Float_t photon_physeta;
      Float_t photon_physphi;
      Float_t photone;
      Float_t pho_seedeta;
      Float_t pho_seedphi;
      Float_t pho_seedet;
      Float_t pho_seeddR;
      Float_t photon_rawE;
      Float_t photon_ecaliso;
      Float_t photon_ecaliso03;
      Float_t photon_hcaliso;
      Float_t photon_hcaliso03;
      Float_t photon_trkisoptsol;
      Float_t photon_trkisoptsol03;
      Float_t photon_trkisopthol;
      Float_t photon_trkisopthol03;
      Int_t photon_ntrksol;
      Int_t photon_ntrksol03;
      Int_t photon_ntrkhol;
      Int_t photon_ntrkhol03;
      Float_t photon_etawid;
      Float_t photon_r9;
      Float_t photon_hadoverem;
      Int_t photon_ebgap;
      Int_t photon_isEB;
      Int_t photon_ebeegap;
      Int_t photon_isloosepho;
      Int_t photon_istightpho;
      Int_t photon_isConv;
      Int_t photon_hasPixSeed;
      Float_t photon_VtxX;
      Float_t photon_VtxY;
      Float_t photon_VtxZ;
      Float_t photon_phiwid;
      Float_t photon_drminjet;
      Float_t photon_drMCpho;

      Int_t nATrackerStripHits;
      Float_t ATrackerStripPosX[10000];
      Float_t ATrackerStripPosY[10000];
      Float_t ATrackerStripPosZ[10000];
      Float_t ATrackerStripPosXErr[10000];
      Float_t ATrackerStripPosYErr[10000];
      Float_t ATrackerStripPosZErr[10000];
      Float_t ATrackerStripDetXPos[10000];
      Float_t ATrackerStripDetYPos[10000];
      Float_t ATrackerStripDetZPos[10000];

      Int_t nCTrackerStripHits;
      Float_t CTrackerStripPosX[10000];
      Float_t CTrackerStripPosY[10000];
      Float_t CTrackerStripPosZ[10000];
      Float_t CTrackerStripPosXErr[10000];
      Float_t CTrackerStripPosYErr[10000];
      Float_t CTrackerStripPosZErr[10000];
      Float_t CTrackerStripDetXPos[10000];
      Float_t CTrackerStripDetYPos[10000];
      Float_t CTrackerStripDetZPos[10000];

      Int_t nATrackerPixelHits;
      Float_t ATrackerPixelPosX[10000];
      Float_t ATrackerPixelPosY[10000];
      Float_t ATrackerPixelPosZ[10000];
      Float_t ATrackerPixelPosXErr[10000];
      Float_t ATrackerPixelPosYErr[10000];
      Float_t ATrackerPixelPosZErr[10000];
      Float_t ATrackerPixelDetXPos[10000];
      Float_t ATrackerPixelDetYPos[10000];
      Float_t ATrackerPixelDetZPos[10000];

      Int_t nCTrackerPixelHits;
      Float_t CTrackerPixelPosX[10000];
      Float_t CTrackerPixelPosY[10000];
      Float_t CTrackerPixelPosZ[10000];
      Float_t CTrackerPixelPosXErr[10000];
      Float_t CTrackerPixelPosYErr[10000];
      Float_t CTrackerPixelPosZErr[10000];
      Float_t CTrackerPixelDetXPos[10000];
      Float_t CTrackerPixelDetYPos[10000];
      Float_t CTrackerPixelDetZPos[10000];
      

      
      

};
#endif
