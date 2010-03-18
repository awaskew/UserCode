#ifndef ClusterAnalyzer_H
#define ClusterAnalyzer_H

/**\class ClusterAnalyzer

 Description: Analyzer to demonstrate getting certain info.

 Implementation:
     \\\author: A. Askew Mar 2010
*/
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
class ClusterAnalyzer : public edm::EDAnalyzer {
   public:
      explicit ClusterAnalyzer( const edm::ParameterSet& );
      ~ClusterAnalyzer();

      virtual void analyze( const edm::Event&, const edm::EventSetup& );
      virtual void beginJob(edm::EventSetup const&);
      virtual void endJob();
 private:

      std::string outputFile_;   // output file
      bool _hasRaw;
      // root file to store histograms
      TFile*  rootFile_;

      TTree* ClusterTree;
      
      Int_t prun;
      Int_t pevt;
      Int_t bx;
      Int_t LumiBlk;
      Int_t nSC;
      Float_t SCX[100];
      Float_t SCY[100];
      Float_t SCZ[100];
      Float_t SCeta[100];
      Float_t SCphi[100];
      Float_t SCet[100];
      Float_t SCe[100];
      Float_t SCrawE[100];

      Int_t nAllCells;
      Float_t AllCellsEta[30000];
      Float_t AllCellsPhi[30000];
      Int_t AllCellsIEta[30000];
      Int_t AllCellsIPhi[30000];
      Float_t AllCellsE[30000];
      Float_t AllCellsEt[30000];
      Int_t AllClustered[30000];
      Int_t AllFlag[30000];
      Float_t AllTime[30000];
};
#endif
