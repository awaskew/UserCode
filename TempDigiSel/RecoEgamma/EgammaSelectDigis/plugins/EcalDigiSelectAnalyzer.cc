///////////////////////////////////////////////////////////////////////
//                    header file for this analyzer                  //
///////////////////////////////////////////////////////////////////////
#include "RecoEgamma/EgammaSelectDigis/plugins/EcalDigiSelectAnalyzer.h"

///////////////////////////////////////////////////////////////////////
//                        CMSSW includes                             //
///////////////////////////////////////////////////////////////////////
#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"
#include "Geometry/CaloTopology/interface/EcalBarrelTopology.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
#include "DataFormats/EcalDigi/interface/EcalDataFrame.h"
#include "DataFormats/EcalDigi/interface/EEDataFrame.h"
#include "DataFormats/EcalDigi/interface/EBDataFrame.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHit.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"


///////////////////////////////////////////////////////////////////////
//                      Root include files                           //
///////////////////////////////////////////////////////////////////////
#include "TH1F.h"
#include "TProfile.h"
#include "TFile.h"
#include "TMath.h"
#include "TTree.h"
 
using namespace std;
 
///////////////////////////////////////////////////////////////////////
//                           Constructor                             //
///////////////////////////////////////////////////////////////////////
EcalDigiSelectAnalyzer::EcalDigiSelectAnalyzer(const edm::ParameterSet& ps)
{
 
// initialize output file
  outputFile_   = ps.getParameter<std::string>("outputFile");
// open output file to store histograms
  rootFile_ = TFile::Open(outputFile_.c_str(),"RECREATE");
}
 
///////////////////////////////////////////////////////////////////////
//                            Destructor                             //
///////////////////////////////////////////////////////////////////////
EcalDigiSelectAnalyzer::~EcalDigiSelectAnalyzer()
{
 
// do anything here that needs to be done at desctruction time
// (e.g. close files, deallocate resources etc.)
  delete rootFile_;
 
}
 
 
///////////////////////////////////////////////////////////////////////
//    method called once each job just before starting event loop    //
///////////////////////////////////////////////////////////////////////
void 
EcalDigiSelectAnalyzer::beginJob(edm::EventSetup const&)
{
 
// go to *OUR* rootfile
  rootFile_->cd();
// book histograms
  ebh_AmplFill_ = new TH1F("ebh_AmplFill_","Amplitudes filled",10,0,10);
  ebh_AmpProf_ = new TProfile("ebh_AmpProf_","Profile of amplitudes",10,0,10,0,10000);
  ebh_nDigis = new TH1F("ebh_nDigis","Number of Digi saved",1000,0,1000);

  eeh_AmplFill_ = new TH1F("eeh_AmplFill_","Amplitudes filled",10,0,10);
  eeh_AmpProf_ = new TProfile("eeh_AmpProf_","Profile of amplitudes",10,0,10,0,10000);
  eeh_nDigis = new TH1F("eeh_nDigis","Number of Digi saved",1000,0,1000);

}
///////////////////////////////////////////////////////////////////////
//                method called to for each event                    //
///////////////////////////////////////////////////////////////////////
void
EcalDigiSelectAnalyzer::analyze(const edm::Event& evt, const edm::EventSetup& es)
{
   
  using namespace std;
  using namespace edm;
  //get barrel digi collection
  edm::Handle<EBDigiCollection> ebpdigis;
  const EBDigiCollection* ebdigis=0;
  evt.getByLabel("selectDigi","selectedEcalEBDigiCollection",ebpdigis);
  ebdigis = ebpdigis.product(); // get a ptr to the product


  //get endcap digi collection
  edm::Handle<EEDigiCollection> eepdigis;
  const EEDigiCollection* eedigis=0;
  evt.getByLabel("selectDigi","selectedEcalEEDigiCollection",eepdigis);
  eedigis = eepdigis.product(); // get a ptr to the product

  ebh_nDigis->Fill(ebdigis->size());
  eeh_nDigis->Fill(ebdigis->size());
      
//   edm::Handle<EBDigiCollection> ebpdigis;
//   const EBDigiCollection* ebdigis=0;
//   evt.getByLabel(edm::InputTag("ecalDigis:ebDigis"),ebpdigis);
//   ebdigis = ebpdigis.product(); // get a ptr to the product


//   //get endcap digi collection
//   edm::Handle<EEDigiCollection> eepdigis;
//   const EEDigiCollection* eedigis=0;
//   evt.getByLabel(edm::InputTag("ecalDigis:eeDigis"),eepdigis);
//   eedigis = eepdigis.product(); // get a ptr to the product


//   edm::Handle<EcalRecHitCollection> ecalhitsCollHb;
//   evt.getByLabel(edm::InputTag("ecalRecHit:EcalRecHitsEB"), ecalhitsCollHb);
//   const EcalRecHitCollection* rechitsCollectionB = ecalhitsCollHb.product();
      
  
//   edm::Handle<EcalRecHitCollection> ecalhitsCollHe;
//   evt.getByLabel(edm::InputTag("ecalRecHit:EcalRecHitsEE"), ecalhitsCollHe);
//   const EcalRecHitCollection* rechitsCollectione = ecalhitsCollHe.product();
      
//   for (EcalRecHitCollection::const_iterator eb = rechitsCollectionB->begin();
//        eb !=rechitsCollectionB->end();++eb){
//     if (eb->energy() > 2.0){
//       EBDetId det = eb->id();
//       EBDigiCollection::const_iterator it = ebdigis->find(det);
//       if (it!=ebdigis->end()){
// 	EBDataFrame myDigi = (*it);
// 	for (int iq =0;iq<myDigi.size();++iq){
// 	  ebh_AmplFill_->Fill(iq,myDigi.sample(iq).adc());
// 	  ebh_AmpProf_->Fill(iq,myDigi.sample(iq).adc());
	  
// 	}
//       }
//     }
//   }

  for (EBDigiCollection::const_iterator myDg = ebdigis->begin();
       myDg!=ebdigis->end();myDg++){
    EBDataFrame myDigi = (*myDg);
    for (int i=0;i<10;++i){
      ebh_AmplFill_->Fill(i,myDigi.sample(i).adc());
      ebh_AmpProf_->Fill(i,myDigi.sample(i).adc());
    }
  }

  for (EEDigiCollection::const_iterator myDg = eedigis->begin();
       myDg!=eedigis->end();myDg++){
    EEDataFrame myDigi = (*myDg);
    for (int i=0;i<10;++i){
      eeh_AmplFill_->Fill(i,myDigi.sample(i).adc());
      eeh_AmpProf_->Fill(i,myDigi.sample(i).adc());
    }
  }


}

///////////////////////////////////////////////////////////////////////
//    method called once each job just after ending the event loop   //
///////////////////////////////////////////////////////////////////////
void 
EcalDigiSelectAnalyzer::endJob()
{
 
// go to *OUR* root file and store histograms
  rootFile_->cd();
 
  ebh_AmplFill_->Write();
  ebh_AmpProf_->Write();
  ebh_nDigis->Write();

  eeh_AmplFill_->Write();
  eeh_AmpProf_->Write();
  eeh_nDigis->Write();


  rootFile_->Close();
 
}

