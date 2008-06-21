#ifndef RecoEgamma_EcalSelectDigis_EcalDigiSelectAnalyzer_H
#define RecoEgamma_EcalSelectDigis_EcalDigiSelectAnalyzer_H
 
#include <memory>
 
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
 
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
 
#include "FWCore/ParameterSet/interface/ParameterSet.h"
 
#include <string>
#include "TH1.h"
#include "TProfile.h"
 
class TFile;
class TProfile;

//
// class declaration
//
class EcalDigiSelectAnalyzer : public edm::EDAnalyzer {
  public:
    explicit EcalDigiSelectAnalyzer( const edm::ParameterSet& );
    ~EcalDigiSelectAnalyzer();
 
 
    virtual void analyze( const edm::Event&, const edm::EventSetup& );
    virtual void beginJob(edm::EventSetup const&);
    virtual void endJob();
  private:
 
    std::string outputFile_; // output file
 
    // root file to store histograms
    TFile*  rootFile_;
 
    TH1F* ebh_AmplFill_;
    TProfile* ebh_AmpProf_;
    TH1F *ebh_nDigis;
    
    TH1F* eeh_AmplFill_;
    TProfile* eeh_AmpProf_;
    TH1F *eeh_nDigis;
	
 
 
};
#endif
