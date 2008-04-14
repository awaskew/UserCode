#include "AnalysisDataFormats/Egamma/interface/PhotonID.h"

using namespace reco;


PhotonID::PhotonID(){
  cutBasedDecision_=false;
  isolationECal_=999;
  isolationSolidTrkCone4_=999;
  isolationHollow1TrkCone4_=999;
  nTrkSolidCone4_=999;
  nTrkHollow1TrkCone4_=999;
  isEBPho_=false;
  isEEPho_=false;
  isEBGap_=false;
  isEEGap_=false;
  isEBEEGap_=false;
  isAlsoElectron_=false;
}

PhotonID::PhotonID(bool Decision, 
		   double BCIso, 
		   double TrkCone4,
		   double HollowCone4, 
		   int nTrkCone4, int nHollow4,
		   bool EBPho, 
		   bool EEPho, 
		   bool EBGap, 
		   bool EEGap, 
		   bool EBEEGap,
		   bool isAlsoElectron){
  cutBasedDecision_=Decision;
  isolationECal_=BCIso;
  isolationSolidTrkCone4_=TrkCone4;
  isolationHollow1TrkCone4_=HollowCone4;
  nTrkSolidCone4_=nTrkCone4;
  nTrkHollow1TrkCone4_=nHollow4;
  isEBPho_=EBPho;
  isEEPho_=EEPho;
  isEBGap_=EBGap;
  isEEGap_=EEGap;
  isEBEEGap_=EBEEGap;
  isAlsoElectron_ = isAlsoElectron;
}


void PhotonID::setFiducialFlags(bool EBPho, 
				bool EEPho, 
				bool EBGap, 
				bool EEGap, 
				bool EBEEGap){
  isEBPho_=EBPho;
  isEEPho_=EEPho;
  isEBGap_=EBGap;
  isEEGap_=EEGap;
  isEBEEGap_=EBEEGap;
}


void PhotonID::setDecision(bool decision){
  cutBasedDecision_ = decision;
}


