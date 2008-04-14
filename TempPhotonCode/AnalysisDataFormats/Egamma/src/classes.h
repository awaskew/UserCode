#include "DataFormats/Common/interface/Wrapper.h"
#include "AnalysisDataFormats/Egamma/interface/PhotonID.h"
#include "AnalysisDataFormats/Egamma/interface/PhotonIDFwd.h"
#include "AnalysisDataFormats/Egamma/interface/PhotonIDAssociation.h"
#include "DataFormats/Common/interface/RefToBase.h"
 
namespace
 {
    namespace
    {
       reco::PhotonIDCollection c1;
       edm::Wrapper<reco::PhotonIDCollection> w1;
       edm::Ref<reco::PhotonIDCollection> r1;
       edm::RefProd<reco::PhotonIDCollection> rp1;
       edm::RefVector<reco::PhotonIDCollection> rv1;
 
       reco::PhotonIDAssociationCollection c2;
       edm::Wrapper<reco::PhotonIDAssociationCollection> w2;
       reco::PhotonIDAssociation va1;
       reco::PhotonIDAssociationRef vr1;
       reco::PhotonIDAssociationRefProd vrp1;
       reco::PhotonIDAssociationRefVector vrv1;
    }
}
