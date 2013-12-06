#ifndef UserCode_CAWZTagFilter_h
#define UserCode_CAWZTagFilter_h

// system include files
#include <memory>
#include <vector>
#include <sstream>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/JetReco/interface/Jet.h"
#include "AnalysisDataFormats/TopObjects/interface/CATopJetTagInfo.h"


// user include files
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/Math/interface/deltaR.h"

#include "DataFormats/Common/interface/View.h"
#include "DataFormats/JetReco/interface/BasicJet.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "PhysicsTools/CandUtils/interface/AddFourMomenta.h"
#include "DataFormats/Candidate/interface/CandMatchMap.h"

#include <Math/VectorUtil.h>
#include <TH1.h>
#include <TH2.h>
#include <TTree.h>

class CAWZJetHelperUser : public std::unary_function<reco::Jet, reco::CATopJetProperties> {
 public:

  CAWZJetHelperUser(double massdropcut) :
  massdropcut_(massdropcut)
  {}

  reco::CATopJetProperties operator()( reco::Jet const & ihardJet ) const;
  
 protected:
  double      massdropcut_;

};



struct GreaterByPtCandPtrUser {
  bool operator()( const edm::Ptr<reco::Candidate> & t1, const edm::Ptr<reco::Candidate> & t2 ) const {
    return t1->pt() > t2->pt();
  }
};



//
// class decleration
//

class CAWZTagFilter : public edm::EDFilter {
   public:
      explicit CAWZTagFilter(const edm::ParameterSet&);
      ~CAWZTagFilter();


   private:
      virtual void beginJob() ;
      virtual bool filter( edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      // ----------member data ---------------------------

  edm::InputTag   src_;

  double      minWMass_;
  double      maxWMass_;
  double      massdropcut_;
  bool        verbose_;

};

reco::CATopJetProperties CAWZJetHelperUser::operator()( reco::Jet const & ihardJet ) const {
  reco::CATopJetProperties properties;
  // Get subjets
  reco::Jet::Constituents subjets = ihardJet.getJetConstituents();
  properties.nSubJets = subjets.size();  // number of subjets
  properties.wMass = 999999.;                  // best W mass
  properties.topMass = 999999.;
  properties.minMass = -1;

  if (properties.nSubJets == 2) {
 
    sort ( subjets.begin(), subjets.end(), GreaterByPtCandPtrUser() );

    reco::Jet::Constituent icandJet = subjets[0];

    reco::Candidate::LorentzVector isubJet = icandJet->p4();
    double imass = isubJet.mass();
    double imw = ihardJet.mass();

    if (imass/imw < massdropcut_) {
     // Get the candidate mass
       properties.wMass = imw;
    }
  }
  
  return properties;
}

#endif
