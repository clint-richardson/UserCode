#ifndef UserCode_CATopTagFilter_h
#define UserCode_CATopTagFilter_h

// system include files
#include <memory>
#include <vector>
#include <sstream>


// user include files
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "AnalysisDataFormats/TopObjects/interface/CATopJetTagInfo.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/JetReco/interface/BasicJet.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "PhysicsTools/CandUtils/interface/AddFourMomenta.h"
#include "DataFormats/Candidate/interface/CandMatchMap.h"
#include "DataFormats/HLTReco/interface/TriggerTypeDefs.h"
#include "HLTrigger/HLTcore/interface/HLTFilter.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include <Math/VectorUtil.h>

class CATopJetHelperUser : public std::unary_function<reco::Jet, reco::CATopJetProperties> {
 public:

  CATopJetHelperUser(double TopMass, double WMass) :
    TopMass_(TopMass), WMass_(WMass)
    {}

    reco::CATopJetProperties operator()( reco::Jet const & ihardJet ) const;
  
 protected:
    double      TopMass_;
    double      WMass_;

};



struct GreaterByPtCandPtrUser {
  bool operator()( const edm::Ptr<reco::Candidate> & t1, const edm::Ptr<reco::Candidate> & t2 ) const {
    return t1->pt() > t2->pt();
  }
};



//
// class decleration
//

class CATopTagFilter : public HLTFilter {
 public:
  explicit CATopTagFilter(const edm::ParameterSet&);
  ~CATopTagFilter();
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  virtual bool hltFilter( edm::Event&, const edm::EventSetup&, trigger::TriggerFilterObjectWithRefs & filterobject);

 private:
  // ----------member data ---------------------------

  edm::InputTag   src_;
  edm::InputTag   pfsrc_;

  double      TopMass_;
  double      minTopMass_;
  double      maxTopMass_;
  double      WMass_;
  double      minWMass_;
  double      maxWMass_;
  double      minMinMass_;
  double      maxMinMass_;

};

reco::CATopJetProperties CATopJetHelperUser::operator()( reco::Jet const & ihardJet ) const {
  reco::CATopJetProperties properties;
  // Get subjets
  reco::Jet::Constituents subjets = ihardJet.getJetConstituents();
  properties.nSubJets = subjets.size();  // number of subjets
  properties.topMass = ihardJet.mass();      // jet mass
  properties.wMass = 999999.;                  // best W mass
  properties.minMass = 999999.;            // minimum mass pairing

  // Require at least three subjets in all cases, if not, untagged
  if ( properties.nSubJets >= 3 ) {

    // Take the highest 3 pt subjets for cuts
    sort ( subjets.begin(), subjets.end(), GreaterByPtCandPtrUser() );
       
    // Now look at the subjets that were formed
    for ( int isub = 0; isub < 2; ++isub ) {

      // Get this subjet
      reco::Jet::Constituent icandJet = subjets[isub];

      // Now look at the "other" subjets than this one, form the minimum invariant mass
      // pairing, as well as the "closest" combination to the W mass
      for ( int jsub = isub + 1; jsub < 3; ++jsub ) {

        // Get the second subjet
	reco::Jet::Constituent jcandJet = subjets[jsub];

	reco::Candidate::LorentzVector wCand = icandJet->p4() + jcandJet->p4();

        // Get the candidate mass
        double imw = wCand.mass();

        // Find the combination closest to the W mass
        if ( fabs( imw - WMass_ ) < fabs(properties.wMass - WMass_) ) {
          properties.wMass = imw;
        }
        // Find the minimum mass pairing. 
        if ( fabs( imw ) < properties.minMass ) {
          properties.minMass = imw;
        }  
      }// end second loop over subjets
    }// end first loop over subjets
  }// endif 3 subjets

  if (properties.minMass == 999999) properties.minMass=-1;

  return properties;
}

#endif
