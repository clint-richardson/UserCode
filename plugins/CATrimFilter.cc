// -*- C++ -*-
//
// Package:    CATrimFilter
// Class:      CATrimFilter
// 
/**\class CATrimFilter CATrimFilter.cc UserCode/CATrimFilter/plugins/CATrimFilter.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  dylan rankin
//         Created:  Wed, 17 Jul 2013 22:11:30 GMT
// $Id$
//
//

// system include files
#include <memory>
#include <vector>
#include <sstream>

// user include files
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "AnalysisDataFormats/TopObjects/interface/CATopJetTagInfo.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/JetReco/interface/BasicJet.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "PhysicsTools/CandUtils/interface/AddFourMomenta.h"
#include "DataFormats/Candidate/interface/CandMatchMap.h"
#include "DataFormats/HLTReco/interface/TriggerTypeDefs.h"
#include "HLTrigger/HLTcore/interface/HLTFilter.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include <Math/VectorUtil.h>

//
// class decleration
//

class CATrimFilter : public HLTFilter {
 public:
  explicit CATrimFilter(const edm::ParameterSet&);
  ~CATrimFilter();
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  virtual bool hltFilter( edm::Event&, const edm::EventSetup&, trigger::TriggerFilterObjectWithRefs & filterobject);


 private:
  // ----------member data ---------------------------

  edm::InputTag   src_;
  edm::InputTag   pfsrc_;
  double      minJetMass_;
  double      minJetpT_;

};

using namespace std;
using namespace reco;
using namespace edm;

//
// constructors and destructor
//
CATrimFilter::CATrimFilter(const edm::ParameterSet& iConfig): HLTFilter(iConfig){

  src_ = iConfig.getParameter<InputTag>("src");
  pfsrc_ = iConfig.getParameter<InputTag>("pfsrc");
  if ( iConfig.exists("minJetMass") ) minJetMass_ = iConfig.getParameter<double>("minJetMass");
  else minJetMass_ = -1;
  if ( iConfig.exists("minJetpT") ) minJetpT_ = iConfig.getParameter<double>("minJetpT");
  else minJetpT_ = -1;

}


CATrimFilter::~CATrimFilter()
{
}

void CATrimFilter::fillDescriptions(edm::ConfigurationDescriptions& descriptions){
  edm::ParameterSetDescription desc;
  makeHLTFilterDescription(desc);
  desc.add<double>("minJetMass",-1.);
  desc.add<double>("minJetpT",-1.);
  desc.add<edm::InputTag>("src",edm::InputTag("hltParticleFlow"));
  desc.add<edm::InputTag>("pfsrc",edm::InputTag("selectedPFJets"));
  desc.add<int>("triggerType",trigger::TriggerJet);
  descriptions.add("hltTrimFilter",desc);
}

// ------------ method called to for each event  ------------
bool CATrimFilter::hltFilter( edm::Event& iEvent, const edm::EventSetup& iSetup, trigger::TriggerFilterObjectWithRefs & filterobject)
{

  // Get the input list of basic jets corresponding to the hard jets
  Handle<reco::BasicJetCollection> pBasicJets;
  iEvent.getByLabel(src_, pBasicJets);

 //get corresponding pf jets
  Handle<reco::PFJetCollection> pfJets;
  iEvent.getByLabel(pfsrc_, pfJets);


  //add filter object
  if(saveTags()){
    filterobject.addCollectionTag(pfsrc_);
  }

  // Now loop over the hard jets and do kinematic cuts
  reco::BasicJetCollection::const_iterator ihardJet = pBasicJets->begin(),
    ihardJetEnd = pBasicJets->end();
  reco::PFJetCollection::const_iterator ipfJet = pfJets->begin();
  bool pass = false;

  for ( ; ihardJet != ihardJetEnd; ++ihardJet, ++ipfJet ) {

    if (ihardJet->pt() > minJetpT_ && ihardJet->mass() > minJetMass_) {
      // Get a ref to the hard jet
      reco::PFJetRef ref = reco::PFJetRef(pfJets,distance(pfJets->begin(),ipfJet));
      //add ref to event
      filterobject.addObject(trigger::TriggerJet,ref);
      pass = true;
    }

  }// end loop over hard jets

  return pass;
}

 

//define this as a plug-in
DEFINE_FWK_MODULE(CATrimFilter);
