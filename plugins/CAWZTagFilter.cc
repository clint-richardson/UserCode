// -*- C++ -*-
//
// Package:    CAWZTagFilter
// Class:      CAWZTagFilter
// 
/**\class CAWZTagFilter CAWZTagFilter.cc UserCode/CAWZTagFilter/plugins/CAWZTagFilter.cc

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

#include "../interface/CAWZTagFilter.h"

using namespace std;
using namespace reco;
using namespace edm;

//
// constructors and destructor
//
CAWZTagFilter::CAWZTagFilter(const edm::ParameterSet& iConfig): HLTFilter(iConfig){

  src_ = iConfig.getParameter<InputTag>("src");
  pfsrc_ = iConfig.getParameter<InputTag>("pfsrc");
  if ( iConfig.exists("minWMass") ) minWMass_ = iConfig.getParameter<double>("minWMass");
  else minWMass_ = -1;
  if ( iConfig.exists("maxWMass") ) maxWMass_ = iConfig.getParameter<double>("maxWMass");
  else maxWMass_ = 999999;
  if ( iConfig.exists("massdropcut") ) massdropcut_ = iConfig.getParameter<double>("massdropcut");
  else massdropcut_ = 1;
  if ( iConfig.exists("verbose") ) verbose_ = iConfig.getParameter<bool>("verbose");
  else verbose_ = false;
}


CAWZTagFilter::~CAWZTagFilter()
{
}

void CAWZTagFilter::fillDescriptions(edm::ConfigurationDescriptions& descriptions){
  edm::ParameterSetDescription desc;
  makeHLTFilterDescription(desc);
  desc.add<double>("maxWMass",230.);
  desc.add<double>("minWMass",140.);
  desc.add<double>("massdropcut",230.);
  desc.add<edm::InputTag>("src",edm::InputTag("hltParticleFlow"));
  desc.add<edm::InputTag>("pfsrc",edm::InputTag("selectedPFJets"));
  desc.add<bool>("verbose",false);
  desc.add<int>("triggerType",trigger::TriggerJet);
  descriptions.add("hltCA8WZTagFilter",desc);
}

// ------------ method called to for each event  ------------
bool CAWZTagFilter::hltFilter( edm::Event& iEvent, const edm::EventSetup& iSetup, trigger::TriggerFilterObjectWithRefs & filterobject)
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

  //initialize the properties
  CAWZJetHelperUser helper( massdropcut_ );
  CATopJetProperties properties;

  // Now loop over the hard jets and do kinematic cuts
  reco::BasicJetCollection::const_iterator ihardJet = pBasicJets->begin(),
    ihardJetEnd = pBasicJets->end();
  reco::PFJetCollection::const_iterator ipfJet = pfJets->begin();
  bool pass = false;

  for ( ; ihardJet != ihardJetEnd; ++ihardJet, ++ipfJet ) {

    if (ihardJet->pt() < 150) continue;

//     if ( verbose_ ) cout << "Processing ihardJet with pt = " << ihardJet->pt() << " , and mass = " << ihardJet->mass() << endl;

    // Get properties
    properties = helper( (reco::Jet&) *ihardJet );

    if (verbose_){
      cout<<"Found high-pt jet while W-tagging"<<endl;
      cout<<"nSubJets: "<<(properties.nSubJets)<<endl;
      cout<<"Mass:     "<<(properties.wMass)<<endl;
      cout<<"nSubjets: "<<ihardJet->numberOfDaughters()<<endl;
      cout<<endl;
    }

    if (properties.wMass < minWMass_ || properties.wMass > maxWMass_) continue;
    else {
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
DEFINE_FWK_MODULE(CAWZTagFilter);
