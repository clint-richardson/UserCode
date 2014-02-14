// -*- C++ -*-
//
// Package:    CATopTagFilter
// Class:      CATopTagFilter
// 
/**\class CATopTagFilter CATopTagFilter.cc UserCode/CATopTagFilter/plugins/CATopTagFilter.cc

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

#include "../interface/CATopTagFilter.h"

using namespace std;
using namespace reco;
using namespace edm;

//
// constructors and destructor
//
CATopTagFilter::CATopTagFilter(const edm::ParameterSet& iConfig) : HLTFilter(iConfig),
								   src_(iConfig.getParameter<InputTag>("src")),
								   inputTag_(iConfig.getParameter<InputTag>("inputTag"))
{

  if ( iConfig.exists("TopMass") ) TopMass_ = iConfig.getParameter<double>("TopMass");
  else TopMass_ = 171.;
  if ( iConfig.exists("minTopMass") ) minTopMass_ = iConfig.getParameter<double>("minTopMass");
  else minTopMass_ = -1;
  if ( iConfig.exists("maxTopMass") ) maxTopMass_ = iConfig.getParameter<double>("maxTopMass");
  else maxTopMass_ = 999999;
  if ( iConfig.exists("WMass") ) WMass_ = iConfig.getParameter<double>("WMass");
  else WMass_ = 80.4;
  if ( iConfig.exists("minWMass") ) minWMass_ = iConfig.getParameter<double>("minWMass");
  else minWMass_ = -1;
  if ( iConfig.exists("maxWMass") ) maxWMass_ = iConfig.getParameter<double>("maxWMass");
  else maxWMass_ = 999999;
  if ( iConfig.exists("minMinMass") ) minMinMass_ = iConfig.getParameter<double>("minMinMass");
  else minMinMass_ = -1;
  if ( iConfig.exists("maxMinMass") ) maxMinMass_ = iConfig.getParameter<double>("maxMinMass");
  else maxMinMass_ = 999999;
  if ( iConfig.exists("verbose") ) verbose_ = iConfig.getParameter<bool>("verbose");
  else verbose_ = false;
}


CATopTagFilter::~CATopTagFilter()
{
}


// ------------ method called to for each event  ------------
bool CATopTagFilter::hltFilter( edm::Event& iEvent, const edm::EventSetup& iSetup, trigger::TriggerFilterObjectWithRefs & filterobject)
{

  // Get the input list of basic jets corresponding to the hard jets
  Handle<View<Jet> > pBasicJets;
  iEvent.getByLabel(src_, pBasicJets);

  // Get a convenient handle
  View<Jet> const & hardJets = *pBasicJets;

  CATopJetHelperUser helper( TopMass_, WMass_ );
  CATopJetProperties properties;

  // Now loop over the hard jets and do kinematic cuts
  View<Jet>::const_iterator ihardJet = hardJets.begin(),
    ihardJetEnd = hardJets.end();
  size_t iihardJet = 0;
  bool pass = false;
  for ( ; ihardJet != ihardJetEnd; ++ihardJet, ++iihardJet ) {

    if (ihardJet->pt() < 350) continue;

//     if ( verbose_ ) cout << "Processing ihardJet with pt = " << ihardJet->pt() << endl;

    // Initialize output variabless
    // Get a ref to the hard jet
    RefToBase<Jet> ref( pBasicJets, iihardJet );    
    // Get properties
    properties = helper( *ihardJet );

    if (verbose_){
      cout<<"Found high-pt jet while top-tagging;"<<endl;
      cout<<"wMass:    "<<(properties.wMass)<<endl;
      cout<<"minMass:  "<<(properties.minMass)<<endl;
      cout<<"topMass:  "<<(properties.topMass)<<endl;
      cout<<"nSubJets: "<<(properties.nSubJets)<<endl; 
    }

    if (properties.minMass < minMinMass_ || properties.minMass > maxMinMass_ || properties.wMass < minWMass_ || properties.wMass > maxWMass_ || properties.topMass < minTopMass_ || properties.topMass > maxTopMass_) continue;
    else {
      pass = true;
      break;
    }

    filterobject.addCollectionTag(inputTag_);


  }// end loop over hard jets

  return pass;
}
// ------------ method called once each job just before starting event loop  ------------
void 
CATopTagFilter::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
CATopTagFilter::endJob() {
}

 

//define this as a plug-in
DEFINE_FWK_MODULE(CATopTagFilter);
