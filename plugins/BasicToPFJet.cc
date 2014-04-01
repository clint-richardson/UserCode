#include "../interface/BasicToPFJet.h"

BasicToPFJet::BasicToPFJet(const edm::ParameterSet& PSet) :
  src_ (PSet.getParameter<edm::InputTag>("src")),
  inputToken_ (consumes<reco::BasicJetCollection>(PSet.getParameter<edm::InputTag>("src")))
{
  produces<reco::PFJetCollection>();
}

void BasicToPFJet::beginJob(){}
void BasicToPFJet::endJob(){}
void BasicToPFJet::produce( edm::Event& Event, const edm::EventSetup& EventSetup){

  //first get the basic jet collection
  edm::Handle<reco::BasicJetCollection> BasicJetColl;
  Event.getByToken(inputToken_, BasicJetColl);

  //now make the new pf jet collection
  reco::PFJetCollection* PFJetColl = new reco::PFJetCollection;
  //make the 'specific'
  reco::PFJet::Specific specific;

  //now get iterator
  reco::BasicJetCollection::const_iterator i = BasicJetColl->begin();

  //loop over basic jets and convert them to pfjets
  for(; i!=BasicJetColl->end(); i++){
    reco::PFJet pfjet(i->p4(),i->vertex(),specific);
    PFJetColl->push_back(pfjet);
  }

  std::auto_ptr<reco::PFJetCollection> selectedPFJets(PFJetColl);
  Event.put(selectedPFJets);
}


//define as plug-in for the framework
//DEFINE_FWK_MODULE(BasicToPFJet);
