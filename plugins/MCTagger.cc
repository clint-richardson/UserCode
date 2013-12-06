
#include <memory>
#include <vector>
#include <sstream>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/JetReco/interface/Jet.h"
#include "AnalysisDataFormats/TopObjects/interface/CATopJetTagInfo.h"

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
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include <Math/VectorUtil.h>
#include <TH1.h>
#include <TH2.h>
#include <TTree.h>


using std::cout;
using std::endl;

class MCTagger : public edm::EDProducer {
  
public:

  MCTagger(const edm::ParameterSet& Pset);
  virtual ~MCTagger(){}

  
private:
  virtual void BeginJob();
  virtual void produce(edm::Event & event, const edm::EventSetup & EventSetup);
  virtual void EndJob(){};


  edm::InputTag             genParticles_it;
  edm::InputTag             pvCollection_it;
  std::vector<reco::Vertex> goodPVs;

  int findMatch(const reco::GenParticleCollection & genParticles, int idToMatch, double eta, double phi);
  double mdeltaR(double eta1, double phi1, double eta2, double phi2);
};

MCTagger::MCTagger(const edm::ParameterSet& Pset){
  //load options
  if (Pset.exists("pvCollection")) pvCollection_it = Pset.getParameter<edm::InputTag>("pvCollection");
  else                              pvCollection_it = edm::InputTag("goodOfflinePrimaryVertices");
  
  if (Pset.exists("genParticles")) genParticles_it = Pset.getParameter<edm::InputTag>("genParticles");
  else                              genParticles_it = edm::InputTag("prunedGenParticles");  
 

  produces<std::vector<int> >( "BoostedTop" ).setBranchAlias( "BoostedTop" );

}
void MCTagger::BeginJob(){ 

}

void MCTagger::produce(edm::Event& iEvent,const edm::EventSetup& iEventSetup){

  //read in particle source
  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByLabel(genParticles_it, genParticles);
  
  //make collection to run over
  std::auto_ptr<reco::GenParticleCollection > GenColl( new reco::GenParticleCollection (*genParticles));
  //make vector of ints to add to event
  std::auto_ptr<std::vector<int> > BTop;

  //First, check to see if the event has a boostedtop event, which is defined by the daughters having a deltaR<=0.8 (which is the size of our 
  //jet algorithm

  for(size_t i = 0; i < GenColl->size();i++){
    
    const reco::GenParticle & p = (*GenColl)[i];
    //    int id = p.pdgId();
    //only use status 3 particles
    if(p.status()==3){
      reco::Candidate* mother = (reco::Candidate*) p.mother();
      if(not mother) continue;
      int Mid = mother->pdgId();
      //make sure the decay is from a top (anti-top)
      if(abs(Mid)==8 ){
	//now make sure the top decays only into W+/- and b(bbar)
	if(mother->numberOfDaughters()==2 && abs( (mother->daughter(0))->pdgId() + (mother->daughter(1))->pdgId())==31){
	  //now I want to check the deltaR between the two daughters, it should be less than 0.8 for the jet to be boosted for our definition
	  double daughtDeltaR = mdeltaR( (mother->daughter(0))->eta(), (mother->daughter(0))->phi(), (mother->daughter(1))->eta(),(mother->daughter(1))->phi());
	  if(daughtDeltaR<=0.8){
	    BTop->push_back(1);
	  }
	  else{
	    BTop->push_back(0);
	  }
	}

      }
    }

 }
  //now add the collection to the event
  iEvent.put( BTop, "BoostedTop");


}


int MCTagger::findMatch(const reco::GenParticleCollection & genParticles, int idToMatch, double eta, double phi){
      
  float dRtmp = 1000;
  float closestDR = 10000.;
  int closestGenPart = -1;
  
  for(size_t j = 0; j < genParticles.size(); ++ j) {
    const reco::GenParticle & p = (genParticles).at(j);
    dRtmp = mdeltaR( eta, phi, p.eta(), p.phi());
    if ( dRtmp < closestDR && abs(p.pdgId()) == idToMatch){// && dRtmp < 0.3) {
      closestDR = dRtmp;
      closestGenPart = j;
    }//end of requirement for matching
  }//end of gen particle loop 
  return closestGenPart;
}

double MCTagger::mdeltaR(double eta1, double phi1, double eta2, double phi2) {
  return std::sqrt( pow( eta1-eta2 , 2) + pow( phi1 -phi2 , 2) );
}



//define this as a plug-in
DEFINE_FWK_MODULE(MCTagger);
