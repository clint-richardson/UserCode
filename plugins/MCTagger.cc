
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
#include "TH1F.h"
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
  TH1F* hist;
  TH1F* passhist;
  TH1F* njethist; 

  edm::InputTag             genParticles_it;
  //  edm::InputTag             pvCollection_it;
  //std::vector<reco::Vertex> goodPVs;

  int findMatch(const reco::GenParticleCollection & genParticles, int idToMatch, double eta, double phi);
  double mdeltaR(double eta1, double phi1, double eta2, double phi2);
};

MCTagger::MCTagger(const edm::ParameterSet& Pset){
  //load options
  //  if (Pset.exists("pvCollection")) pvCollection_it = Pset.getParameter<edm::InputTag>("pvCollection");
  // else                              pvCollection_it = edm::InputTag("goodOfflinePrimaryVertices");
  
  if (Pset.exists("genParticles")) genParticles_it = Pset.getParameter<edm::InputTag>("genParticles");
  else                              genParticles_it = edm::InputTag("prunedGenParticles");  
 

  produces<std::vector<int> >( "BoostedTop" ).setBranchAlias( "BoostedTop" );
  //attempt to add histogram
  edm::Service<TFileService> fs;
  hist = fs->make<TH1F>("Tags", "tags", 4, 0 , 2);
  passhist = fs->make<TH1F>("Passes", "passes", 3, 0 , 2);
  njethist = fs->make<TH1F>("NJets", "njets", 6, 0,3);
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
  std::auto_ptr<std::vector<int> > BTop( new std::vector<int> );



  //First, check to see if the event has a boostedtop event, which is defined by the daughters having a deltaR<=0.8 (which is the size of our 
  //jet algorithm
  bool topchk = false;
  int count = 0;
  for(size_t i = 0; i < GenColl->size();i++){

    const reco::GenParticle & p = (*GenColl)[i];
    int id = p.pdgId();
    //only use status 3 particles
    if(p.status()==3){
      //reco::Candidate* mother = (reco::Candidate*) p.mother();
      //if(not mother) continue;
      // int Mid = mother->pdgId();
      //make sure the decay is from a top (anti-top)
      if(abs(id)==6 ){
	float PT = p.pt();
	//now make sure the top decays only into W+/- and b(bbar)
	if((abs( (p.daughter(0))->pdgId() + (p.daughter(1))->pdgId())==29)){
	  //now I want to check the deltaR between the three daughter quarks, so first for check that the W decays hadronically
	  // also will get index of w
	  bool hadron = false;
	  int wInd = -1;
	  for(unsigned int j=0;j<p.numberOfDaughters();j++){
	    if(abs(p.daughter(j)->pdgId())==24){
	      wInd=j;
	      if(abs(((p.daughter(j))->daughter(0))->pdgId())<11){
		hadron = true;
		}
	    }
	  }
	  //then check for b quark index
	  int bInd = -1;
	  for(unsigned int j=0;j<p.numberOfDaughters();j++){
	    if(abs(p.daughter(j)->pdgId())==5){
	      bInd=j;
	    }
	  }
	  //now get the pair-wise delta-R's of the three quarks
	  float dr1 = mdeltaR( ((p.daughter(wInd))->daughter(0))->eta(), ((p.daughter(wInd))->daughter(0))->phi(), ((p.daughter(wInd))->daughter(1))->eta(), ((p.daughter(wInd))->daughter(1))->phi()); 

	  float dr2 = mdeltaR( (p.daughter(bInd))->eta(), (p.daughter(wInd))->phi(), ((p.daughter(wInd))->daughter(1))->eta(), ((p.daughter(wInd))->daughter(1))->phi()); 
	  
	  float dr3 = mdeltaR( (p.daughter(bInd))->eta(), (p.daughter(wInd))->phi(), ((p.daughter(wInd))->daughter(0))->eta(), ((p.daughter(wInd))->daughter(0))->phi()); 
	  //find the max delta r between the three quarks
	  float maxdR = 0.0;
	  if(dr1>dr2 && dr1>dr3) maxdR = dr1;
	  if(dr2>dr1 && dr2>dr3) maxdR = dr2;
	  if(dr3>dr2 && dr3>dr1) maxdR = dr3;
			 

	  //double daughtDeltaR = mdeltaR( (p.daughter(0))->eta(), (p.daughter(0))->phi(), (p.daughter(1))->eta(),(p.daughter(1))->phi());
	  if(maxdR<=0.8 && hadron && PT >280.0){
	    BTop->push_back(1);
	    count = count+1;
	    hist->Fill(1.0);
	    
	    topchk = true;

	    cout<<"max delta R = "<<maxdR<<" and top pt is "<<PT<<endl;
	  }
	  else{
	    BTop->push_back(0);
	    hist->Fill(0.0);
	  }
	}
	else{
	  for(unsigned int j=0;j<p.numberOfDaughters();j++){
	    cout<<"Daughter "<<j<<"'s id is "<<p.daughter(j)->pdgId()<<endl;
	  }
	}

      }
    }

 }
  //now add the collection to the event
  iEvent.put( BTop, "BoostedTop");

  cout<<"Number of boosted jets is: "<<count<<endl;
  njethist->Fill(count);
  if (topchk) {
    passhist->Fill(1.0);
  }
  else {
    passhist->Fill(0.0);
  }

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
