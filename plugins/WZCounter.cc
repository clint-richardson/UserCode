#include <memory>
#include <vector>
#include <sstream>
#include <string>

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
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/JetCollection.h"
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

struct GreaterByPtCandPtrUser {
  bool operator()( const edm::Ptr<reco::Candidate> & t1, const edm::Ptr<reco::Candidate> & t2 ) const {
    return t1->pt() > t2->pt();
  }
};

class WZCounter : public edm::EDProducer {
  
public:

  WZCounter(const edm::ParameterSet& Pset);
  virtual ~WZCounter(){}

  
private:
  virtual void BeginJob();
  virtual void produce(edm::Event & event, const edm::EventSetup & EventSetup);
  virtual void EndJob(){};




  double minWMass;
  double maxWMass;
  double massdropcut;
  TH1F* hist;
  edm::InputTag             jetColl;

};

WZCounter::WZCounter(const edm::ParameterSet& Pset){
  //load options
  //  if (Pset.exists("pvCollection")) pvCollection_it = Pset.getParameter<edm::InputTag>("pvCollection");
  // else                              pvCollection_it = edm::InputTag("goodOfflinePrimaryVertices");
  
  if (Pset.exists("jetCollection")) jetColl = Pset.getParameter<edm::InputTag>("jetCollection");
  else                              jetColl = edm::InputTag("ca8PFJetsCHS");  
  if ( Pset.exists("minWMass") ) minWMass = Pset.getParameter<double>("minWMass");
  else minWMass = 60;
  if ( Pset.exists("maxWMass") ) maxWMass = Pset.getParameter<double>("maxWMass");
  else maxWMass = 130;
  if ( Pset.exists("massdropcut") ) massdropcut = Pset.getParameter<double>("massdropcut");
  else massdropcut = 0.4;

  //attempt to add histogram
  edm::Service<TFileService> fs;

  hist = fs->make<TH1F>("WZCounterHist", "1=pass 0=fail", 2,0,2);

}

void WZCounter::BeginJob(){
}

void WZCounter::produce(edm::Event& iEvent,const edm::EventSetup& iEventSetup){
  //read in jets
  edm::Handle<reco::BasicJetCollection> Jets;
  iEvent.getByLabel(jetColl, Jets);

  //make collection to run over
  std::auto_ptr<reco::BasicJetCollection > JetColl( new reco::BasicJetCollection (*Jets));
  bool pass = false;
  for(size_t i = 0; i<JetColl->size();i++){
    //get the hard jet
    const reco::BasicJet& ijet = (*JetColl)[i];
    float massdrop=10;
    //get the subjets
    reco::Jet::Constituents subjets = ijet.getJetConstituents();
    //find the mass drop first cut on number of subjets >=3
    if(subjets.size()== 2){
      // Take the highest 3 pt subjets for cuts
      sort ( subjets.begin(), subjets.end(), GreaterByPtCandPtrUser() );
      //get the highest pt subjet
      reco::Jet::Constituent icandJet = subjets[0];

      reco::Candidate::LorentzVector isubJet = icandJet->p4();
      double imass = isubJet.mass();
      double imw = ijet.mass();
      massdrop=imass/imw;

    }// endif 2 subjets

    //make the jet pass all the top tag cuts
    if(subjets.size()==2 && massdrop<massdropcut && ijet.mass()>minWMass && ijet.mass()<maxWMass){
      pass=true;
    }
  }
  if(pass){
    hist->Fill(1);
  }
  else hist->Fill(0);
}

//define this as a plug-in
DEFINE_FWK_MODULE(WZCounter);
