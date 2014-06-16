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

class TrimCounter : public edm::EDProducer {
  
public:

  TrimCounter(const edm::ParameterSet& Pset);
  virtual ~TrimCounter(){}

  
private:
  virtual void BeginJob();
  virtual void produce(edm::Event & event, const edm::EventSetup & EventSetup);
  virtual void EndJob(){};




  double minJetMass;
  double minJetpT;
  TH1F* hist;
  edm::InputTag             jetColl;

};

TrimCounter::TrimCounter(const edm::ParameterSet& Pset){
  //load options
  //  if (Pset.exists("pvCollection")) pvCollection_it = Pset.getParameter<edm::InputTag>("pvCollection");
  // else                              pvCollection_it = edm::InputTag("goodOfflinePrimaryVertices");
  
  if (Pset.exists("jetCollection")) jetColl = Pset.getParameter<edm::InputTag>("jetCollection");
  else                              jetColl = edm::InputTag("ca8PFJetsCHS");  
  if ( Pset.exists("minJetMass") ) minJetMass = Pset.getParameter<double>("minJetMass");
  else minJetMass = -1;
  if ( Pset.exists("minJetpT") ) minJetpT = Pset.getParameter<double>("minJetpT");
  else minJetpT = -1;

  //attempt to add histogram
  edm::Service<TFileService> fs;

  hist = fs->make<TH1F>("TrimCounterHist", "1=pass 0=fail", 2,0,2);

}

void TrimCounter::BeginJob(){
}

void TrimCounter::produce(edm::Event& iEvent,const edm::EventSetup& iEventSetup){
  //read in jets
  edm::Handle<reco::BasicJetCollection> Jets;
  iEvent.getByLabel(jetColl, Jets);

  //make collection to run over
  std::auto_ptr<reco::BasicJetCollection > JetColl( new reco::BasicJetCollection (*Jets));
  bool pass = false;
  for(size_t i = 0; i<JetColl->size();i++){
    //get the hard jet
    const reco::BasicJet& ijet = (*JetColl)[i];
    if(ijet.mass()>minJetMass && ijet.pt() > minJetpT){
      pass=true;
      continue;
    }
  }
  if(pass){
    hist->Fill(1);
  }
  else hist->Fill(0);
}

//define this as a plug-in
DEFINE_FWK_MODULE(TrimCounter);
