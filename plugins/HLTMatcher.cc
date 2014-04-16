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

class HLTMatcher : public edm::EDProducer {
  
public:

  HLTMatcher(const edm::ParameterSet& Pset);
  virtual ~HLTMatcher(){}

  
private:
  virtual void BeginJob();
  virtual void produce(edm::Event & event, const edm::EventSetup & EventSetup);
  virtual void EndJob(){};




  TH1F* pthist;
  std::string histname_;
  edm::InputTag   OfflinejetColl;
  edm::InnputTag  OnlinejetColl;
};

HLTMatcher::HLTMatcher(const edm::ParameterSet& Pset){
  //load options
  //  if (Pset.exists("pvCollection")) pvCollection_it = Pset.getParameter<edm::InputTag>("pvCollection");
  // else                              pvCollection_it = edm::InputTag("goodOfflinePrimaryVertices");
  
  if (Pset.exists("OfflinejetCollection")) OfflinejetColl = Pset.getParameter<edm::InputTag>("OfflinejetCollection");
  else                              OfflinejetColl = edm::InputTag("ca8PFJetsCHS");  
  if (Pset.exists("OnlinejetCollection")) OnlinejetColl = Pset.getParameter<edm::InputTag>("OnlinejetCollection");
  else                              OfflinejetColl = edm::InputTag("ca8PFJetsCHS");  
  histname_ = Pset.getParameter<std::string>("histname");
  //attempt to add histogram
  edm::Service<TFileService> fs;
  float ptbins[18] = { 0,300,320,340,360,380,400,420,440,460,480,500,550,600,650,700,800,1000};
  pthist = fs->make<TH1F>(histname_.c_str(), "JetPT", 17, ptbins);
}

void HLTMatcher::BeginJob(){
}

void HLTMatcher::produce(edm::Event& iEvent,const edm::EventSetup& iEventSetup){
  //read in jets
  edm::Handle<reco::BasicJetCollection> recoJets;
  iEvent.getByLabel(OfflinejetColl, recoJets);
  edm::Handle<reco::BasicJetCollection> hltJets;
  iEvent.getByLabel(OnlinejetColl, hltJets);

  //make collection to run over
  std::auto_ptr<reco::BasicJetCollection > RecoJetColl( new reco::BasicJetCollection (*recoJets));
  std::auto_ptr<reco::BasicJetCollection > hltJetColl( new reco::BasicJetCollection (*hltJets));

  for(size_t i = 0; i<RecoJetColl->size();i++){
    //get the offlinehard jet
    const reco::BasicJet& irecojet = (*RecoJetColl)[i];
    //match the reco jet with hlt jet by deltaR and get iterator
    int ihlt = GetIteratorByDeltaR(irecojet,hltJetColl);
    if(ihlt==-1) continue;
    const reco::BasicJet& ihltjet = (*hltJetColl)[ihlt];

    //get jet masses;
    float recomass = irecojet.mass();
    float hltmass  = ihltjet.mass();

    //get the subjets
    reco::Jet::Constituents recosubjets = irecojet.getJetConstituents();
    reco::Jet::Constituents hltsubjets = ihltjet.getJetConstituents();
    bool reco_pass = toptag(recosubjets,recomass);
    bool hlt_pass = toptag(hltsubjets,hltmass);


    //make the both jets pass toptag for hist to be filled

    if(reco_pass && hlt_pass){
      float pt = irecojet.pt();
      cout<<"Jet pt is: "<<pt<<endl;
      pthist->Fill(pt);
    }
  }  
}    

int GetIteratorByDeltaR(const reco::BasicJet& recojet,std::auto_ptr<reco::BasicJetCollection > hltJetColl){

  float eta1 = recojet.eta();
  float ph1 = recojet.phi();

  int iter=0;
  float mindR=1000;

  for(int i=0; i<hltJetColl->size(),i++){
    const reco::BasicJet& ihltjet = (*hltJetColl)[i];
    float eta2 = ihltjet.eta();
    float phi2 = ihltjet.phi();

    float dR = pow( pow(eta1-eta2,2) + pow(phi1-ph2,2), 0.5);
    if(dR<mindR){
      mindR=dR;
      iter=i;
    }
  }
  if(midR>0.1) return -1;
  return iter;
}

bool toptag(reco::Jet::Constituents subjets,float jetmass){
  //instantiate min mass
  float minmass = 99999.;

  //find the mass drop first cut on number of subjets >=3
  if(subjets.size()>= 3){
    // Take the highest 3 pt subjets for cuts
    sort ( subjets.begin(), subjets.end(), GreaterByPtCandPtrUser() );
    // Now look at the subjets that were formed
    for(int isub=0; isub<2; ++isub){
      // Get this subjet
      reco::Jet::Constituent icandJet = subjets[isub];
      // Now look at the "other" subjets than this one, form the minimum invariant mass pairing
      for ( int jsub = isub + 1; jsub < 3; ++jsub ) {
	// Get the second subjet
	reco::Jet::Constituent jcandJet = subjets[jsub];
	//Combine the 4 vectors
	reco::Candidate::LorentzVector Cand = icandJet->p4() + jcandJet->p4();
	// Find the minimum mass pairing. 
	if ( fabs( Cand.mass() ) < minmass ) {
	  minmass = Cand.mass();
	  }  
      }// end second loop over subjets
    }// end first loop over subjets
  }// endif 3 subjets

  if(subjets.size()>=3 && minmass>=50 && jetmass()>140 && jetmass()<230){
    return true;
  }
  else return false;
}

//define this as a plug-in
DEFINE_FWK_MODULE(HLTMatcher);
