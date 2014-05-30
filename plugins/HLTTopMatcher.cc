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


bool toptag(reco::Jet::Constituents subjets,float jetmass, TH1F* sjhist, TH1F* mmhist1, TH1F* mmhist2, double minTM, double maxTM, double minMM);

struct GreaterByPtCandPtrUser {
  bool operator()( const edm::Ptr<reco::Candidate> & t1, const edm::Ptr<reco::Candidate> & t2 ) const {
    return t1->pt() > t2->pt();
  }
};

class HLTTopMatcher : public edm::EDProducer {
  
public:

  HLTTopMatcher(const edm::ParameterSet& Pset);
  virtual ~HLTTopMatcher(){}

  
private:
  virtual void BeginJob();
  virtual void produce(edm::Event & event, const edm::EventSetup & EventSetup);
  virtual void EndJob(){};




  TH1F* hltpthist;
  TH1F* recopthist;
  TH1F* postpthist;
  TH1F* deltaRhist;
  TH1F* recomasshist;
  TH1F* hltmasshist;
  TH1F* ratiomasshist_nc;
  TH1F* ratiomasshist_c;
  TH1F* ratiomasshist_f;
  TH1F* recommhist_nc;
  TH1F* hltmmhist_nc;
  TH1F* recommhist_c;
  TH1F* hltmmhist_c;
  TH1F* recosjhist;
  TH1F* hltsjhist;
  TH1F* ratiosjhist_nc;
  TH1F* ratiosjhist_c;

  std::string histname_;
  edm::InputTag   OfflinejetColl;
  edm::InputTag  OnlinejetColl;

  double minTopMass;
  double maxTopMass;
  double minMinMass;
};

HLTTopMatcher::HLTTopMatcher(const edm::ParameterSet& Pset){
  //load options
  //  if (Pset.exists("pvCollection")) pvCollection_it = Pset.getParameter<edm::InputTag>("pvCollection");
  // else                              pvCollection_it = edm::InputTag("goodOfflinePrimaryVertices");
  
  if (Pset.exists("OfflinejetCollection")) OfflinejetColl = Pset.getParameter<edm::InputTag>("OfflinejetCollection");
  else                              OfflinejetColl = edm::InputTag("ca8PFJetsCHS");  
  if (Pset.exists("OnlinejetCollection")) OnlinejetColl = Pset.getParameter<edm::InputTag>("OnlinejetCollection");
  else                              OnlinejetColl = edm::InputTag("hltCA8TopJets");  
  if (Pset.exists("minTopMass")) minTopMass = Pset.getParameter<double>("minTopMass");
  else                              minTopMass = 140;
  if (Pset.exists("maxTopMass")) maxTopMass = Pset.getParameter<double>("maxTopMass");
  else                              minTopMass = 230;
  if (Pset.exists("minMinMass")) minMinMass = Pset.getParameter<double>("minMinMass");
  else                              minMinMass = 50;
  histname_ = Pset.getParameter<std::string>("histname");
  //attempt to add histogram
  edm::Service<TFileService> fs;
  recopthist = fs->make<TH1F>("recopthist", "RECO JetPT", 100, 0.,2000.);
  hltpthist = fs->make<TH1F>("hltpthist", "HLT JetPT", 100, 0.,2000.);

  float ptbins[18] = { 0,300,320,340,360,380,400,420,440,460,480,500,550,600,650,700,800,1000};
  postpthist = fs->make<TH1F>("PostPTHist", "JetPT", 17, ptbins);

  deltaRhist = fs->make<TH1F>("deltaRhist","#Delta_{R} between offline and hlt matched jets",50,0,0.5);
  recomasshist = fs->make<TH1F>("recomasshist","Reco Jet Mass",100,0,300);
  hltmasshist = fs->make<TH1F>("hltmasshist","HLT Jet Mass",100,0,300);
  ratiomasshist_nc = fs->make<TH1F>("ratiomasshist_nc","HLT/RECO mass no cuts",50,0.,2);
  ratiomasshist_c = fs->make<TH1F>("ratiomasshist_c","HLT/RECO mass w/ HLT pass",50,0.,2);
  ratiomasshist_f = fs->make<TH1F>("ratiomasshist_f","HLT/RECO mass w/ HLT fail",50,0.,2);
  recosjhist = fs->make<TH1F>("recosjhist","reco nsubjets",5,0,5);
  hltsjhist = fs->make<TH1F>("hltsjhist","hlt nsubjets",5,0,5);
  ratiosjhist_nc = fs->make<TH1F>("ratiosjhist_nc","hlt/reco nsubjets no cuts",25,0,2.5);
  ratiosjhist_c = fs->make<TH1F>("ratiosjhist_c","hlt/reco nsubjets with hlt pass",25,0,2.5);
  recommhist_nc = fs->make<TH1F>("recommhist_nc", "reco min mass no cuts",100,0,200);
  hltmmhist_nc = fs->make<TH1F>("hltmmhist_nc", "hlt min mass no cuts", 100,0,200);
  recommhist_c = fs->make<TH1F>("recommhist_c", "reco min mass cuts",100,0,200);
  hltmmhist_c = fs->make<TH1F>("hltmmhist_c", "hlt min mass cuts", 100,0,200);
}

void HLTTopMatcher::BeginJob(){
}

void HLTTopMatcher::produce(edm::Event& iEvent,const edm::EventSetup& iEventSetup){
  //read in jets
  edm::Handle<reco::BasicJetCollection> recoJets;
  iEvent.getByLabel(OfflinejetColl, recoJets);
  edm::Handle<reco::BasicJetCollection> hltJets;
  iEvent.getByLabel(OnlinejetColl, hltJets);
  std::auto_ptr<reco::BasicJetCollection > hltJetColl( new reco::BasicJetCollection (*hltJets));
  reco::BasicJetCollection::const_iterator irecojet = recoJets->begin();
  int trial=0;
  for(; irecojet!=recoJets->end();irecojet++){
    trial++;
    //cout<<"trial: "<<trial<<endl;

    //match the reco jet with hlt jet by deltaR and get iterator
    float eta1 = irecojet->eta();
    float phi1 = irecojet->phi();

    size_t iter=0;
    float mindR=1000;
    size_t itemp=0;
    for(reco::BasicJetCollection::const_iterator i=hltJets->begin(); i!=hltJets->end();i++){
      //    cout<<"looping over "<<i<<"th jet"<<endl;
      float eta2 = i->eta();
      float phi2 = i->phi();
      
      float dR = pow( pow(eta1-eta2,2) + pow(phi1-phi2,2), 0.5);
      if(dR<mindR){
	mindR=dR;
	iter=itemp;
      }
      itemp+=1;
    }  
  
    if(mindR>1) continue;
    //cout<<"Got iterator and it's: "<<iter<<" with dR of "<<mindR<<endl;
    reco::BasicJetCollection::const_iterator ihltjet = hltJets->begin();
    for(size_t i=0;i<iter;i++){
      ihltjet++;
    }
    //fill deltaR hist
    deltaRhist->Fill(mindR);
    //cout<<"Got hltjet"<<endl;
    //get jet masses;
    float recomass = irecojet->mass();
    float hltmass  = ihltjet->mass();
    recomasshist->Fill(recomass);
    hltmasshist->Fill(hltmass);
    ratiomasshist_nc->Fill(hltmass/recomass);
    //cout<<"got jet masses"<<endl;
    //get the subjets
    reco::Jet::Constituents recosubjets = irecojet->getJetConstituents();
    //cout<<"got reco subjets"<<endl;
    reco::Jet::Constituents hltsubjets = ihltjet->getJetConstituents();
    //cout<<"got hlt subjets"<<endl;
    ratiosjhist_nc->Fill(hltsubjets.size()/recosubjets.size());
    bool reco_pass = toptag(recosubjets,recomass,recosjhist,recommhist_nc,recommhist_c,minTopMass,maxTopMass,minMinMass);
    //if(reco_pass) cout<<"got reco pass"<<endl;
    //if(!reco_pass) cout<<"got reco fail"<<endl;
    bool hlt_pass = toptag(hltsubjets,hltmass,hltsjhist,hltmmhist_nc,hltmmhist_c,minTopMass,maxTopMass,minMinMass);
    if(hlt_pass){
      //cout<<"got hlt pass"<<endl;
      ratiosjhist_c->Fill(hltsubjets.size()/recosubjets.size());
      ratiomasshist_c->Fill(hltmass/recomass);
    }
    if(!hlt_pass){
      //cout<<"got hlt fail"<<endl;
      ratiomasshist_f->Fill(hltmass/recomass);
    }
    //make the both jets pass toptag for hist to be filled
    //if(reco_pass && !hlt_pass) cout<<"I'm handling things fine"<<endl;
    if(reco_pass && hlt_pass){
      float pt1 = irecojet->pt();
      float pt2 = ihltjet->pt();
      //cout<<"hlt Jet pt is: "<<pt2<<" reco pt is: "<<pt1<<endl;
      recopthist->Fill(pt1);
      postpthist->Fill(pt1);
      //cout<<"filled reco pt hist"<<endl;
      hltpthist->Fill(pt2);
      //cout<<"filled hlt pt hist"<<endl;
    }
  
  }  

}  

bool toptag(reco::Jet::Constituents subjets,float jetmass, TH1F* sjhist,TH1F* mmhist1,TH1F* mmhist2, double minTM, double maxTM, double minMM){
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
  
  sjhist->Fill(subjets.size());
  mmhist1->Fill(minmass);

  if(subjets.size()>=3 && jetmass>minTM && jetmass<maxTM){
      mmhist2->Fill(minmass);
      if(minmass>minMM) return true;
      else return false;
  }
  else return false;
}

//define this as a plug-in
DEFINE_FWK_MODULE(HLTTopMatcher);
