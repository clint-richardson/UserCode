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

class hltTreeMaker : public edm::EDProducer {
  
public:

  hltTreeMaker(const edm::ParameterSet& Pset);
  virtual ~hltTreeMaker(){}

  
private:
  virtual void BeginJob();
  virtual void produce(edm::Event & event, const edm::EventSetup & EventSetup);
  virtual void EndJob(){};


  edm::InputTag   OfflinejetColl;
  edm::InputTag  OnlinejetColl;

  double minTopMass;
  double maxTopMass;
  double minMinMass;
};

hltTreeMaker::hltTreeMaker(const edm::ParameterSet& Pset){
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

  //set branch aliases
  produces<std::vector<float> >( "hltJetMass" ).setBranchAlias( "hltJetMass" );
  produces<std::vector<float> >( "recoJetMass" ).setBranchAlias( "recoJetMass" );
  produces<std::vector<float> >( "hltPT" ).setBranchAlias( "hltPT" );
  produces<std::vector<float> >( "recoPT" ).setBranchAlias( "recoPT" );
  produces<std::vector<int> >( "N_hltSubjets" ).setBranchAlias( "N_hltSubjets" );
  produces<std::vector<int> >( "N_recoSubjets" ).setBranchAlias( "N_recoSubjets" );
  produces<std::vector<int> >( "hltTopPass" ).setBranchAlias( "hltTopPass" );
  produces<std::vector<int> >( "recoTopPass" ).setBranchAlias( "recoTopPass" );
  produces<std::vector<float> >("hltMinMass").setBranchAlias( "hltMinMass");
  
}

void hltTreeMaker::BeginJob(){
}

void hltTreeMaker::produce(edm::Event& iEvent,const edm::EventSetup& iEventSetup){
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
DEFINE_FWK_MODULE(hltTreeMaker);
