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



/*
This module should store variables of interest in the form of a root tree
for the leading pT jet in an event. Later hopefully I can figure out how to
make it loop over all jets. It should store variables corresponding to hlt
and reco quantities for both top and wz tagging, as well as hlt quantities
from the trimmed jet producer

*/


using std::cout;
using std::endl;


bool toptag(reco::Jet::Constituents subjets,float jetmass, double minTM, double maxTM, double minMM);

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


  edm::InputTag   OfflineTopjetColl;
  edm::InputTag  OnlineTopjetColl;

  edm::InputTag   OfflineWZjetColl;
  edm::InputTag  OnlineWZjetColl;

  edm::InputTag   OfflineTrimjetColl;
  edm::InputTag  OnlineTrimjetColl;

  double minTopMass;
  double maxTopMass;
  double minMinMass;
};

hltTreeMaker::hltTreeMaker(const edm::ParameterSet& Pset){
  //load options
  
  if (Pset.exists("OfflineTopjetCollection")) OfflinejetColl = Pset.getParameter<edm::InputTag>("OfflineTopjetCollection");
  else                              OfflinejetColl = edm::InputTag("ca8PFJetsCHS");  
  if (Pset.exists("OnlineTopjetCollection")) OnlinejetColl = Pset.getParameter<edm::InputTag>("OnlineTopjetCollection");
  else                              OnlinejetColl = edm::InputTag("hltCA8TopJets"); 

  if (Pset.exists("OfflineWZjetCollection")) OfflinejetColl = Pset.getParameter<edm::InputTag>("OfflineWZjetCollection");
  else                              OfflinejetColl = edm::InputTag("ca8PFJetsCHS");  
  if (Pset.exists("OnlineWZjetCollection")) OnlinejetColl = Pset.getParameter<edm::InputTag>("OnlineWZjetCollection");
  else                              OnlinejetColl = edm::InputTag("hltCA8WZJets");  

  if (Pset.exists("OfflineTrimjetCollection")) OfflinejetColl = Pset.getParameter<edm::InputTag>("OfflineTrimjetCollection");
  else                              OfflinejetColl = edm::InputTag("ca8PFJetsCHS");  
  if (Pset.exists("OnlineTrimjetCollection")) OnlinejetColl = Pset.getParameter<edm::InputTag>("OnlineTrimjetCollection");
  else                              OnlinejetColl = edm::InputTag("hltCA8TopJets");
  
  if (Pset.exists("minTopMass")) minTopMass = Pset.getParameter<double>("minTopMass");
  else                              minTopMass = 140;
  if (Pset.exists("maxTopMass")) maxTopMass = Pset.getParameter<double>("maxTopMass");
  else                              minTopMass = 230;
  if (Pset.exists("minMinMass")) minMinMass = Pset.getParameter<double>("minMinMass");
  else                              minMinMass = 50;

  //set branch aliases

  //top jet variables first
  produces<std::vector<float> >( "hltTopJetMass" ).setBranchAlias( "hltTopJetMass" );
  produces<std::vector<float> >( "recoTopJetMass" ).setBranchAlias( "recoTopJetMass" );
  produces<std::vector<float> >( "hltTopPT" ).setBranchAlias( "hltTopPT" );
  produces<std::vector<float> >( "recoTopPT" ).setBranchAlias( "recoTopPT" );
  produces<std::vector<int> >( "N_hltTopSubjets" ).setBranchAlias( "N_hltTopSubjets" );
  produces<std::vector<int> >( "N_recoTopSubjets" ).setBranchAlias( "N_recoTopSubjets" );
  produces<std::vector<int> >( "hltTopPass" ).setBranchAlias( "hltTopPass" );
  produces<std::vector<int> >( "recoTopPass" ).setBranchAlias( "recoTopPass" );
  produces<std::vector<float> >("hltMinMass").setBranchAlias( "hltMinMass");
  produces<std::vector<float> >("recoMinMass").setBranchAlias( "recoMinMass");

  //now wz jet variables
  produces<std::vector<float> >( "hltWZJetMass" ).setBranchAlias( "hltWZJetMass" );
  produces<std::vector<float> >( "recoWZJetMass" ).setBranchAlias( "recoWZJetMass" );
  produces<std::vector<float> >( "hltWZPT" ).setBranchAlias( "hltWZPT" );
  produces<std::vector<float> >( "recoWZPT" ).setBranchAlias( "recoWZPT" );
  produces<std::vector<int> >( "N_hltWZSubjets" ).setBranchAlias( "N_hltWZSubjets" );
  produces<std::vector<int> >( "N_recoWZSubjets" ).setBranchAlias( "N_recoWZSubjets" );
  produces<std::vector<int> >( "hltWZPass" ).setBranchAlias( "hltWZPass" );
  produces<std::vector<int> >( "recoWZPass" ).setBranchAlias( "recoWZPass" );
  produces<std::vector<float> >("hltMassDrop").setBranchAlias( "hltMassDrop");
  produces<std::vector<float> >("recoMassDrop").setBranchAlias( "recoMassDrop");


  //now trimmed jet variables
  produces<std::vector<float> >( "hltTrimJetMass" ).setBranchAlias( "hltTrimJetMass" );
  produces<std::vector<float> >( "recoTrimJetMass" ).setBranchAlias( "recoTrimJetMass" );
  produces<std::vector<float> >( "hltTrimPT" ).setBranchAlias( "hltTrimPT" );
  produces<std::vector<int> >( "N_hltTrimSubjets" ).setBranchAlias( "N_hltTrimSubjets" );

  /*note that I'm not including a trim pass since we want to investigate moving this
  variable around anyway and it's easy enough to plot things with cuts using trees
  I am including hlt top and wz pass since those are more fixed cuts right now though
  with the stored information one should be able to redefine them as needed afterwards
  */
  
}

void hltTreeMaker::BeginJob(){
}

void hltTreeMaker::produce(edm::Event& iEvent,const edm::EventSetup& iEventSetup){

  //read in top jets
  edm::Handle<reco::BasicJetCollection> recoTopJets;
  iEvent.getByLabel(OfflineTopjetColl, recoTopJets);
  edm::Handle<reco::BasicJetCollection> hltTopJets;
  iEvent.getByLabel(OnlineTopjetColl, hltTopJets);
  std::auto_ptr<reco::BasicJetCollection > hltTopJetColl( new reco::BasicJetCollection (*hltTopJets));

  //read in WZ jets
  edm::Handle<reco::BasicJetCollection> recoWZJets;
  iEvent.getByLabel(OfflineWZjetColl, recoWZJets);
  edm::Handle<reco::BasicJetCollection> hltWZJets;
  iEvent.getByLabel(OnlineWZjetColl, hltWZJets);
  std::auto_ptr<reco::BasicJetCollection > hltWZJetColl( new reco::BasicJetCollection (*hltWZJets));

  //read in Trim jets
  edm::Handle<reco::BasicJetCollection> hltTrimJets;
  iEvent.getByLabel(OnlineTrimjetColl, hltTrimJets);
  std::auto_ptr<reco::BasicJetCollection > hltTrimJetColl( new reco::BasicJetCollection (*hltTrimJets));


  //now sort them by pT, let's use recojets to sort since they should have slightly better pT values
  sort ( recoTopJets.begin(), recoTopJets.end(), GreaterByPtCandPtrUser() );

  //and get the highest pT jet
  reco::BasicJetCollection::const_iterator irecoTopjet = recoTopJets->begin();

  //match the reco jet with hlt jet by deltaR and get iterator
  float eta1 = irecoTopjet->eta();
  float phi1 = irecoTopjet->phi();
  
  size_t iter=0;
  float mindR=1000;
  size_t itemp=0;

  for(reco::BasicJetCollection::const_iterator i=hltTopJets->begin(); i!=hltTopJets->end();i++){
    float eta2 = i->eta();
    float phi2 = i->phi();
    
    float dR = pow( pow(eta1-eta2,2) + pow(phi1-phi2,2), 0.5);
    if(dR<mindR){
      mindR=dR;
      iter=itemp;
    }
    itemp+=1;
  }  
 
  //iterate to the correct hlt jet
  reco::BasicJetCollection::const_iterator ihltTopjet = hltTopJets->begin();
  for(size_t i=0;i<iter;i++){
    ihltTopjet++;
  }

  /*Ok, I'm going to get the WZ and Trimmed jets here, but I'm going to get the jet that most
   closely matches in deltaR the highest pT top jet. this is to avoid amibuities about having
  different jets be the highest pT jet in different collections - now though I need to match
  both offline and online WZ jets to the offline top jet....*/ 

  //reset iteration variables
  iter=0;
  mindR=1000;
  itemp=0;

  for(reco::BasicJetCollection::const_iterator i=hltWZJets->begin(); i!=hltWZJets->end();i++){
    float eta2 = i->eta();
    float phi2 = i->phi();
    
    float dR = pow( pow(eta1-eta2,2) + pow(phi1-phi2,2), 0.5);
    if(dR<mindR){
      mindR=dR;
      iter=itemp;
    }
    itemp+=1;
  }  
 
  //iterate to the correct hlt jet
  reco::BasicJetCollection::const_iterator ihltWZjet = hltWZJets->begin();
  for(size_t i=0;i<iter;i++){
    ihltWZjet++;
  }


  //WZ offline jets
  //reset iteration variables
  iter=0;
  mindR=1000;
  itemp=0;


  for(reco::BasicJetCollection::const_iterator i=recoWZJets->begin(); i!=recoWZJets->end();i++){
    float eta2 = i->eta();
    float phi2 = i->phi();
    
    float dR = pow( pow(eta1-eta2,2) + pow(phi1-phi2,2), 0.5);
    if(dR<mindR){
      mindR=dR;
      iter=itemp;
    }
    itemp+=1;
  }  


  //iterate to the correct reco jet
  reco::BasicJetCollection::const_iterator irecoWZjet = recoWZJets->begin();
  for(size_t i=0;i<iter;i++){
    irecoWZjet++;
  }


  //now trimmed jets
  //reset iteration variables
  iter=0;
  mindR=1000;
  itemp=0;

  for(reco::BasicJetCollection::const_iterator i=hltTrimJets->begin(); i!=hltTrimJets->end();i++){
    float eta2 = i->eta();
    float phi2 = i->phi();
    
    float dR = pow( pow(eta1-eta2,2) + pow(phi1-phi2,2), 0.5);
    if(dR<mindR){
      mindR=dR;
      iter=itemp;
    }
    itemp+=1;
  }  
 
  //iterate to the correct hlt jet
  reco::BasicJetCollection::const_iterator ihltTrimjet = hltTrimJets->begin();
  for(size_t i=0;i<iter;i++){
    ihltTrimjet++;
  }


  /*now I should have 5 jets:
   *ihltTopjet
   *irecoTopjet
   *ihltWZjet
   *recoWZjet
   *ihltTrimjet
   */
  
  
  //get jet masses;
  float recoTopmass = irecoTopjet->mass();
  float hltTopmass  = ihltTopjet->mass();
  float hltWZmass  = ihltWZjet->mass();
  float recoWZmass = irecoWZjet->mass();
  float hltTrimmass  = ihltTrimjet->mass();


  //get the subjets
  reco::Jet::Constituents recoTopsubjets = irecoTopjet->getJetConstituents();
  reco::Jet::Constituents hltTopsubjets = ihltTopjet->getJetConstituents();
  reco::Jet::Constituents recoWZsubjets = irecoWZjet->getJetConstituents();
  reco::Jet::Constituents hltWZsubjets = ihltWZjet->getJetConstituents();
  reco::Jet::Constituents hltTrimsubjets = ihltTrimjet->getJetConstituents();

  bool recoTopPass = toptag(recosubjets,recomass,minTopMass,maxTopMass,minMinMass);
  bool hltTopPass = toptag(hltsubjets,hltmass,minTopMass,maxTopMass,minMinMass);

  bool recoWZPass = WZtag(recosubjets,recomass,recosjhist,recommhist_nc,recommhist_c,minWZMass,maxWZMass,minMinMass);
  bool hltWZPass = WZtag(hltsubjets,hltmass,hltsjhist,hltmmhist_nc,hltmmhist_c,minWZMass,maxWZMass,minMinMass);

  
}    

bool toptag(reco::Jet::Constituents subjets,float jetmass, double minTM, double maxTM, double minMM){
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
  
  if(subjets.size()>=3 && jetmass>minTM && jetmass<maxTM){
      if(minmass>minMM) return true;
      else return false;
  }
  else return false;
}

//define this as a plug-in
DEFINE_FWK_MODULE(hltTreeMaker);
