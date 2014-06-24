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

using namespace std;

//declare needed functions
int Toptag(reco::Jet::Constituents subjets,float jetmass, double minTM, double maxTM, double minMM);
int WZtag(reco::Jet::Constituents subjets,float jetmass, double minWM, double maxWM, double mdc);

float minmass(reco::Jet::Constituents subjets);
float massdrop(float jetmass,reco::Jet::Constituents subjets);


struct GreaterByPtCandPtrUser {
  bool operator()( const edm::Ptr<reco::Candidate> & t1, const edm::Ptr<reco::Candidate> & t2 ) const {
    return t1->pt() > t2->pt();
  }
};



//plugin declaration
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


  edm::InputTag  OnlineTrimjetColl;
  //  edm::InputTag  OnlineCA8jetColl;

  double minTopMass;
  double maxTopMass;
  double minMinMass;

  double minWZMass;
  double maxWZMass;
  double MassDropCut;

  TTree* tree;

  //all the tree variables
  float hltTopJetMass_;
  float recoTopJetMass_;
  float hltWZJetMass_;
  float recoWZJetMass_;
  float hltTrimJetMass_;
  //  float hltCA8JetMass_;

  
  float hltTopPT_;
  float recoTopPT_;
  float hltWZPT_;
  float recoWZPT_;
  float hltTrimPT_;
  // float hltCA8PT_;

  float hltMinMass_;
  float recoMinMass_;
  float hltMassDrop_;
  float recoMassDrop_;

  int NhltTopSubjets_;
  int NrecoTopSubjets_;
  int NhltWZSubjets_;
  int NrecoWZSubjets_;
  int NhltTrimSubjets_;

  int hltTopPass_;
  int recoTopPass_;
  int hltWZPass_;
  int recoWZPass_;

  float hltTopEta_;
  float hltTopPhi_;
  float recoTopEta_;
  float recoTopPhi_;
  float hltWZEta_;
  float hltWZPhi_;
  float recoWZEta_;
  float recoWZPhi_;
  float hltTrimEta_;
  float hltTrimPhi_;

};

hltTreeMaker::hltTreeMaker(const edm::ParameterSet& Pset){
  //load options
  
  if (Pset.exists("OfflineTopjetCollection")) OfflineTopjetColl = Pset.getParameter<edm::InputTag>("OfflineTopjetCollection");
  else                              OfflineTopjetColl = edm::InputTag("ca8PFJetsCHS");  
  if (Pset.exists("OnlineTopjetCollection")) OnlineTopjetColl = Pset.getParameter<edm::InputTag>("OnlineTopjetCollection");
  else                              OnlineTopjetColl = edm::InputTag("hltCA8TopJets"); 

  if (Pset.exists("OfflineWZjetCollection")) OfflineWZjetColl = Pset.getParameter<edm::InputTag>("OfflineWZjetCollection");
  else                              OfflineWZjetColl = edm::InputTag("ca8PFJetsCHS");  
  if (Pset.exists("OnlineWZjetCollection")) OnlineWZjetColl = Pset.getParameter<edm::InputTag>("OnlineWZjetCollection");
  else                              OnlineWZjetColl = edm::InputTag("hltCA8WZJets");  


  if (Pset.exists("OnlineTrimjetCollection")) OnlineTrimjetColl = Pset.getParameter<edm::InputTag>("OnlineTrimjetCollection");
  else                              OnlineTrimjetColl = edm::InputTag("hltCA8TopJets");

  /*  if (Pset.exists("OnlineCA8jetCollection")) OnlineTrimjetColl = Pset.getParameter<edm::InputTag>("OnlineCA8jetCollection");
  else                              OnlineTrimjetColl = edm::InputTag("hltCambridgeAachen8PFJets");
  */
  
  if (Pset.exists("minTopMass")) minTopMass = Pset.getParameter<double>("minTopMass");
  else                              minTopMass = 140;
  if (Pset.exists("maxTopMass")) maxTopMass = Pset.getParameter<double>("maxTopMass");
  else                              minTopMass = 230;
  if (Pset.exists("minMinMass")) minMinMass = Pset.getParameter<double>("minMinMass");
  else                              minMinMass = 50;


  if (Pset.exists("minWZMass")) minWZMass = Pset.getParameter<double>("minWZMass");
  else                              minWZMass = 90;
  if (Pset.exists("maxWZMass")) maxWZMass = Pset.getParameter<double>("maxWZMass");
  else                              minWZMass = 170;
  if (Pset.exists("MassDropCut")) minMinMass = Pset.getParameter<double>("MassDropCut");
  else                              minMinMass = 0.4;



  //set branch aliases
  edm::Service<TFileService> fs;
  tree = fs->make<TTree>("JetVariables", "JetVariables");

  tree->Branch("hltTopJetMass", &hltTopJetMass_,"hltTopJetMass/F");
  tree->Branch("recoTopJetMass",&recoTopJetMass_,"recoTopJetMass/F");
  tree->Branch("hltTopPT",&hltTopPT_,"hltTopPT/F");
  tree->Branch("recoTopPT",&recoTopPT_,"recoTopPT/F");
  tree->Branch("NhltTopSubjets",&NhltTopSubjets_,"NhltTopSubjets/i");
  tree->Branch("NrecoTopSubjets",&NrecoTopSubjets_,"NrecoTopSubjets/i");
  tree->Branch("hltTopPass",&hltTopPass_,"hltTopPass/i");
  tree->Branch("recoTopPass",&recoTopPass_,"recoTopPass/i");
  tree->Branch("hltMinMass",&hltMinMass_,"hltMinMass/F");
  tree->Branch("recoMinMass",&recoMinMass_,"recoMinMass/F");
  tree->Branch("hltWZJetMass",&hltWZJetMass_,"hltWZJetMass/F");
  tree->Branch("recoWZJetMass",&recoWZJetMass_,"recoWZJetMass/F");
  tree->Branch("hltWZPT",&hltWZPT_,"hltWZPT/F");
  tree->Branch("recoWZPT",&recoWZPT_,"recoWZPT/F");
  tree->Branch("NhltWZSubjets",&NhltWZSubjets_,"NhltWZSubjets/i");
  tree->Branch("NrecoWZSubjets",&NrecoWZSubjets_,"NrecoWZSubjets/i");
  tree->Branch("hltWZPass",&hltWZPass_,"hltWZPass/i");
  tree->Branch("recoWZPass",&recoWZPass_,"recoWZPass/i");
  tree->Branch("hltMassDrop",&hltMassDrop_,"hltMassDrop/F");
  tree->Branch("recoMassDrop",&recoMassDrop_,"recoMassDrop/F");
  tree->Branch("hltTrimJetMass",&hltTrimJetMass_,"hltTrimJetMass/F");
  tree->Branch("hltTrimPT",&hltTrimPT_,"hltTrimPT/F");
  tree->Branch("NhltTrimSubjets",&NhltTrimSubjets_,"NhltTrimSubjets/i");
  tree->Branch("hltTopEta",&hltTopEta_,"hltTopEta/F");
  tree->Branch("hltTopPhi",&hltTopPhi_,"hltTopPhi/F");
  tree->Branch("recoTopEta",&recoTopEta_,"recoTopEta/F");
  tree->Branch("recoTopPhi",&recoTopPhi_,"recoTopPhi/F");
  tree->Branch("hltWZEta",&hltWZEta_,"hltWZEta/F");
  tree->Branch("hltWZPhi",&hltWZPhi_,"hltWZPhi/F");
  tree->Branch("recoWZEta",&recoWZEta_,"recoWZEta/F");
  tree->Branch("recoWZPhi",&recoWZPhi_,"recoWZPhi/F");
  tree->Branch("hltTrimEta",&hltTrimEta_,"hltTrimEta/F");
  tree->Branch("hltTrimPhi",&hltTrimPhi_,"hltTrimPhi/F");
  /*tree->Branch("hltCA8JetMass",&hltCA8JetMass_,"hltCA8JetMass/F");
  tree->Branch("hltCA8PT",&hltCA8PT_,"hltCA8PT/F");
  */
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

  /*  //read in Trim jets
  edm::Handle<reco::BasicJetCollection> hltCA8Jets;
  iEvent.getByLabel(OnlineCA8jetColl, hltCA8Jets);
  std::auto_ptr<reco::BasicJetCollection > hltCA8JetColl( new reco::BasicJetCollection (*hltCA8Jets));
  */

  //protection against there not being a reco top jet

  if(recoTopJets->size()==0 || hltTopJets->size()==0 || recoWZJets->size()==0 || hltWZJets->size()==0 || hltTrimJets->size()==0){ //|| hltCA8Jets->size()==0){
    recoTopJetMass_=-9999;
    recoTopPass_=-1;
    NrecoTopSubjets_=-1;
    recoTopPT_=-9999;
    recoMinMass_=-9999;
    hltTopJetMass_=-9999;
    hltTopPass_=-1;
    NhltTopSubjets_=-1;
    hltTopPT_=-9999;
    hltMinMass_=-9999;
    recoWZJetMass_=-9999;
    recoWZPass_=-1;
    NrecoWZSubjets_=-1;
    recoWZPT_=-9999;
    recoMassDrop_=9999;
    hltWZJetMass_=-9999;
    hltWZPass_=-1;
    NhltWZSubjets_=-1;
    hltWZPT_=-9999;
    hltMassDrop_=9999;
    hltTrimJetMass_=-9999;
    hltTrimPT_=-9999;
    NhltTrimSubjets_=-1;
    hltTopEta_=9999;
    hltTopPhi_=9999;
    recoTopPhi_=9999;
    recoTopEta_=9999;
    hltWZEta_=9999;
    hltWZPhi_=9999;
    recoWZEta_=9999;
    recoWZPhi_=9999;
    hltTrimEta_=9999;
    hltTrimPhi_=9999;
    //hltCA8JetMass_=-9999;
    //hltCA8PT_=-9999;
  }
  else{

    //get the top jet
    reco::BasicJetCollection::const_iterator irecoTopjet;
    
    //sort the top jets by pt
    float dummyPT=0;
    
    for(reco::BasicJetCollection::const_iterator i=recoTopJets->begin(); i!=recoTopJets->end();i++){
      if(i->pt()>dummyPT){
	dummyPT=i->pt();
	irecoTopjet=i;
      }
      
    }
    
    
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
    /*
    //now CA8 ungroomed jets
    //reset iteration variables
    iter=0;
    mindR=1000;
    itemp=0;
    
    for(reco::BasicJetCollection::const_iterator i=hltCA8Jets->begin(); i!=hltCA8Jets->end();i++){
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
    reco::BasicJetCollection::const_iterator ihltCA8jet = hltCA8Jets->begin();
    for(size_t i=0;i<iter;i++){
    ihltCA8jet++;
    }
    */
    
    
    /*now I should have 6 jets:
     *ihltTopjet
     *irecoTopjet
     *ihltWZjet
     *recoWZjet
     *ihltTrimjet
     *ihltCA8jet
     */
    
    
    
    //get jet masses;
    recoTopJetMass_ =irecoTopjet->mass();
    hltTopJetMass_ =ihltTopjet->mass();
    hltWZJetMass_ = ihltWZjet->mass();
    recoWZJetMass_ = irecoWZjet->mass();
    hltTrimJetMass_ = ihltTrimjet->mass();
    //hltCA8JetMass_ = ihltCA8jet->mass();

    //get the subjets
    reco::Jet::Constituents recoTopsubjets = irecoTopjet->getJetConstituents();
    reco::Jet::Constituents hltTopsubjets = ihltTopjet->getJetConstituents();
    reco::Jet::Constituents recoWZsubjets = irecoWZjet->getJetConstituents();
    reco::Jet::Constituents hltWZsubjets = ihltWZjet->getJetConstituents();
    reco::Jet::Constituents hltTrimsubjets = ihltTrimjet->getJetConstituents();
     
    
    recoTopPass_ =  Toptag(recoTopsubjets,irecoTopjet->mass(),minTopMass,maxTopMass,minMinMass);
    hltTopPass_ = Toptag(hltTopsubjets,ihltTopjet->mass(),minTopMass,maxTopMass,minMinMass);
    recoWZPass_ = WZtag(recoWZsubjets,irecoWZjet->mass(),minWZMass,maxWZMass,MassDropCut); 
    hltWZPass_ = WZtag(hltWZsubjets,ihltWZjet->mass(),minWZMass,maxWZMass,MassDropCut);
    
    
    hltMinMass_ = minmass(hltTopsubjets);
    recoMinMass_ =minmass(recoTopsubjets);
    hltMassDrop_ = massdrop(ihltWZjet->mass(),hltWZsubjets);  
    recoMassDrop_ = massdrop(irecoWZjet->mass(),recoWZsubjets);
    
    NhltTopSubjets_ = hltTopsubjets.size();  
    NrecoTopSubjets_ = recoTopsubjets.size();
    NhltWZSubjets_ = hltWZsubjets.size();
    NrecoWZSubjets_ = recoWZsubjets.size();
    NhltTrimSubjets_ = hltTrimsubjets.size();
    
    hltTopPT_ = ihltTopjet->pt();
    recoTopPT_ = irecoTopjet->pt();
    hltWZPT_ = ihltWZjet->pt();
    recoWZPT_ = irecoWZjet->pt();
    hltTrimPT_ = ihltTrimjet->pt();
    //hltCA8PT_ = ihltCA8jet->pt();

    hltTopEta_=ihltTopjet->eta();
    hltTopPhi_=ihltTopjet->phi();
    recoTopEta_=irecoTopjet->eta();
    recoTopPhi_=irecoTopjet->phi();
    hltWZEta_=ihltWZjet->eta();
    hltWZPhi_=ihltWZjet->phi();
    recoWZEta_=irecoWZjet->eta();
    recoWZPhi_=irecoWZjet->phi();
    hltTrimEta_=ihltTrimjet->eta();
    hltTrimPhi_=ihltTrimjet->phi();

  }
    
  //fill the tree
  tree->Fill();
  
}    

int Toptag(reco::Jet::Constituents subjets,float jetmass, double minTM, double maxTM, double minMM){
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
      if(minmass>minMM) return 1;
      else return 0;
  }
  else return 0;
}



int WZtag(reco::Jet::Constituents subjets,float jetmass, double minWM, double maxWM, double mdc){
  //instantiate min mass
  float massdrop = 99999.;

 //find the mass drop first cut on number of subjets >=3
    if(subjets.size()== 2){
      // Take the highest 3 pt subjets for cuts
      sort ( subjets.begin(), subjets.end(), GreaterByPtCandPtrUser() );
      //get the highest pt subjet
      reco::Jet::Constituent icandJet = subjets[0];

      reco::Candidate::LorentzVector isubJet = icandJet->p4();
      double imass = isubJet.mass();
      massdrop=imass/jetmass;
    }


    if(subjets.size()==2 && jetmass>minWM && jetmass<maxWM){
      if(massdrop<mdc) return 1;
      else return 0;
    }
    else return 0;
}

float minmass(reco::Jet::Constituents subjets){
  //instantiate min mass
  float minmass = 99999.;

  //find the mass drop first cut on number of subjets >=3
  if(subjets.size()>= 3){
    // Take the highest 3 pt subjets for cuts
    sort ( subjets.begin(), subjets.end(), GreaterByPtCandPtrUser() );
    // Now look at the subjets that were formed
    for(unsigned int isub=0; isub<subjets.size(); ++isub){
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
  
  return minmass;
}


float massdrop(float jetmass,reco::Jet::Constituents subjets){
  //instantiate min mass
  float massdrop = 99999.;

 //find the mass drop first cut on number of subjets >=3
    if(subjets.size()== 2){
      // Take the highest 3 pt subjets for cuts
      sort ( subjets.begin(), subjets.end(), GreaterByPtCandPtrUser() );
      //get the highest pt subjet
      reco::Jet::Constituent icandJet = subjets[0];

      reco::Candidate::LorentzVector isubJet = icandJet->p4();
      double imass = isubJet.mass();
      massdrop=imass/jetmass;
    }

    return massdrop;
}


//define this as a plug-in
DEFINE_FWK_MODULE(hltTreeMaker);
