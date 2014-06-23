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
int toptag(reco::Jet::Constituents subjets,float jetmass, double minTM, double maxTM, double minMM);
int wztag(reco::Jet::Constituents subjets,float jetmass, double minWM, double maxWM, double mdc);

float MinMass(reco::Jet::Constituents subjets);
float MassDrop(float jetmass,reco::Jet::Constituents subjets);


struct GreaterByPtCandPtrUser {
  bool operator()( const edm::Ptr<reco::Candidate> & t1, const edm::Ptr<reco::Candidate> & t2 ) const {
    return t1->pt() > t2->pt();
  }
};



//plugin declaration
class hltEDMTreeMaker : public edm::EDProducer {
  
public:

  hltEDMTreeMaker(const edm::ParameterSet& Pset);
  virtual ~hltEDMTreeMaker(){}

  
private:
  virtual void BeginJob();
  virtual void produce(edm::Event & event, const edm::EventSetup & EventSetup);
  virtual void EndJob(){};


  edm::InputTag   OfflineTopjetColl;
  edm::InputTag  OnlineTopjetColl;

  edm::InputTag   OfflineWZjetColl;
  edm::InputTag  OnlineWZjetColl;


  edm::InputTag  OnlineTrimjetColl;

  double minTopMass;
  double maxTopMass;
  double minMinMass;

  double minWZMass;
  double maxWZMass;
  double MassDropCut;



};

hltEDMTreeMaker::hltEDMTreeMaker(const edm::ParameterSet& Pset){
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


  //top jet variables first
  produces<std::vector<float> >( "hltTopJetMass" ).setBranchAlias( "hltTopJetMass" );
  produces<std::vector<float> >( "recoTopJetMass" ).setBranchAlias( "recoTopJetMass" );
  produces<std::vector<float> >( "hltTopPT" ).setBranchAlias( "hltTopPT" );
  produces<std::vector<float> >( "recoTopPT" ).setBranchAlias( "recoTopPT" );
  produces<std::vector<int> >( "NhltTopSubjets" ).setBranchAlias( "NhltTopSubjets" );
  produces<std::vector<int> >( "NrecoTopSubjets" ).setBranchAlias( "NrecoTopSubjets" );
  produces<std::vector<int> >( "hltTopPass" ).setBranchAlias( "hltTopPass" );
  produces<std::vector<int> >( "recoTopPass" ).setBranchAlias( "recoTopPass" );
  produces<std::vector<float> >("hltMinMass").setBranchAlias( "hltMinMass");
  produces<std::vector<float> >("recoMinMass").setBranchAlias( "recoMinMass");

  //now wz jet variables
  produces<std::vector<float> >( "hltWZJetMass" ).setBranchAlias( "hltWZJetMass" );
  produces<std::vector<float> >( "recoWZJetMass" ).setBranchAlias( "recoWZJetMass" );
  produces<std::vector<float> >( "hltWZPT" ).setBranchAlias( "hltWZPT" );
  produces<std::vector<float> >( "recoWZPT" ).setBranchAlias( "recoWZPT" );
  produces<std::vector<int> >( "NhltWZSubjets" ).setBranchAlias( "NhltWZSubjets" );
  produces<std::vector<int> >( "NrecoWZSubjets" ).setBranchAlias( "NrecoWZSubjets" );
  produces<std::vector<int> >( "hltWZPass" ).setBranchAlias( "hltWZPass" );
  produces<std::vector<int> >( "recoWZPass" ).setBranchAlias( "recoWZPass" );
  produces<std::vector<float> >("hltMassDrop").setBranchAlias( "hltMassDrop");
  produces<std::vector<float> >("recoMassDrop").setBranchAlias( "recoMassDrop");


  //now trimmed jet variables
  produces<std::vector<float> >( "hltTrimJetMass" ).setBranchAlias( "hltTrimJetMass" );
  produces<std::vector<float> >( "recoTrimJetMass" ).setBranchAlias( "recoTrimJetMass" );
  produces<std::vector<float> >( "hltTrimPT" ).setBranchAlias( "hltTrimPT" );
  produces<std::vector<int> >( "NhltTrimSubjets" ).setBranchAlias( "NhltTrimSubjets" );
  
  /*note that I'm not including a trim pass since we want to investigate moving this
  variable around anyway and it's easy enough to plot things with cuts using trees
  I am including hlt top and wz pass since those are more fixed cuts right now though
  with the stored information one should be able to redefine them as needed afterwards
  */
  
}

void hltEDMTreeMaker::BeginJob(){
}

void hltEDMTreeMaker::produce(edm::Event& iEvent,const edm::EventSetup& iEventSetup){

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
  
  cout<<"got jet collections"<<endl;

  //get the higest pt jet
  reco::BasicJetCollection::const_iterator irecoTopjet = recoTopJets->begin();
  for(reco::BasicJetCollection::const_iterator dummyjet=recoTopJets->begin(); dummyjet!=recoTopJets->end();dummyjet++){

    if(dummyjet->pt()>irecoTopjet->pt()){
      irecoTopjet=dummyjet;
    }

  }

  


  //match the reco jet with hlt jet by deltaR and get iterator
  float eta1 = irecoTopjet->eta();
  float phi1 = irecoTopjet->phi();
  cout<<"reco top eta: "<<eta1<<" and phi: "<<phi1<<endl;
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
  cout<<"found matched hlttopjet"<<endl;

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

  cout<<"found matched hlt wz jet"<<endl;

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

  cout<<"found matche reco wz jet"<<endl;

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

  cout<<"found matched hlt trim jet"<<endl;

  /*now I should have 5 jets:
   *ihltTopjet
   *irecoTopjet
   *ihltWZjet
   *recoWZjet
   *ihltTrimjet
   */
  
  
  //get jet masses;
  std::auto_ptr<std::vector<float> > recotopmass(new vector<float>);
  cout<<"made std auto ptr and the jet mass is"<<irecoTopjet->mass()<<endl;
  recotopmass->push_back(irecoTopjet->mass());
  cout<<"reco top mass "<<irecoTopjet->mass()<<endl;
  std::auto_ptr<std::vector<float> > hlttopmass(new vector<float>);
  hlttopmass->push_back( ihltTopjet->mass());
  cout<<"hlt top mass "<<ihltTopjet->mass();
  std::auto_ptr<std::vector<float> > hltwzmass(new vector<float>);
  hltwzmass->push_back(ihltWZjet->mass());
  std::auto_ptr<std::vector<float> > recowzmass(new vector<float>);
  recowzmass->push_back(irecoWZjet->mass());
  std::auto_ptr<std::vector<float> > hlttrimmass(new vector<float>);
  hlttrimmass->push_back( ihltTrimjet->mass());

  cout<<"got jet masses"<<endl;

  //get the subjets
  reco::Jet::Constituents recoTopsubjets = irecoTopjet->getJetConstituents();
  reco::Jet::Constituents hltTopsubjets = ihltTopjet->getJetConstituents();
  reco::Jet::Constituents recoWZsubjets = irecoWZjet->getJetConstituents();
  reco::Jet::Constituents hltWZsubjets = ihltWZjet->getJetConstituents();
  reco::Jet::Constituents hltTrimsubjets = ihltTrimjet->getJetConstituents();

  cout<<"got subjets"<<endl;

  std::auto_ptr<std::vector<int> > recotoppass(new vector<int>);
  recotoppass->push_back(toptag(recoTopsubjets,irecoTopjet->mass(),minTopMass,maxTopMass,minMinMass));
  std::auto_ptr<std::vector<int> > hlttoppass(new vector<int>);
  hlttoppass->push_back(toptag(hltTopsubjets,ihltTopjet->mass(),minTopMass,maxTopMass,minMinMass));

  std::auto_ptr<std::vector<int> > recowzpass(new vector<int>);
  recowzpass->push_back(wztag(recoWZsubjets,irecoWZjet->mass(),minWZMass,maxWZMass,MassDropCut));
  std::auto_ptr<std::vector<int> > hltwzpass(new vector<int>);
  hltwzpass->push_back(wztag(hltWZsubjets,ihltWZjet->mass(),minWZMass,maxWZMass,MassDropCut));

  std::auto_ptr<std::vector<float> > hltminmass(new vector<float>);
  hltminmass->push_back(MinMass(hltTopsubjets));
  std::auto_ptr<std::vector<float> > recominmass(new vector<float>);
  recominmass->push_back(MinMass(recoTopsubjets));
  std::auto_ptr<std::vector<float> > hltmassdrop(new vector<float>);
  hltmassdrop->push_back(MassDrop(ihltWZjet->mass(),hltWZsubjets));  
  std::auto_ptr<std::vector<float> > recomassdrop(new vector<float>);
  recomassdrop->push_back(MassDrop(irecoWZjet->mass(),recoWZsubjets));

  std::auto_ptr<std::vector<int> > Nhlttopsubjets(new vector<int>);
  Nhlttopsubjets->push_back( hltTopsubjets.size());
  std::auto_ptr<std::vector<int> > Nrecotopsubjets(new vector<int>);
  Nrecotopsubjets->push_back(recoTopsubjets.size());
  std::auto_ptr<std::vector<int> > Nhltwzsubjets(new vector<int>);
  Nhltwzsubjets->push_back(hltWZsubjets.size());
  std::auto_ptr<std::vector<int> > Nrecowzsubjets(new vector<int>);
  Nrecowzsubjets->push_back(recoWZsubjets.size());

  std::auto_ptr<std::vector<int> > Nhlttrimsubjets(new vector<int>);
  Nhlttrimsubjets->push_back(hltTrimsubjets.size());

  std::auto_ptr<std::vector<float> > hlttopPT(new vector<float>);
  hlttopPT->push_back(ihltTopjet->pt());
  std::auto_ptr<std::vector<float> > recotopPT(new vector<float>);
  recotopPT->push_back( irecoTopjet->pt());
  std::auto_ptr<std::vector<float> > hltwzPT(new vector<float>);
  hltwzPT->push_back(ihltWZjet->pt());
  std::auto_ptr<std::vector<float> > recowzPT(new vector<float>);
  recowzPT->push_back(irecoWZjet->pt());
  std::auto_ptr<std::vector<float> > hlttrimPT(new vector<float>);
  hlttrimPT->push_back(ihltTrimjet->pt());

  cout<<"got the variables"<<endl;

  //put pt info in event
  iEvent.put(hlttopPT, "hltTopPT");
  iEvent.put(recotopPT, "recoTopPT");
  iEvent.put(hltwzPT, "hltWZPT");
  iEvent.put(recowzPT, "recoWZPT");
  iEvent.put(hlttrimPT, "hltTrimPT");
  //put n subjets in event
  iEvent.put(Nhlttopsubjets,"NhltTopSubjets");
  iEvent.put(Nrecotopsubjets,"NrecoTopSubjets");
  iEvent.put(Nhltwzsubjets,"NhltWZSubjets");
  iEvent.put(Nrecowzsubjets,"NrecoWZSubjets");
  iEvent.put(Nhlttrimsubjets, "NhltTrimSubjets");
  //put jet masses in event
  iEvent.put(hlttopmass, "hltTopJetMass");
  iEvent.put(recotopmass,"recoTopJetMass");
  iEvent.put(hltwzmass, "hltWZJetMass");
  iEvent.put(recowzmass,"recoWZJetMass");
  iEvent.put(hlttrimmass,"hltTrimJetMass");
  //put the min masses in event
  iEvent.put(hltminmass,"hltMinMass");
  iEvent.put(recominmass,"recoMinMass");
  //put the mass drops in event
  iEvent.put(hltmassdrop,"hltMassDrop");
  iEvent.put(recomassdrop,"recoMassDrop");
  //put the tagging decisions in event
  iEvent.put(hlttoppass,"hltTopPass");
  iEvent.put(recotoppass,"recoTopPass");
  iEvent.put(hltwzpass,"hltWZPass");
  iEvent.put(recowzpass,"recoWZPass");

  cout<<"put variables in the event"<<endl;
  
}    

int toptag(reco::Jet::Constituents subjets,float jetmass, double minTM, double maxTM, double minMM){
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



int wztag(reco::Jet::Constituents subjets,float jetmass, double minWM, double maxWM, double mdc){
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

float MinMass(reco::Jet::Constituents subjets){
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


float MassDrop(float jetmass,reco::Jet::Constituents subjets){
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
DEFINE_FWK_MODULE(hltEDMTreeMaker);
