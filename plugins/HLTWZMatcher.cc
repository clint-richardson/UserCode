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


bool wztag(reco::Jet::Constituents subjets,float jetmass, TH1F* sjhist, TH1F* mdhist1, TH1F* mdhist2,int mode);

struct GreaterByPtCandPtrUser {
  bool operator()( const edm::Ptr<reco::Candidate> & t1, const edm::Ptr<reco::Candidate> & t2 ) const {
    return t1->pt() > t2->pt();
  }
};

class HLTWZMatcher : public edm::EDProducer {
  
public:

  HLTWZMatcher(const edm::ParameterSet& Pset);
  virtual ~HLTWZMatcher(){}

  
private:
  virtual void BeginJob();
  virtual void produce(edm::Event & event, const edm::EventSetup & EventSetup);
  virtual void EndJob(){};




  TH1F* pthist;
  TH1F* deltaRhist;
  TH1F* recomasshist;
  TH1F* hltmasshist;
  TH1F* ratiomasshist_nc;
  TH1F* ratiomasshist_c;
  TH1F* ratiomasshist_f;
  TH1F* recomdhist_nc;
  TH1F* hltmdhist_nc;
  TH1F* recomdhist_c;
  TH1F* hltmdhist_c;
  TH1F* recosjhist;
  TH1F* hltsjhist;
  TH1F* ratiosjhist_nc;
  TH1F* ratiosjhist_c;

  std::string histname_;
  edm::InputTag   OfflinejetColl;
  edm::InputTag  OnlinejetColl;
};

HLTWZMatcher::HLTWZMatcher(const edm::ParameterSet& Pset){
  //load options
  //  if (Pset.exists("pvCollection")) pvCollection_it = Pset.getParameter<edm::InputTag>("pvCollection");
  // else                              pvCollection_it = edm::InputTag("goodOfflinePrimaryVertices");
  
  if (Pset.exists("OfflinejetCollection")) OfflinejetColl = Pset.getParameter<edm::InputTag>("OfflinejetCollection");
  else                              OfflinejetColl = edm::InputTag("ca8PFJetsCHS");  
  if (Pset.exists("OnlinejetCollection")) OnlinejetColl = Pset.getParameter<edm::InputTag>("OnlinejetCollection");
  else                              OnlinejetColl = edm::InputTag("hltCA8TopJets");  
  histname_ = Pset.getParameter<std::string>("histname");
  //attempt to add histogram
  edm::Service<TFileService> fs;
  float ptbins[25] = { 0,160,180,200,220,240,260,280,300,320,340,360,380,400,420,440,460,480,500,550,600,650,700,800,1000};
  pthist = fs->make<TH1F>(histname_.c_str(), "JetPT", 24, ptbins);
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
  recomdhist_nc = fs->make<TH1F>("recomdhist_nc", "reco massdrop no cuts",50,0,1);
  hltmdhist_nc = fs->make<TH1F>("hltmdhist_nc", "hlt massdrop no cuts", 50,0,1);
  recomdhist_c = fs->make<TH1F>("recomdhist_c", "reco massdrop cuts",50,0,1);
  hltmdhist_c = fs->make<TH1F>("hltmdhist_c", "hlt massdrop cuts", 50,0,1);
}

void HLTWZMatcher::BeginJob(){
}

void HLTWZMatcher::produce(edm::Event& iEvent,const edm::EventSetup& iEventSetup){
  //read in jets
  edm::Handle<reco::BasicJetCollection> recoJets;
  iEvent.getByLabel(OfflinejetColl, recoJets);
  edm::Handle<reco::BasicJetCollection> hltJets;
  iEvent.getByLabel(OnlinejetColl, hltJets);
  std::auto_ptr<reco::BasicJetCollection > hltJetColl( new reco::BasicJetCollection (*hltJets));
  reco::BasicJetCollection::const_iterator irecojet = recoJets->begin();
  int trial=0;
  int njets=0;
  for(; irecojet!=recoJets->end();irecojet++){
    trial++;
    njets++;
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
    cout<<"mindR is: "<<mindR<<endl;
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
    bool reco_pass = wztag(recosubjets,recomass,recosjhist,recomdhist_nc,recomdhist_c,0);
    //if(reco_pass) //cout<<"got reco pass"<<endl;
      //if(!reco_pass) //cout<<"got reco fail"<<endl;
    bool hlt_pass = wztag(hltsubjets,hltmass,hltsjhist,hltmdhist_nc,hltmdhist_c,1);
    if(hlt_pass){
      //cout<<"got hlt pass"<<endl;
      ratiosjhist_c->Fill(hltsubjets.size()/recosubjets.size());
      ratiomasshist_c->Fill(hltmass/recomass);
    }
    if(!hlt_pass){
      //cout<<"got hlt fail"<<endl;
      ratiomasshist_f->Fill(hltmass/recomass);
    }
    //make the both jets pass wztag for hist to be filled
    //if(reco_pass && !hlt_pass) //cout<<"I'm handling things fine"<<endl;
    if(reco_pass && hlt_pass){
      float pt = irecojet->pt();
      //cout<<"Jet pt is: "<<pt<<endl;
      pthist->Fill(pt);
    }
  
  }
  //cout<<"Number of offline jets in event is: "<<njets<<endl;
}  

bool wztag(reco::Jet::Constituents subjets,float jetmass, TH1F* sjhist,TH1F* mdhist1,TH1F* mdhist2,int mode){
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

    sjhist->Fill(subjets.size());
    mdhist1->Fill(massdrop);
    if(mode==0){
      cout<<"number of subjets: "<<subjets.size()<<"; jet mass: "<<jetmass<<"; mass drop: "<<massdrop<<endl;
    }
    if(subjets.size()==2 && jetmass>60 && jetmass<130){
      mdhist2->Fill(massdrop);
      if(massdrop<0.4) return true;
      else return false;
    }
    else return false;
}

//define this as a plug-in
DEFINE_FWK_MODULE(HLTWZMatcher);
