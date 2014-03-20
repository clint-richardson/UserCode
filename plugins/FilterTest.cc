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
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEventWithRefs.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include <Math/VectorUtil.h>
#include "TH1F.h"
#include <TH2.h>
#include <TTree.h>



using std::cout;
using std::endl;

class FilterTest : public edm::EDProducer {
  
public:

  FilterTest(const edm::ParameterSet& Pset);
  virtual ~FilterTest(){}

  
private:
  virtual void BeginJob();
  virtual void produce(edm::Event & event, const edm::EventSetup & EventSetup);
  virtual void EndJob(){};
  

  HLTConfigProvider         hltConfig_;

};

FilterTest::FilterTest(const edm::ParameterSet& Pset){
  

}

void FilterTest::BeginJob(){
}

void FilterTest::produce(edm::Event& iEvent,const edm::EventSetup& iEventSetup){

  //Trigger Filters/Objects
  edm::Handle<trigger::TriggerEvent> aodTriggerEvent;
  iEvent.getByLabel("hltTriggerSummaryAOD",aodTriggerEvent);


  //get trigger objects
  trigger::TriggerObjectCollection allObjects = aodTriggerEvent->getObjects();

  //  for(int j=0; j<aodTriggerEvent->sizeFilters();j++){
  //cout<<"filter is: "<<aodTriggerEvent->filterTag(j).label()<<endl;
  //}

  //check to see if trigger was fired
  for(int j = 0; j<aodTriggerEvent->sizeFilters(); j++){
    if(aodTriggerEvent->filterTag(j).label()!="hltCATopTagFilter") continue;
    //get keys
    trigger::Keys keys = aodTriggerEvent->filterKeys(j);
    for(size_t n=0; n<keys.size();n++){
      float pt = allObjects[keys[n]].pt();
      cout<<"Transverse momentum is: "<<pt<<endl;
      cout<<"Filter is: "<<aodTriggerEvent->filterTag(j).label()<<endl;
    }
  }
  /*
  //2nd attempt at accessing trigger results following strategy of HLTInfo.cc
  edm::Handle<edm::TriggerResults> hltresults;
  edm::InputTag trigRes = edm::InputTag("TriggerResults");
  iEvent.getByLabel(trigRes,hltresults);
  if(hltresults.isValid()){
    int nTrig = hltresults->size();
    edm::TriggerNames const& triggerNames = iEvent.triggerNames(*hltresults);
    for( int iTrig =0; iTrig<nTrig; iTrig++){
      std::string trigName = triggerNames.triggerName(iTrig);
      cout<<"Trigger name is: "<<trigName<<endl;
    }
    }*/
  
  HLTConfigProvider hltConfig;
  for(unsigned int iHLT=0; iHLT< hltConfig.size(); iHLT++){
    cout<<"Path Name is: "<<hltConfig.triggerName(iHLT);
  }
  

}


//define this as a plug-in
DEFINE_FWK_MODULE(FilterTest);
