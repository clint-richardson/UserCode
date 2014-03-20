#include "../interface/hltJetMassAnalyzer.h"



hltJetMassAnalyzer::hltJetMassAnalyzer(const edm::ParameterSet & PSet){
  
  filter = PSet.getParameter<std::string>("hltFilter");

  //initialize the file service
  edm::Service<TFileService> fs;

  //make the histogram titles
  std::string masstitle,ptTitle;
  if(filter=="hltCATopTagFilter"){
    masstitle="Top Tagged Jet Mass Distribution; Mass (GeV); N_{Events}";
    ptTitle="Top Tagged Jet Pt Distribution; Pt (GeV); N_{Events}";
  }
  else{
    masstitle="W Tagged Jet Mass Distribution; Mass (GeV); N_{Events}";
    ptTitle="W Tagged Jet Pt Distribution; Pt (GeV); N_{Events}";
  }

  //initialize the histograms
  masshist = fs->make<TH1F>(masstitle.c_str(),"masshist",100,0.,300);
  pthist = fs->make<TH1F>(ptTitle.c_str(),"pthist",100,0.,2000);

};

void hltJetMassAnalyzer::beginJob(){};

void hltJetMassAnalyzer::analyze(const edm::Event &Event, const edm::EventSetup & EventSetup){

  //get the trigger summary
  edm::Handle<trigger::TriggerEvent> aodTriggerEvent;
  Event.getByLabel("hltTriggerSummaryAOD",aodTriggerEvent);

  //get the trigger objects
  trigger::TriggerObjectCollection allObjects = aodTriggerEvent->getObjects();

  for(int j=0;j<aodTriggerEvent->sizeFilters();j++){
    //don't go through unless it's the filter we want
    if(aodTriggerEvent->filterTag(j).label()!=filter) continue;
    //get trigger keys
    trigger::Keys keys = aodTriggerEvent->filterKeys(j);

    for(size_t n=0; n<keys.size();n++){
      float pt = allObjects[keys[n]].pt();
      float mass = allObjects[keys[n]].mass();
      masshist->Fill(mass);
      pthist->Fill(pt);
      //      std::cout<<"Mass is: "<<mass<<std::endl;
    }
  }
}

//define as framework plugin
DEFINE_FWK_MODULE(hltJetMassAnalyzer);
