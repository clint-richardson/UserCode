#ifndef UserCode_hltJetMassAnalyzer_h
#define UserCode_hltJetMassAnalyzer_h

// system include files
#include <memory>
#include <vector>
#include <sstream>
#include <string>
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/JetReco/interface/CATopJetTagInfo.h"


// user include files
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/Math/interface/deltaR.h"

#include "DataFormats/Common/interface/View.h"
#include "DataFormats/JetReco/interface/BasicJet.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "PhysicsTools/CandUtils/interface/AddFourMomenta.h"
#include "DataFormats/Candidate/interface/CandMatchMap.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEventWithRefs.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerTypeDefs.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include <Math/VectorUtil.h>
#include <TH1.h>
#include <TH2.h>
#include <TTree.h>



class hltJetMassAnalyzer : public edm::EDAnalyzer{
 public:

  explicit hltJetMassAnalyzer(const edm::ParameterSet& PSet);
  ~hltJetMassAnalyzer(){};

 private:
  TH1F* masshist;
  TH1F* pthist;
  virtual void beginJob();
  virtual void analyze(const edm::Event & event, const edm::EventSetup & EventSetup);
  virtual void endJob(){};
  std::string filter;
};


#endif
