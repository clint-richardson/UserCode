#ifndef HadronTagger_h
#define HadronTagger_h
// -*- C++ -*-
//
// Package:    HadronTagger
// Class:      HadronTagger
// 
/*
 Description: edmFilter to match to top tag filter

 Implementation:
                    
*/
//
// Original Author:  Michael Maes
//         Created:  Wed Dec  3 12:07:13 CET 2009
// $Id: HadronTagger.h,v 1.4 2010/07/21 04:23:24 wmtan Exp $
//
// Adapted for TopTagMCMatching: Dylan Rankin


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "DataFormats/Common/interface/View.h"


#include <map>
#include <vector>

#include <string>

//
// class declaration
//

class HadronTagger : public edm::EDFilter {
   public:
      explicit HadronTagger(const edm::ParameterSet&);
      ~HadronTagger();

      virtual bool filter(edm::Event&, const edm::EventSetup&);

   private:

      std::string label_;

     

      // ----------member data ---------------------------
		
};
#endif
