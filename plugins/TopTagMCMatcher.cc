#include <vector>
#include <iostream>
#include <cmath>
#include <sstream>
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "../interface/TopTagMCMatcher.h"


//
// constants, enums and typedefs
//

//
// static data member definitions
//
using namespace edm;
using namespace reco;
using namespace std;

//
// constructors and destructor
//

TopTagMCMatcher::TopTagMCMatcher(const edm::ParameterSet& iConfig) : 
  label_(iConfig.getUntrackedParameter("moduleLabel",std::string("genParticles"))),
  ptmin_(iConfig.getUntrackedParameter("pTmin",0.))
{
}


TopTagMCMatcher::~TopTagMCMatcher()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
TopTagMCMatcher::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  bool accept=false;

  bool tWb = false;
  bool tbarWb = false;

  Handle<reco::GenParticleCollection> genParticles ;
  iEvent.getByLabel( label_, genParticles ) ;

  for ( size_t i= 0;i< genParticles->size(); ++i ) {

    const GenParticle & p = (*genParticles)[i];

    int pdgID = p.pdgId();
    
    int status = p.status();

    double pT = p.pt();

    bool bjet = false;
    double beta = 5;
    double bphi = 3;
    bool W1pjet = false;
    double W1peta = 5;
    double W1pphi = 3;
    bool W2pjet = false;
    double W2peta = 5;
    double W2pphi = 3;
    int Wpjetsnum = 0;
    bool bbjet = false;
    double bbeta = 5;
    double bbphi = 3;
    bool W1mjet = false;
    double W1meta = 5;
    double W1mphi = 3;
    bool W2mjet = false;
    double W2meta = 5;
    double W2mphi = 3;
    int Wmjetsnum = 0;


    if (pT > ptmin_ && status == 3 && pdgID == 6) {
      //std::cout << "found top" << std::endl;

      for (size_t j = 0; j<p.numberOfDaughters(); j++){
        const Candidate * d = p.daughter(j);
        int dauId = d->pdgId();
        //std::cout << "top decayed to a " << dauId << std::endl;
        if (dauId == 5) {
          //std::cout << "found bjet" << std::endl;
          bjet = true;
          beta = d->eta();
          bphi = d->phi();
        }
        if (abs(dauId) == 24) {
          //std::cout << "found W" << std::endl;


          for (size_t Wj = 0; Wj<d->numberOfDaughters(); Wj++){
            const Candidate * Wd = d->daughter(Wj);
            int dauId = Wd->pdgId();
            if (abs(dauId) < 10 && Wpjetsnum == 0) {
              //std::cout << "found first hardonic W decay" << std::endl;
              //std::cout << "particle 1 from W+: " << dauId << std::endl;
              W1pjet = true;
              W1peta = Wd->eta();
              W1pphi = Wd->phi();
              Wpjetsnum = 1;
            }
            if (abs(dauId) < 10 && Wpjetsnum == 1) {
              //std::cout << "found second hardonic W decay" << std::endl;
              //std::cout << "particle 2 from W+: " << dauId << std::endl;
              W2pjet = true;
              W2peta = Wd->eta();
              W2pphi = Wd->phi();
            }
          }
        }
      }
    }// daughters 

    double bW1pDR = deltaR (beta, bphi, W1peta, W1pphi);
    double bW2pDR = deltaR (beta, bphi, W2peta, W2pphi);
    double W1pW2pDR = deltaR (W1peta, W1pphi, W2peta, W2pphi);
    double tjDR = TMath::Max(bW1pDR,bW2pDR);
    tjDR = TMath::Max(tjDR,W1pW2pDR);

    if (bjet && W1pjet && W2pjet && tjDR < 0.8) {
    //if (bjet && W1pjet && W2pjet) {
      //std::cout << "maxDeltaR = " << tjDR << endl;
      //std::cout << pT << " pT boosted top decaying hadronically found in:" << std::endl;
      //std::cout << "Run " << iEvent.id().run() << std::endl;
      //std::cout << "Lumi " << iEvent.id().luminosityBlock() << std::endl;
      //std::cout << "Event " << iEvent.id().event() << std::endl;
      tWb = true;
    }
    if (pT > ptmin_ && status == 3 && pdgID == -6) {
      //std::cout << "found anti-top" << std::endl;

      for (size_t j = 0; j<p.numberOfDaughters(); j++){
        const Candidate * d = p.daughter(j);
        int dauId = d->pdgId();
        //std::cout << "anti-top decayed to a " << dauId << std::endl;
        if (dauId == -5) {
          //std::cout << "found bbarjet" << std::endl;
          bbjet = true;
          bbeta = d->eta();
          bbphi = d->phi();
        }
        if (abs(dauId) == 24) {
          //std::cout << "found W" << std::endl;

          for (size_t Wj = 0; Wj<d->numberOfDaughters(); Wj++){
            const Candidate * Wd = d->daughter(Wj);
            int dauId = Wd->pdgId();
            if (abs(dauId) < 10 && Wmjetsnum == 0) {
              //std::cout << "found first hardonic W decay" << std::endl;
              //std::cout << "particle 1 from W-: " << dauId << std::endl;
              W1mjet = true;
              W1meta = Wd->eta();
              W1mphi = Wd->phi();
              Wmjetsnum = 1;
            }
            if (abs(dauId) < 10 && Wmjetsnum == 1) {
              //std::cout << "found second hardonic W decay" << std::endl;
              //std::cout << "particle 2 from W-: " << dauId << std::endl;
              W2mjet = true;
              W2meta = Wd->eta();
              W2mphi = Wd->phi();
            }
          }
        }
      }
    }// daughters 

    double bbW1mDR = deltaR (bbeta, bbphi, W1meta, W1mphi);
    double bbW2mDR = deltaR (bbeta, bbphi, W2meta, W2mphi);
    double W1mW2mDR = deltaR (W1meta, W1mphi, W2meta, W2mphi);
    double tbjDR = TMath::Max(bbW1mDR,bbW2mDR);
    tbjDR = TMath::Max(tbjDR,W1mW2mDR);

    if (bbjet && W1mjet && W2mjet && tbjDR < 0.8) {
    //if (bbjet && W1mjet && W2mjet) {
      //std::cout << "maxDeltaR = " << tbjDR << endl;
      //std::cout << pT << " pT boosted anti-top decaying hadronically found in:" << std::endl;
      //std::cout << "Run " << iEvent.id().run() << std::endl;
      //std::cout << "Lumi " << iEvent.id().luminosityBlock() << std::endl;
      //std::cout << "Event " << iEvent.id().event() << std::endl;
      tbarWb = true;
    }
  }
  if (tWb && tbarWb) {
    std::cout << "Found hadronic top and anti-top event" << std::endl;
  }
  else if (tWb && !tbarWb) {
    std::cout << "Found hadronic top only event" << std::endl;
  }
  else if (!tWb && tbarWb) {
    std::cout << "Found hadronic anti-top only event" << std::endl;
  }
  else {
    std::cout << "Found no hadronic top event" << std::endl;
  }

  if (tWb || tbarWb) {
    accept = true;
  }
 
  return accept;
}

DEFINE_FWK_MODULE(TopTagMCMatcher);
