import FWCore.ParameterSet.Config as cms

topmcmatch = cms.EDFilter("TopTagMCMatcher",
	pTmin = cms.untracked.double(0.),
									  )
