import FWCore.ParameterSet.Config as cms

process = cms.Process("TEST")

## Load Gen Filter
#process.load("UserCode.HLTJetTagging.TopTagMCMatcher_cfi")
process.load("UserCode.HLTJetTagging.HadronTagger_cfi")

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))

process.source = cms.Source("PoolSource",
							fileNames = cms.untracked.vstring(
	'file:./outputA.root'
        
	#'file:///uscms_data/d2/avetisya/HLT/CMSSW_5_3_9/src/reco_DIGI_L1_DIGI2RAW.root'
	)
)

# Output definition
process.output = cms.OutputModule("PoolOutputModule",
								  #fastCloning = cms.untracked.bool(False),
								  outputCommands = process.RECOSIMEventContent.outputCommands,
								  fileName = cms.untracked.string('TopTagMCMatcher.root')
								  )

process.p = cms.Path(process.hadrontag)
process.outpath = cms.EndPath(process.output)
