import FWCore.ParameterSet.Config as cms

process = cms.Process("JetMassHist")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32 ( -1 ) )

#load in source file
process.source = cms.Source( "PoolSource",
                            fileNames = cms.untracked.vstring('file:./outputA.root'),
#                           fileNames = cms.untracked.vstring('file:./ZprimeTottbar_RECO_withTriggerBits_v2.root'),
                            inputCommands = cms.untracked.vstring('keep *'),
                            )


#run the analyzer
process.JetMassHist = cms.EDAnalyzer("hltJetMassAnalyzer",
                                     hltFilter = cms.string("hltCATopTagFilter"),
                                )

#run the file service to save histograms
process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string('hlt_topjet_hists.root')
)

process.p = cms.Path(process.JetMassHist)
