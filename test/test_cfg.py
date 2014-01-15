import FWCore.ParameterSet.Config as cms

process = cms.Process("Tagging")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32 ( -1 ) )

#load in source file
process.source = cms.Source( "PoolSource",
                            fileNames = cms.untracked.vstring('file:./outputA.root'),
                            inputCommands = cms.untracked.vstring('keep *'),
                            )


#run the producer
process.Tagger = cms.EDProducer("MCTagger",
                                genParticles = cms.InputTag("genParticles"),
                                )

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string('histodemo.root')
                                   )

process.out = cms.OutputModule("PoolOutputModule",
                                       fileName = cms.untracked.string('myOutputFile.root') ,
                                       outputCommands = cms.untracked.vstring("keep *"),
                                       )

process.p = cms.Path(process.Tagger)

process.e = cms.EndPath(process.out)
