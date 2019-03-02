import FWCore.ParameterSet.Config as cms


process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("FULLSIMReader.FULLSIMReader.FULLSIMReader_cff")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'file:/eos/cms/store/user/pablom/GluGluHH/output_25.root'
    )
)

process.p = cms.Path(process.fullsimreader)


