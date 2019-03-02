import FWCore.ParameterSet.Config as cms


fullsimreader = cms.EDAnalyzer('FULLSIMReader',
    outputFileName = cms.string('output.root'),
    bDiscriminatorName = cms.string('pfCombinedInclusiveSecondaryVertexV2BJetTags'),
    ElectronCollection = cms.InputTag("slimmedElectrons"),
    PuppiJetCollection = cms.InputTag("slimmedJetsPuppi"),
    MuonCollection = cms.InputTag("slimmedMuons"),
    PhotonCollection = cms.InputTag("slimmedPhotons"),
    IsoTrackCollection = cms.InputTag("isolatedTracks"),
    PrimaryVertexCollection = cms.InputTag("offlineSlimmedPrimaryVertices")
    
)


