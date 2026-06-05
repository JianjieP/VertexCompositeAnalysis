import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras
process = cms.Process('ANASKIM', eras.Run3_2023_UPC)

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Reconstruction_Data_cff')

# Limit the output messages
process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 5000
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))

# Define the input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring("root://xrootd-cms.infn.it//store/hidata/HIRun2023A/HIForward0/MINIAOD/14Feb2025-v1/100000/04b121ec-419e-4278-81b1-b24f6c9f4159.root"),
)
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(60000))

# Set the global tag
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = cms.string('141X_dataRun3_v6')

# Add the Particle producer
process.load("VertexCompositeAnalysis.VertexCompositeProducer.generalParticles_cff")

# Switch the AOD-style track inputs to the MiniAOD unpacked collections.
process.load("VertexCompositeAnalysis.VertexCompositeProducer.unpackedTracksAndVertices_cfi")
from Configuration.Applications.ConfigBuilder import MassReplaceInputTag
process = MassReplaceInputTag(process, "offlinePrimaryVertices", "offlineSlimmedPrimaryVertices")
process = MassReplaceInputTag(process, "generalTracks", "unpackedTracksAndVertices")
process = MassReplaceInputTag(process, "particleFlow", "packedPFCandidates")

# FourPi selection
pionSelection = cms.string("(pt > 0.0 && abs(eta) < 3.0) && quality(\"highPurity\")")
pionFinalSelection = cms.string("") #"abs(userFloat(\"dzSig\"))<3.0 && abs(userFloat(\"dxySig\"))<3.0")
fourPiSelection = cms.string("charge==0")
process.fourPi = process.generalParticles.clone(
#    mass = cms.double(2.0),
#    pdgId = cms.uint32(113),
#    width = cms.double(3.1),
    preSelection = fourPiSelection,
    # daughter information
    daughterInfo = cms.VPSet([
        cms.PSet(pdgId = cms.uint32(211), charge = cms.int32(+1), selection = pionSelection, finalSelection = pionFinalSelection),
        cms.PSet(pdgId = cms.uint32(211), charge = cms.int32(+1), selection = pionSelection, finalSelection = pionFinalSelection),
        cms.PSet(pdgId = cms.uint32(211), charge = cms.int32(-1), selection = pionSelection, finalSelection = pionFinalSelection),
        cms.PSet(pdgId = cms.uint32(211), charge = cms.int32(-1), selection = pionSelection, finalSelection = pionFinalSelection),
    ]),
)
process.oneFourPi = cms.EDFilter("CandViewCountFilter", src = cms.InputTag("fourPi"), minNumber = cms.uint32(1))

# Add fourPi event selection
process.fourTracks = cms.EDFilter("TrackCountFilter", src = cms.InputTag("generalTracks"), minNumber = cms.uint32(4))
process.hpTracks = cms.EDFilter("TrackSelector", src = cms.InputTag("generalTracks"), cut = cms.string("quality(\"highPurity\")"))
process.hpCands = cms.EDProducer("ChargedCandidateProducer", src = cms.InputTag("hpTracks"), particleType = cms.string('pi+'))
process.maxFourHPCands = cms.EDFilter("PATCandViewCountFilter", src = cms.InputTag("hpCands"), minNumber = cms.uint32(0), maxNumber = cms.uint32(4))
process.goodTracks = cms.EDFilter("TrackSelector",
            src = cms.InputTag("generalTracks"),
            cut = pionSelection,
            )
process.fourGoodTracks = cms.EDFilter("TrackCountFilter", src = cms.InputTag("goodTracks"), minNumber = cms.uint32(4))
#process.goodPions = cms.EDProducer("ChargedCandidateProducer",
#            src = cms.InputTag("goodTracks"),
 #           particleType = cms.string('pi+')
 #           )
#process.goodFourPions = cms.EDProducer("CandViewShallowCloneCombiner",
 #           cut = fourPiSelection,
 #           checkCharge = cms.bool(False),
 #           decay = cms.string('goodPions@+ goodPions@+ goodPions@- goodPions@-')
 #           )
#process.oneGoodFourPi = cms.EDFilter("CandViewCountFilter", src = cms.InputTag("goodFourPions"), minNumber = cms.uint32(1))
process.fourPiEvtSel = cms.Sequence(process.fourTracks * process.hpTracks * process.hpCands * process.maxFourHPCands * process.goodTracks * process.fourGoodTracks  )

# Add trigger selection
import HLTrigger.HLTfilters.hltHighLevel_cfi
process.hltFilter = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.hltFilter.andOr = cms.bool(True)
process.hltFilter.throw = cms.bool(False)
process.hltFilter.HLTPaths = [
    # UPC zero bias triggers
    'HLT_HIUPC_ZeroBias_SinglePixelTrack_MaxPixelTrack_v*',
    'HLT_HIUPC_ZeroBias_SinglePixelTrackLowPt_MaxPixelCluster400_v*',
]

# Add PbPb collision event selection
process.load('VertexCompositeAnalysis.VertexCompositeProducer.collisionEventSelection_cff')
process.load('VertexCompositeAnalysis.VertexCompositeProducer.hfCoincFilter_cff')
process.primaryVertexFilter.src = cms.InputTag("offlineSlimmedPrimaryVertices")
process.colEvtSel = cms.Sequence(process.primaryVertexFilter)
process.clustercomp = cms.Sequence(process.clusterCompatibilityFilter)

# Define the event selection sequence
process.eventFilter = cms.Sequence(
    process.hltFilter *
    process.fourPiEvtSel
)
process.eventFilter_step = cms.Path( process.eventFilter )

# Define the analysis steps
process.fourPi_rereco_step = cms.Path(process.eventFilter * process.fourPi * process.oneFourPi)

# tree producer
from VertexCompositeAnalysis.VertexCompositeAnalyzer.particle_tree_cff import particleAna
process.fourPiAna = particleAna.clone(
  recoParticles = cms.InputTag("fourPi"),
  triggerInfo = cms.untracked.VPSet([
    # UPC zero bias triggers
    cms.PSet(path = cms.string('HLT_HIUPC_ZeroBias_SinglePixelTrack_MaxPixelTrack_v*')),
    cms.PSet(path = cms.string('HLT_HIUPC_ZeroBias_SinglePixelTrackLowPt_MaxPixelCluster400_v*')),
  ]),
  selectEvents = cms.string("fourPi_rereco_step"),
  eventFilterNames = cms.untracked.vstring(
      'Flag_colEvtSel',
      'Flag_clustercomp',
      'Flag_primaryVertexFilter',
      'Flag_hfPosFilterNTh9p2',
      'Flag_hfNegFilterNTh8p6',
  )
)
process.generalana_step = cms.EndPath( process.fourPiAna )

# Define the output
process.TFileService = cms.Service("TFileService", fileName = cms.string('fourpi_ana.root'))

# Define the process schedule
process.schedule = cms.Schedule(
    process.eventFilter_step,
    process.fourPi_rereco_step,
    process.generalana_step
)

# Add the event selection filters
process.Flag_colEvtSel = cms.Path(process.eventFilter * process.colEvtSel)
process.Flag_clustercomp = cms.Path(process.clustercomp)
#process.Flag_hfCoincFilter2Th4 = cms.Path(process.eventFilter * process.hfCoincFilter2Th4)
process.Flag_primaryVertexFilter = cms.Path(process.eventFilter * process.primaryVertexFilter)
process.Flag_hfPosFilterNTh9p2 = cms.Path(process.eventFilter * process.hfPosFilterNTh9p2_seq)
process.Flag_hfNegFilterNTh8p6 = cms.Path(process.eventFilter * process.hfNegFilterNTh8p6_seq)

eventFilterPaths = [ process.Flag_colEvtSel , process.Flag_clustercomp , process.Flag_primaryVertexFilter, process.Flag_hfPosFilterNTh9p2, process.Flag_hfNegFilterNTh8p6 ]

process.eventFilter = cms.Sequence(
    process.hltFilter *
    process.primaryVertexFilter *
    process.fourPiEvtSel
)

for P in eventFilterPaths:
    process.schedule.insert(0, P)
