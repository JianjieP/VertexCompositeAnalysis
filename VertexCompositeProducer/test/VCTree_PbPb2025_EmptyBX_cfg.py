import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras
process = cms.Process('ANASKIM', eras.Run3_2025_UPC)

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Reconstruction_Data_cff')

# Limit the output messages
process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 200
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))

# Define the input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring("root://xrootd-cms.infn.it//store/hidata/HIRun2025A/HIEmptyBX/AOD/PromptReco-v1/000/398/954/00000/e169ad6f-9085-4919-af93-55acb7d9893e.root"),
    )
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(2000))

# Set the global tag
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = cms.string('151X_dataRun3_Prompt_v1')

# Add trigger selection
import HLTrigger.HLTfilters.hltHighLevel_cfi
process.hltFilter = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.hltFilter.andOr = cms.bool(True)
process.hltFilter.throw = cms.bool(False)
process.hltFilter.HLTPaths = [
    # Empty BX triggers
    'HLT_HIL1NotBptxOR_v*', #0
    'HLT_HIL1UnpairedBunchBptxMinus_v*', #1
    'HLT_HIL1UnpairedBunchBptxPlus_v*', #2
]
process.eventFilterHLT = cms.Sequence(process.hltFilter)
process.eventFilterHLT_step = cms.Path(process.eventFilterHLT)

# Add the Particle tree
from VertexCompositeAnalysis.VertexCompositeAnalyzer.particle_tree_cff import particleAna

process.eventAna = particleAna.clone(
  recoParticles = cms.InputTag(""),
  selectEvents = cms.string("eventFilterHLT_step"),
  eventFilterNames = cms.untracked.vstring(
      'Flag_primaryVertexFilter',
      'Flag_hiClusterCompatibility',
  ),
  triggerInfo = cms.untracked.VPSet([
    # Empty BX triggers
    cms.PSet(path = cms.string('HLT_HIL1NotBptxOR_v*')), #0
    cms.PSet(path = cms.string('HLT_HIL1UnpairedBunchBptxMinus_v*')), #1
    cms.PSet(path = cms.string('HLT_HIL1UnpairedBunchBptxPlus_v*')), #2
  ]),
#   dataset = cms.untracked.string("Dataset_HIForward*"),
)

# Define the output
process.TFileService = cms.Service("TFileService", fileName = cms.string('ana_EmptyBX.root'))
process.p = cms.EndPath(process.eventAna)

# Define the process schedule
process.schedule = cms.Schedule(
    process.eventFilterHLT_step,
    process.p,
)

# Add the event selection filters
process.load('VertexCompositeAnalysis.VertexCompositeProducer.collisionEventSelection_cff')
process.Flag_primaryVertexFilter = cms.Path(process.primaryVertexFilter)
process.Flag_hiClusterCompatibility = cms.Path(process.hiClusterCompatibility)

eventFilterPaths = [ process.Flag_primaryVertexFilter , process.Flag_hiClusterCompatibility ]

for P in eventFilterPaths:
    process.schedule.insert(0, P)