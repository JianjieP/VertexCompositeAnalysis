import FWCore.ParameterSet.Config as cms

dimuana = cms.EDAnalyzer('PATCompositeTreeProducer',
  doRecoNtuple = cms.untracked.bool(True),
  doMuonNtuple = cms.untracked.bool(False),
  doTrackNtuple = cms.untracked.bool(False),
  doGenNtuple = cms.untracked.bool(False),
  doGenMatching = cms.untracked.bool(False),
  doGenMatchingTOF = cms.untracked.bool(False),
  decayInGen = cms.untracked.bool(False),
  twoLayerDecay = cms.untracked.bool(False),
  threeProngDecay = cms.untracked.bool(False),
  isGammaGamma = cms.untracked.bool(False),

  #PID used only for GEN and/or GEN match
  PID_dau = cms.untracked.vint32(13, 13),
  beamSpotSrc = cms.untracked.InputTag("offlineBeamSpot"),
  VertexCollection = cms.untracked.InputTag("offlinePrimaryVertices"),
  VertexCompositeCollection = cms.untracked.InputTag("generalDiMuCandidates:DiMu"),
  MuonCollection = cms.untracked.InputTag("patMuonsWithTrigger"),
  TrackCollection = cms.untracked.InputTag("generalTracks"),
  PFCandidateCollection = cms.untracked.InputTag("particleFlow"),
  GenParticleCollection = cms.untracked.InputTag("genMuons"),
  doMuon = cms.untracked.bool(True),
  doMuonFull = cms.untracked.bool(True),

  #Trigger info
  TriggerResultCollection = cms.untracked.InputTag("TriggerResults::HLT"),
  triggerPathNames = cms.untracked.vstring(
    "HLT_HIL2DoubleMu0_M1p5to6_Open_v",
    "HLT_HIL2DoubleMu0_M7to15_Open_v",
    "HLT_HIL3DoubleMu0_M2to4p5_Open_v",
    "HLT_HIL3DoubleMu0_M7to15_Open_v",
    "HLT_HIL3DoubleMu0_Quarkonia_Open_v",
    "HLT_HIL3SingleMu12_v",
    "HLT_HIUPC_SingleMuCosmic_NotMBHF2AND_v",
    "HLT_HIUPC_SingleMuCosmic_NotMBHF2OR_v",
    "HLT_HIUPC_SingleMuOpen_NotMBHF2AND_v",
    "HLT_HIUPC_SingleMuOpen_NotMBHF2OR_v",
    "HLT_HIMinimumBiasHF1ANDZDC1nOR_v",
    "HLT_HIMinimumBiasHF1AND_v",
    "HLT_HIUPC_ZeroBias_SinglePixelTrack_MaxPixelTrack_v",
    "HLT_HIUPC_ZeroBias_SinglePixelTrackLowPt_MaxPixelCluster400_v",
    "HLT_HIUPC_ZeroBias_MinPixelCluster400_MaxPixelCluster10000_v",
    "HLT_HIZeroBias_HighRate_v",
  ),
  triggerFilterNames = cms.untracked.vstring(
    'hltL2fL1fL1sDoubleMuOpenL2Filtered0',
    'hltL2DoubleMuOpenMassFiltered7to15',
    'hltL3DoubleMuOpenMassFiltered2to4p5',
    'hltL3DoubleMuOpenMassFiltered7to15',
    'hltL3DoubleMuOpenMassFiltered2to4p5OR7to15',
    'hltL3fL1fL1sSingleMu7L3Filtered12',
    'hltL1sSingleMuCosmicNotMBHF2AND',
    'hltL1sSingleMuCosmicNotMBHF2OR',
    'hltL1sSingleMuOpenNotMBHF2AND',
    'hltL1sSingleMuOpenNotMBHF2OR',
    '',
    '',
    '',
    '',
    '',
    '',
  ),
  stageL1Trigger = cms.uint32(2),

  #Filter info
  FilterResultCollection = cms.untracked.InputTag("TriggerResults::ANASKIM"),
  eventFilterNames = cms.untracked.vstring(
      'Flag_colEvtSel',
      'Flag_hfCoincFilter2Th4',
      'Flag_primaryVertexFilter',
      'Flag_clusterCompatibilityFilter',
  ),
  selectEvents = cms.untracked.string(""),

  isCentrality = cms.bool(True),
  centralityBinLabel = cms.InputTag("centralityBin","HFtowers"),
  centralitySrc = cms.InputTag("hiCentrality"),

  isEventPlane = cms.bool(True),
  eventplaneSrc = cms.InputTag("hiEvtPlaneFlat"),

  saveTree = cms.untracked.bool(True),
  saveHistogram = cms.untracked.bool(False),
  saveAllHistogram = cms.untracked.bool(False),
  massHistPeak = cms.untracked.double(60.0),
  massHistWidth = cms.untracked.double(60.0),
  massHistBins = cms.untracked.int32(1200),

  pTBins = cms.untracked.vdouble(0.0,0.2,1.8,3.0,4.5,6.0,8.0,10.,20.),
  yBins = cms.untracked.vdouble(-2.4,-1.4,0,1.4,2.4),

  useAnyMVA = cms.bool(False),
  isSkimMVA = cms.untracked.bool(False),
  MVACollection = cms.InputTag("generalDiMuCandidates:MVAValues")
)

dimuana_mc = dimuana.clone(
  doGenNtuple = cms.untracked.bool(True),
  doGenMatching = cms.untracked.bool(True),
  doGenMatchingTOF = cms.untracked.bool(True),
  decayInGen = cms.untracked.bool(True),
)
