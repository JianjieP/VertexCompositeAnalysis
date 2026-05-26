
from CRABClient.UserUtilities import config
from CRABAPI.RawCommand import crabCommand
#from CRABClient.UserUtilities import getUsername

config = config()

config.section_("General")
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = False

config.section_('JobType')
config.JobType.pluginName = 'Analysis'
config.JobType.allowUndistributedCMSSW = True
#config.JobType.psetName = 'PbPbSkimTree2023_DiMu_Analyzer_cfg.py'
config.JobType.psetName = 'PbPbSkimAndTree2025_FourPi_cfg.py'
#config.JobType.inputFiles = ['emap_2023_newZDC_v3.txt','CentralityTable_HFtowers200_DataPbPb_periHYDJETshape_run3v1302x04_offline_374289.db']
config.JobType.maxMemoryMB = 2500
#config.JobType.maxJobRuntimeMin = 2100
config.JobType.maxJobRuntimeMin = 1100

config.section_('Data')


config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased'
#config.Data.lumiMask = 'Cert_Collisions2023HI_374288_375659_Muon.json'
config.Data.lumiMask = 'Cert_Collisions2025_HI_399465_400426_Golden.json'
#config.Data.lumiMask = 'Cert_Collisions2023HI_374288_375823_Golden.json'
config.Data.publication = False
config.Data.allowNonValidInputDataset = True

config.section_('Site')
config.Site.storageSite = 'T3_CH_CERNBOX'

dataMap = {
    0: '/HIForward0/HIRun2025A-PromptReco-v1/AOD',
    1: '/HIForward1/HIRun2025A-PromptReco-v1/AOD',
    2: '/HIForward2/HIRun2025A-PromptReco-v1/AOD',
    3: '/HIForward3/HIRun2025A-PromptReco-v1/AOD',
    4: '/HIForward4/HIRun2025A-PromptReco-v1/AOD',
    5: '/HIForward5/HIRun2025A-PromptReco-v1/AOD',
    6: '/HIForward6/HIRun2025A-PromptReco-v1/AOD',
    7: '/HIForward7/HIRun2025A-PromptReco-v1/AOD',
    8: '/HIForward8/HIRun2025A-PromptReco-v1/AOD',
    9: '/HIForward9/HIRun2025A-PromptReco-v1/AOD',
    10: '/HIForward10/HIRun2025A-PromptReco-v1/AOD',
    11: '/HIForward11/HIRun2025A-PromptReco-v1/AOD',
    12: '/HIForward12/HIRun2025A-PromptReco-v1/AOD',
    13: '/HIForward13/HIRun2025A-PromptReco-v1/AOD',
    14: '/HIForward14/HIRun2025A-PromptReco-v1/AOD',
    15: '/HIForward15/HIRun2025A-PromptReco-v1/AOD',
    16: '/HIForward16/HIRun2025A-PromptReco-v1/AOD',
    17: '/HIForward17/HIRun2025A-PromptReco-v1/AOD',
    18: '/HIForward18/HIRun2025A-PromptReco-v1/AOD',
    19: '/HIForward19/HIRun2025A-PromptReco-v1/AOD',
    20: '/HIForward20/HIRun2025A-PromptReco-v1/AOD',
    21: '/HIForward21/HIRun2025A-PromptReco-v1/AOD',
    22: '/HIForward22/HIRun2025A-PromptReco-v1/AOD',
    23: '/HIForward23/HIRun2025A-PromptReco-v1/AOD'
}

for i in range(0,24):

    ## Submit the pion PDs
    config.General.requestName = f'HIForward{i}_HIRun2025A_FourPi_0410'
    config.Data.inputDataset = dataMap[i]
    config.Data.unitsPerJob = 20
    config.Data.outputDatasetTag = config.General.requestName
    config.Data.outLFNDirBase = f'/store/user/jianjie/fourpi/HIForward{i}/'

    crabCommand('submit', config=config)