
from CRABClient.UserUtilities import config
from CRABAPI.RawCommand import crabCommand
#from CRABClient.UserUtilities import getUsername

config = config()

config.section_("General")
config.General.workArea = 'crab_projects_0520'
config.General.transferOutputs = True
config.General.transferLogs = False

config.section_('JobType')
config.JobType.pluginName = 'Analysis'
config.JobType.allowUndistributedCMSSW = True
#config.JobType.psetName = 'PbPbSkimTree2023_DiMu_Analyzer_cfg.py'
config.JobType.psetName = 'PbPbSkimAndTree2023_FourPi_cfg.py'
#config.JobType.inputFiles = ['emap_2023_newZDC_v3.txt','CentralityTable_HFtowers200_DataPbPb_periHYDJETshape_run3v1302x04_offline_374289.db']
config.JobType.maxMemoryMB = 2500
#config.JobType.maxJobRuntimeMin = 2100
config.JobType.maxJobRuntimeMin = 1100

config.section_('Data')


config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased'
#config.Data.lumiMask = 'Cert_Collisions2023HI_374288_375659_Muon.json'
config.Data.lumiMask = 'Cert_Collisions2023HI_374288_375823_Golden.json'
#config.Data.lumiMask = 'Cert_Collisions2023HI_374288_375823_Golden.json'
config.Data.publication = False
config.Data.allowNonValidInputDataset = True

config.section_('Site')
config.Site.storageSite = 'T3_CH_CERNBOX'

dataMap = {
    0: '/HIForward0/HIRun2023A-14Feb2025-v1/MINIAOD',
    1: '/HIForward1/HIRun2023A-14Feb2025-v1/MINIAOD',
    2: '/HIForward2/HIRun2023A-14Feb2025-v1/MINIAOD',
    3: '/HIForward3/HIRun2023A-14Feb2025-v1/MINIAOD',
    4: '/HIForward4/HIRun2023A-14Feb2025-v1/MINIAOD',
    5: '/HIForward5/HIRun2023A-14Feb2025-v1/MINIAOD',
    6: '/HIForward6/HIRun2023A-14Feb2025-v1/MINIAOD',
    7: '/HIForward7/HIRun2023A-14Feb2025-v1/MINIAOD',
    8: '/HIForward8/HIRun2023A-14Feb2025-v1/MINIAOD',
    9: '/HIForward9/HIRun2023A-14Feb2025-v1/MINIAOD',
    10: '/HIForward10/HIRun2023A-14Feb2025-v1/MINIAOD',
    11: '/HIForward11/HIRun2023A-14Feb2025-v1/MINIAOD',
    12: '/HIForward12/HIRun2023A-14Feb2025-v1/MINIAOD',
    13: '/HIForward13/HIRun2023A-14Feb2025-v1/MINIAOD',
    14: '/HIForward14/HIRun2023A-14Feb2025-v1/MINIAOD',
    15: '/HIForward15/HIRun2023A-14Feb2025-v1/MINIAOD',
    16: '/HIForward16/HIRun2023A-14Feb2025-v1/MINIAOD',
    17: '/HIForward17/HIRun2023A-14Feb2025-v1/MINIAOD',
    18: '/HIForward18/HIRun2023A-14Feb2025-v1/MINIAOD',
    19: '/HIForward19/HIRun2023A-14Feb2025-v1/MINIAOD',
}

for i in range(0,20):

    ## Submit the pion PDs
    config.General.requestName = f'HIForward{i}_HIRun2023A-14Feb2025-v1_FourPi_0520'
    config.Data.inputDataset = dataMap[i]
    config.Data.unitsPerJob = 20
    config.Data.outputDatasetTag = config.General.requestName
    config.Data.outLFNDirBase = f'/store/user/jianjie/fourpi/HIForward{i}/'

    crabCommand('submit', config=config)