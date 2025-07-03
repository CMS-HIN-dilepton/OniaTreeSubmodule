from CRABAPI.RawCommand import crabCommand
from CRABClient.ClientExceptions import ClientException
from http.client import HTTPException

# We want to put all the CRAB project directories from the tasks we submit here into one common directory.
# That's why we need to set this parameter (here or above in the configuration file, it does not matter, we will not overwrite it).
from CRABClient.UserUtilities import config
config = config()

config.section_("General")
#config.General.requestName = "HIPhysicsRawPrime5_PromptReco_v2"
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = False

config.section_("JobType")
config.JobType.pluginName = "Analysis"
config.JobType.psetName = "hioniaanalyzer_pO_privateRecoMINIAOD_cfg.py"
config.JobType.maxMemoryMB = 2000         # request high memory machines.
config.JobType.numCores = 4
config.JobType.allowUndistributedCMSSW = True #Problems with slc7
config.JobType.maxJobRuntimeMin = 200 # max = 2750

config.section_("Data")
config.Data.inputDBS = 'global'
#config.Data.totalUnits = -1
config.Data.splitting = "FileBased"
config.Data.unitsPerJob = 50

config.Data.allowNonValidInputDataset = True
config.Data.publication = False
config.Data.runRange = '393952-394007'
config.Data.lumiMask = 'goodLSfromDCS.json' #temporary JSON file (local!!)

config.section_("Site")
config.Site.storageSite = "T3_CH_CERNBOX"
#config.Site.whitelist = ["T2_US_*","T2_CH_CERN","T1_US_*"]

# Multi crab part

def submit(config):
    try:
        crabCommand('submit', config = config, dryrun=False)
    except HTTPException as hte:
        print("Failed submitting task: %s" % (hte.headers))
    except ClientException as cle:
        print("Failed submitting task: %s" % (cle))

# Submit the jobs: 20 HIForward PDs, ~140k files each, average of 100k events/file

config.Data.outLFNDirBase = '/store/user/fdamas/LightIon2025/pO/'


for i in range(60):

    config.General.requestName = f'DimuonSkim_{i}'
    config.Data.inputDataset = f"/IonPhysics{i}/pORun2025-IonDimuon-PromptReco-v1/USER"
    config.Data.outputDatasetTag = config.General.requestName


    print("Submitting CRAB job for: "+ config.Data.inputDataset)
    submit(config)