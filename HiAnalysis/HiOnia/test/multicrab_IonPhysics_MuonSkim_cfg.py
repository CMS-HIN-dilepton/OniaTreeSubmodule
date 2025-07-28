from CRABAPI.RawCommand import crabCommand
from CRABClient.ClientExceptions import ClientException
from http.client import HTTPException

# We want to put all the CRAB project directories from the tasks we submit here into one common directory.
# That's why we need to set this parameter (here or above in the configuration file, it does not matter, we will not overwrite it).
from CRABClient.UserUtilities import config
config = config()

from CRABClient.UserUtilities import getUsername
username = getUsername()

##########################

collisionSystem = 'OO' # pO, OO, or NeNe

muonSkim = 'IonDimuon' # IonDimuon or IonHighPtMuon


config.section_("General")
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = False

config.section_("JobType")
config.JobType.pluginName = "Analysis"
config.JobType.psetName = "hioniaanalyzer_LightIon2025_DATA_cfg.py"

config.JobType.maxMemoryMB = 2000         # request high memory machines.
#config.JobType.numCores = 4
config.JobType.allowUndistributedCMSSW = True
config.JobType.maxJobRuntimeMin = 200 # max = 2750

config.section_("Data")
config.Data.inputDBS = 'global'
#config.Data.totalUnits = -1
config.Data.splitting = "FileBased"
config.Data.unitsPerJob = 50

config.Data.allowNonValidInputDataset = True
config.Data.publication = False

config.Data.outLFNDirBase = '/store/user/' + username + '/LightIon2025/' + collisionSystem


config.section_("Site")
config.Site.storageSite = "T3_CH_CERNBOX"
#config.Site.whitelist = ["T2_US_*","T2_CH_CERN","T1_US_*"]

# settings based on collision system name
if collisionSystem == 'pO':
    config.Data.runRange = '393952-394007'
    config.Data.lumiMask = 'https://cms-service-dqmdc.web.cern.ch/CAF/certification/Collisions25pO/pO_muon.json' # final muon json

elif collisionSystem == 'OO':
    config.Data.runRange = '394153-394217'
    config.Data.lumiMask = 'https://cms-service-dqmdc.web.cern.ch/CAF/certification/Collisions25OO/Cert_Collisions2025OO_394153_394217_muon.json'

elif collisionSystem == 'NeNe':
    config.Data.runRange = '394269-394272'
    config.Data.lumiMask = 'https://cms-service-dqmdc.web.cern.ch/CAF/certification/Collisions25NeNe/Cert_Collisions2025Nene_394269_394272_muon.json'

else:
    print("This config script does not support CRAB job submission for collision name: %s. Check the settings!" % (collisionSystem))

# Multi crab part

def submit(config):
    try:
        crabCommand('submit', config = config, dryrun=False)
    except HTTPException as hte:
        print("Failed submitting task: %s" % (hte.headers))
    except ClientException as cle:
        print("Failed submitting task: %s" % (cle))

# Submit the jobs: 60 IonPhysics PDs

for i in range(60):

    config.General.requestName = f'{collisionSystem}_{muonSkim}_{i}'
    config.Data.inputDataset = f"/IonPhysics{i}/{collisionSystem}Run2025-{muonSkim}-PromptReco-v1/USER"
    config.Data.outputDatasetTag = config.General.requestName

    print("Submitting CRAB job for: "+ config.Data.inputDataset)
    submit(config)