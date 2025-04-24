from WMCore.Configuration import Configuration

config = Configuration()

config.section_("General")
config.General.requestName = "Upsilon1S"
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = False

config.section_("JobType")
config.JobType.pluginName = "Analysis"
config.JobType.psetName = "hioniaanalyzer_PbPbPrompt_13_2_X_MC_cfg.py"
config.JobType.maxMemoryMB = 5000         # request high memory machines.
config.JobType.numCores = 4
config.JobType.allowUndistributedCMSSW = True #Problems with slc7
config.JobType.maxJobRuntimeMin = 1000 #2750    # request longer runtime, ~48 hours.


config.section_("Data")
config.Data.inputDataset = '/Upsilon1SToMuMu_Pthat2_TuneCP5_HydjetDrumMB_5p36TeV_pythia8/HINPbPbSpring23MiniAOD-132X_mcRun3_2023_realistic_HI_v9-v2/MINIAODSIM'
config.Data.inputDBS = 'global'
config.Data.unitsPerJob = 20
#config.Data.totalUnits = -1
config.Data.splitting = "FileBased"
config.Data.allowNonValidInputDataset = True
config.Data.outputDatasetTag = config.General.requestName

config.Data.outLFNDirBase = '/store/user/fdamas/PbPb2023/Upsilon/EmbeddedMC'
config.Data.publication = False

config.section_("Site")
config.Site.storageSite = "T3_CH_CERNBOX"
#config.Site.whitelist = ["T2_US_*","T1_US_*","T2_CH_CERN","T2_FR_*"]