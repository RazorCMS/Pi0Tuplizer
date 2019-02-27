from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
#config.General.requestName = 'CMSSW_8_3_0_QCD_Pt-15to20_EMEnriched_TuneCUETP8M1_13TeV_pythia8_withNxtalCut_test500' 
config.General.requestName = 'CMSSW_9_4_0_Gun_MultiPion_FlatPt-1To15_mc2017_27Feb2019_V1' 
config.General.workArea = 'crab'

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'python/Pi0Tuplizer_94X_GEN_SIM_RECO.py'
config.JobType.outputFiles = ['pi0Ntuple.root']
config.JobType.allowUndistributedCMSSW = True

config.section_("Data")

#config.Data.inputDataset = '/QCD_Pt-15to20_EMEnriched_TuneCUETP8M1_13TeV_pythia8/PhaseIFall16DR-FlatPU28to62HcalNZSRAW_90X_upgrade2017_realistic_v6_C1-v1/GEN-SIM-RAW'
config.Data.inputDataset = '/MultiPion_FlatPt-1To15_PhotonPtFilter/RunIIFall17DRPremix-94X_mc2017_realistic_v10-v1/GEN-SIM-RECO'
#config.Data.lumiMask = 'data/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt'

config.Data.inputDBS = 'global' #change this according to the DBS instance (usually 'global') of the target dataset
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
#config.Data.totalUnits = 500
config.Data.publication = False

config.section_("Site")
config.Site.storageSite = 'T2_US_Caltech'
config.Data.outLFNDirBase = '/store/group/phys_susy/razor/Pi0Calib2019/'
