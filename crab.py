from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
#config.General.requestName = 'CMSSW_8_3_0_QCD_Pt-15to20_EMEnriched_TuneCUETP8M1_13TeV_pythia8_withNxtalCut_test500' 
config.General.requestName = 'CMSSW_8_3_0_Gun_SinglePion_FlatPt-1To15_07Aug2017_V1' 
config.General.workArea = 'crab'

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'python/Pi0Tuplizer_830_GEN_SIM_RECO.py'
config.JobType.outputFiles = ['pi0Ntuple.root']

config.section_("Data")

#config.Data.inputDataset = '/QCD_Pt-15to20_EMEnriched_TuneCUETP8M1_13TeV_pythia8/PhaseIFall16DR-FlatPU28to62HcalNZSRAW_90X_upgrade2017_realistic_v6_C1-v1/GEN-SIM-RAW'
config.Data.inputDataset = '/SinglePion_FlatPt-1To15/PhaseIFall16DR-FlatPU10to50RECO_90X_upgrade2017_realistic_v6_C1-v1/GEN-SIM-RECO'
#config.Data.lumiMask = 'data/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt'

config.Data.inputDBS = 'global' #change this according to the DBS instance (usually 'global') of the target dataset
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 3
#config.Data.totalUnits = 500
config.Data.publication = False
config.Data.ignoreLocality = True #enable AAA

config.section_("Site")
config.Site.storageSite = 'T2_CH_CERN'
config.Data.outLFNDirBase = '/store/group/dpg_ecal/alca_ecalcalib/piZero2017/zhicaiz/Gun_SinglePion_FlatPt-1To15'
