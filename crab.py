from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = 'CMSSW_8_0_25_ALlCaP0_Run2016B_v2_RAW_Pi0Ntuple_PhoNtuple_test500' 
config.General.workArea = 'crab'

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'python/Pi0Tuplizer_80X_AlCaP0_RAW.py'
config.JobType.outputFiles = ['pi0Ntuple.root']

config.section_("Data")

config.Data.inputDataset = '/AlCaP0/Run2016B-v2/RAW'
config.Data.lumiMask = 'data/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt'

config.Data.inputDBS = 'global' #change this according to the DBS instance (usually 'global') of the target dataset
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
config.Data.totalUnits = 500
config.Data.publication = False
config.Data.ignoreLocality = True #enable AAA

config.section_("Site")
config.Site.storageSite = 'T2_CH_CERN'
config.Data.outLFNDirBase = '/store/group/dpg_ecal/alca_ecalcalib/piZero2016/zhicaiz/Pi0Tuplizer_2016B'
