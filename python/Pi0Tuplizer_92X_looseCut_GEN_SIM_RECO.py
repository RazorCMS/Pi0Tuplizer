#########################options##############################
isMC_ = True
MCAssoc_ = True
useDRcutPair_ = False
FillL1SeedFinalDecision_ = False
FillDiPhotonNtuple_ = True
FillPhotonNtuple_ = True
#########################options##############################
import FWCore.ParameterSet.Config as cms
import RecoLocalCalo.EcalRecProducers.ecalRecalibRecHit_cfi
import os, sys, imp, re


CMSSW_VERSION=os.getenv("CMSSW_VERSION")
process = cms.Process("Pi0Tuplizer")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
process.load("Configuration.StandardSequences.Services_cff")
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load("Configuration.EventContent.EventContent_cff")
process.load("SimGeneral.MixingModule.mixNoPU_cfi")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.RawToDigi_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("CommonTools.ParticleFlow.EITopPAG_cff")
process.load("DQMOffline.Configuration.DQMOfflineMC_cff")

process.GlobalTag.globaltag = '92X_upgrade2017_realistic_v10'

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.MessageLogger.cerr.FwkReport.reportEvery = 10000
process.options = cms.untracked.PSet(
   wantSummary = cms.untracked.bool(True),
   SkipEvent = cms.untracked.vstring('ProductNotFound','CrystalIDError')
)

#load input file
process.source = cms.Source('PoolSource',
    fileNames = cms.untracked.vstring(
        'file:/afs/cern.ch/work/z/zhicaiz/public/release/PiZero_git/CMSSW_9_2_8/src/PiZero/Pi0Tuplizer/python/step2_1_pi0.root',
        'file:/afs/cern.ch/work/z/zhicaiz/public/release/PiZero_git/CMSSW_9_2_8/src/PiZero/Pi0Tuplizer/python/step2_107_pi0.root'
        #'file:/afs/cern.ch/work/z/zhicaiz/public/release/PiZero_git/CMSSW_9_2_8/src/PiZero/Pi0Tuplizer/python/step2_10_eta.root'
	#'file:/afs/cern.ch/work/z/zhicaiz/public/release/McM/pi0/CMSSW_9_2_8/src/step2_eta_gg.root'
	#'file:/afs/cern.ch/work/z/zhicaiz/public/release/McM/pi0/CMSSW_9_2_8/src/step2_eta_ppp.root'
    )
)

#define output file
process.TFileService = cms.Service("TFileService",
    fileName = cms.string("pi0Ntuple.root"),
    closeFileFast = cms.untracked.bool(True)
)

#provide input parameters
process.ntuples = cms.EDAnalyzer('Pi0Tuplizer',
FillL1SeedFinalDecision = cms.untracked.bool(True),
uGtAlgInputTag = cms.untracked.InputTag('hltGtStage2Digis'),
EBRecHitCollectionTag_Pi0 = cms.untracked.InputTag("ecalRecHit","EcalRecHitsEB","RECO"),
EERecHitCollectionTag_Pi0 = cms.untracked.InputTag("ecalRecHit","EcalRecHitsEE","RECO"),
ESRecHitCollectionTag_Pi0 = cms.untracked.InputTag("ecalPreshowerRecHit","EcalRecHitsES","RECO"),
EBRecHitCollectionTag_Eta = cms.untracked.InputTag("ecalRecHit","EcalRecHitsEB","RECO"),
EERecHitCollectionTag_Eta = cms.untracked.InputTag("ecalRecHit","EcalRecHitsEE","RECO"),
ESRecHitCollectionTag_Eta = cms.untracked.InputTag("ecalPreshowerRecHit","EcalRecHitsES","RECO"),
genParticles = cms.InputTag("genParticles"),
PhotonOrderOption = cms.untracked.string("SeedEBased"),# "SeedEBased" (g1 is the one with larger seed Energy) or "PhoPtBased"(g1 is the one with larger pt)
MC_Asssoc_DeltaR = cms.untracked.double(0.1),
EB_Seed_E_Pi0_ = cms.untracked.double(0.5),
EE_Seed_E_Pi0_ = cms.untracked.double(0.5),
pairPtCut_barrel1_Pi0_ = cms.untracked.double(1.0),
pairPtCut_barrel2_Pi0_ = cms.untracked.double(1.0),
pairPtCut_endcap1_Pi0_ = cms.untracked.double(1.0),
pairPtCut_endcap2_Pi0_ = cms.untracked.double(1.0),
gPtCut_barrel1_Pi0_ = cms.untracked.double(1.00),
gPtCut_barrel2_Pi0_ = cms.untracked.double(1.00),
gPtCut_endcap1_Pi0_ = cms.untracked.double(1.00),
gPtCut_endcap2_Pi0_ = cms.untracked.double(1.00),
s4s9Cut_barrel1_Pi0_ = cms.untracked.double(0.75),
s4s9Cut_barrel2_Pi0_ = cms.untracked.double(0.75),
s4s9Cut_endcap1_Pi0_ = cms.untracked.double(0.75),
s4s9Cut_endcap2_Pi0_ = cms.untracked.double(0.75),
nxtal1Cut_barrel1_Pi0_ = cms.untracked.double(6.5),
nxtal1Cut_barrel2_Pi0_ = cms.untracked.double(0.5),
nxtal1Cut_endcap1_Pi0_ = cms.untracked.double(6.5),
nxtal1Cut_endcap2_Pi0_ = cms.untracked.double(6.5),
nxtal2Cut_barrel1_Pi0_ = cms.untracked.double(6.5),
nxtal2Cut_barrel2_Pi0_ = cms.untracked.double(6.5),
nxtal2Cut_endcap1_Pi0_ = cms.untracked.double(6.5),
nxtal2Cut_endcap2_Pi0_ = cms.untracked.double(6.5),
EB_Seed_E_Eta_ = cms.untracked.double(0.5),
EE_Seed_E_Eta_ = cms.untracked.double(0.5),
pairPtCut_barrel1_Eta_ = cms.untracked.double(1.0),
pairPtCut_barrel2_Eta_ = cms.untracked.double(1.0),
pairPtCut_endcap1_Eta_ = cms.untracked.double(1.0),
pairPtCut_endcap2_Eta_ = cms.untracked.double(1.0),
gPtCut_barrel1_Eta_ = cms.untracked.double(1.00),
gPtCut_barrel2_Eta_ = cms.untracked.double(1.00),
gPtCut_endcap1_Eta_ = cms.untracked.double(1.00),
gPtCut_endcap2_Eta_ = cms.untracked.double(1.00),
s4s9Cut_barrel1_Eta_ = cms.untracked.double(0.75),
s4s9Cut_barrel2_Eta_ = cms.untracked.double(0.75),
s4s9Cut_endcap1_Eta_ = cms.untracked.double(0.75),
s4s9Cut_endcap2_Eta_ = cms.untracked.double(0.75),
nxtal1Cut_barrel1_Eta_ = cms.untracked.double(6.5),
nxtal1Cut_barrel2_Eta_ = cms.untracked.double(6.5),
nxtal1Cut_endcap1_Eta_ = cms.untracked.double(6.5),
nxtal1Cut_endcap2_Eta_ = cms.untracked.double(6.5),
nxtal2Cut_barrel1_Eta_ = cms.untracked.double(6.5),
nxtal2Cut_barrel2_Eta_ = cms.untracked.double(6.5),
nxtal2Cut_endcap1_Eta_ = cms.untracked.double(6.5),
nxtal2Cut_endcap2_Eta_ = cms.untracked.double(6.5),
DRcutPair_Pi0_	= cms.untracked.double(0.30),
DRcutPair_Eta_	= cms.untracked.double(0.30),
isoGammaBeltdR_Zone_Pi0_        = cms.untracked.double(0.2),
isoGammaBeltdR_Zone_Eta_        = cms.untracked.double(0.2),
isoPairBeltdR_Zone_Pi0_         = cms.untracked.double(0.2),
isoPairBeltdR_Zone_Eta_         = cms.untracked.double(0.3),
isoGammaBeltdEta_Zone_Pi0_      = cms.untracked.double(0.05),
isoGammaBeltdEta_Zone_Eta_      = cms.untracked.double(0.1),
isoPairBeltdEta_Zone_Pi0_       = cms.untracked.double(0.05),
isoPairBeltdEta_Zone_Eta_       = cms.untracked.double(0.1),
isoPairCut_                     = cms.untracked.double(999.9),
isoGammaCut_                    = cms.untracked.double(999.9)
)

#########################paratmeters for the tuplizer##############################
if isMC_:
        process.ntuples.isMC = cms.untracked.bool(True)
else:
        process.ntuples.isMC = cms.untracked.bool(False)
if MCAssoc_:
        process.ntuples.MCAssoc = cms.untracked.bool(True)
else:
        process.ntuples.MCAssoc = cms.untracked.bool(False)
if useDRcutPair_:
	process.ntuples.useDRcutPair = cms.untracked.bool(True)
else:
	process.ntuples.useDRcutPair = cms.untracked.bool(False)
if FillL1SeedFinalDecision_:
        process.ntuples.FillL1SeedFinalDecision = cms.untracked.bool(True)
else:
        process.ntuples.FillL1SeedFinalDecision = cms.untracked.bool(False)
if FillDiPhotonNtuple_:
        process.ntuples.FillDiPhotonNtuple = cms.untracked.bool(True)
else:
        process.ntuples.FillDiPhotonNtuple = cms.untracked.bool(False)
if FillPhotonNtuple_:
        process.ntuples.FillPhotonNtuple = cms.untracked.bool(True)
else:
        process.ntuples.FillPhotonNtuple = cms.untracked.bool(False)

#define path
process.p = cms.Path()
process.p *= process.ntuples
