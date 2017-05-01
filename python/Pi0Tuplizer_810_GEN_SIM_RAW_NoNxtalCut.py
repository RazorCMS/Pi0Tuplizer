#########################options##############################
isMC_ = True
FillL1SeedFinalDecision_ = True
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

process.GlobalTag.globaltag = '81X_upgrade2017_realistic_v26'

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.options = cms.untracked.PSet(
   wantSummary = cms.untracked.bool(True),
   SkipEvent = cms.untracked.vstring('ProductNotFound','CrystalIDError')
)

#load input file
process.source = cms.Source('PoolSource',
    fileNames = cms.untracked.vstring(
        #'root://cmsxrootd.gridka.de//store/mc/PhaseIFall16DR/QCD_Pt-15to20_EMEnriched_TuneCUETP8M1_13TeV_pythia8/GEN-SIM-RAW/FlatPU28to62HcalNZSRAW_81X_upgrade2017_realistic_v26-v1/100000/002CF4F7-A4F0-E611-8C1F-141877411EA2.root'
        'root://cms-xrd-global.cern.ch//store/mc/PhaseIFall16DR/QCD_Pt-15to20_EMEnriched_TuneCUETP8M1_13TeV_pythia8/GEN-SIM-RAW/FlatPU28to62HcalNZSRAW_81X_upgrade2017_realistic_v26-v1/100000/002CF4F7-A4F0-E611-8C1F-141877411EA2.root'
	#'root://cms-xrd-global.cern.ch//store/mc/PhaseIFall16DR/QCD_Pt-20to30_EMEnriched_TuneCUETP8M1_13TeV_pythia8/GEN-SIM-RAW/FlatPU28to62HcalNZSRAW_81X_upgrade2017_realistic_v26-v1/100000/00013679-37EE-E611-A51C-E0071B7A45D0.root'
	#'root://cms-xrd-global.cern.ch//store/mc/PhaseIFall16DR/QCD_Pt-30to50_EMEnriched_TuneCUETP8M1_13TeV_pythia8/GEN-SIM-RAW/FlatPU28to62HcalNZSRAW_81X_upgrade2017_realistic_v26-v1/100000/000A5CBA-FCEF-E611-9FAE-0CC47A4D760C.root'
	#'root://cms-xrd-global.cern.ch//store/mc/PhaseIFall16DR/QCD_Pt-50to80_EMEnriched_TuneCUETP8M1_13TeV_pythia8/GEN-SIM-RAW/FlatPU28to62HcalNZSRAW_81X_upgrade2017_realistic_v26-v1/100000/0000B688-8DF0-E611-819F-02163E01446C.root'
	#'root://cms-xrd-global.cern.ch//store/mc/PhaseIFall16DR/QCD_Pt-80to120_EMEnriched_TuneCUETP8M1_13TeV_pythia8/GEN-SIM-RAW/FlatPU28to62HcalNZSRAW_81X_upgrade2017_realistic_v26-v1/110000/00288A9D-91EF-E611-8E6B-0025905B85CA.root'
	#'root://cms-xrd-global.cern.ch//store/mc/PhaseIFall16DR/QCD_Pt-120to170_EMEnriched_TuneCUETP8M1_13TeV_pythia8/GEN-SIM-RAW/FlatPU28to62HcalNZSRAW_81X_upgrade2017_realistic_v26-v1/100000/003C4BE7-4EF0-E611-9459-001EC94BFFEB.root'
    )
)

#define output file
process.TFileService = cms.Service("TFileService",
    fileName = cms.string("pi0Ntuple_withoutNxtalCut.root"),
    closeFileFast = cms.untracked.bool(True)
)

#provide input parameters
process.ntuples = cms.EDAnalyzer('Pi0Tuplizer',
FillL1SeedFinalDecision = cms.untracked.bool(True),
uGtAlgInputTag = cms.untracked.InputTag('hltGtStage2Digis'),
EBRecHitCollectionTag_Eta = cms.untracked.InputTag("ecalRecHit","EcalRecHitsEB","Pi0Tuplizer"),
EERecHitCollectionTag_Eta = cms.untracked.InputTag("ecalRecHit","EcalRecHitsEE","Pi0Tuplizer"),
ESRecHitCollectionTag_Eta = cms.untracked.InputTag("ecalPreshowerRecHit","EcalRecHitsES","Pi0Tuplizer"),
EBRecHitCollectionTag_Pi0 = cms.untracked.InputTag("ecalRecHit","EcalRecHitsEB","Pi0Tuplizer"),
EERecHitCollectionTag_Pi0 = cms.untracked.InputTag("ecalRecHit","EcalRecHitsEE","Pi0Tuplizer"),
ESRecHitCollectionTag_Pi0 = cms.untracked.InputTag("ecalPreshowerRecHit","EcalRecHitsES","Pi0Tuplizer"),
PhotonOrderOption = cms.untracked.string("SeedEBased"),# "SeedEBased" (g1 is the one with larger seed Energy) or "PhoPtBased"(g1 is the one with larger pt)
EB_Seed_E_Pi0_ = cms.untracked.double(0.5),
EE_Seed_E_Pi0_ = cms.untracked.double(0.5),
pairPtCut_barrel1_Pi0_ = cms.untracked.double(2.0),
pairPtCut_barrel2_Pi0_ = cms.untracked.double(3.75),
pairPtCut_endcap1_Pi0_ = cms.untracked.double(2.0),
pairPtCut_endcap2_Pi0_ = cms.untracked.double(2.0),
gPtCut_barrel1_Pi0_ = cms.untracked.double(0.65),
gPtCut_barrel2_Pi0_ = cms.untracked.double(1.1),
gPtCut_endcap1_Pi0_ = cms.untracked.double(0.95),
gPtCut_endcap2_Pi0_ = cms.untracked.double(0.95),
s4s9Cut_barrel1_Pi0_ = cms.untracked.double(0.88),
s4s9Cut_barrel2_Pi0_ = cms.untracked.double(0.90),
s4s9Cut_endcap1_Pi0_ = cms.untracked.double(0.85),
s4s9Cut_endcap2_Pi0_ = cms.untracked.double(0.92),
nxtal1Cut_barrel1_Pi0_ = cms.untracked.double(0.0),
nxtal1Cut_barrel2_Pi0_ = cms.untracked.double(0.0),
nxtal1Cut_endcap1_Pi0_ = cms.untracked.double(0.0),
nxtal1Cut_endcap2_Pi0_ = cms.untracked.double(0.0),
nxtal2Cut_barrel1_Pi0_ = cms.untracked.double(0.0),
nxtal2Cut_barrel2_Pi0_ = cms.untracked.double(0.0),
nxtal2Cut_endcap1_Pi0_ = cms.untracked.double(0.0),
nxtal2Cut_endcap2_Pi0_ = cms.untracked.double(0.0),
EB_Seed_E_Eta_ = cms.untracked.double(0.5),
EE_Seed_E_Eta_ = cms.untracked.double(0.5),
pairPtCut_barrel1_Eta_ = cms.untracked.double(3.0),
pairPtCut_barrel2_Eta_ = cms.untracked.double(3.0),
pairPtCut_endcap1_Eta_ = cms.untracked.double(3.0),
pairPtCut_endcap2_Eta_ = cms.untracked.double(3.0),
gPtCut_barrel1_Eta_ = cms.untracked.double(0.65),
gPtCut_barrel2_Eta_ = cms.untracked.double(1.4),
gPtCut_endcap1_Eta_ = cms.untracked.double(1.0),
gPtCut_endcap2_Eta_ = cms.untracked.double(1.0),
s4s9Cut_barrel1_Eta_ = cms.untracked.double(0.9),
s4s9Cut_barrel2_Eta_ = cms.untracked.double(0.9),
s4s9Cut_endcap1_Eta_ = cms.untracked.double(0.9),
s4s9Cut_endcap2_Eta_ = cms.untracked.double(0.9),
nxtal1Cut_barrel1_Eta_ = cms.untracked.double(0.0),
nxtal1Cut_barrel2_Eta_ = cms.untracked.double(0.0),
nxtal1Cut_endcap1_Eta_ = cms.untracked.double(0.0),
nxtal1Cut_endcap2_Eta_ = cms.untracked.double(0.0),
nxtal2Cut_barrel1_Eta_ = cms.untracked.double(0.0),
nxtal2Cut_barrel2_Eta_ = cms.untracked.double(0.0),
nxtal2Cut_endcap1_Eta_ = cms.untracked.double(0.0),
nxtal2Cut_endcap2_Eta_ = cms.untracked.double(0.0),
isoGammaBeltdR_Zone_Pi0_        = cms.untracked.double(0.2),
isoGammaBeltdR_Zone_Eta_        = cms.untracked.double(0.2),
isoPairBeltdR_Zone_Pi0_         = cms.untracked.double(0.2),
isoPairBeltdR_Zone_Eta_         = cms.untracked.double(0.3),
isoGammaBeltdEta_Zone_Pi0_      = cms.untracked.double(0.05),
isoGammaBeltdEta_Zone_Eta_      = cms.untracked.double(0.1),
isoPairBeltdEta_Zone_Pi0_       = cms.untracked.double(0.05),
isoPairBeltdEta_Zone_Eta_       = cms.untracked.double(0.1),
isoPairCut_                     = cms.untracked.double(999.),
isoGammaCut_                    = cms.untracked.double(999.)
)

#########################paratmeters for the tuplizer##############################
if isMC_:
        process.ntuples.isMC = cms.untracked.bool(True)
else:
        process.ntuples.isMC = cms.untracked.bool(False)
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
process.p *= process.ecalDigis*process.ecalPreshowerDigis
process.p *= process.bunchSpacingProducer*process.ecalLocalRecoSequence
process.p *= process.ntuples
