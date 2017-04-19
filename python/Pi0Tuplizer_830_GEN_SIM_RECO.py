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

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )

process.MessageLogger.cerr.FwkReport.reportEvery = 10000
process.options = cms.untracked.PSet(
   wantSummary = cms.untracked.bool(True),
   SkipEvent = cms.untracked.vstring('ProductNotFound','CrystalIDError')
)

#load input file
process.source = cms.Source('PoolSource',
    fileNames = cms.untracked.vstring(
        '/store/mc/PhaseIFall16DR/SinglePion_FlatPt-1To15/GEN-SIM-RECO/FlatPU10to50RECO_90X_upgrade2017_realistic_v6_C1-v1/510000/42C9D4EC-A413-E711-8E3B-FA163E66742D.root'
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
PhotonOrderOption = cms.untracked.string("SeedEBased"),# "SeedEBased" (g1 is the one with larger seed Energy) or "PhoPtBased"(g1 is the one with larger pt)
EB_Seed_E_Pi0_ = cms.untracked.double(0.5),
EE_Seed_E_Pi0_ = cms.untracked.double(0.5),
pairPtCut_barrel1_Pi0_ = cms.untracked.double(2.6),
pairPtCut_barrel2_Pi0_ = cms.untracked.double(2.6),
pairPtCut_endcap1_Pi0_ = cms.untracked.double(3.0),
pairPtCut_endcap2_Pi0_ = cms.untracked.double(1.5),
gPtCut_barrel1_Pi0_ = cms.untracked.double(1.3),
gPtCut_barrel2_Pi0_ = cms.untracked.double(1.3),
gPtCut_endcap1_Pi0_ = cms.untracked.double(0.95),
gPtCut_endcap2_Pi0_ = cms.untracked.double(0.65),
s4s9Cut_barrel1_Pi0_ = cms.untracked.double(0.83),
s4s9Cut_barrel2_Pi0_ = cms.untracked.double(0.83),
s4s9Cut_endcap1_Pi0_ = cms.untracked.double(0.95),
s4s9Cut_endcap2_Pi0_ = cms.untracked.double(0.95),
nxtal1Cut_barrel1_Pi0_ = cms.untracked.double(0.),
nxtal1Cut_barrel2_Pi0_ = cms.untracked.double(0.),
nxtal1Cut_endcap1_Pi0_ = cms.untracked.double(0.),
nxtal1Cut_endcap2_Pi0_ = cms.untracked.double(0.),
nxtal2Cut_barrel1_Pi0_ = cms.untracked.double(0.),
nxtal2Cut_barrel2_Pi0_ = cms.untracked.double(0.),
nxtal2Cut_endcap1_Pi0_ = cms.untracked.double(0.),
nxtal2Cut_endcap2_Pi0_ = cms.untracked.double(0.),
EB_Seed_E_Eta_ = cms.untracked.double(0.5),
EE_Seed_E_Eta_ = cms.untracked.double(0.5),
pairPtCut_barrel1_Eta_ = cms.untracked.double(4.0),
pairPtCut_barrel2_Eta_ = cms.untracked.double(4.0),
pairPtCut_endcap1_Eta_ = cms.untracked.double(3.0),
pairPtCut_endcap2_Eta_ = cms.untracked.double(3.0),
gPtCut_barrel1_Eta_ = cms.untracked.double(1.2),
gPtCut_barrel2_Eta_ = cms.untracked.double(1.2),
gPtCut_endcap1_Eta_ = cms.untracked.double(1.0),
gPtCut_endcap2_Eta_ = cms.untracked.double(0.70),
s4s9Cut_barrel1_Eta_ = cms.untracked.double(0.87),
s4s9Cut_barrel2_Eta_ = cms.untracked.double(0.87),
s4s9Cut_endcap1_Eta_ = cms.untracked.double(0.90),
s4s9Cut_endcap2_Eta_ = cms.untracked.double(0.90),
nxtal1Cut_barrel1_Eta_ = cms.untracked.double(0.),
nxtal1Cut_barrel2_Eta_ = cms.untracked.double(0.),
nxtal1Cut_endcap1_Eta_ = cms.untracked.double(0.),
nxtal1Cut_endcap2_Eta_ = cms.untracked.double(0.),
nxtal2Cut_barrel1_Eta_ = cms.untracked.double(0.),
nxtal2Cut_barrel2_Eta_ = cms.untracked.double(0.),
nxtal2Cut_endcap1_Eta_ = cms.untracked.double(0.),
nxtal2Cut_endcap2_Eta_ = cms.untracked.double(0.)
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
process.p *= process.ntuples
