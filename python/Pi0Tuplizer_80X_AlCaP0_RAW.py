#########################options##############################
isMC_ = False
isPi0_ = True
FillL1SeedFinalDecision_ = True
FillDiPhotonNtuple_ = True
FillPhotonNtuple_ = True
#########################options##############################
import FWCore.ParameterSet.Config as cms
import RecoLocalCalo.EcalRecProducers.ecalRecalibRecHit_cfi
import os, sys, imp, re
import RecoLocalCalo.EcalRecProducers.ecalMultiFitUncalibRecHit_cfi
from RecoLocalCalo.EcalRecProducers.ecalRecHit_cfi import *


CMSSW_VERSION=os.getenv("CMSSW_VERSION")
process = cms.Process("Pi0Tuplizer")


#load modules
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
process.load('Configuration.StandardSequences.Reconstruction_cff')

process.GlobalTag.globaltag = '80X_dataRun2_2016LegacyRepro_Candidate_v2'
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.MessageLogger.cerr.FwkReport.reportEvery = 10000
process.options = cms.untracked.PSet(
   wantSummary = cms.untracked.bool(True),
   SkipEvent = cms.untracked.vstring('ProductNotFound','CrystalIDError')
)

#load input file
process.source = cms.Source('PoolSource',
    fileNames = cms.untracked.vstring(
        '/store/data/Run2016B/AlCaP0/RAW/v2/000/273/158/00000/06BD06B1-1618-E611-ADA7-02163E014286.root'
    )
)

#process.Out = cms.OutputModule("PoolOutputModule",
#         fileName = cms.untracked.string ("MyOutputFile.root")
#)

#define output file
process.TFileService = cms.Service("TFileService", 
    fileName = cms.string("pi0Ntuple.root"),
    closeFileFast = cms.untracked.bool(True)
)

#########################paratmeters for the tuplizer##############################
process.ntuples = cms.EDAnalyzer('Pi0Tuplizer',
uGtAlgInputTag = cms.untracked.InputTag('hltGtStage2Digis'),
EBRecHitCollectionTag = cms.untracked.InputTag("ecalRecHit","EcalRecHitsEB","Pi0Tuplizer"),
EERecHitCollectionTag = cms.untracked.InputTag("ecalRecHit","EcalRecHitsEE","Pi0Tuplizer"),
ESRecHitCollectionTag = cms.untracked.InputTag("hltAlCaPi0RecHitsFilterEEonlyRegional", "pi0EcalRecHitsES"),
PhotonOrderOption = cms.untracked.string("SeedEBased"),# "SeedEBased" (g1 is the one with larger seed Energy) or "PhoPtBased"(g1 is the one with larger pt)
EB_Seed_E = cms.untracked.double(0.5),
EE_Seed_E = cms.untracked.double(0.5),
pi0PtCut_barrel1 = cms.untracked.double(2.6),
pi0PtCut_barrel2 = cms.untracked.double(2.6),
pi0PtCut_endcap1 = cms.untracked.double(3.0),
pi0PtCut_endcap2 = cms.untracked.double(1.5),
gPtCut_barrel1 = cms.untracked.double(1.3),
gPtCut_barrel2 = cms.untracked.double(1.3),
gPtCut_endcap1 = cms.untracked.double(0.95),
gPtCut_endcap2 = cms.untracked.double(0.65),
s4s9Cut_barrel1 = cms.untracked.double(0.83),
s4s9Cut_barrel2 = cms.untracked.double(0.83),
s4s9Cut_endcap1 = cms.untracked.double(0.95),
s4s9Cut_endcap2 = cms.untracked.double(0.95),
nxtal1Cut_barrel1 = cms.untracked.double(0.),
nxtal1Cut_barrel2 = cms.untracked.double(0.),
nxtal1Cut_endcap1 = cms.untracked.double(0.),
nxtal1Cut_endcap2 = cms.untracked.double(0.),
nxtal2Cut_barrel1 = cms.untracked.double(0.),
nxtal2Cut_barrel2 = cms.untracked.double(0.),
nxtal2Cut_endcap1 = cms.untracked.double(0.),
nxtal2Cut_endcap2 = cms.untracked.double(0.)
)
#########################paratmeters for the tuplizer##############################
if isMC_:
	process.ntuples.isMC = cms.untracked.bool(True)
else:
	process.ntuples.isMC = cms.untracked.bool(False)
if isPi0_:
	process.ntuples.isPi0 = cms.untracked.bool(True)
else:
	process.ntuples.isPi0 = cms.untracked.bool(False)
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

#DUMMY RECHIT
process.dummyHits = cms.EDProducer('DummyRechitDigis',
                                     doDigi = cms.untracked.bool(True),
                                     # rechits
                                     barrelHitProducer      = cms.InputTag('hltAlCaPi0EBUncalibrator','pi0EcalRecHitsEB'),
                                     endcapHitProducer      = cms.InputTag('hltAlCaPi0EEUncalibrator','pi0EcalRecHitsEE'),
                                     barrelRecHitCollection = cms.untracked.string('dummyBarrelRechits'),
                                     endcapRecHitCollection = cms.untracked.string('dummyEndcapRechits'),
                                     # digis
                                     barrelDigis            = cms.InputTag("hltAlCaPi0EBRechitsToDigis","pi0EBDigis"),
                                     endcapDigis            = cms.InputTag("hltAlCaPi0EERechitsToDigis","pi0EEDigis"),
                                     barrelDigiCollection   = cms.untracked.string('dummyBarrelDigis'),
                                     endcapDigiCollection   = cms.untracked.string('dummyEndcapDigis'))

#DIGI to UNCALIB
import RecoLocalCalo.EcalRecProducers.ecalMultiFitUncalibRecHit_cfi
process.ecalMultiFitUncalibRecHit = RecoLocalCalo.EcalRecProducers.ecalMultiFitUncalibRecHit_cfi.ecalMultiFitUncalibRecHit.clone()
process.ecalMultiFitUncalibRecHit.EBdigiCollection = cms.InputTag('dummyHits','dummyBarrelDigis','Pi0Tuplizer')
process.ecalMultiFitUncalibRecHit.EEdigiCollection = cms.InputTag('dummyHits','dummyEndcapDigis','Pi0Tuplizer')
process.ecalMultiFitUncalibRecHit.algoPSet.useLumiInfoRunHeader = cms.bool( False )

#UNCALIB to CALIB
from RecoLocalCalo.EcalRecProducers.ecalRecHit_cfi import *
process.ecalDetIdToBeRecovered =  RecoLocalCalo.EcalRecProducers.ecalDetIdToBeRecovered_cfi.ecalDetIdToBeRecovered.clone()
process.ecalRecHit.killDeadChannels = cms.bool( False )
process.ecalRecHit.recoverEBVFE = cms.bool( False )
process.ecalRecHit.recoverEEVFE = cms.bool( False )
process.ecalRecHit.recoverEBFE = cms.bool( False )
process.ecalRecHit.recoverEEFE = cms.bool( False )
process.ecalRecHit.recoverEEIsolatedChannels = cms.bool( False )
process.ecalRecHit.recoverEBIsolatedChannels = cms.bool( False )
process.ecalLocalRecoSequence = cms.Sequence(ecalRecHit)

#define path
process.p = cms.Path()
process.p *= process.dummyHits
process.p *= process.ecalMultiFitUncalibRecHit
process.p *= process.ecalLocalRecoSequence
process.p *= process.ntuples

#process.end = cms.EndPath(process.Out)
