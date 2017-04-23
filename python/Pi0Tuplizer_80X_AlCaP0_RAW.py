#########################options##############################
isMC_ = False
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
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )
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
EBRecHitCollectionTag_Eta = cms.untracked.InputTag("ecalRecHitEta","EcalRecHitsEB","Pi0Tuplizer"),
EERecHitCollectionTag_Eta = cms.untracked.InputTag("ecalRecHitEta","EcalRecHitsEE","Pi0Tuplizer"),
ESRecHitCollectionTag_Eta = cms.untracked.InputTag('dummyHitsEta','dummyPreshowerRechits',"Pi0Tuplizer"),
EBRecHitCollectionTag_Pi0 = cms.untracked.InputTag("ecalRecHitPi0","EcalRecHitsEB","Pi0Tuplizer"),
EERecHitCollectionTag_Pi0 = cms.untracked.InputTag("ecalRecHitPi0","EcalRecHitsEE","Pi0Tuplizer"),
ESRecHitCollectionTag_Pi0 = cms.untracked.InputTag('dummyHitsPi0','dummyPreshowerRechits',"Pi0Tuplizer"),
PhotonOrderOption = cms.untracked.string("SeedEBased"),# "SeedEBased" (g1 is the one with larger seed Energy) or "PhoPtBased"(g1 is the one with larger pt)
EB_Seed_E_Pi0_ = cms.untracked.double(0.0),#(0.5),
EE_Seed_E_Pi0_ = cms.untracked.double(0.0),#(0.5),
pairPtCut_barrel1_Pi0_ = cms.untracked.double(2.6),
pairPtCut_barrel2_Pi0_ = cms.untracked.double(2.6),
pairPtCut_endcap1_Pi0_ = cms.untracked.double(3.0),
pairPtCut_endcap2_Pi0_ = cms.untracked.double(1.5),
gPtCut_barrel1_Pi0_ = cms.untracked.double(0.0),#(1.3),
gPtCut_barrel2_Pi0_ = cms.untracked.double(0.0),#(1.3),
gPtCut_endcap1_Pi0_ = cms.untracked.double(0.0),#(0.95),
gPtCut_endcap2_Pi0_ = cms.untracked.double(0.0),#(0.65),
s4s9Cut_barrel1_Pi0_ = cms.untracked.double(0.0),#(0.83),
s4s9Cut_barrel2_Pi0_ = cms.untracked.double(0.0),#(0.83),
s4s9Cut_endcap1_Pi0_ = cms.untracked.double(0.0),#(0.95),
s4s9Cut_endcap2_Pi0_ = cms.untracked.double(0.0),#(0.95),
nxtal1Cut_barrel1_Pi0_ = cms.untracked.double(0.),
nxtal1Cut_barrel2_Pi0_ = cms.untracked.double(0.),
nxtal1Cut_endcap1_Pi0_ = cms.untracked.double(0.),
nxtal1Cut_endcap2_Pi0_ = cms.untracked.double(0.),
nxtal2Cut_barrel1_Pi0_ = cms.untracked.double(0.),
nxtal2Cut_barrel2_Pi0_ = cms.untracked.double(0.),
nxtal2Cut_endcap1_Pi0_ = cms.untracked.double(0.),
nxtal2Cut_endcap2_Pi0_ = cms.untracked.double(0.),
EB_Seed_E_Eta_ = cms.untracked.double(0.0),#(0.5),
EE_Seed_E_Eta_ = cms.untracked.double(0.0),#(0.5),
pairPtCut_barrel1_Eta_ = cms.untracked.double(4.0),
pairPtCut_barrel2_Eta_ = cms.untracked.double(4.0),
pairPtCut_endcap1_Eta_ = cms.untracked.double(3.0),
pairPtCut_endcap2_Eta_ = cms.untracked.double(3.0),
gPtCut_barrel1_Eta_ = cms.untracked.double(0.0),#(1.2),
gPtCut_barrel2_Eta_ = cms.untracked.double(0.0),#(1.2),
gPtCut_endcap1_Eta_ = cms.untracked.double(0.0),#(1.0),
gPtCut_endcap2_Eta_ = cms.untracked.double(0.0),#(0.70),
s4s9Cut_barrel1_Eta_ = cms.untracked.double(0.0),#(0.87),
s4s9Cut_barrel2_Eta_ = cms.untracked.double(0.0),#(0.87),
s4s9Cut_endcap1_Eta_ = cms.untracked.double(0.0),#(0.90),
s4s9Cut_endcap2_Eta_ = cms.untracked.double(0.0),#(0.90),
nxtal1Cut_barrel1_Eta_ = cms.untracked.double(0.),
nxtal1Cut_barrel2_Eta_ = cms.untracked.double(0.),
nxtal1Cut_endcap1_Eta_ = cms.untracked.double(0.),
nxtal1Cut_endcap2_Eta_ = cms.untracked.double(0.),
nxtal2Cut_barrel1_Eta_ = cms.untracked.double(0.),
nxtal2Cut_barrel2_Eta_ = cms.untracked.double(0.),
nxtal2Cut_endcap1_Eta_ = cms.untracked.double(0.),
nxtal2Cut_endcap2_Eta_ = cms.untracked.double(0.),
isoGammaBeltdR_Zone_Pi0_	= cms.untracked.double(0.2),
isoGammaBeltdR_Zone_Eta_ 	= cms.untracked.double(0.2),
isoPairBeltdR_Zone_Pi0_		= cms.untracked.double(0.2),
isoPairBeltdR_Zone_Eta_		= cms.untracked.double(0.3),
isoGammaBeltdEta_Zone_Pi0_ 	= cms.untracked.double(0.05),
isoGammaBeltdEta_Zone_Eta_ 	= cms.untracked.double(0.1),
isoPairBeltdEta_Zone_Pi0_  	= cms.untracked.double(0.05),
isoPairBeltdEta_Zone_Eta_  	= cms.untracked.double(0.1)
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

#DUMMY RECHIT
process.dummyHitsPi0 = cms.EDProducer('DummyRechitDigis',
                                     doDigi = cms.untracked.bool(False),
                                     # rechits
                                     barrelHitProducer      = cms.InputTag('',''),
                                     endcapHitProducer      = cms.InputTag('',''),
                                     preshowerHitProducer      = cms.InputTag("hltAlCaPi0RecHitsFilterEEonlyRegional", "pi0EcalRecHitsES"),
                                     barrelRecHitCollection = cms.untracked.string('dummyBarrelRechits'),
                                     endcapRecHitCollection = cms.untracked.string('dummyEndcapRechits'),
                                     preshowerRecHitCollection = cms.untracked.string('dummyPreshowerRechits'),
                                     # digis
                                     barrelDigis            = cms.InputTag("hltAlCaPi0EBRechitsToDigis","pi0EBDigis"),
                                     endcapDigis            = cms.InputTag("hltAlCaPi0EERechitsToDigis","pi0EEDigis"),
                                     barrelDigiCollection   = cms.untracked.string('dummyBarrelDigis'),
                                     endcapDigiCollection   = cms.untracked.string('dummyEndcapDigis'))

process.dummyHitsEta = cms.EDProducer('DummyRechitDigis',
                                     doDigi = cms.untracked.bool(False),
                                     # rechits
                                     barrelHitProducer      = cms.InputTag('',''),
                                     endcapHitProducer      = cms.InputTag('',''),
                                     preshowerHitProducer      = cms.InputTag("hltAlCaEtaRecHitsFilterEEonlyRegional", "etaEcalRecHitsES"),
                                     barrelRecHitCollection = cms.untracked.string('dummyBarrelRechits'),
                                     endcapRecHitCollection = cms.untracked.string('dummyEndcapRechits'),
                                     preshowerRecHitCollection = cms.untracked.string('dummyPreshowerRechits'),
                                     # digis
                                     barrelDigis            = cms.InputTag("hltAlCaEtaEBRechitsToDigis","etaEBDigis"),
                                     endcapDigis            = cms.InputTag("hltAlCaEtaEERechitsToDigis","etaEEDigis"),
                                     barrelDigiCollection   = cms.untracked.string('dummyBarrelDigis'),
                                     endcapDigiCollection   = cms.untracked.string('dummyEndcapDigis'))


#DUMMY DIGIS
process.dummyDigisPi0 = cms.EDProducer('DummyRechitDigis',
                                     doDigi = cms.untracked.bool(True),
                                     # rechits
                                     barrelHitProducer      = cms.InputTag('hltAlCaPi0EBUncalibrator','pi0EcalRecHitsEB'),
                                     endcapHitProducer      = cms.InputTag('hltAlCaPi0EEUncalibrator','pi0EcalRecHitsEE'),
				     preshowerHitProducer      = cms.InputTag('',''),
                                     barrelRecHitCollection = cms.untracked.string('dummyBarrelRechits'),
                                     endcapRecHitCollection = cms.untracked.string('dummyEndcapRechits'),
				     preshowerRecHitCollection = cms.untracked.string('dummyPreshowerRechits'),
                                     # digis
                                     barrelDigis            = cms.InputTag("hltAlCaPi0EBRechitsToDigis","pi0EBDigis"),
                                     endcapDigis            = cms.InputTag("hltAlCaPi0EERechitsToDigis","pi0EEDigis"),
                                     barrelDigiCollection   = cms.untracked.string('dummyBarrelDigis'),
                                     endcapDigiCollection   = cms.untracked.string('dummyEndcapDigis'))
process.dummyDigisEta = cms.EDProducer('DummyRechitDigis',
                                     doDigi = cms.untracked.bool(True),
                                     # rechits
                                     barrelHitProducer      = cms.InputTag('hltAlCaEtaEBUncalibrator','etaEcalRecHitsEB'),
                                     endcapHitProducer      = cms.InputTag('hltAlCaEtaEEUncalibrator','etaEcalRecHitsEE'),
				     preshowerHitProducer      = cms.InputTag('',''),
                                     barrelRecHitCollection = cms.untracked.string('dummyBarrelRechits'),
                                     endcapRecHitCollection = cms.untracked.string('dummyEndcapRechits'),
				     preshowerRecHitCollection = cms.untracked.string('dummyPreshowerRechits'),
                                     # digis
                                     barrelDigis            = cms.InputTag("hltAlCaEtaEBRechitsToDigis","etaEBDigis"),
                                     endcapDigis            = cms.InputTag("hltAlCaEtaEERechitsToDigis","etaEEDigis"),
                                     barrelDigiCollection   = cms.untracked.string('dummyBarrelDigis'),
                                     endcapDigiCollection   = cms.untracked.string('dummyEndcapDigis'))


#DIGI to UNCALIB
import RecoLocalCalo.EcalRecProducers.ecalMultiFitUncalibRecHit_cfi

import hltEcalUncalibRecHit_cfi

#process.ecalMultiFitUncalibRecHitPi0 = RecoLocalCalo.EcalRecProducers.ecalMultiFitUncalibRecHit_cfi.ecalMultiFitUncalibRecHit.clone()
process.ecalMultiFitUncalibRecHitPi0 = hltEcalUncalibRecHit_cfi.hltEcalUncalibRecHit.clone()
process.ecalMultiFitUncalibRecHitPi0.EBdigiCollection = cms.InputTag('dummyDigisPi0','dummyBarrelDigis','Pi0Tuplizer')
process.ecalMultiFitUncalibRecHitPi0.EEdigiCollection = cms.InputTag('dummyDigisPi0','dummyEndcapDigis','Pi0Tuplizer')
process.ecalMultiFitUncalibRecHitPi0.algoPSet.useLumiInfoRunHeader = cms.bool( False )

#process.ecalMultiFitUncalibRecHitEta = RecoLocalCalo.EcalRecProducers.ecalMultiFitUncalibRecHit_cfi.ecalMultiFitUncalibRecHit.clone()
process.ecalMultiFitUncalibRecHitEta = hltEcalUncalibRecHit_cfi.hltEcalUncalibRecHit.clone()
process.ecalMultiFitUncalibRecHitEta.EBdigiCollection = cms.InputTag('dummyDigisEta','dummyBarrelDigis','Pi0Tuplizer')
process.ecalMultiFitUncalibRecHitEta.EEdigiCollection = cms.InputTag('dummyDigisEta','dummyEndcapDigis','Pi0Tuplizer')
process.ecalMultiFitUncalibRecHitEta.algoPSet.useLumiInfoRunHeader = cms.bool( False )


#UNCALIB to CALIB
from RecoLocalCalo.EcalRecProducers.ecalRecHit_cfi import *
process.ecalDetIdToBeRecovered =  RecoLocalCalo.EcalRecProducers.ecalDetIdToBeRecovered_cfi.ecalDetIdToBeRecovered.clone()

process.ecalRecHitPi0 = cms.EDProducer("EcalRecHitProducer",
    				EEuncalibRecHitCollection = cms.InputTag("ecalMultiFitUncalibRecHitPi0","EcalUncalibRecHitsEE"),
    				EBuncalibRecHitCollection = cms.InputTag("ecalMultiFitUncalibRecHitPi0","EcalUncalibRecHitsEB"),
    				EBrechitCollection = cms.string('EcalRecHitsEB'),
				EErechitCollection = cms.string('EcalRecHitsEE'),
				killDeadChannels = cms.bool( False ),
				recoverEBVFE = cms.bool( False ),
				recoverEEVFE = cms.bool( False ),
				recoverEBFE = cms.bool( False ),
				recoverEEFE = cms.bool( False ),
				recoverEEIsolatedChannels = cms.bool( False ),
				recoverEBIsolatedChannels = cms.bool( False )
				)

process.ecalRecHitEta = cms.EDProducer("EcalRecHitProducer",
    				EEuncalibRecHitCollection = cms.InputTag("ecalMultiFitUncalibRecHitEta","EcalUncalibRecHitsEE"),
    				EBuncalibRecHitCollection = cms.InputTag("ecalMultiFitUncalibRecHitEta","EcalUncalibRecHitsEB"),
    				EBrechitCollection = cms.string('EcalRecHitsEB'),
				EErechitCollection = cms.string('EcalRecHitsEE'),
				killDeadChannels = cms.bool( False ),
				recoverEBVFE = cms.bool( False ),
				recoverEEVFE = cms.bool( False ),
				recoverEBFE = cms.bool( False ),
				recoverEEFE = cms.bool( False ),
				recoverEEIsolatedChannels = cms.bool( False ),
				recoverEBIsolatedChannels = cms.bool( False )
				)



process.ecalLocalRecoSequence = cms.Sequence(process.ecalRecHitPi0+process.ecalRecHitEta)

#define path
process.p = cms.Path()
process.p *= process.dummyDigisPi0
process.p *= process.dummyDigisEta
process.p *= process.ecalMultiFitUncalibRecHitPi0
process.p *= process.ecalMultiFitUncalibRecHitEta
process.p *= process.ecalLocalRecoSequence
process.p *= process.dummyHitsPi0
process.p *= process.dummyHitsEta
process.p *= process.ntuples

#process.end = cms.EndPath(process.Out)
