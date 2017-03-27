
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

process.MessageLogger.cerr.FwkReport.reportEvery = 10
process.options = cms.untracked.PSet(
   wantSummary = cms.untracked.bool(True),
   SkipEvent = cms.untracked.vstring('ProductNotFound','CrystalIDError')
)

#load input file
process.source = cms.Source('PoolSource',
    fileNames = cms.untracked.vstring(
        'root://cmsxrootd.fnal.gov//store/mc/PhaseIFall16DR/QCD_Pt-15to20_EMEnriched_TuneCUETP8M1_13TeV_pythia8/GEN-SIM-RAW/FlatPU28to62HcalNZSRAW_81X_upgrade2017_realistic_v26-v1/100000/002CF4F7-A4F0-E611-8C1F-141877411EA2.root'
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
uGtAlgInputTag = cms.untracked.InputTag('hltGtStage2Digis')
)

#define path
process.p = cms.Path()
process.p *= process.ecalDigis*process.ecalPreshowerDigis
process.p *= process.bunchSpacingProducer*process.ecalLocalRecoSequence
process.p *= process.ntuples
