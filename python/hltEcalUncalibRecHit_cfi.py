import FWCore.ParameterSet.Config as cms

from RecoLocalCalo.EcalRecProducers.ecalPulseShapeParameters_cff import *

hltEcalUncalibRecHit = cms.EDProducer( "EcalUncalibRecHitProducer",
  EEdigiCollection = cms.InputTag( "hltEcalDigis:eeDigis" ),
  EBdigiCollection = cms.InputTag( "hltEcalDigis:ebDigis" ),
  EEhitCollection = cms.string( "EcalUncalibRecHitsEE" ),
  EBhitCollection = cms.string( "EcalUncalibRecHitsEB" ),
  algo = cms.string( "EcalUncalibRecHitWorkerMultiFit" ),
  algoPSet = cms.PSet(
    outOfTimeThresholdGain61pEB = cms.double( 5.0 ),
    EBtimeFitParameters = cms.vdouble( -2.015452, 3.130702, -12.3473, 41.88921, -82.83944, 91.01147, -50.35761, 11.05621 ),
    activeBXs = cms.vint32( -5, -4, -3, -2, -1, 0, 1, 2 ),
    amplitudeThresholdEE = cms.double( 10.0 ),
    EBtimeConstantTerm = cms.double( 0.6 ), 
    EEtimeFitLimits_Lower = cms.double( 0.2 ),
    outOfTimeThresholdGain61pEE = cms.double( 1000.0 ),
    ebSpikeThreshold = cms.double( 1.042 ),
    EBtimeNconst = cms.double( 28.5 ),
    ampErrorCalculation = cms.bool( False ),
    kPoorRecoFlagEB = cms.bool( True ), 
    EBtimeFitLimits_Lower = cms.double( 0.2 ),
    kPoorRecoFlagEE = cms.bool( False ),
    chi2ThreshEB_ = cms.double( 65.0 ),
    EEtimeFitParameters = cms.vdouble( -2.390548, 3.553628, -17.62341, 67.67538, -133.213, 140.7432, -75.41106, 16.20277 ),
    useLumiInfoRunHeader = cms.bool( False ),
    outOfTimeThresholdGain12mEE = cms.double( 1000.0 ),
    outOfTimeThresholdGain12mEB = cms.double( 5.0 ),
    EEtimeFitLimits_Upper = cms.double( 1.4 ),
    prefitMaxChiSqEB = cms.double( 15.0 ),
    EEamplitudeFitParameters = cms.vdouble( 1.89, 1.4 ),
    prefitMaxChiSqEE = cms.double( 10.0 ),
    EBamplitudeFitParameters = cms.vdouble( 1.138, 1.652 ),
    EBtimeFitLimits_Upper = cms.double( 1.4 ),
    timealgo = cms.string( "None" ),
    amplitudeThresholdEB = cms.double( 10.0 ),
    outOfTimeThresholdGain12pEE = cms.double( 1000.0 ),
    outOfTimeThresholdGain12pEB = cms.double( 5.0 ),
    EEtimeNconst = cms.double( 31.8 ),
    outOfTimeThresholdGain61mEB = cms.double( 5.0 ),
    outOfTimeThresholdGain61mEE = cms.double( 1000.0 ),
    EEtimeConstantTerm = cms.double( 1.0 ),
    chi2ThreshEE_ = cms.double( 50.0 ),
    doPrefitEE = cms.bool( True ),
    doPrefitEB = cms.bool( True )
  )
)

