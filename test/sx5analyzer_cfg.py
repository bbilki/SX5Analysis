# -*- coding: utf-8 -*-
import FWCore.ParameterSet.Config as cms

process = cms.Process("SX5Analysis")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

import sys
runNumber = sys.argv[2]
runType = int(sys.argv[3])

process.source = cms.Source("HcalTBSource",
    fileNames = cms.untracked.vstring(
        'file:./Data/SX5_'+runNumber+'.root'
    )
)

process.options = cms.untracked.PSet(
        wantSummary = cms.untracked.bool(False)
        )

process.tbunpack = cms.EDProducer("HcalTBObjectUnpacker",
        IncludeUnmatchedHits = cms.untracked.bool(False),
	ConfigurationFile=cms.untracked.string('HFcommissioning/SX5Analysis/test/configQADCTDC.txt'),
        HcalSlowDataFED = cms.untracked.int32(14),
        HcalTriggerFED = cms.untracked.int32(1),
        HcalTDCFED = cms.untracked.int32(8),
        HcalQADCFED = cms.untracked.int32(8),
        fedRawDataCollectionTag = cms.InputTag('source')
)

process.hcalDigis = cms.EDProducer("HcalRawToDigi",
                                   #       UnpackHF = cms.untracked.bool(True),
                                   ### Falg to enable unpacking of TTP channels(default = false)
                                   ### UnpackTTP = cms.untracked.bool(True),
                                   FilterDataQuality = cms.bool(False),
                                   InputLabel = cms.InputTag('source'),
                                   #HcalFirstFED = cms.untracked.int32(928),
                                   HcalFirstFED = cms.untracked.int32(938),
                                   #ComplainEmptyData = cms.untracked.bool(False),
                                   ComplainEmptyData = cms.untracked.bool(True),
                                   #       UnpackCalib = cms.untracked.bool(True),
                                   #FEDs = cms.untracked.vint32(928,930,932),
                                   #FEDs = cms.untracked.vint32(928,930,932,938),
                                   #FEDs = cms.untracked.vint32(928),
                                   FEDs = cms.untracked.vint32(938),
                                   firstSample = cms.int32(0),
                                   lastSample = cms.int32(10)
                                   )


process.hcalAnalyzer = cms.EDAnalyzer('SX5Analysis',
        OutFileName = cms.untracked.string('N_'+runNumber+'.root'),
        RunType = cms.int32(runType),
	histoFED =  cms.int32(74)
)

process.load('Configuration.Geometry.GeometryIdeal_cff')

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.autoCond import autoCond
#from CondCore.DBCommon.CondDBSetup_cfi import *

from CondCore.CondDB.CondDB_cfi import *

process.GlobalTag.globaltag = autoCond['startup'] 

#   EMAP Needed for H2 DATA
process.es_ascii = cms.ESSource('HcalTextCalibrations',
        input = cms.VPSet(
               cms.PSet(
                object = cms.string('ElectronicsMap'),
		file = cms.FileInPath('HFcommissioning/SX5Analysis/test/emap_SX5_v3.txt')
               )
        )
)

process.es_prefer = cms.ESPrefer('HcalTextCalibrations', 'es_ascii')
process.p = cms.Path(process.hcalDigis*process.hcalAnalyzer)


