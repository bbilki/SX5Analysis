// -*- C++ -*-
//
// Package:    HFcommissioning/SX5Analysis
// Class:      SX5Analysis
// 
/**\class SX5Analysis SX5Analysis.cc HFcommissioning/SX5Analysis/plugins/SX5Analysis.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Burak Bilki
//         Created:  Sun, 04 Dec 2016 12:10:07 GMT
//
//

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "EventFilter/HcalRawToDigi/interface/HcalHTRData.h"
#include "EventFilter/HcalRawToDigi/interface/HcalDCCHeader.h"
#include "EventFilter/HcalRawToDigi/interface/HcalUnpacker.h"
#include "DataFormats/HcalDetId/interface/HcalOtherDetId.h"
#include "DataFormats/HcalDigi/interface/HcalQIESample.h"
#include "DataFormats/HcalDigi/interface/QIE10DataFrame.h"
#include "DataFormats/HcalDetId/interface/HcalSubdetector.h"
#include "DataFormats/HcalDetId/interface/HcalCalibDetId.h"
#include "EventFilter/HcalRawToDigi/interface/AMC13Header.h"
#include "EventFilter/HcalRawToDigi/interface/HcalUHTRData.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/FEDRawData/interface/FEDRawDataCollection.h"
#include "DataFormats/FEDRawData/interface/FEDHeader.h"
#include "DataFormats/FEDRawData/interface/FEDTrailer.h"
#include "DataFormats/FEDRawData/interface/FEDNumbering.h"
#include "DataFormats/FEDRawData/interface/FEDRawData.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"
#include "CalibFormats/HcalObjects/interface/HcalDbService.h"
#include "CalibFormats/HcalObjects/interface/HcalDbRecord.h"
#include "CalibFormats/HcalObjects/interface/HcalCalibrations.h"
#include "CalibFormats/HcalObjects/interface/HcalCoderDb.h"

#include "TH1D.h"
#include "TH2D.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TProfile.h"
#include "TFile.h"
#include "TSystem.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TStyle.h"

#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;

/// Per Event Header Structure
struct eventHeader
{
	uint32_t cdf0;
	uint32_t cdf1;
	uint32_t cdf2;
	uint32_t cdf3;
	uint32_t h0;
	uint32_t h1;
	uint32_t h2;
	uint32_t h3;
};

struct edata
{
	int ieta[144];
	int iphi[144];
	int depth[144];
	int pulse[144][10];
	int histo[144][60];
};

//station - PMT - 1+2|3+4 - ieta|iphi|depth (or fiber|channel|uHTR)
int MAP2Ch[3][24][2][3]={{{{7,2,1},{10,3,1}},{{10,2,1},{7,3,1}},{{6,0,1},{9,1,1}},{{9,0,1},{6,1,1}},{{8,0,1},{11,1,1}},{{11,0,1},{8,1,1}},{{6,2,1},{9,3,1}},{{9,2,1},{6,3,1}},{{8,2,1},{11,3,1}},{{11,2,1},{8,3,1}},{{7,0,1},{10,1,1}},{{10,0,1},{7,1,1}},{{1,2,1},{4,3,1}},{{4,2,1},{1,3,1}},{{0,0,1},{3,1,1}},{{3,0,1},{0,1,1}},{{2,0,1},{5,1,1}},{{5,0,1},{2,1,1}},{{0,2,1},{3,3,1}},{{3,2,1},{0,3,1}},{{2,2,1},{5,3,1}},{{5,2,1},{2,3,1}},{{1,0,1},{4,1,1}},{{4,0,1},{1,1,1}}},{{{1,2,2},{4,3,2}},{{4,2,2},{1,3,2}},{{0,0,2},{3,1,2}},{{3,0,2},{0,1,2}},{{2,0,2},{5,1,2}},{{5,0,2},{2,1,2}},{{0,2,2},{3,3,2}},{{3,2,2},{0,3,2}},{{2,2,2},{5,3,2}},{{5,2,2},{2,3,2}},{{1,0,2},{4,1,2}},{{4,0,2},{1,1,2}},{{7,2,2},{10,3,2}},{{10,2,2},{7,3,2}},{{6,0,2},{9,1,2}},{{9,0,2},{6,1,2}},{{8,0,2},{11,1,2}},{{11,0,2},{8,1,2}},{{6,2,2},{9,3,2}},{{9,2,2},{6,3,2}},{{8,2,2},{11,3,2}},{{11,2,2},{8,3,2}},{{7,0,2},{10,1,2}},{{10,0,2},{7,1,2}}},{{{19,2,2},{22,3,2}},{{22,2,2},{19,3,2}},{{18,0,2},{21,1,2}},{{21,0,2},{18,1,2}},{{20,0,2},{23,1,2}},{{23,0,2},{20,1,2}},{{18,2,2},{21,3,2}},{{21,2,2},{18,3,2}},{{20,2,2},{23,3,2}},{{23,2,2},{20,3,2}},{{19,0,2},{22,1,2}},{{22,0,2},{19,1,2}},{{13,2,2},{16,3,2}},{{16,2,2},{13,3,2}},{{12,0,2},{15,1,2}},{{15,0,2},{12,1,2}},{{14,0,2},{17,1,2}},{{17,0,2},{14,1,2}},{{12,2,2},{15,3,2}},{{15,2,2},{12,3,2}},{{14,2,2},{17,3,2}},{{17,2,2},{14,3,2}},{{13,0,2},{16,1,2}},{{16,0,2},{13,1,2}}}};

class SX5Analysis : public edm::one::EDAnalyzer<edm::one::SharedResources>
{
	public:
		explicit SX5Analysis(const edm::ParameterSet&);
		~SX5Analysis();
		static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
	private:
		virtual void beginJob() override;
		virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
		virtual void endJob() override;
		TFile *_file;
		TTree *_tree;
		int histoFED;
		int EID;
		int numChannels;
		string outFileName;
		int runType;
		string digiCollection;
		edata ed;
		
		edm::EDGetTokenT<HcalDataFrameContainer<QIE10DataFrame> > tok_QIE10DigiCollection_;
		edm::EDGetTokenT<HFDigiCollection> hf_token;
		edm::EDGetTokenT<FEDRawDataCollection> raw_token;  
		edm::Handle<QIE10DigiCollection> qie10DigiCollection;
		edm::Handle<FEDRawDataCollection> raw_collection;  
};

SX5Analysis::SX5Analysis(const edm::ParameterSet& iConfig)
{
	tok_QIE10DigiCollection_ = consumes<HcalDataFrameContainer<QIE10DataFrame> >(edm::InputTag("hcalDigis"));
	hf_token = consumes<HFDigiCollection>(edm::InputTag("hcalDigis"));
	raw_token = consumes<FEDRawDataCollection>(edm::InputTag("source"));
	runType = iConfig.getParameter<int>("RunType");//1:pedestal 2:LED 3:histogram
	outFileName=iConfig.getUntrackedParameter<string>("OutFileName");
	histoFED = iConfig.getParameter<int>("histoFED");
	_file = new TFile(outFileName.c_str(), "recreate");
	_tree = new TTree("Events", "Events");
	_tree->Branch("ieta", ed.ieta, "ieta[144]/I");
	_tree->Branch("iphi", ed.iphi, "iphi[144]/I");
	_tree->Branch("depth", ed.depth, "depth[144]/I");
	if(runType==1 || runType==2)//pedestal or LED
	{
		_tree->Branch("pulse", ed.pulse, "pulse[144][10]/I");
	}
	else if(runType==3)//histogram
	{
		_tree->Branch("histo", ed.histo, "histo[144][60]/I");
	}
	numChannels=0;
	EID=0;
}

SX5Analysis::~SX5Analysis()
{
	_file->cd();
	_file->Write();
	_file->Close();
}

void SX5Analysis::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	using namespace edm;
	using namespace hcal;
	
	iEvent.getByToken(tok_QIE10DigiCollection_,qie10DigiCollection);
	const QIE10DigiCollection& qie10dc=*(qie10DigiCollection);
	iEvent.getByToken(raw_token,raw_collection);  
	
	edm::ESHandle<HcalElectronicsMap> item;
	edm::ESHandle<HcalDbService> pSetup;
	iSetup.get<HcalDbRecord>().get(pSetup);
	iSetup.get<HcalElectronicsMapRcd>().get(item);
	const FEDRawData& raw = raw_collection->FEDData(histoFED);
	
	edata eda;
	if(runType==1 || runType==2)//pedestal or LED
	{
		for (unsigned int j=0; j < qie10dc.size(); j++)
		{
			QIE10DataFrame qie10df = static_cast<QIE10DataFrame>(qie10dc[j]);
			DetId detid = qie10df.detid();
			HcalDetId hcaldetid = HcalDetId(detid);
			eda.ieta[j] = hcaldetid.ieta();
			eda.iphi[j] = hcaldetid.iphi();
			eda.depth[j] = hcaldetid.depth();
			int nTS = qie10df.samples();
			for(int i=0; i<nTS; ++i)
			{
				eda.pulse[j][i]=qie10df[i].adc();
			}
		}
	}
	else if(runType==3)//histogram
	{
		//the histos
		const struct eventHeader* eh =(const struct eventHeader*)(raw.data());
		const uint32_t* pData = (const uint32_t*) raw.data(); 
		uint32_t numHistos  = ((eh->h3)>>16)&0xFFFF;
		uint32_t numBins    = ((eh->h3)>>1)&0x0000FFFE; //includes overflow and header word
		uint32_t fiber   = 0;
		uint32_t channel = 0;
		pData+=8;
		int iz=0;
		int nH=-1;
		
		for (unsigned int iHist = 0; iHist<numHistos; iHist++)
		{
			if(iHist==96) iz++;
			fiber   = (*pData>>7)&0x1F;
			channel = (*pData>>2)&0x1F;
			
			bool found=false;
			for(int is1=0;is1<3;is1++)
			{
				for(int is2=0;is2<24;is2++)
				{
					for(int is3=0;is3<2;is3++)
					{
						if(MAP2Ch[is1][is2][is3][0]==((int)fiber) && MAP2Ch[is1][is2][is3][1]==((int)channel) && MAP2Ch[is1][is2][is3][2]==(iz+1)) {found=true;break;}
					}
					if(found){break;}
				}
				if(found){break;}
			}
			if(found)
			{
				nH++;
				eda.ieta[nH] = ((int)fiber);
				eda.iphi[nH] = ((int)channel);
				eda.depth[nH] = (iz+1);
			}
			for(unsigned int iBin = 0; iBin<numBins+1; iBin++)
			{
				if(iBin<60 && found)
				{
					eda.histo[nH][iBin]=pData[iBin+1];
				}
			}
			if(iHist<(numHistos-1))
			{
				pData+=(numBins+2);
			}
		}
	}
	ed=eda;
	_tree->Fill();
	EID++;
}


// ------------ method called once each job just before starting event loop  ------------
void 
SX5Analysis::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
SX5Analysis::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
SX5Analysis::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(SX5Analysis);
