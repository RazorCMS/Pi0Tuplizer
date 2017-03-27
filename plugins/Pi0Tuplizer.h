

#ifndef PI0TUPLIZER_H
#define PI0TUPLIZER_H


/***********************c++*********************/
#include <memory>
#include <string>
#include <vector>
#include <tuple>

/***********************root********************/
#include "TFile.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TLorentzVector.h"
#include "TTree.h"

/***********************cmssw*******************/
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/ESProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "CondFormats/EcalObjects/interface/EcalChannelStatus.h"
#include "CondFormats/DataRecord/interface/EcalChannelStatusRcd.h"
#include "DataFormats/L1Trigger/interface/EGamma.h"
#include "DataFormats/L1Trigger/interface/Muon.h"
#include "DataFormats/L1Trigger/interface/Tau.h"
#include "DataFormats/L1Trigger/interface/Jet.h"
#include "DataFormats/L1Trigger/interface/EtSum.h"
#include "DataFormats/L1TGlobal/interface/GlobalAlgBlk.h"
#include "DataFormats/L1TGlobal/interface/GlobalExtBlk.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "DataFormats/CaloRecHit/interface/CaloID.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
#include "DataFormats/Math/interface/LorentzVector.h"

#include "RecoEcal/EgammaCoreTools/interface/EcalRecHitLess.h"
#include "RecoEcal/EgammaCoreTools/interface/PositionCalc.h"

#include "L1Trigger/L1TGlobal/interface/L1TGlobalUtil.h"

#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/CaloTopology/interface/EcalBarrelHardcodedTopology.h"
#include "Geometry/CaloTopology/interface/EcalEndcapHardcodedTopology.h"
#include "Geometry/CaloTopology/interface/EcalPreshowerTopology.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloGeometry/interface/TruncatedPyramid.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/EcalAlgo/interface/EcalPreshowerGeometry.h"


using namespace std;
using namespace reco;
using namespace edm;

#define NL1SEED 300
#define NPI0MAX 300

struct PosCalcParams {
     float  param_LogWeighted_;
     float  param_T0_barl_    ;
     float  param_T0_endc_    ;
     float  param_T0_endcES_  ;
     float  param_W0_         ;
     float  param_X0_         ;
};

class Pi0Tuplizer : public edm::EDAnalyzer {
public:
 //analyzer constructor and destructor
  	explicit Pi0Tuplizer(const edm::ParameterSet&);
  	 ~Pi0Tuplizer();

	static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
   
private:
      	virtual void beginJob() ;
      	virtual void analyze(const edm::Event&, const edm::EventSetup&);
      	virtual void endJob() ;

      	virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      	virtual void endRun(edm::Run const&, edm::EventSetup const&);
   
 //specific functions 
   	void loadEvent(const edm::Event& iEvent, const edm::EventSetup& iSetup); //call at the beginning of each event to get input handles from the python config
	virtual void resetBranches(); // clear all variables
   	virtual void setBranches(); // set branch of ntuple
	void GetL1SeedBit();
	void recoPhoCluster_EB();
	void recoPhoCluster_EE();

	void recoDiPhoEvents_EB();
	void recoDiPhoEvents_EE();
	
	int diff_neta_s(int neta1, int neta2);	
	int diff_nphi_s(int nphi1, int nphi2);	
	
	float GetDeltaR(float eta1, float eta2, float phi1, float phi2);
	float DeltaPhi(float phi1, float phi2);
 //output TTree and file
      	TTree *Pi0Events;
      	TTree *PhoEvents;

 //variables to be saved in the pi0 ntuple
      	uint    runNum;
      	uint    lumiNum;
      	uint    eventNum;
	bool    allL1SeedFinalDecision[NL1SEED];
	vector<int> * L1SeedBitFinalDecision;

	int 	N_Pi0_rec;
	float 	mPi0_rec[NPI0MAX];	
	float 	ptPi0_rec[NPI0MAX];	
	float 	etaPi0_rec[NPI0MAX];	
	float 	phiPi0_rec[NPI0MAX];	
	float 	enG1_rec[NPI0MAX];	
	float 	enG2_rec[NPI0MAX];	
	float 	etaG1_rec[NPI0MAX];	
	float 	etaG2_rec[NPI0MAX];	
	float 	phiG1_rec[NPI0MAX];	
	float 	phiG2_rec[NPI0MAX];	
	float 	ptG1_rec[NPI0MAX];	
	float 	ptG2_rec[NPI0MAX];	
	int 	iEtaG1_rec[NPI0MAX];	
	int 	iEtaG2_rec[NPI0MAX];	
	int 	iPhiG1_rec[NPI0MAX];	
	int 	iPhiG2_rec[NPI0MAX];	
	float 	deltaRG1G2_rec[NPI0MAX];	
	int 	nxtalG1_rec[NPI0MAX];	
	int 	nxtalG2_rec[NPI0MAX];	
	float 	seedTimeG1_rec[NPI0MAX];	
	float 	seedTimeG2_rec[NPI0MAX];	
	float 	s4s9G1_rec[NPI0MAX];	
	float 	s4s9G2_rec[NPI0MAX];	
	float 	s2s9G1_rec[NPI0MAX];	
	float 	s2s9G2_rec[NPI0MAX];	
	float 	s1s9G1_rec[NPI0MAX];	
	float 	s1s9G2_rec[NPI0MAX];	

 //variables to be saved in the photon ntuple
	float pho_E;	
	float pho_Eta;	
	int pho_iEta;	
	float pho_Phi;	
	int pho_iPhi;	
	float pho_Pt;	
	float pho_SeedTime;	
	float pho_ClusterTime;	
	float pho_S4S9;	
	float pho_S2S9;	
	float pho_S1S9;	
	int pho_Nxtal;	
		
 //reconstructed variables shared for each event
 	std::vector< CaloCluster > ebclusters;
	std::vector< CaloCluster > eeclusters;
	
	CaloTopology *ebtopology_;
      	CaloTopology *eetopology_;
      	CaloSubdetectorTopology *estopology_;
	const CaloGeometry* geometry;

	vector<float>  ebSeedTime;
	vector<float>  eeSeedTime;
	std::vector<int> ebNxtal;
	std::vector<int> eeNxtal;
	vector<float>  ebS4S9;
	vector<float>  eeS4S9;
	vector<float>  ebS2S9;
	vector<float>  eeS2S9;
	vector<float>  ebS1S9;
	vector<float>  eeS1S9;

	PosCalcParams PCparams_ = {
				true,//param_LogWeighted_
				7.4,//param_T0_barl_
				3.1,//param_T0_endc_
				1.2,//param_T0_endcES_
				4.2,//param_W0_ 
				0.89//param_X0_
				}; //shower shape paramters
	 
 //Collections and related variables read from input file
 	edm::EDGetTokenT<BXVector<GlobalAlgBlk>> uGtAlgToken_;
	edm::Handle<BXVector<GlobalAlgBlk>> uGtAlg;

	edm::EDGetTokenT<EBRecHitCollection> EBRecHitCollectionToken_;
      	edm::EDGetTokenT<EERecHitCollection> EERecHitCollectionToken_;
      	edm::EDGetTokenT<ESRecHitCollection> ESRecHitCollectionToken_;
	edm::Handle<EBRecHitCollection> ebRecHit;
      	edm::Handle<EBRecHitCollection> eeRecHit;
      	edm::Handle<ESRecHitCollection> esRecHit;

//cuts and options read from cfg file	
 	bool isMC_;
	bool isPi0_;
	bool FillL1SeedFinalDecision_;
	bool FillDiPhotonNtuple_;
	bool FillPhotonNtuple_;
	std::string PhotonOrderOption_;
	double EB_Seed_E_; //seed energy threshold setting for EB
	double EE_Seed_E_; //seed energy threshold setting for EE
	double pi0PtCut_barrel1;//0.0 < abs(eta) < 1.0
	double pi0PtCut_barrel2;//1.0 < abs(eta) < 1.5
	double pi0PtCut_endcap1;//1.5 < abs(eta) < 1.8
	double pi0PtCut_endcap2;//1.8 < abs(eta)
	double gPtCut_barrel1;
	double gPtCut_barrel2;
	double gPtCut_endcap1;
	double gPtCut_endcap2;
	double s4s9Cut_barrel1;
	double s4s9Cut_barrel2;
	double s4s9Cut_endcap1;
	double s4s9Cut_endcap2;
	double nxtal1Cut_barrel1;
	double nxtal1Cut_barrel2;
	double nxtal1Cut_endcap1;
	double nxtal1Cut_endcap2;
	double nxtal2Cut_barrel1;
	double nxtal2Cut_barrel2;
	double nxtal2Cut_endcap1;
	double nxtal2Cut_endcap2;
		
};

#endif
