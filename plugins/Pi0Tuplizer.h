

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
#include "DataFormats/EcalDetId/interface/ESDetId.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Provenance/interface/Timestamp.h"
#include "DataFormats/Provenance/interface/EventAuxiliary.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/EgammaReco/interface/PreshowerCluster.h"

#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"


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

#include "RecoCaloTools/Navigation/interface/EcalPreshowerNavigator.h"


using namespace std;
using namespace reco;
using namespace edm;

typedef std::map<DetId, EcalRecHit> RecHitsMap;
typedef std::set<DetId> HitsID;

#define NL1SEED 300
#define NPI0MAX 100

struct PosCalcParams {
     float  param_LogWeighted_;
     float  param_T0_barl_    ;
     float  param_T0_endc_    ;
     float  param_T0_endcES_  ;
     float  param_W0_         ;
     float  param_X0_         ;
};

struct cluster3x3map {
	float ehit[9];
	int iEtaiX[9];
	int iPhiiY[9];
	bool isUsedByOthers[9];

	void initAll()
	{
		for(int i=0; i<9; i++)
		{
			ehit[i] = 0;
			iEtaiX[i] = -999;
			iPhiiY[i] = -999;
			isUsedByOthers[i] = false;
		}
	}
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
   	void loadEvent_Pi0(const edm::Event& iEvent, const edm::EventSetup& iSetup); //call at the beginning of each event to get input handles from the python config
   	void loadEvent_Eta(const edm::Event& iEvent, const edm::EventSetup& iSetup); //call at the beginning of each event to get input handles from the python config
	void loadCut_Pi0(const edm::ParameterSet& iConfig);
	void loadCut_Eta(const edm::ParameterSet& iConfig);
	virtual void resetBranches(); // clear all variables
	virtual void resetPhoBranches(); // clear all variables in PhoEvents tree
   	virtual void setBranches(); // set branch of ntuple
	void GetL1SeedBit();
	void GetMCTruth();
	void MCTruthAssoc(bool isPi0, double deltaR);
	void recoPhoCluster_EB(bool isPi0_);
	void recoPhoCluster_EE(bool isPi0_);

	void recoDiPhoEvents_EB(bool isPi0_);
	void recoDiPhoEvents_EE(bool isPi0_);
	
	int diff_neta_s(int neta1, int neta2);	
	int diff_nphi_s(int nphi1, int nphi2);	
	
	float GetDeltaR(float eta1, float eta2, float phi1, float phi2);
	float DeltaPhi(float phi1, float phi2);
	void setmapforCNN(cluster3x3map g1map, cluster3x3map g2map, int index_Pi0, bool isEB);
	float getGapParameter(int ieta, int iphi);
 //output TTree and file
      	TTree *Pi0Events;
      	TTree *PhoEvents;

 //variables to be saved in the pi0 ntuple
      	uint    runNum;
      	uint    lumiNum;
      	uint    eventNum;
	uint 	eventTime;//in second, since 1970

	bool    allL1SeedFinalDecision[NL1SEED];
	vector<int> * L1SeedBitFinalDecision;

	int 	N_Pair_rec;
	int 	N_ebPair_rec;
	int 	N_eePair_rec;
	int 	N_Pho_rec;
	int 	N_ebPho_rec;
	int 	N_eePho_rec;
	int 	N_ebRecHit;
	int 	N_eeRecHit;
	int 	N_esRecHit;

	int 	N_Pi0_rec;
	int 	N_Pi0_gen;
	int 	N_Pi0_genall;
	int 	N_Pi0_match;
	int 	N_ebPi0_rec;
	int 	N_eePi0_rec;
	int 	N_Pho_rec_Pi0_;
	int 	N_ebPho_rec_Pi0_;
	int 	N_eePho_rec_Pi0_;
	int 	N_ebRecHit_Pi0_;
	int 	N_eeRecHit_Pi0_;
	int 	N_esRecHit_Pi0_;
	
	int 	N_Eta_rec;
	int 	N_Eta_gen;
	int 	N_Eta_genall;
	int 	N_Eta_match;
	int 	N_ebEta_rec;
	int 	N_eeEta_rec;
	int 	N_Pho_rec_Eta_;
	int 	N_ebPho_rec_Eta_;
	int 	N_eePho_rec_Eta_;
	int 	N_ebRecHit_Eta_;
	int 	N_eeRecHit_Eta_;
	int 	N_esRecHit_Eta_;

	int 	nIsoGamma0p3Pi0_genall[NPI0MAX];	
	int 	nIsoGamma0p2Pi0_genall[NPI0MAX];	
	int 	nIsoGamma0p1Pi0_genall[NPI0MAX];	
	float 	ptPi0_genall[NPI0MAX];	
	float 	etaPi0_genall[NPI0MAX];	
	float 	phiPi0_genall[NPI0MAX];	

	int 	nIsoGamma0p3Pi0_gen[NPI0MAX];	
	int 	nIsoGamma0p2Pi0_gen[NPI0MAX];	
	int 	nIsoGamma0p1Pi0_gen[NPI0MAX];	
	float 	ptPi0_gen[NPI0MAX];	
	float 	etaPi0_gen[NPI0MAX];	
	float 	phiPi0_gen[NPI0MAX];	
	float 	deltaRG1G2Pi0_gen[NPI0MAX];	
	float 	enG1_Pi0_gen[NPI0MAX];	
	float 	ptG1_Pi0_gen[NPI0MAX];	
	float 	etaG1_Pi0_gen[NPI0MAX];	
	float 	phiG1_Pi0_gen[NPI0MAX];	
	float 	enG2_Pi0_gen[NPI0MAX];	
	float 	ptG2_Pi0_gen[NPI0MAX];	
	float 	etaG2_Pi0_gen[NPI0MAX];	
	float 	phiG2_Pi0_gen[NPI0MAX];	

	int 	nIsoGamma0p3Eta_genall[NPI0MAX];	
	int 	nIsoGamma0p2Eta_genall[NPI0MAX];	
	int 	nIsoGamma0p1Eta_genall[NPI0MAX];	
	float 	ptEta_genall[NPI0MAX];	
	float 	etaEta_genall[NPI0MAX];	
	float 	phiEta_genall[NPI0MAX];	

	int 	nIsoGamma0p3Eta_gen[NPI0MAX];	
	int 	nIsoGamma0p2Eta_gen[NPI0MAX];	
	int 	nIsoGamma0p1Eta_gen[NPI0MAX];	
	float 	ptEta_gen[NPI0MAX];	
	float 	etaEta_gen[NPI0MAX];	
	float 	phiEta_gen[NPI0MAX];	
	float 	deltaRG1G2Eta_gen[NPI0MAX];	
	float 	enG1_Eta_gen[NPI0MAX];	
	float 	ptG1_Eta_gen[NPI0MAX];	
	float 	etaG1_Eta_gen[NPI0MAX];	
	float 	phiG1_Eta_gen[NPI0MAX];	
	float 	enG2_Eta_gen[NPI0MAX];	
	float 	ptG2_Eta_gen[NPI0MAX];	
	float 	etaG2_Eta_gen[NPI0MAX];	
	float 	phiG2_Eta_gen[NPI0MAX];	

	bool 	fromPi0[NPI0MAX];	
	float 	mPi0_rec[NPI0MAX];	
	float 	ptPi0_rec[NPI0MAX];	
	float 	isoPi0_rec[NPI0MAX];	
	float 	etaPi0_rec[NPI0MAX];	
	float 	phiPi0_rec[NPI0MAX];	
	float 	isoG1_rec[NPI0MAX];	
	float 	isoG2_rec[NPI0MAX];	
	float 	enG1_rec[NPI0MAX];	
	float 	enG1_true[NPI0MAX];	
	float 	dRG1_withtrue[NPI0MAX];	
	float 	enG2_rec[NPI0MAX];	
	float 	enG2_true[NPI0MAX];	
	float 	dRG2_withtrue[NPI0MAX];	
	float 	etaG1_rec[NPI0MAX];	
	float 	etaG2_rec[NPI0MAX];	
	float 	phiG1_rec[NPI0MAX];	
	float 	phiG2_rec[NPI0MAX];	
	float 	ptG1_rec[NPI0MAX];	
	float 	ptG2_rec[NPI0MAX];	
	int 	iEtaG1_rec[NPI0MAX];	
	int 	iXG1_rec[NPI0MAX];	
	int 	iEtaG2_rec[NPI0MAX];	
	int 	iXG2_rec[NPI0MAX];	
	int 	iPhiG1_rec[NPI0MAX];	
	int 	iYG1_rec[NPI0MAX];	
	int 	iPhiG2_rec[NPI0MAX];	
	int 	iYG2_rec[NPI0MAX];	
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
	float ehitG1_rec[NPI0MAX][9];
	float overlapG1_rec[NPI0MAX][9];
	float gapG1_rec[NPI0MAX][9];
	float ehitG2_rec[NPI0MAX][9];
	float overlapG2_rec[NPI0MAX][9];
	float gapG2_rec[NPI0MAX][9];

 //variables to be saved in the photon ntuple
	float pho_E;	
	float pho_seedE;	
	float pho_Eta;	
	int pho_iEta;	
	int pho_iX;	
	float pho_Phi;	
	int pho_iPhi;	
	int pho_iY;	
	float pho_Pt;	
	float pho_SeedTime;	
	float pho_ClusterTime;	
	float pho_S4S9;	
	float pho_S2S9;	
	float pho_S1S9;	
	int pho_Nxtal;	
	float pho_x;
	float pho_y;
	float pho_z;
		
 //reconstructed variables shared for each event
 	std::vector< CaloCluster > ebclusters_Pi0_;
 	std::vector< CaloCluster > ebclusters_Eta_;
	std::vector< CaloCluster > eeclusters_Pi0_;
	std::vector< CaloCluster > eeclusters_Eta_;
 	std::vector< cluster3x3map > ebclusters_Pi0_hitmap;
 	std::vector< cluster3x3map > ebclusters_Eta_hitmap;
 	std::vector< cluster3x3map > eeclusters_Pi0_hitmap;
 	std::vector< cluster3x3map > eeclusters_Eta_hitmap;

	vector<int> ebclusters_Pi0_MC1_index;
	vector<int> ebclusters_Eta_MC1_index;
	vector<int> eeclusters_Pi0_MC1_index;
	vector<int> eeclusters_Eta_MC1_index;
	vector<int> ebclusters_Pi0_MC2_index;
	vector<int> ebclusters_Eta_MC2_index;
	vector<int> eeclusters_Pi0_MC2_index;
	vector<int> eeclusters_Eta_MC2_index;
	
	
	CaloTopology *ebtopology_;
      	CaloTopology *eetopology_;
      	CaloSubdetectorTopology *estopology_;
	const EcalPreshowerGeometry *esGeometry_;
	const CaloGeometry* geometry;

	vector<float>  ebSeedTime_Pi0_;
	vector<float>  eeSeedTime_Pi0_;
	std::vector<int> ebNxtal_Pi0_;
	std::vector<int> eeNxtal_Pi0_;
	vector<float>  ebS4S9_Pi0_;
	vector<float>  eeS4S9_Pi0_;
	vector<float>  ebS2S9_Pi0_;
	vector<float>  eeS2S9_Pi0_;
	vector<float>  ebS1S9_Pi0_;
	vector<float>  eeS1S9_Pi0_;

	vector<float>  ebSeedTime_Eta_;
	vector<float>  eeSeedTime_Eta_;
	std::vector<int> ebNxtal_Eta_;
	std::vector<int> eeNxtal_Eta_;
	vector<float>  ebS4S9_Eta_;
	vector<float>  eeS4S9_Eta_;
	vector<float>  ebS2S9_Eta_;
	vector<float>  eeS2S9_Eta_;
	vector<float>  ebS1S9_Eta_;
	vector<float>  eeS1S9_Eta_;

	vector<int> GammaIndex_Pi0_;
	vector<int> GammaIndex_Eta_;
	vector<TLorentzVector> Gamma1MC_Pi0_;
	vector<TLorentzVector> Gamma1MC_Eta_;
	vector<TLorentzVector> Gamma2MC_Pi0_;
	vector<TLorentzVector> Gamma2MC_Eta_;
	
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

	edm::EDGetTokenT<EBRecHitCollection> EBRecHitCollectionToken_Pi0_;
      	edm::EDGetTokenT<EERecHitCollection> EERecHitCollectionToken_Pi0_;
      	edm::EDGetTokenT<ESRecHitCollection> ESRecHitCollectionToken_Pi0_;
	
	edm::EDGetTokenT<EBRecHitCollection> EBRecHitCollectionToken_Eta_;
      	edm::EDGetTokenT<EERecHitCollection> EERecHitCollectionToken_Eta_;
      	edm::EDGetTokenT<ESRecHitCollection> ESRecHitCollectionToken_Eta_;

	edm::EDGetTokenT<reco::GenParticleCollection> genParticlesToken_;

	edm::Handle<EBRecHitCollection> ebRecHit;
      	edm::Handle<EBRecHitCollection> eeRecHit;
      	edm::Handle<ESRecHitCollection> esRecHit;
	
	edm::Handle<reco::GenParticleCollection> genParticles;
	
	edm::EDGetTokenT<edm::SimTrackContainer>  g4_simTk_Token_;
      	edm::EDGetTokenT<edm::SimVertexContainer> g4_simVtx_Token_;

        edm::Handle<SimTrackContainer> simTracks_h;
	edm::Handle<SimVertexContainer> simVert_h;

	bool foundEB;
	bool foundEE;
	bool foundES;

//cuts and options read from cfg file	
 	bool isMC_;
 	bool MCAssoc_;
	bool useDRcutPair_;
	double MC_Asssoc_DeltaR;
	//bool isPi0_;
	bool FillL1SeedFinalDecision_;
	bool FillDiPhotonNtuple_;
	bool FillPhotonNtuple_;
	std::string PhotonOrderOption_;
	double EB_Seed_E_Pi0_; //seed energy threshold setting for EB
	double EE_Seed_E_Pi0_; //seed energy threshold setting for EE
	double pairPtCut_barrel1_Pi0_;//0.0 < abs(eta) < 1.0
	double pairPtCut_barrel2_Pi0_;//1.0 < abs(eta) < 1.5
	double pairPtCut_endcap1_Pi0_;//1.5 < abs(eta) < 1.8
	double pairPtCut_endcap2_Pi0_;//1.8 < abs(eta)
	double gPtCut_barrel1_Pi0_;
	double gPtCut_barrel2_Pi0_;
	double gPtCut_endcap1_Pi0_;
	double gPtCut_endcap2_Pi0_;
	double s4s9Cut_barrel1_Pi0_;
	double s4s9Cut_barrel2_Pi0_;
	double s4s9Cut_endcap1_Pi0_;
	double s4s9Cut_endcap2_Pi0_;
	double nxtal1Cut_barrel1_Pi0_;
	double nxtal1Cut_barrel2_Pi0_;
	double nxtal1Cut_endcap1_Pi0_;
	double nxtal1Cut_endcap2_Pi0_;
	double nxtal2Cut_barrel1_Pi0_;
	double nxtal2Cut_barrel2_Pi0_;
	double nxtal2Cut_endcap1_Pi0_;
	double nxtal2Cut_endcap2_Pi0_;
	double DRcutPair_Pi0_;

	double EB_Seed_E_Eta_; //seed energy threshold setting for EB
	double EE_Seed_E_Eta_; //seed energy threshold setting for EE
	double pairPtCut_barrel1_Eta_;//0.0 < abs(eta) < 1.0
	double pairPtCut_barrel2_Eta_;//1.0 < abs(eta) < 1.5
	double pairPtCut_endcap1_Eta_;//1.5 < abs(eta) < 1.8
	double pairPtCut_endcap2_Eta_;//1.8 < abs(eta)
	double gPtCut_barrel1_Eta_;
	double gPtCut_barrel2_Eta_;
	double gPtCut_endcap1_Eta_;
	double gPtCut_endcap2_Eta_;
	double s4s9Cut_barrel1_Eta_;
	double s4s9Cut_barrel2_Eta_;
	double s4s9Cut_endcap1_Eta_;
	double s4s9Cut_endcap2_Eta_;
	double nxtal1Cut_barrel1_Eta_;
	double nxtal1Cut_barrel2_Eta_;
	double nxtal1Cut_endcap1_Eta_;
	double nxtal1Cut_endcap2_Eta_;
	double nxtal2Cut_barrel1_Eta_;
	double nxtal2Cut_barrel2_Eta_;
	double nxtal2Cut_endcap1_Eta_;
	double nxtal2Cut_endcap2_Eta_;
	double DRcutPair_Eta_;

	double isoGammaBeltdR_Zone_Pi0_;
	double isoGammaBeltdEta_Zone_Pi0_;
	double isoPairBeltdR_Zone_Pi0_;
	double isoPairBeltdEta_Zone_Pi0_;
	double isoGammaBeltdR_Zone_Eta_;
	double isoGammaBeltdEta_Zone_Eta_;
	double isoPairBeltdR_Zone_Eta_;
	double isoPairBeltdEta_Zone_Eta_;
	double isoPairCut_;
	double isoGammaCut_;
};

class PreshowerTools{
    public:
	PreshowerTools(const CaloGeometry* extGeom, CaloSubdetectorTopology* topology_p,  edm::Handle< ESRecHitCollection > & esHandle);
      	PreshowerCluster makeOnePreshowerCluster(int stripwindow, ESDetId *strip);
	static const double mip_;
      	static const double gamma_;
      	static const double calib_planeX_;
      	static const double calib_planeY_;
      	static const int clusterwindowsize_;


    private:
      	const CaloGeometry* geom_;
      	CaloSubdetectorTopology *estopology_;
      	std::vector<ESDetId> esroad_2d;
	HitsID  used_strips;
	RecHitsMap  rechits_map;
	void findESRoad(int stripwindow, ESDetId strip, EcalPreshowerNavigator theESNav, int plane);
      	bool goodStrip(RecHitsMap::iterator candidate_it);
};

#endif
