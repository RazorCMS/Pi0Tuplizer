// -*- C++ -*-
/*
 *   Description: reconstruction of pi0/eta and save the ntuple
*/
// Author: Zhicai Zhang
// Created: Fri Mar 24 11:01:09 CET 2017


#include "Pi0Tuplizer.h"

using namespace std;

#define DEBUG


const double PreshowerTools::mip_ = 8.108e-05;
const double PreshowerTools::gamma_ = 0.024;
const double PreshowerTools::calib_planeX_ = 1.0;
const double PreshowerTools::calib_planeY_ = 0.7;
const int    PreshowerTools::clusterwindowsize_ = 15;


class ecalRecHitLess : public std::binary_function<EcalRecHit, EcalRecHit, bool>
{
public:
  bool operator()(EcalRecHit x, EcalRecHit y)
  {
    return (x.energy() > y.energy());
  }
};



Pi0Tuplizer::Pi0Tuplizer(const edm::ParameterSet& iConfig)
{
	//get parameters from iConfig
	isMC_	= iConfig.getUntrackedParameter<bool>("isMC",false);
	MCAssoc_	= iConfig.getUntrackedParameter<bool>("MCAssoc",false);
	MC_Asssoc_DeltaR                   = iConfig.getUntrackedParameter<double>("MC_Asssoc_DeltaR",0.1);
	//isPi0_	= iConfig.getUntrackedParameter<bool>("isPi0",false);
	FillL1SeedFinalDecision_	= iConfig.getUntrackedParameter<bool>("FillL1SeedFinalDecision",false);
	FillDiPhotonNtuple_	= iConfig.getUntrackedParameter<bool>("FillDiPhotonNtuple",false);
	FillPhotonNtuple_	= iConfig.getUntrackedParameter<bool>("FillPhotonNtuple",false);

        if(FillL1SeedFinalDecision_)   uGtAlgToken_ = consumes<BXVector<GlobalAlgBlk>>(iConfig.getUntrackedParameter<edm::InputTag>("uGtAlgInputTag"));	
	if(MCAssoc_) 	
	{
		g4_simTk_Token_  = consumes<edm::SimTrackContainer>(edm::InputTag("g4SimHits"));
	    	g4_simVtx_Token_ = consumes<edm::SimVertexContainer>(edm::InputTag("g4SimHits"));
		genParticlesToken_      = consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genParticles"));
	}
	
	edm::InputTag tag_EBRecHit_Pi0_ = iConfig.getUntrackedParameter<edm::InputTag>("EBRecHitCollectionTag_Pi0");	
	edm::InputTag tag_EERecHit_Pi0_ = iConfig.getUntrackedParameter<edm::InputTag>("EERecHitCollectionTag_Pi0");	
	edm::InputTag tag_ESRecHit_Pi0_ = iConfig.getUntrackedParameter<edm::InputTag>("ESRecHitCollectionTag_Pi0");	

	EBRecHitCollectionToken_Pi0_    = consumes<EBRecHitCollection>(tag_EBRecHit_Pi0_);
	EERecHitCollectionToken_Pi0_    = consumes<EERecHitCollection>(tag_EERecHit_Pi0_);
	ESRecHitCollectionToken_Pi0_    = consumes<ESRecHitCollection>(tag_ESRecHit_Pi0_);

	edm::InputTag tag_EBRecHit_Eta_ = iConfig.getUntrackedParameter<edm::InputTag>("EBRecHitCollectionTag_Eta");	
	edm::InputTag tag_EERecHit_Eta_ = iConfig.getUntrackedParameter<edm::InputTag>("EERecHitCollectionTag_Eta");	
	edm::InputTag tag_ESRecHit_Eta_ = iConfig.getUntrackedParameter<edm::InputTag>("ESRecHitCollectionTag_Eta");	

	EBRecHitCollectionToken_Eta_    = consumes<EBRecHitCollection>(tag_EBRecHit_Eta_);
	EERecHitCollectionToken_Eta_    = consumes<EERecHitCollection>(tag_EERecHit_Eta_);
	ESRecHitCollectionToken_Eta_    = consumes<ESRecHitCollection>(tag_ESRecHit_Eta_);

	PhotonOrderOption_		= iConfig.getUntrackedParameter<std::string>("PhotonOrderOption");

	loadCut_Pi0(iConfig);
	loadCut_Eta(iConfig);
		
	ebtopology_ = new CaloTopology();
    	EcalBarrelHardcodedTopology* ebHTopology = new EcalBarrelHardcodedTopology();
    	ebtopology_->setSubdetTopology(DetId::Ecal,EcalBarrel,ebHTopology);

	eetopology_ = new CaloTopology();
    	EcalEndcapHardcodedTopology* eeHTopology=new EcalEndcapHardcodedTopology();
    	eetopology_->setSubdetTopology(DetId::Ecal,EcalEndcap,eeHTopology);

#ifdef DEBUG
	cout<<"InputTag -- "<<endl;
	cout<<"tag_EBRecHit_Pi0_ =  "<<tag_EBRecHit_Pi0_.label()<<", "<<tag_EBRecHit_Pi0_.instance()<<", "<<tag_EBRecHit_Pi0_.process()<<endl;
	cout<<"tag_EERecHit_Pi0_ =  "<<tag_EERecHit_Pi0_.label()<<", "<<tag_EERecHit_Pi0_.instance()<<", "<<tag_EERecHit_Pi0_.process()<<endl;
	cout<<"tag_ESRecHit_Pi0_ =  "<<tag_ESRecHit_Pi0_.label()<<", "<<tag_ESRecHit_Pi0_.instance()<<", "<<tag_ESRecHit_Pi0_.process()<<endl;
	cout<<"tag_EBRecHit_Eta_ =  "<<tag_EBRecHit_Eta_.label()<<", "<<tag_EBRecHit_Eta_.instance()<<", "<<tag_EBRecHit_Eta_.process()<<endl;
	cout<<"tag_EERecHit_Eta_ =  "<<tag_EERecHit_Eta_.label()<<", "<<tag_EERecHit_Eta_.instance()<<", "<<tag_EERecHit_Eta_.process()<<endl;
	cout<<"tag_ESRecHit_Eta_ =  "<<tag_ESRecHit_Eta_.label()<<", "<<tag_ESRecHit_Eta_.instance()<<", "<<tag_ESRecHit_Eta_.process()<<endl;

#endif
	
	//create output file and tree	
	edm::Service<TFileService> fs;
	if(FillDiPhotonNtuple_) 
	{
		Pi0Events = fs->make<TTree>("Pi0Events", "reconstructed pi0/eta ntuple");	
	}
	if(FillPhotonNtuple_) PhoEvents = fs->make<TTree>("PhoEvents", "reconstructed photon ntuple");	

}


Pi0Tuplizer::~Pi0Tuplizer()
{

}


//------ Method called for each event ------//
void Pi0Tuplizer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

	resetBranches();
	loadEvent(iEvent, iSetup);

//fill event info
  	runNum = iEvent.id().run();
  	lumiNum = iEvent.luminosityBlock();
  	eventNum = iEvent.id().event();
	eventTime = iEvent.eventAuxiliary().time().unixTime();

//call specific functions
	if(FillL1SeedFinalDecision_) GetL1SeedBit();

	if(MCAssoc_) GetMCTruth();

	loadEvent_Pi0(iEvent, iSetup); 
	resetPhoBranches();
	if(foundEB) recoPhoCluster_EB(true);//reconstruct photon clusters in EB
	resetPhoBranches();
	if(foundEE) recoPhoCluster_EE(true);//reconstruct photon clusters in EE
	

	if(MCAssoc_) MCTruthAssoc(true, MC_Asssoc_DeltaR);
	

	if(FillDiPhotonNtuple_) recoDiPhoEvents_EB(true);//reconstruct pi0/eta events from the photon clusters in EB
	N_ebPair_rec += N_Pair_rec;
	if(FillDiPhotonNtuple_) recoDiPhoEvents_EE(true);//reconstruct pi0/eta events from the photon clusters in EE
	N_eePair_rec += N_Pair_rec - N_ebPair_rec;
	
	N_Pi0_rec = N_Pair_rec;
	N_ebPi0_rec = N_ebPair_rec;
	N_eePi0_rec = N_eePair_rec;
		

	loadEvent_Eta(iEvent, iSetup); 
	resetPhoBranches();
	if(foundEB) recoPhoCluster_EB(false);//reconstruct photon clusters in EB
	resetPhoBranches();
	if(foundEE) recoPhoCluster_EE(false);//reconstruct photon clusters in EE

 	if(MCAssoc_) MCTruthAssoc(false, MC_Asssoc_DeltaR);
	

	if(FillDiPhotonNtuple_) recoDiPhoEvents_EB(false);//reconstruct pi0/eta events from the photon clusters in EB
	N_ebPair_rec = N_Pair_rec - N_eePair_rec;
	if(FillDiPhotonNtuple_) recoDiPhoEvents_EE(false);//reconstruct pi0/eta events from the photon clusters in EE
	N_eePair_rec = N_Pair_rec - N_ebPair_rec;
	
	N_Eta_rec = N_Pair_rec - N_Pi0_rec;
	N_ebEta_rec = N_ebPair_rec - N_ebPi0_rec;
	N_eeEta_rec = N_eePair_rec - N_eePi0_rec;
//fill ntuple
#ifdef DEBUG
	cout<<"N_Pho: "<<N_Pho_rec<<"  N_ebRecHit: "<<N_ebRecHit<<"   N_eeRecHit:  "<<N_eeRecHit<<endl;
#endif
	if(FillDiPhotonNtuple_ && (N_Pair_rec > 0 || (MCAssoc_ && ((N_Pi0_gen+N_Eta_gen) > 0) )) ) Pi0Events->Fill();
	//Pi0Events->Fill();

}


//
void Pi0Tuplizer::recoPhoCluster_EB(bool isPi0_)
{
	std::vector<EcalRecHit> ebseeds;

	for(EBRecHitCollection::const_iterator itb= ebRecHit->begin(); itb != ebRecHit->end(); ++itb)
	{
		if(itb->energy() > (isPi0_?EB_Seed_E_Pi0_:EB_Seed_E_Eta_))  ebseeds.push_back( *itb );	
	}
	
	//if(ebseeds.size() < 2) return;
	
	
	sort(ebseeds.begin(), ebseeds.end(), ecalRecHitLess());	
	//make photon clusters
	std::set<EBDetId> isUsed;
	for (std::vector<EcalRecHit>::iterator itseed=ebseeds.begin(); itseed!=ebseeds.end(); itseed++)
	{
		EBDetId seed_id( itseed->id() );
		if(isUsed.count(seed_id)!=0) continue;

		int seed_ieta = seed_id.ieta();
    		int seed_iphi = seed_id.iphi();
		
		std::vector<DetId> clus_v = ebtopology_->getWindow(seed_id,3,3);
		vector<const EcalRecHit*> RecHitsInWindow;
			
		float simple_energy = 0.0;
		float posTotalEnergy = 0.0;
		for (std::vector<DetId>::const_iterator det=clus_v.begin(); det!=clus_v.end(); det++)
		{
			EBDetId thisId( *det );
			if(isUsed.count(thisId)!=0) continue; 			
			EBRecHitCollection::const_iterator ixtal = ebRecHit->find( thisId );
			if( ixtal == ebRecHit->end() ) continue; 
			RecHitsInWindow.push_back( &(*ixtal) );
			simple_energy +=  ixtal->energy();
			if(ixtal->energy()>0.) posTotalEnergy += ixtal->energy();	
		}
	
		if(simple_energy <= 0) continue;

		//variables to be saved in the 3x3 cluster
		float e3x3 = 0.0;
		float e2x2 = 0.0;
		float s4s9_tmp[4]={0.,0.,0.,0.};

		float maxEne = itseed->energy();
		float maxEne2 = -999.;

		std::vector<std::pair<DetId,float> > enFracs;

		float total_weight = 0.0;
		float xclu(0.), yclu(0.), zclu(0.);
		
		for(unsigned int j=0; j<RecHitsInWindow.size();j++)
    		{
			EBDetId det(RecHitsInWindow[j]->id());	
			int ieta = det.ieta();
	        	int iphi = det.iphi();
			int deta = diff_neta_s(seed_ieta,ieta);
        		int dphi = diff_nphi_s(seed_iphi,iphi);
			float en = RecHitsInWindow[j]->energy();
	
			if(abs(deta)<=1 && abs(dphi)<=1)
        		{
				e3x3 +=  en;
				if(en<maxEne && en>maxEne2) maxEne2 = en;

				if(deta <= 0 && dphi <=0){ s4s9_tmp[0] += en; }
          			if(deta >= 0 && dphi <=0){ s4s9_tmp[1] += en; }
          			if(deta <= 0 && dphi >=0){ s4s9_tmp[2] += en; }
          			if(deta >= 0 && dphi >=0){ s4s9_tmp[3] += en; }

				isUsed.insert(RecHitsInWindow[j]->id());
		
				enFracs.push_back( std::make_pair( RecHitsInWindow[j]->id(), en ) );
	
				float weight = std::max( float(0.), PCparams_.param_W0_ + log(en/posTotalEnergy) );

				float maxDepth = PCparams_.param_X0_ * ( PCparams_.param_T0_barl_ + log( posTotalEnergy ) );
				float maxToFront;
				float pos_geo;
				
				const CaloCellGeometry* cell=geometry->getGeometry(det);
            			GlobalPoint posit = ( dynamic_cast<const TruncatedPyramid*>(cell) )->getPosition( 0. );
            			pos_geo = posit.mag();
			
				const CaloCellGeometry* cell_seed=geometry->getGeometry( seed_id );
        			GlobalPoint posit_seed = ( dynamic_cast<const TruncatedPyramid*>(cell_seed) )->getPosition( 0. );
        			maxToFront = posit_seed.mag();
		
				float depth = maxDepth + maxToFront - pos_geo;
				GlobalPoint posThis = ( dynamic_cast<const TruncatedPyramid*>(cell) )->getPosition( depth );
				xclu += weight*posThis.x();
          			yclu += weight*posThis.y();
          			zclu += weight*posThis.z();
          			total_weight += weight;
			}	
		}//end loop of xtals in 3x3 clusters

		e2x2 = *max_element( s4s9_tmp,s4s9_tmp+4);
	
		math::XYZPoint clusPos( xclu/total_weight, yclu/total_weight, zclu/total_weight );	
		
		//apply pt and s4s9 cut on photon cluster
		float s4s9 = e2x2/e3x3;
		float s2s9 =maxEne/e3x3;
		if(maxEne2 > -900.0 ) s2s9 = (maxEne+maxEne2)/e3x3;

		float ptClus = e3x3/cosh(clusPos.eta());
		if(fabs( clusPos.eta() ) < 1.0 ) 
		{
			if(ptClus<(isPi0_?gPtCut_barrel1_Pi0_:gPtCut_barrel1_Eta_)) continue;
			if(s4s9<(isPi0_?s4s9Cut_barrel1_Pi0_:s4s9Cut_barrel1_Eta_)) continue;
		}
		else if (fabs( clusPos.eta() ) < 1.5 ) 
		{
			if(ptClus<(isPi0_?gPtCut_barrel2_Pi0_:gPtCut_barrel2_Eta_)) continue;
			if(s4s9<(isPi0_?s4s9Cut_barrel2_Pi0_:s4s9Cut_barrel2_Eta_)) continue;
		}
		else 
		{
			continue;
		}
	
		if(isPi0_)
		{	
		ebclusters_Pi0_.push_back( CaloCluster( e3x3, clusPos, CaloID(CaloID::DET_ECAL_BARREL), enFracs, CaloCluster::undefined, seed_id ) );
		ebSeedTime_Pi0_.push_back( itseed->time() );	
		ebNxtal_Pi0_.push_back(RecHitsInWindow.size());
		ebS4S9_Pi0_.push_back(s4s9);		
		ebS2S9_Pi0_.push_back(s2s9);		
		ebS1S9_Pi0_.push_back((maxEne)/e3x3);		
		}
	
		else
		{
		ebclusters_Eta_.push_back( CaloCluster( e3x3, clusPos, CaloID(CaloID::DET_ECAL_BARREL), enFracs, CaloCluster::undefined, seed_id ) );
		ebSeedTime_Eta_.push_back( itseed->time() );	
		ebNxtal_Eta_.push_back(RecHitsInWindow.size());
		ebS4S9_Eta_.push_back(s4s9);		
		ebS2S9_Eta_.push_back(s2s9);		
		ebS1S9_Eta_.push_back((maxEne)/e3x3);		
		}
		//fill photon cluster
		pho_E = e3x3;
		pho_seedE = maxEne;
		pho_Eta = clusPos.eta();
		pho_iEta = seed_id.ieta();
		pho_Phi = clusPos.phi();
		pho_iPhi = seed_id.iphi();
		pho_Pt = ptClus;
		pho_SeedTime = itseed->time();
		pho_ClusterTime = -999.9;//to be constructed
		pho_S4S9 = s4s9;
		pho_S2S9 = s2s9;
		pho_S1S9 = (maxEne)/e3x3; 
		pho_Nxtal = RecHitsInWindow.size();
		pho_x = clusPos.x();
		pho_y = clusPos.y();
		pho_z = clusPos.z();
		if(FillPhotonNtuple_ ) PhoEvents->Fill();
		N_Pho_rec ++;
		N_ebPho_rec ++;
		if(isPi0_)
		{
		N_Pho_rec_Pi0_ ++;
		N_ebPho_rec_Pi0_ ++;
		}
		else
		{
		N_Pho_rec_Eta_ ++;
		N_ebPho_rec_Eta_ ++;
		}
		
	}//end loop of all seed xtals	
	
}

void Pi0Tuplizer::recoPhoCluster_EE(bool isPi0_)
{

	PreshowerTools esClusteringAlgo(geometry, estopology_, esRecHit);
	

	std::vector<EcalRecHit> eeseends;

	for(EERecHitCollection::const_iterator itb= eeRecHit->begin(); itb != eeRecHit->end(); ++itb)
	{
		if(itb->energy() > (isPi0_?EE_Seed_E_Pi0_:EE_Seed_E_Eta_))  eeseends.push_back( *itb );	
	}
	
	//if(eeseends.size() < 2) return;
	
	
	sort(eeseends.begin(), eeseends.end(), ecalRecHitLess());	
	//make photon clusters
	std::set<EEDetId> isUsed;
	for (std::vector<EcalRecHit>::iterator itseed=eeseends.begin(); itseed!=eeseends.end(); itseed++)
	{
		EEDetId seed_id( itseed->id() );
		if(isUsed.count(seed_id)!=0) continue;

		int seed_ix = seed_id.ix();
    		int seed_iy = seed_id.iy();
		
		std::vector<DetId> clus_v = eetopology_->getWindow(seed_id,3,3);
		vector<const EcalRecHit*> RecHitsInWindow;
			
		float simple_energy = 0.0;
		float posTotalEnergy = 0.0;
		for (std::vector<DetId>::const_iterator det=clus_v.begin(); det!=clus_v.end(); det++)
		{
			EEDetId thisId( *det );
			if(isUsed.count(thisId)!=0) continue; 			
			EERecHitCollection::const_iterator ixtal = eeRecHit->find( thisId );
			if( ixtal == eeRecHit->end() ) continue; 
			RecHitsInWindow.push_back( &(*ixtal) );
			simple_energy +=  ixtal->energy();
			if(ixtal->energy()>0.) posTotalEnergy += ixtal->energy();	
		}
	
		if(simple_energy <= 0) continue;

		//variables to be saved in the 3x3 cluster
		float e3x3 = 0.0;
		float e2x2 = 0.0;
		float s4s9_tmp[4]={0.,0.,0.,0.};

		float maxEne = itseed->energy();
		float maxEne2 = -999.;

		std::vector<std::pair<DetId,float> > enFracs;

		float total_weight = 0.0;
		float xclu(0.), yclu(0.), zclu(0.);
		
		for(unsigned int j=0; j<RecHitsInWindow.size();j++)
    		{
			EEDetId det(RecHitsInWindow[j]->id());	
			int ix = det.ix();
	        	int iy = det.iy();
			int dix = seed_ix-ix;
        		int diy = seed_iy-iy;
			float en = RecHitsInWindow[j]->energy();
	
			if(abs(dix)<=1 && abs(diy)<=1)
        		{
				e3x3 +=  en;
				if(en<maxEne && en>maxEne2) maxEne2 = en;

				if(dix <= 0 && diy <=0){ s4s9_tmp[0] += en; }
          			if(dix >= 0 && diy <=0){ s4s9_tmp[1] += en; }
          			if(dix <= 0 && diy >=0){ s4s9_tmp[2] += en; }
          			if(dix >= 0 && diy >=0){ s4s9_tmp[3] += en; }

				isUsed.insert(RecHitsInWindow[j]->id());
		
				enFracs.push_back( std::make_pair( RecHitsInWindow[j]->id(), en ) );
	
				float weight = std::max( float(0.), PCparams_.param_W0_ + log(en/posTotalEnergy) );

				float maxDepth = PCparams_.param_X0_ * ( PCparams_.param_T0_barl_ + log( posTotalEnergy ) );
				float maxToFront;
				float pos_geo;
				
				const CaloCellGeometry* cell=geometry->getGeometry(det);
            			GlobalPoint posit = ( dynamic_cast<const TruncatedPyramid*>(cell) )->getPosition( 0. );
            			pos_geo = posit.mag();
			
				const CaloCellGeometry* cell_seed=geometry->getGeometry( seed_id );
        			GlobalPoint posit_seed = ( dynamic_cast<const TruncatedPyramid*>(cell_seed) )->getPosition( 0. );
        			maxToFront = posit_seed.mag();
		
				float depth = maxDepth + maxToFront - pos_geo;
				GlobalPoint posThis = ( dynamic_cast<const TruncatedPyramid*>(cell) )->getPosition( depth );
				xclu += weight*posThis.x();
          			yclu += weight*posThis.y();
          			zclu += weight*posThis.z();
          			total_weight += weight;
			}	
		}//end loop of xtals in 3x3 clusters

		e2x2 = *max_element( s4s9_tmp,s4s9_tmp+4);
	
		math::XYZPoint clusPos( xclu/total_weight, yclu/total_weight, zclu/total_weight );	
		
		//apply pt and s4s9 cut on photon cluster
		float s4s9 = e2x2/e3x3;
		float s2s9 =maxEne/e3x3;
		if(maxEne2 > -900.0 ) s2s9 = (maxEne+maxEne2)/e3x3;

		float ptClus = e3x3/cosh(clusPos.eta());
		if(fabs( clusPos.eta() ) < 1.5 ) continue;
		else if(fabs( clusPos.eta() ) < 1.8 ) 
		{
			if(ptClus<(isPi0_?gPtCut_endcap1_Pi0_:gPtCut_endcap1_Eta_)) continue;
			if(s4s9<(isPi0_?s4s9Cut_endcap1_Pi0_:s4s9Cut_endcap1_Eta_)) continue;
		}
		else if (fabs( clusPos.eta() ) < 3.0 ) 
		{
			if(ptClus<(isPi0_?gPtCut_endcap2_Pi0_:gPtCut_endcap2_Eta_)) continue;
			if(s4s9<(isPi0_?s4s9Cut_endcap2_Pi0_:s4s9Cut_endcap2_Eta_)) continue;
		}
		else 
		{
			continue;
		}

// add preshower energy
		double deltaE = 0.0;
		if(fabs(clusPos.eta())>1.7  &&  fabs(clusPos.eta())<2.55)
		{
			double X = clusPos.x();
			double Y = clusPos.y();
			double Z = clusPos.z();
			const GlobalPoint point(X,Y,Z);

			DetId tmp1 = esGeometry_->getClosestCellInPlane(point,1);
        		DetId tmp2 = esGeometry_->getClosestCellInPlane(point,2);

			if ((tmp1.rawId()!=0) && (tmp2.rawId()!=0))
        		{
				ESDetId tmp1_conversion (tmp1);
          			ESDetId tmp2_conversion (tmp2);
				float es_clusterwindowsize = 5;
          			PreshowerCluster preshowerclusterp1 = esClusteringAlgo.makeOnePreshowerCluster( es_clusterwindowsize, &tmp1_conversion);
          			PreshowerCluster preshowerclusterp2 = esClusteringAlgo.makeOnePreshowerCluster( es_clusterwindowsize, &tmp2_conversion);
				double e1 = preshowerclusterp1.energy();
          			double e2 = preshowerclusterp2.energy();
//		cout<<"DEBUG deltaE  e1 = "<<e1<<"   e2 = "<<e2<<endl;
          			// GeV to #MIPs
          			e1 = e1 / PreshowerTools::mip_;
          			e2 = e2 / PreshowerTools::mip_;
				if(e1 > 2.0 && e2 > 2.0)
				{
					deltaE = PreshowerTools::gamma_*(PreshowerTools::calib_planeX_*e1 + PreshowerTools::calib_planeY_*e2);
				}
			}
//		cout<<" e3x3 = "<<e3x3<<"   deltaE = "<<deltaE<<endl;
		}

//
		double tempEnergy = e3x3 + deltaE;
	
		if(isPi0_)
		{	
		eeclusters_Pi0_.push_back( CaloCluster( tempEnergy, clusPos, CaloID(CaloID::DET_ECAL_ENDCAP), enFracs, CaloCluster::undefined, seed_id ) );
		eeSeedTime_Pi0_.push_back( itseed->time() );	
		eeNxtal_Pi0_.push_back(RecHitsInWindow.size());
		eeS4S9_Pi0_.push_back(s4s9);		
		eeS2S9_Pi0_.push_back(s2s9);		
		eeS1S9_Pi0_.push_back((maxEne)/e3x3);		
		}
		else
		{	
		eeclusters_Eta_.push_back( CaloCluster( tempEnergy, clusPos, CaloID(CaloID::DET_ECAL_ENDCAP), enFracs, CaloCluster::undefined, seed_id ) );
		eeSeedTime_Eta_.push_back( itseed->time() );	
		eeNxtal_Eta_.push_back(RecHitsInWindow.size());
		eeS4S9_Eta_.push_back(s4s9);		
		eeS2S9_Eta_.push_back(s2s9);		
		eeS1S9_Eta_.push_back((maxEne)/e3x3);		
		}

		//fill photon cluster
		pho_E = tempEnergy;
		pho_seedE = maxEne;
		pho_Eta = clusPos.eta();
		pho_iX = seed_id.ix();
		pho_Phi = clusPos.phi();
		pho_iY = seed_id.iy();
		pho_Pt = tempEnergy/cosh(clusPos.eta());
		pho_SeedTime = itseed->time();
		pho_ClusterTime = -999.9;//to be constructed
		pho_S4S9 = s4s9;
		pho_S2S9 = s2s9;
		pho_S1S9 = (maxEne)/e3x3; 
		pho_Nxtal = RecHitsInWindow.size();
		pho_x = clusPos.x();
		pho_y = clusPos.y();
		pho_z = clusPos.z();
		if(FillPhotonNtuple_ ) PhoEvents->Fill();
		N_Pho_rec ++;
		N_eePho_rec ++;
	
		if(isPi0_)
		{
		N_Pho_rec_Pi0_ ++;
		N_eePho_rec_Pi0_ ++;
		}
		else
		{
		N_Pho_rec_Eta_ ++;
		N_eePho_rec_Eta_ ++;
		}
	
	}//end loop of all seed xtals	
	

}


void Pi0Tuplizer::recoDiPhoEvents_EB(bool isPi0_)
{
	//N_Pair_rec = 0;

	int i=0;
	
	cout<<"DEBUG recoDiphoEvents 000  _EB  "<<isPi0_<<endl;
	if(MCAssoc_ && (N_Pi0_match==0) && (N_Eta_match==0)) return;	
	cout<<"DEBUG recoDiphoEvents 001  _EB  "<<isPi0_<<endl;

	for(unsigned int i=0;i<ebclusters_Pi0_MC1_index.size();i++)
	{
		if(isPi0_ && ebclusters_Pi0_MC1_index[i]>=0) cout<<"eb pi0 cluster "<<i<<" matched to GEN pho1 of index "<<ebclusters_Pi0_MC1_index[i]<<endl;
		if(isPi0_ && ebclusters_Pi0_MC2_index[i]>=0) cout<<"eb pi0 cluster "<<i<<" matched to GEN pho2 of index "<<ebclusters_Pi0_MC2_index[i]<<endl;
	}

	for(std::vector<CaloCluster>::const_iterator g1  = (isPi0_? ebclusters_Pi0_.begin() : ebclusters_Eta_.begin() ); g1 != (isPi0_? ebclusters_Pi0_.end() : ebclusters_Eta_.end()); ++g1, ++i)
  	{
		//cout<<"DEBUG recoDiphoEvents 002-1   -  ebcluster "<<i<<endl;
		if(MCAssoc_ && ( (isPi0_? ebclusters_Pi0_MC1_index[i] : ebclusters_Eta_MC1_index[i]) <0) && ( (isPi0_ ? ebclusters_Pi0_MC2_index[i] : ebclusters_Eta_MC2_index[i]) <0) ) continue;
		int MC_index_i = isPi0_? ( (ebclusters_Pi0_MC1_index[i]>=0) ? ebclusters_Pi0_MC1_index[i]: ebclusters_Pi0_MC2_index[i]): ( (ebclusters_Eta_MC1_index[i]>=0) ? ebclusters_Eta_MC1_index[i]: ebclusters_Eta_MC2_index[i]);
		
		cout<<"DEBUG recoDiphoEvents 002   -  ebcluster "<<i<<" GammaMC "<<MC_index_i<<endl;

		if(g1->seed().subdetId()!=1) continue;

		int j=i+1;
		for(std::vector<CaloCluster>::const_iterator g2 = g1+1; g2 != (isPi0_? ebclusters_Pi0_.end() : ebclusters_Eta_.end()); ++g2, ++j ) 
		{
			if(MCAssoc_ && isPi0_ && !(ebclusters_Pi0_MC1_index[j]==MC_index_i || ebclusters_Pi0_MC2_index[j]==MC_index_i ) ) continue;
			if(MCAssoc_ && (!isPi0_) && !(ebclusters_Eta_MC1_index[j]==MC_index_i || ebclusters_Eta_MC2_index[j]==MC_index_i )) continue;

			cout<<"DEBUG recoDiphoEvents 003   -  ebcluster "<<j<<endl;

			if(g2->seed().subdetId()!=1) continue;
			int ind1 = i;//
        		int ind2 = j;//
        		bool Inverted=false;// if pt(i)>pt(j), Inverted = false; else true
			bool MCInverted=false;
			if((isPi0_? ebclusters_Pi0_MC2_index[i] : ebclusters_Eta_MC2_index[i]) < 0) MCInverted = true;
 
 			//if keep this part, then the leading photon (g1) is the one with higher pt; 
 			//otherwise the leading photon is the one with higher seed energy
			if( PhotonOrderOption_ == "PhoPtBased" && g1->energy()/cosh(g1->eta()) < g2->energy()/cosh(g2->eta()) ) 
			{
				Inverted=true;
				ind1 = j;
				ind2 = i;
			}
		
			float Corr1 = 1.0;
			float Corr2 = 1.0;// whatever corrections to be added in the future
	
		 	math::PtEtaPhiMLorentzVector g1P4( (Corr1*g1->energy())/cosh(g1->eta()), g1->eta(), g1->phi(), 0. );
        		math::PtEtaPhiMLorentzVector g2P4( (Corr2*g2->energy())/cosh(g2->eta()), g2->eta(), g2->phi(), 0. );
			math::PtEtaPhiMLorentzVector pi0P4 = g1P4 + g2P4;
			math::PtEtaPhiMLorentzVector g1P4_nocor( (g1->energy())/cosh(g1->eta()), g1->eta(), g1->phi(), 0. );
        		math::PtEtaPhiMLorentzVector g2P4_nocor( (g2->energy())/cosh(g2->eta()), g2->eta(), g2->phi(), 0. );
        		math::PtEtaPhiMLorentzVector pi0P4_nocor = g1P4_nocor + g2P4_nocor;

			if( pi0P4_nocor.mass()<0.05 && pi0P4.mass() < 0.05 ) continue;
			cout<<"DEBUG recoDiphoEvents 004"<<endl;
			//apply kinamatics cut on diphoton and nxtal cut			
			if(fabs( pi0P4.eta() ) < 1.0 ) 
			{
				if(pi0P4.Pt()<(isPi0_?pairPtCut_barrel1_Pi0_:pairPtCut_barrel1_Eta_)) continue;
				if((isPi0_? ebNxtal_Pi0_[ind1] : ebNxtal_Eta_[ind1])<(isPi0_?nxtal1Cut_barrel1_Pi0_:nxtal1Cut_barrel1_Eta_)) continue;
				if((isPi0_? ebNxtal_Pi0_[ind2] : ebNxtal_Eta_[ind2])<(isPi0_?nxtal2Cut_barrel1_Pi0_:nxtal2Cut_barrel1_Eta_)) continue;
			}
			else if (fabs( pi0P4.eta() ) < 1.5 ) 
			{
				if(pi0P4.Pt()<(isPi0_?pairPtCut_barrel2_Pi0_:pairPtCut_barrel2_Eta_)) continue;
				if((isPi0_? ebNxtal_Pi0_[ind1] : ebNxtal_Eta_[ind1])<(isPi0_?nxtal1Cut_barrel2_Pi0_:nxtal1Cut_barrel2_Eta_)) continue;
				if((isPi0_? ebNxtal_Pi0_[ind2] : ebNxtal_Eta_[ind2])<(isPi0_?nxtal2Cut_barrel2_Pi0_:nxtal2Cut_barrel2_Eta_)) continue;
			}
			else 
			{
				continue;
			}

			if( g1P4.eta() == g2P4.eta() && g1P4.phi() == g2P4.phi() ) continue;
			cout<<"DEBUG recoDiphoEvents 005"<<endl;
			
			//calculate diphoton isolation
			double isoPi0_temp = 0.0;
			double isoG1_temp = 0.0;
			double isoG2_temp = 0.0;

			for(std::vector<CaloCluster>::const_iterator giso  = (isPi0_? ebclusters_Pi0_.begin() : ebclusters_Eta_.begin() ); giso != (isPi0_? ebclusters_Pi0_.end() : ebclusters_Eta_.end()); ++giso)
			{
				
				double dR_iso_Pair = GetDeltaR(pi0P4.Eta(), giso->eta(), pi0P4.Phi(), giso->phi());	
				double dR_iso_G1 = GetDeltaR(g1P4.Eta(), giso->eta(), g1P4.Phi(), giso->phi());	
				double dR_iso_G2 = GetDeltaR(g2P4.Eta(), giso->eta(), g2P4.Phi(), giso->phi());	
				double dEta_iso_Pair = fabs(pi0P4.Eta() - giso->eta());
				double dEta_iso_G1 = fabs(g1P4.Eta() - giso->eta());
				double dEta_iso_G2 = fabs(g2P4.Eta() - giso->eta());
				double pTiso = giso->energy()/cosh(giso->eta());

				if(giso->seed()!=g1->seed() && 
					dR_iso_G1 <= (isPi0_?isoGammaBeltdR_Zone_Pi0_:isoGammaBeltdR_Zone_Eta_) && 
					dEta_iso_G1 <= (isPi0_?isoGammaBeltdEta_Zone_Pi0_:isoGammaBeltdEta_Zone_Eta_)
				  )
				{
					isoG1_temp += pTiso;
				}
			
				if(giso->seed()!=g2->seed() && 
					dR_iso_G2 <= (isPi0_?isoGammaBeltdR_Zone_Pi0_:isoGammaBeltdR_Zone_Eta_) && 
					dEta_iso_G2 <= (isPi0_?isoGammaBeltdEta_Zone_Pi0_:isoGammaBeltdEta_Zone_Eta_)
				  )
				{
					isoG2_temp += pTiso;
				}
			

				if(giso->seed()!=g1->seed() && 
					giso->seed()!=g2->seed() &&
					dR_iso_Pair <= (isPi0_?isoPairBeltdR_Zone_Pi0_:isoPairBeltdR_Zone_Eta_) && 
					dEta_iso_Pair <= (isPi0_?isoPairBeltdEta_Zone_Pi0_:isoPairBeltdEta_Zone_Eta_)
				  )
				{
					isoPi0_temp += pTiso;
				}
			
			}	

			//apply isolation cut
			if(isoPi0_temp > isoPairCut_) continue;
			if(isoG1_temp > isoGammaCut_) continue;
			if(isoG2_temp > isoGammaCut_) continue;
			cout<<"DEBUG recoDiphoEvents 006"<<endl;

			//fill pi0/eta ntuple
			if(N_Pair_rec >= NPI0MAX-1) break; // too many pi0s
			if( FillDiPhotonNtuple_ && pi0P4.mass() > ((isPi0_)?0.05:0.1) && pi0P4.mass() < ((isPi0_)?0.35:0.9) )
			{
				fromPi0[N_Pair_rec]  =  isPi0_;
				mPi0_rec[N_Pair_rec]  =  pi0P4.mass();
				ptPi0_rec[N_Pair_rec] =  pi0P4.Pt();
				isoPi0_rec[N_Pair_rec] =  isoPi0_temp;
				etaPi0_rec[N_Pair_rec] =  pi0P4.Eta();
				phiPi0_rec[N_Pair_rec] =  pi0P4.Phi();

				if(Inverted)
				{
					isoG1_rec[N_Pair_rec] = isoG2_temp;
					isoG2_rec[N_Pair_rec] = isoG1_temp;
				}
				else
				{
					isoG1_rec[N_Pair_rec] = isoG1_temp;
					isoG2_rec[N_Pair_rec] = isoG2_temp;
				}	
				deltaRG1G2_rec[N_Pair_rec] = GetDeltaR( g1P4.eta(), g2P4.eta(), g1P4.phi(), g2P4.phi() );

				EBDetId  id_1(g1->seed());
              			EBDetId  id_2(g2->seed());
              			if(Inverted)
                		{
                  			id_1 = g2->seed();
                  			id_2 = g1->seed();
                		}

				if(!Inverted)
				{
					enG1_rec[N_Pair_rec] =  g1P4.E();
					enG2_rec[N_Pair_rec] =  g2P4.E();
					etaG1_rec[N_Pair_rec] =  g1P4.Eta();
					etaG2_rec[N_Pair_rec] =  g2P4.Eta();
					phiG1_rec[N_Pair_rec] =  g1P4.Phi();
					phiG2_rec[N_Pair_rec] =  g2P4.Phi();
					ptG1_rec[N_Pair_rec] =  g1P4.Pt();
					ptG2_rec[N_Pair_rec] =  g2P4.Pt();
				}
				else
				{
					enG1_rec[N_Pair_rec] =  g2P4.E();
					enG2_rec[N_Pair_rec] =  g1P4.E();
					etaG1_rec[N_Pair_rec] =  g2P4.Eta();
					etaG2_rec[N_Pair_rec] =  g1P4.Eta();
					phiG1_rec[N_Pair_rec] =  g2P4.Phi();
					phiG2_rec[N_Pair_rec] =  g1P4.Phi();
					ptG1_rec[N_Pair_rec] =  g2P4.Pt();
					ptG2_rec[N_Pair_rec] =  g1P4.Pt();
				}
				
				iEtaG1_rec[N_Pair_rec] =  id_1.ieta();
				iEtaG2_rec[N_Pair_rec] =  id_2.ieta();
				iPhiG1_rec[N_Pair_rec] =  id_1.iphi();
				iPhiG1_rec[N_Pair_rec] =  id_2.iphi();
				
				cout<<"DEBUG recoDiphoEvents 007"<<endl;
				if((Inverted && MCInverted ) || ( (!Inverted) && (!MCInverted)))
				{
					if(isPi0_ && ebclusters_Pi0_MC1_index[ind1]>=0) 
					{
						enG1_true[N_Pair_rec] = Gamma1MC_Pi0_[ebclusters_Pi0_MC1_index[ind1]].E();
						enG2_true[N_Pair_rec] = Gamma2MC_Pi0_[ebclusters_Pi0_MC1_index[ind1]].E();
					}
					if((!isPi0_) && ebclusters_Eta_MC1_index[ind1]>=0) 
					{
						enG1_true[N_Pair_rec] = Gamma1MC_Eta_[ebclusters_Eta_MC1_index[ind1]].E();
						enG2_true[N_Pair_rec] = Gamma2MC_Eta_[ebclusters_Eta_MC1_index[ind1]].E();
					}
					if(isPi0_ && ebclusters_Pi0_MC2_index[ind1]>=0) 
					{
						enG1_true[N_Pair_rec] = Gamma1MC_Pi0_[ebclusters_Pi0_MC2_index[ind1]].E();
						enG2_true[N_Pair_rec] = Gamma2MC_Pi0_[ebclusters_Pi0_MC2_index[ind1]].E();
					}
					if((!isPi0_) && ebclusters_Eta_MC2_index[ind1]>=0) 
					{
						enG1_true[N_Pair_rec] = Gamma1MC_Eta_[ebclusters_Eta_MC2_index[ind1]].E();
						enG2_true[N_Pair_rec] = Gamma2MC_Eta_[ebclusters_Eta_MC2_index[ind1]].E();
					}
				}
				else
				{
					if(isPi0_ && ebclusters_Pi0_MC1_index[ind1]>=0) 
					{
						enG1_true[N_Pair_rec] = Gamma2MC_Pi0_[ebclusters_Pi0_MC1_index[ind1]].E();
						enG2_true[N_Pair_rec] = Gamma1MC_Pi0_[ebclusters_Pi0_MC1_index[ind1]].E();
					}
					if((!isPi0_) && ebclusters_Eta_MC1_index[ind1]>=0) 
					{
						enG1_true[N_Pair_rec] = Gamma2MC_Eta_[ebclusters_Eta_MC1_index[ind1]].E();
						enG2_true[N_Pair_rec] = Gamma1MC_Eta_[ebclusters_Eta_MC1_index[ind1]].E();
					}
					if(isPi0_ && ebclusters_Pi0_MC2_index[ind1]>=0) 
					{
						enG1_true[N_Pair_rec] = Gamma2MC_Pi0_[ebclusters_Pi0_MC2_index[ind1]].E();
						enG2_true[N_Pair_rec] = Gamma1MC_Pi0_[ebclusters_Pi0_MC2_index[ind1]].E();
					}
					if((!isPi0_) && ebclusters_Eta_MC2_index[ind1]>=0) 
					{
						enG1_true[N_Pair_rec] = Gamma2MC_Eta_[ebclusters_Eta_MC2_index[ind1]].E();
						enG2_true[N_Pair_rec] = Gamma1MC_Eta_[ebclusters_Eta_MC2_index[ind1]].E();
					}

				}

				if(isPi0_)
				{
				seedTimeG1_rec[N_Pair_rec] = ebSeedTime_Pi0_[ind1];	
				seedTimeG2_rec[N_Pair_rec] = ebSeedTime_Pi0_[ind2];	
				s4s9G1_rec[N_Pair_rec] = ebS4S9_Pi0_[ind1];	
				s4s9G2_rec[N_Pair_rec] = ebS4S9_Pi0_[ind2];	
				s2s9G1_rec[N_Pair_rec] = ebS2S9_Pi0_[ind1];	
				s2s9G2_rec[N_Pair_rec] = ebS2S9_Pi0_[ind2];	
				s1s9G1_rec[N_Pair_rec] = ebS1S9_Pi0_[ind1];	
				s1s9G2_rec[N_Pair_rec] = ebS1S9_Pi0_[ind2];	
				nxtalG1_rec[N_Pair_rec] = ebNxtal_Pi0_[ind1];	
				nxtalG2_rec[N_Pair_rec] = ebNxtal_Pi0_[ind2];	
				}
				else
				{
				seedTimeG1_rec[N_Pair_rec] = ebSeedTime_Eta_[ind1];	
				seedTimeG2_rec[N_Pair_rec] = ebSeedTime_Eta_[ind2];	
				s4s9G1_rec[N_Pair_rec] = ebS4S9_Eta_[ind1];	
				s4s9G2_rec[N_Pair_rec] = ebS4S9_Eta_[ind2];	
				s2s9G1_rec[N_Pair_rec] = ebS2S9_Eta_[ind1];	
				s2s9G2_rec[N_Pair_rec] = ebS2S9_Eta_[ind2];	
				s1s9G1_rec[N_Pair_rec] = ebS1S9_Eta_[ind1];	
				s1s9G2_rec[N_Pair_rec] = ebS1S9_Eta_[ind2];	
				nxtalG1_rec[N_Pair_rec] = ebNxtal_Eta_[ind1];	
				nxtalG2_rec[N_Pair_rec] = ebNxtal_Eta_[ind2];	
				}
				N_Pair_rec ++;			
			}
			if(N_Pair_rec >= NPI0MAX-1) break; // too many pi0s
		
		}//end loop of g2	
	}//end loop of g1
}

void Pi0Tuplizer::recoDiPhoEvents_EE(bool isPi0_)
{


	int i=0;
	
	if(MCAssoc_ && (N_Pi0_match==0) && (N_Eta_match==0)) return;	

	cout<<"DEBUG recoDiphoEvents 001  _EE"<<endl;
/*
	for(unsigned int i=0;i<eeclusters_Pi0_MC1_index.size();i++)
	{
		if(eeclusters_Pi0_MC1_index[i]>=0) cout<<"eb pi0 cluster "<<i<<" matched to GEN pho1 of index "<<eeclusters_Pi0_MC1_index[i]<<endl;
		if(eeclusters_Pi0_MC2_index[i]>=0) cout<<"eb pi0 cluster "<<i<<" matched to GEN pho2 of index "<<eeclusters_Pi0_MC2_index[i]<<endl;
	}
*/
	for(std::vector<CaloCluster>::const_iterator g1  = (isPi0_? eeclusters_Pi0_.begin() : eeclusters_Eta_.begin() ); g1 != (isPi0_? eeclusters_Pi0_.end() : eeclusters_Eta_.end()); ++g1, ++i)
  	{
		if(MCAssoc_ && ((isPi0_? eeclusters_Pi0_MC1_index[i] : eeclusters_Eta_MC1_index[i]) <0) && ((isPi0_ ? eeclusters_Pi0_MC2_index[i] : eeclusters_Eta_MC2_index[i]) <0) ) continue;
		if(g1->seed().subdetId()!=2) continue;
		int MC_index_i = isPi0_? ((eeclusters_Pi0_MC1_index[i]>=0) ? eeclusters_Pi0_MC1_index[i]: eeclusters_Pi0_MC2_index[i]): ((eeclusters_Eta_MC1_index[i]>=0) ? eeclusters_Eta_MC1_index[i]: eeclusters_Eta_MC2_index[i]);

		cout<<"DEBUG recoDiphoEvents 002   -  eecluster "<<i<<" GammaMC "<<MC_index_i<<endl;

		int j=i+1;
		for(std::vector<CaloCluster>::const_iterator g2 = g1+1; g2 != (isPi0_? eeclusters_Pi0_.end() : eeclusters_Eta_.end()); ++g2, ++j ) 
		{

			if(MCAssoc_ && isPi0_ && !(eeclusters_Pi0_MC1_index[j]==MC_index_i || eeclusters_Pi0_MC2_index[j]==MC_index_i ) ) continue;
                        if(MCAssoc_ && (!isPi0_) && !(eeclusters_Eta_MC1_index[j]==MC_index_i || eeclusters_Eta_MC2_index[j]==MC_index_i )) continue;
			cout<<"DEBUG recoDiphoEvents 003 "<<endl;


			if(g2->seed().subdetId()!=2) continue;
			int ind1 = i;//
        		int ind2 = j;//
        		bool Inverted=false;// if pt(i)>pt(j), Inverted = false; else true
			 
 			//if keep this part, then the leading photon (g1) is the one with higher pt; 
 			//otherwise the leading photon is the one with higher seed energy
			if( PhotonOrderOption_ == "PhoPtBased" && g1->energy()/cosh(g1->eta()) < g2->energy()/cosh(g2->eta()) ) 
			{
				Inverted=true;
				ind1 = j;
				ind2 = i;
			}
		
			float Corr1 = 1.0;
			float Corr2 = 1.0;// whatever corrections to be added in the future
	
		 	math::PtEtaPhiMLorentzVector g1P4( (Corr1*g1->energy())/cosh(g1->eta()), g1->eta(), g1->phi(), 0. );
        		math::PtEtaPhiMLorentzVector g2P4( (Corr2*g2->energy())/cosh(g2->eta()), g2->eta(), g2->phi(), 0. );
			math::PtEtaPhiMLorentzVector pi0P4 = g1P4 + g2P4;
			math::PtEtaPhiMLorentzVector g1P4_nocor( (g1->energy())/cosh(g1->eta()), g1->eta(), g1->phi(), 0. );
        		math::PtEtaPhiMLorentzVector g2P4_nocor( (g2->energy())/cosh(g2->eta()), g2->eta(), g2->phi(), 0. );
        		math::PtEtaPhiMLorentzVector pi0P4_nocor = g1P4_nocor + g2P4_nocor;

			if( pi0P4_nocor.mass()<0.05 && pi0P4.mass() < 0.05 ) continue;
			cout<<"DEBUG recoDiphoEvents 004"<<endl;
			//apply kinamatics cut on diphoton and nxtal cut			
			if(fabs( pi0P4.eta() ) < 1.5 ) continue;
			else if(fabs( pi0P4.eta() ) < 1.8 ) 
			{
				if(pi0P4.Pt()<(isPi0_?pairPtCut_endcap1_Pi0_:pairPtCut_endcap1_Eta_)) continue;
				if((isPi0_? eeNxtal_Pi0_[ind1] : eeNxtal_Eta_[ind1])<(isPi0_?nxtal1Cut_endcap1_Pi0_:nxtal1Cut_endcap1_Eta_)) continue;
				if((isPi0_? eeNxtal_Pi0_[ind2] : eeNxtal_Eta_[ind2])<(isPi0_?nxtal2Cut_endcap1_Pi0_:nxtal2Cut_endcap1_Eta_)) continue;
		
			}
			else if (fabs( pi0P4.eta() ) < 3.0 ) 
			{
				if(pi0P4.Pt()<(isPi0_?pairPtCut_endcap2_Pi0_:pairPtCut_endcap2_Eta_)) continue;
				if((isPi0_? eeNxtal_Pi0_[ind1] : eeNxtal_Eta_[ind1])<(isPi0_?nxtal1Cut_endcap2_Pi0_:nxtal1Cut_endcap2_Eta_)) continue;
				if((isPi0_? eeNxtal_Pi0_[ind2] : eeNxtal_Eta_[ind2])<(isPi0_?nxtal2Cut_endcap2_Pi0_:nxtal2Cut_endcap2_Eta_)) continue;
			}
			else 
			{
				continue;
			}

			if( g1P4.eta() == g2P4.eta() && g1P4.phi() == g2P4.phi() ) continue;
			cout<<"DEBUG recoDiphoEvents 005"<<endl;
			
			//calculate diphoton isolation
			double isoPi0_temp = 0.0;
			double isoG1_temp = 0.0;
			double isoG2_temp = 0.0;

			for(std::vector<CaloCluster>::const_iterator giso  = (isPi0_? eeclusters_Pi0_.begin() : eeclusters_Eta_.begin() ); giso != (isPi0_? eeclusters_Pi0_.end() : eeclusters_Eta_.end()); ++giso)
			{
				
				double dR_iso_Pair = GetDeltaR(pi0P4.Eta(), giso->eta(), pi0P4.Phi(), giso->phi());	
				double dR_iso_G1 = GetDeltaR(g1P4.Eta(), giso->eta(), g1P4.Phi(), giso->phi());	
				double dR_iso_G2 = GetDeltaR(g2P4.Eta(), giso->eta(), g2P4.Phi(), giso->phi());	
				double dEta_iso_Pair = fabs(pi0P4.Eta() - giso->eta());
				double dEta_iso_G1 = fabs(g1P4.Eta() - giso->eta());
				double dEta_iso_G2 = fabs(g2P4.Eta() - giso->eta());
				double pTiso = giso->energy()/cosh(giso->eta());

				if(giso->seed()!=g1->seed() && 
					dR_iso_G1 <= (isPi0_?isoGammaBeltdR_Zone_Pi0_:isoGammaBeltdR_Zone_Eta_) && 
					dEta_iso_G1 <= (isPi0_?isoGammaBeltdEta_Zone_Pi0_:isoGammaBeltdEta_Zone_Eta_)
				  )
				{
					isoG1_temp += pTiso;
				}
			
				if(giso->seed()!=g2->seed() && 
					dR_iso_G2 <= (isPi0_?isoGammaBeltdR_Zone_Pi0_:isoGammaBeltdR_Zone_Eta_) && 
					dEta_iso_G2 <= (isPi0_?isoGammaBeltdEta_Zone_Pi0_:isoGammaBeltdEta_Zone_Eta_)
				  )
				{
					isoG2_temp += pTiso;
				}
			

				if(giso->seed()!=g1->seed() && 
					giso->seed()!=g2->seed() &&
					dR_iso_Pair <= (isPi0_?isoPairBeltdR_Zone_Pi0_:isoPairBeltdR_Zone_Eta_) && 
					dEta_iso_Pair <= (isPi0_?isoPairBeltdEta_Zone_Pi0_:isoPairBeltdEta_Zone_Eta_)
				  )
				{
					isoPi0_temp += pTiso;
				}
			
			}	

			//apply isolation cut
			if(isoPi0_temp > isoPairCut_) continue;
			if(isoG1_temp > isoGammaCut_) continue;
			if(isoG2_temp > isoGammaCut_) continue;

			cout<<"DEBUG recoDiphoEvents 006"<<endl;

			//fill pi0/eta ntuple
			if(N_Pair_rec >= NPI0MAX-1) break; // too many pi0s
			if( FillDiPhotonNtuple_ && pi0P4.mass() > ((isPi0_)?0.05:0.10) && pi0P4.mass() < ((isPi0_)?0.35:0.9) )
			{
				fromPi0[N_Pair_rec]  =  isPi0_;
				mPi0_rec[N_Pair_rec]  =  pi0P4.mass();
				ptPi0_rec[N_Pair_rec] =  pi0P4.Pt();
				isoPi0_rec[N_Pair_rec] =  isoPi0_temp;
                                etaPi0_rec[N_Pair_rec] =  pi0P4.Eta();
                                phiPi0_rec[N_Pair_rec] =  pi0P4.Phi();

                                if(Inverted)
                                {
                                        isoG1_rec[N_Pair_rec] = isoG2_temp;
                                        isoG2_rec[N_Pair_rec] = isoG1_temp;
                                }
                                else
                                {
                                        isoG1_rec[N_Pair_rec] = isoG1_temp;
                                        isoG2_rec[N_Pair_rec] = isoG2_temp;
                                }

				deltaRG1G2_rec[N_Pair_rec] = GetDeltaR( g1P4.eta(), g2P4.eta(), g1P4.phi(), g2P4.phi() );

				EEDetId  id_1(g1->seed());
              			EEDetId  id_2(g2->seed());
              			if(Inverted)
                		{
                  			id_1 = g2->seed();
                  			id_2 = g1->seed();
                		}

				if(!Inverted)
				{
					enG1_rec[N_Pair_rec] =  g1P4.E();
					enG2_rec[N_Pair_rec] =  g2P4.E();
					etaG1_rec[N_Pair_rec] =  g1P4.Eta();
					etaG2_rec[N_Pair_rec] =  g2P4.Eta();
					phiG1_rec[N_Pair_rec] =  g1P4.Phi();
					phiG2_rec[N_Pair_rec] =  g2P4.Phi();
					ptG1_rec[N_Pair_rec] =  g1P4.Pt();
					ptG2_rec[N_Pair_rec] =  g2P4.Pt();
				}
				else
				{
					enG1_rec[N_Pair_rec] =  g2P4.E();
					enG2_rec[N_Pair_rec] =  g1P4.E();
					etaG1_rec[N_Pair_rec] =  g2P4.Eta();
					etaG2_rec[N_Pair_rec] =  g1P4.Eta();
					phiG1_rec[N_Pair_rec] =  g2P4.Phi();
					phiG2_rec[N_Pair_rec] =  g1P4.Phi();
					ptG1_rec[N_Pair_rec] =  g2P4.Pt();
					ptG2_rec[N_Pair_rec] =  g1P4.Pt();
				}
				
				iXG1_rec[N_Pair_rec] =  id_1.ix();
				iXG2_rec[N_Pair_rec] =  id_2.ix();
				iYG1_rec[N_Pair_rec] =  id_1.iy();
				iYG1_rec[N_Pair_rec] =  id_2.iy();
	
				if(isPi0_ && eeclusters_Pi0_MC1_index[ind1]>=0) 
				{
					enG1_true[N_Pair_rec] = Gamma1MC_Pi0_[eeclusters_Pi0_MC1_index[ind1]].E();
					enG2_true[N_Pair_rec] = Gamma2MC_Pi0_[eeclusters_Pi0_MC1_index[ind1]].E();
				}
				if((!isPi0_) && eeclusters_Eta_MC1_index[ind1]>=0) 
				{
					enG1_true[N_Pair_rec] = Gamma1MC_Eta_[eeclusters_Eta_MC1_index[ind1]].E();
					enG2_true[N_Pair_rec] = Gamma2MC_Eta_[eeclusters_Eta_MC1_index[ind1]].E();
				}
				if(isPi0_ && eeclusters_Pi0_MC2_index[ind1]>=0) 
				{
					enG1_true[N_Pair_rec] = Gamma2MC_Pi0_[eeclusters_Pi0_MC2_index[ind1]].E();
					enG2_true[N_Pair_rec] = Gamma1MC_Pi0_[eeclusters_Pi0_MC2_index[ind1]].E();
				}
				if((!isPi0_) && eeclusters_Eta_MC2_index[ind1]>=0) 
				{
					enG1_true[N_Pair_rec] = Gamma2MC_Eta_[eeclusters_Eta_MC2_index[ind1]].E();
					enG2_true[N_Pair_rec] = Gamma1MC_Eta_[eeclusters_Eta_MC2_index[ind1]].E();
				}
			
				if(isPi0_)
				{
				seedTimeG1_rec[N_Pair_rec] = eeSeedTime_Pi0_[ind1];	
				seedTimeG2_rec[N_Pair_rec] = eeSeedTime_Pi0_[ind2];	
				s4s9G1_rec[N_Pair_rec] = eeS4S9_Pi0_[ind1];	
				s4s9G2_rec[N_Pair_rec] = eeS4S9_Pi0_[ind2];	
				s2s9G1_rec[N_Pair_rec] = eeS2S9_Pi0_[ind1];	
				s2s9G2_rec[N_Pair_rec] = eeS2S9_Pi0_[ind2];	
				s1s9G1_rec[N_Pair_rec] = eeS1S9_Pi0_[ind1];	
				s1s9G2_rec[N_Pair_rec] = eeS1S9_Pi0_[ind2];	
				nxtalG1_rec[N_Pair_rec] = eeNxtal_Pi0_[ind1];	
				nxtalG2_rec[N_Pair_rec] = eeNxtal_Pi0_[ind2];	
				}
				else
				{
				seedTimeG1_rec[N_Pair_rec] = eeSeedTime_Eta_[ind1];	
				seedTimeG2_rec[N_Pair_rec] = eeSeedTime_Eta_[ind2];	
				s4s9G1_rec[N_Pair_rec] = eeS4S9_Eta_[ind1];	
				s4s9G2_rec[N_Pair_rec] = eeS4S9_Eta_[ind2];	
				s2s9G1_rec[N_Pair_rec] = eeS2S9_Eta_[ind1];	
				s2s9G2_rec[N_Pair_rec] = eeS2S9_Eta_[ind2];	
				s1s9G1_rec[N_Pair_rec] = eeS1S9_Eta_[ind1];	
				s1s9G2_rec[N_Pair_rec] = eeS1S9_Eta_[ind2];	
				nxtalG1_rec[N_Pair_rec] = eeNxtal_Eta_[ind1];	
				nxtalG2_rec[N_Pair_rec] = eeNxtal_Eta_[ind2];	
				}
				cout<<"DEBUG recoDiphoEvents 007"<<endl;
				N_Pair_rec ++;			
			}
			if(N_Pair_rec >= NPI0MAX-1) break; // too many pi0s
		}//end loop of g2	
	}//end loop of g1

}



void Pi0Tuplizer::GetL1SeedBit()
{
	if(uGtAlg.isValid()) 
	{
        	for(int bx=uGtAlg->getFirstBX(); bx<=uGtAlg->getLastBX(); bx++) 
		{
        		for(std::vector<GlobalAlgBlk>::const_iterator algBlk = uGtAlg->begin(bx); algBlk != uGtAlg->end(bx); ++algBlk) 
			{
                		for(uint i=0;i<NL1SEED;i++)
				{
                			if(algBlk->getAlgoDecisionFinal(i)) 
					{
					allL1SeedFinalDecision[i] = true;
					L1SeedBitFinalDecision->push_back(i);
					}
                		}
              		}
           	}
        }
	
}


void Pi0Tuplizer::GetMCTruth()
{
	for(size_t iG=0; iG<genParticles->size();iG++)
        {
		//cout<<"DEBUG   genP:  "<<iG<<"   id= "<<(*genParticles)[iG].pdgId()<<"  status= "<<(*genParticles)[iG].status()<<endl;
			
		if((*genParticles)[iG].status()!=2) continue;

		unsigned int ndau = (*genParticles)[iG].numberOfDaughters();
		if(((*genParticles)[iG].pdgId() != 111) && ((*genParticles)[iG].pdgId() != 221)) continue;

		if((*genParticles)[iG].pdgId() == 111)
		{
			ptPi0_genall[N_Pi0_genall] = (*genParticles)[iG].pt(); 	
			etaPi0_genall[N_Pi0_genall] = (*genParticles)[iG].p4().Eta(); 	
			phiPi0_genall[N_Pi0_genall] = (*genParticles)[iG].p4().Phi(); 		
			N_Pi0_genall ++;
		}

		if((*genParticles)[iG].pdgId() == 221)
		{
			ptEta_genall[N_Eta_genall] = (*genParticles)[iG].pt(); 	
			etaEta_genall[N_Eta_genall] = (*genParticles)[iG].p4().Eta(); 	
			phiEta_genall[N_Eta_genall] = (*genParticles)[iG].p4().Phi(); 		
			N_Eta_genall ++;
		}

		
		if(ndau != 2 ) continue;	
		bool isDiphoton = true;	
		for (unsigned int jD=0; jD < ndau; ++jD) 
                {
                	const reco::Candidate *dau = (*genParticles)[iG].daughter(jD);
			if(dau->pdgId() != 22) isDiphoton=false;
                }

		if(!isDiphoton) continue;

		//fill GEN pi0
		if((*genParticles)[iG].pdgId() == 111)
		{
			TLorentzVector gamma1_temp, gamma2_temp;

			const reco::Candidate *dau1 = (*genParticles)[iG].daughter(0);
			const reco::Candidate *dau2 = (*genParticles)[iG].daughter(1);
			if(dau1->pt() > dau2->pt()) 
			{
				gamma1_temp.SetPtEtaPhiE(dau1->pt(), dau1->p4().Eta(), dau1->p4().Phi(), dau1->p4().E());
				gamma2_temp.SetPtEtaPhiE(dau2->pt(), dau2->p4().Eta(), dau2->p4().Phi(), dau2->p4().E());
			}
			else	
			{
				gamma2_temp.SetPtEtaPhiE(dau1->pt(), dau1->p4().Eta(), dau1->p4().Phi(), dau1->p4().E());
				gamma1_temp.SetPtEtaPhiE(dau2->pt(), dau2->p4().Eta(), dau2->p4().Phi(), dau2->p4().E());
			}
			Gamma1MC_Pi0_.push_back(gamma1_temp);
			Gamma2MC_Pi0_.push_back(gamma2_temp);
			GammaIndex_Pi0_.push_back(N_Pi0_genall-1);
		}
		//fill GEN eta
		if((*genParticles)[iG].pdgId() == 221)
		{
			TLorentzVector gamma1_temp, gamma2_temp;

			const reco::Candidate *dau1 = (*genParticles)[iG].daughter(0);
			const reco::Candidate *dau2 = (*genParticles)[iG].daughter(1);
			if(dau1->pt() > dau2->pt()) 
			{
				gamma1_temp.SetPtEtaPhiE(dau1->pt(), dau1->p4().Eta(), dau1->p4().Phi(), dau1->p4().E());
				gamma2_temp.SetPtEtaPhiE(dau2->pt(), dau2->p4().Eta(), dau2->p4().Phi(), dau2->p4().E());
			}
			else	
			{
				gamma2_temp.SetPtEtaPhiE(dau1->pt(), dau1->p4().Eta(), dau1->p4().Phi(), dau1->p4().E());
				gamma1_temp.SetPtEtaPhiE(dau2->pt(), dau2->p4().Eta(), dau2->p4().Phi(), dau2->p4().E());
			}
			Gamma1MC_Eta_.push_back(gamma1_temp);
			Gamma2MC_Eta_.push_back(gamma2_temp);
			GammaIndex_Eta_.push_back(N_Eta_genall-1);
		}
	}

	N_Pi0_gen = Gamma2MC_Pi0_.size();
	N_Eta_gen = Gamma2MC_Eta_.size();
	
	cout<<"DEBUG  = number of true pi0s(all, gg): "<<N_Pi0_genall<<"   "<<Gamma2MC_Pi0_.size()<<endl;	
	cout<<"DEBUG  = number of true etas(all, gg): "<<N_Eta_genall<<"   "<<Gamma2MC_Eta_.size()<<endl;	

	for(int i=0;i<N_Pi0_genall && i<NPI0MAX;i++)
	{
		//calculate how many gammas from other pi0/eta are in this pi0 cone		
		for(int j=0;j<N_Pi0_gen;j++)
		{
			if(GammaIndex_Pi0_[j] == i) continue;
		
			double dR1_this = GetDeltaR(Gamma1MC_Pi0_[j].Eta(), etaPi0_genall[i], Gamma1MC_Pi0_[j].Phi(), phiPi0_genall[i]);
			double dR2_this = GetDeltaR(Gamma2MC_Pi0_[j].Eta(), etaPi0_genall[i], Gamma2MC_Pi0_[j].Phi(), phiPi0_genall[i]);
			if(dR1_this > 0.3 && dR2_this > 0.3) continue;
			
			if(dR1_this < 0.1)
			{
				nIsoGamma0p3Pi0_genall[i] += 1;
				nIsoGamma0p2Pi0_genall[i] += 1;
				nIsoGamma0p1Pi0_genall[i] += 1;
			}
			else if(dR1_this < 0.2)
			{
				
				nIsoGamma0p3Pi0_genall[i] += 1;
				nIsoGamma0p2Pi0_genall[i] += 1;
			}
			else if(dR1_this < 0.3)
			{
				nIsoGamma0p3Pi0_genall[i] += 1;
			}

			if(dR2_this < 0.1)
			{
				nIsoGamma0p3Pi0_genall[i] += 1;
				nIsoGamma0p2Pi0_genall[i] += 1;
				nIsoGamma0p1Pi0_genall[i] += 1;
			}
			else if(dR2_this < 0.2)
			{
				
				nIsoGamma0p3Pi0_genall[i] += 1;
				nIsoGamma0p2Pi0_genall[i] += 1;
			}
			else if(dR2_this < 0.3)
			{
				nIsoGamma0p3Pi0_genall[i] += 1;
			}

		}

		for(int j=0;j<N_Eta_gen;j++)
		{
			double dR1_this = GetDeltaR(Gamma1MC_Eta_[j].Eta(), etaPi0_genall[i], Gamma1MC_Eta_[j].Phi(), phiPi0_genall[i]);
			double dR2_this = GetDeltaR(Gamma2MC_Eta_[j].Eta(), etaPi0_genall[i], Gamma2MC_Eta_[j].Phi(), phiPi0_genall[i]);
			if(dR1_this > 0.3 && dR2_this > 0.3) continue;
			
			if(dR1_this < 0.1)
			{
				nIsoGamma0p3Pi0_genall[i] += 1;
				nIsoGamma0p2Pi0_genall[i] += 1;
				nIsoGamma0p1Pi0_genall[i] += 1;
			}
			else if(dR1_this < 0.2)
			{
				
				nIsoGamma0p3Pi0_genall[i] += 1;
				nIsoGamma0p2Pi0_genall[i] += 1;
			}
			else if(dR1_this < 0.3)
			{
				nIsoGamma0p3Pi0_genall[i] += 1;
			}

			if(dR2_this < 0.1)
			{
				nIsoGamma0p3Pi0_genall[i] += 1;
				nIsoGamma0p2Pi0_genall[i] += 1;
				nIsoGamma0p1Pi0_genall[i] += 1;
			}
			else if(dR2_this < 0.2)
			{
				
				nIsoGamma0p3Pi0_genall[i] += 1;
				nIsoGamma0p2Pi0_genall[i] += 1;
			}
			else if(dR2_this < 0.3)
			{
				nIsoGamma0p3Pi0_genall[i] += 1;
			}

		}

	}


	for(int i=0;i<N_Eta_genall && i<NPI0MAX;i++)
	{
		//calculate how many gammas from other pi0/eta are in this pi0 cone		
		for(int j=0;j<N_Eta_gen;j++)
		{
			if(GammaIndex_Eta_[j] == i) continue;
		
			double dR1_this = GetDeltaR(Gamma1MC_Eta_[j].Eta(), etaEta_genall[i], Gamma1MC_Eta_[j].Phi(), phiEta_genall[i]);
			double dR2_this = GetDeltaR(Gamma2MC_Eta_[j].Eta(), etaEta_genall[i], Gamma2MC_Eta_[j].Phi(), phiEta_genall[i]);
			if(dR1_this > 0.3 && dR2_this > 0.3) continue;
			
			if(dR1_this < 0.1)
			{
				nIsoGamma0p3Eta_genall[i] += 1;
				nIsoGamma0p2Eta_genall[i] += 1;
				nIsoGamma0p1Eta_genall[i] += 1;
			}
			else if(dR1_this < 0.2)
			{
				
				nIsoGamma0p3Eta_genall[i] += 1;
				nIsoGamma0p2Eta_genall[i] += 1;
			}
			else if(dR1_this < 0.3)
			{
				nIsoGamma0p3Eta_genall[i] += 1;
			}

			if(dR2_this < 0.1)
			{
				nIsoGamma0p3Eta_genall[i] += 1;
				nIsoGamma0p2Eta_genall[i] += 1;
				nIsoGamma0p1Eta_genall[i] += 1;
			}
			else if(dR2_this < 0.2)
			{
				
				nIsoGamma0p3Eta_genall[i] += 1;
				nIsoGamma0p2Eta_genall[i] += 1;
			}
			else if(dR2_this < 0.3)
			{
				nIsoGamma0p3Eta_genall[i] += 1;
			}

		}

		for(int j=0;j<N_Pi0_gen;j++)
		{
			double dR1_this = GetDeltaR(Gamma1MC_Pi0_[j].Eta(), etaEta_genall[i], Gamma1MC_Pi0_[j].Phi(), phiEta_genall[i]);
			double dR2_this = GetDeltaR(Gamma2MC_Pi0_[j].Eta(), etaEta_genall[i], Gamma2MC_Pi0_[j].Phi(), phiEta_genall[i]);
			if(dR1_this > 0.3 && dR2_this > 0.3) continue;
			
			if(dR1_this < 0.1)
			{
				nIsoGamma0p3Eta_genall[i] += 1;
				nIsoGamma0p2Eta_genall[i] += 1;
				nIsoGamma0p1Eta_genall[i] += 1;
			}
			else if(dR1_this < 0.2)
			{
				
				nIsoGamma0p3Eta_genall[i] += 1;
				nIsoGamma0p2Eta_genall[i] += 1;
			}
			else if(dR1_this < 0.3)
			{
				nIsoGamma0p3Eta_genall[i] += 1;
			}

			if(dR2_this < 0.1)
			{
				nIsoGamma0p3Eta_genall[i] += 1;
				nIsoGamma0p2Eta_genall[i] += 1;
				nIsoGamma0p1Eta_genall[i] += 1;
			}
			else if(dR2_this < 0.2)
			{
				
				nIsoGamma0p3Eta_genall[i] += 1;
				nIsoGamma0p2Eta_genall[i] += 1;
			}
			else if(dR2_this < 0.3)
			{
				nIsoGamma0p3Eta_genall[i] += 1;
			}

		}

	}




	for(int i=0;i<N_Pi0_gen && i<NPI0MAX;i++)
	{
		ptPi0_gen[i] = (Gamma1MC_Pi0_[i]+Gamma2MC_Pi0_[i]).Pt();
		etaPi0_gen[i] = (Gamma1MC_Pi0_[i]+Gamma2MC_Pi0_[i]).Eta();
		phiPi0_gen[i] = (Gamma1MC_Pi0_[i]+Gamma2MC_Pi0_[i]).Phi();
		deltaRG1G2Pi0_gen[i] = Gamma1MC_Pi0_[i].DeltaR(Gamma2MC_Pi0_[i]);
		enG1_Pi0_gen[i] = Gamma1MC_Pi0_[i].E();
		ptG1_Pi0_gen[i] = Gamma1MC_Pi0_[i].Pt();
		etaG1_Pi0_gen[i] = Gamma1MC_Pi0_[i].Eta();
		phiG1_Pi0_gen[i] = Gamma1MC_Pi0_[i].Phi();
		enG2_Pi0_gen[i] = Gamma2MC_Pi0_[i].E();
		ptG2_Pi0_gen[i] = Gamma2MC_Pi0_[i].Pt();
		etaG2_Pi0_gen[i] = Gamma2MC_Pi0_[i].Eta();
		phiG2_Pi0_gen[i] = Gamma2MC_Pi0_[i].Phi();
		//calculate how many gammas from other pi0/eta are in this pi0 cone
		for(int j=0;j<N_Pi0_gen;j++)
		{
			if(j==i) continue;
			double dR1_this = Gamma1MC_Pi0_[j].DeltaR(Gamma1MC_Pi0_[i]+Gamma2MC_Pi0_[i]);
			double dR2_this = Gamma2MC_Pi0_[j].DeltaR(Gamma1MC_Pi0_[i]+Gamma2MC_Pi0_[i]);
			if(dR1_this > 0.3 && dR2_this > 0.3) continue;

			if(dR1_this < 0.1)
			{
				nIsoGamma0p3Pi0_gen[i] += 1;
				nIsoGamma0p2Pi0_gen[i] += 1;
				nIsoGamma0p1Pi0_gen[i] += 1;
			}
			else if(dR1_this < 0.2)
			{
				
				nIsoGamma0p3Pi0_gen[i] += 1;
				nIsoGamma0p2Pi0_gen[i] += 1;
			}
			else if(dR1_this < 0.3)
			{
				nIsoGamma0p3Pi0_gen[i] += 1;
			}

			if(dR2_this < 0.1)
			{
				nIsoGamma0p3Pi0_gen[i] += 1;
				nIsoGamma0p2Pi0_gen[i] += 1;
				nIsoGamma0p1Pi0_gen[i] += 1;
			}
			else if(dR2_this < 0.2)
			{
				
				nIsoGamma0p3Pi0_gen[i] += 1;
				nIsoGamma0p2Pi0_gen[i] += 1;
			}
			else if(dR2_this < 0.3)
			{
				nIsoGamma0p3Pi0_gen[i] += 1;
			}
		}
		
		for(int j=0;j<N_Eta_gen;j++)
		{
			double dR1_this = Gamma1MC_Eta_[j].DeltaR(Gamma1MC_Pi0_[i]+Gamma2MC_Pi0_[i]);
			double dR2_this = Gamma2MC_Eta_[j].DeltaR(Gamma1MC_Pi0_[i]+Gamma2MC_Pi0_[i]);
			if(dR1_this > 0.3 && dR2_this > 0.3) continue;

			if(dR1_this < 0.1)
			{
				nIsoGamma0p3Pi0_gen[i] += 1;
				nIsoGamma0p2Pi0_gen[i] += 1;
				nIsoGamma0p1Pi0_gen[i] += 1;
			}
			else if(dR1_this < 0.2)
			{
				
				nIsoGamma0p3Pi0_gen[i] += 1;
				nIsoGamma0p2Pi0_gen[i] += 1;
			}
			else if(dR1_this < 0.3)
			{
				nIsoGamma0p3Pi0_gen[i] += 1;
			}

			if(dR2_this < 0.1)
			{
				nIsoGamma0p3Pi0_gen[i] += 1;
				nIsoGamma0p2Pi0_gen[i] += 1;
				nIsoGamma0p1Pi0_gen[i] += 1;
			}
			else if(dR2_this < 0.2)
			{
				
				nIsoGamma0p3Pi0_gen[i] += 1;
				nIsoGamma0p2Pi0_gen[i] += 1;
			}
			else if(dR2_this < 0.3)
			{
				nIsoGamma0p3Pi0_gen[i] += 1;
			}
		}


	}

	for(int i=0;i<N_Eta_gen && i<NPI0MAX;i++)
	{
		ptEta_gen[i] = (Gamma1MC_Eta_[i]+Gamma2MC_Eta_[i]).Pt();
		etaEta_gen[i] = (Gamma1MC_Eta_[i]+Gamma2MC_Eta_[i]).Eta();
		phiEta_gen[i] = (Gamma1MC_Eta_[i]+Gamma2MC_Eta_[i]).Phi();
		deltaRG1G2Eta_gen[i] = Gamma1MC_Eta_[i].DeltaR(Gamma2MC_Eta_[i]);
		enG1_Eta_gen[i] = Gamma1MC_Eta_[i].E();
		ptG1_Eta_gen[i] = Gamma1MC_Eta_[i].Pt();
		etaG1_Eta_gen[i] = Gamma1MC_Eta_[i].Eta();
		phiG1_Eta_gen[i] = Gamma1MC_Eta_[i].Phi();
		enG2_Eta_gen[i] = Gamma2MC_Eta_[i].E();
		ptG2_Eta_gen[i] = Gamma2MC_Eta_[i].Pt();
		etaG2_Eta_gen[i] = Gamma2MC_Eta_[i].Eta();
		phiG2_Eta_gen[i] = Gamma2MC_Eta_[i].Phi();

		for(int j=0;j<N_Eta_gen;j++)
		{
			if(j==i) continue;
			double dR1_this = Gamma1MC_Eta_[j].DeltaR(Gamma1MC_Eta_[i]+Gamma2MC_Eta_[i]);
			double dR2_this = Gamma2MC_Eta_[j].DeltaR(Gamma1MC_Eta_[i]+Gamma2MC_Eta_[i]);
			if(dR1_this > 0.3 && dR2_this > 0.3) continue;

			if(dR1_this < 0.1)
			{
				nIsoGamma0p3Eta_gen[i] += 1;
				nIsoGamma0p2Eta_gen[i] += 1;
				nIsoGamma0p1Eta_gen[i] += 1;
			}
			else if(dR1_this < 0.2)
			{
				
				nIsoGamma0p3Eta_gen[i] += 1;
				nIsoGamma0p2Eta_gen[i] += 1;
			}
			else if(dR1_this < 0.3)
			{
				nIsoGamma0p3Eta_gen[i] += 1;
			}

			if(dR2_this < 0.1)
			{
				nIsoGamma0p3Eta_gen[i] += 1;
				nIsoGamma0p2Eta_gen[i] += 1;
				nIsoGamma0p1Eta_gen[i] += 1;
			}
			else if(dR2_this < 0.2)
			{
				
				nIsoGamma0p3Eta_gen[i] += 1;
				nIsoGamma0p2Eta_gen[i] += 1;
			}
			else if(dR2_this < 0.3)
			{
				nIsoGamma0p3Eta_gen[i] += 1;
			}
		}

		for(int j=0;j<N_Pi0_gen;j++)
		{
			double dR1_this = Gamma1MC_Pi0_[j].DeltaR(Gamma1MC_Eta_[i]+Gamma2MC_Eta_[i]);
			double dR2_this = Gamma2MC_Pi0_[j].DeltaR(Gamma1MC_Eta_[i]+Gamma2MC_Eta_[i]);
			if(dR1_this > 0.3 && dR2_this > 0.3) continue;

			if(dR1_this < 0.1)
			{
				nIsoGamma0p3Eta_gen[i] += 1;
				nIsoGamma0p2Eta_gen[i] += 1;
				nIsoGamma0p1Eta_gen[i] += 1;
			}
			else if(dR1_this < 0.2)
			{
				
				nIsoGamma0p3Eta_gen[i] += 1;
				nIsoGamma0p2Eta_gen[i] += 1;
			}
			else if(dR1_this < 0.3)
			{
				nIsoGamma0p3Eta_gen[i] += 1;
			}

			if(dR2_this < 0.1)
			{
				nIsoGamma0p3Eta_gen[i] += 1;
				nIsoGamma0p2Eta_gen[i] += 1;
				nIsoGamma0p1Eta_gen[i] += 1;
			}
			else if(dR2_this < 0.2)
			{
				
				nIsoGamma0p3Eta_gen[i] += 1;
				nIsoGamma0p2Eta_gen[i] += 1;
			}
			else if(dR2_this < 0.3)
			{
				nIsoGamma0p3Eta_gen[i] += 1;
			}
		}



	}

}

void Pi0Tuplizer::MCTruthAssoc(bool isPi0, double deltaR)
{
if(isPi0)
{
	for(unsigned int i=0;i<ebclusters_Pi0_.size();i++)
	{
		ebclusters_Pi0_MC1_index.push_back(-1);
		ebclusters_Pi0_MC2_index.push_back(-1);
	}
	for(unsigned int i=0;i<eeclusters_Pi0_.size();i++)
	{
		eeclusters_Pi0_MC1_index.push_back(-1);
		eeclusters_Pi0_MC2_index.push_back(-1);
	}
	
//	cout<<"DEBUG   entering MCTruthAssoc for pi0, gen size = "<<Gamma1MC_Pi0_.size()<<endl;
	
	for(unsigned int i=0;i<Gamma1MC_Pi0_.size();i++)
	{
//		cout<<"DEBUG  trying to match GEN pi0s with photon energies:  "<<Gamma1MC_Pi0_[i].E()<<"  "<<Gamma2MC_Pi0_[i].E()<<endl;
		double min_deltaR_G1 = 9999.9;
		double min_deltaR_G2 = 9999.9;
		bool G1_found_in_eb = false;
		bool G1_found_in_ee = false;
		bool G2_found_in_eb = false;
		bool G2_found_in_ee = false;
		unsigned int G1_cluster_ind = -1;
		unsigned int G2_cluster_ind = -1;
		double G1_cluster_eta = -9999.99;
		double G2_cluster_eta = -9999.99;

		//match G1
//		cout<<"DEBUG   pi0, trying to match G1...."<<endl;
		for(unsigned int j=0;j<ebclusters_Pi0_.size();j++)
		{
			if(ebclusters_Pi0_MC1_index[j]>=0 || ebclusters_Pi0_MC2_index[j]>=0) continue;
			double deltaR_this = GetDeltaR(Gamma1MC_Pi0_[i].Eta(), ebclusters_Pi0_[j].eta(), Gamma1MC_Pi0_[i].Phi(), ebclusters_Pi0_[j].phi());
			if(deltaR_this < min_deltaR_G1 && deltaR_this < deltaR) 
			{
				G1_found_in_eb	= true;
				G1_cluster_ind = j;
				G1_cluster_eta = ebclusters_Pi0_[j].eta();
				min_deltaR_G1 = deltaR_this;
			}
		}	
		if(!G1_found_in_eb)
		{
			for(unsigned int j=0;j<eeclusters_Pi0_.size();j++)
                	{
				if(eeclusters_Pi0_MC1_index[j]>=0 || eeclusters_Pi0_MC2_index[j]>=0) continue;
                        	double deltaR_this = GetDeltaR(Gamma1MC_Pi0_[i].Eta(), eeclusters_Pi0_[j].eta(), Gamma1MC_Pi0_[i].Phi(), eeclusters_Pi0_[j].phi());
                        	if(deltaR_this < min_deltaR_G1 && deltaR_this < deltaR)
                        	{       
                                	G1_found_in_ee  = true;
                                	G1_cluster_ind = j;
					G1_cluster_eta = eeclusters_Pi0_[j].eta();
					min_deltaR_G1 = deltaR_this;
                        	}
                	} 
		}
		//match G2, if G1 not matched to any reco photon, then don't need to match G2
		if((!G1_found_in_eb) && (!G1_found_in_ee)) continue;
		else
		{
//			cout<<"DEBUG   pi0, G1 matched, now trying to match G2..."<<endl;
			for(unsigned int j=0;j<ebclusters_Pi0_.size() && j!=G1_cluster_ind ;j++)
                        {
				if(ebclusters_Pi0_MC1_index[j]>=0 || ebclusters_Pi0_MC2_index[j]>=0) continue;
                                double deltaR_this = GetDeltaR(Gamma2MC_Pi0_[i].Eta(), ebclusters_Pi0_[j].eta(), Gamma2MC_Pi0_[i].Phi(), ebclusters_Pi0_[j].phi());
                                if(deltaR_this < min_deltaR_G2 && deltaR_this < deltaR)
                                {
                                        G2_found_in_eb  = true;
                                        G2_cluster_ind = j;
					G2_cluster_eta = ebclusters_Pi0_[j].eta();
					min_deltaR_G2 = deltaR_this;
                                }
                        }
			if(!G2_found_in_eb)
			{
				for(unsigned int j=0;j<eeclusters_Pi0_.size() && j!=G1_cluster_ind;j++)
                        	{
					if(eeclusters_Pi0_MC1_index[j]>=0 || eeclusters_Pi0_MC2_index[j]>=0) continue;
                                	double deltaR_this = GetDeltaR(Gamma2MC_Pi0_[i].Eta(), eeclusters_Pi0_[j].eta(), Gamma2MC_Pi0_[i].Phi(), eeclusters_Pi0_[j].phi());
                                	if(deltaR_this < min_deltaR_G2 && deltaR_this < deltaR)
                                	{
                                        	G2_found_in_ee  = true;
                                        	G2_cluster_ind = j;
						G2_cluster_eta = eeclusters_Pi0_[j].eta();
                                        	min_deltaR_G2 = deltaR_this;
                                	}
                        	}

			}	
		}
		if((G1_found_in_eb || G1_found_in_ee) && (G2_found_in_eb || G2_found_in_ee))
		{
			if(G1_found_in_eb) ebclusters_Pi0_MC1_index[G1_cluster_ind] = i;
			if(G1_found_in_ee) eeclusters_Pi0_MC1_index[G1_cluster_ind] = i;
			if(G2_found_in_eb) ebclusters_Pi0_MC2_index[G2_cluster_ind] = i;
			if(G2_found_in_ee) eeclusters_Pi0_MC2_index[G2_cluster_ind] = i;
#ifdef DEBUG
			cout<<"DEBUG  GEN pi0 "<<i<<" matched to reco pi0 (g1, g2) = .... ("<<G1_cluster_ind<<(G1_found_in_eb ? " EB " : " EE ")<<", "<<G2_cluster_ind<<(G2_found_in_eb ? " EB " : " EE ")<<"); \n E1/2 = "<<Gamma1MC_Pi0_[i].E() <<" / "<<Gamma2MC_Pi0_[i].E()<<"  Eta1/2(gen) = "<<Gamma1MC_Pi0_[i].Eta() <<" / "<<Gamma2MC_Pi0_[i].Eta() <<" Eta1/2(reco) = "<<G1_cluster_eta<<" / "<<G2_cluster_eta<< " deltaR1/2 = "<<min_deltaR_G1<<" / "<<min_deltaR_G2<<endl;
#endif
			N_Pi0_match ++;
		}
	}
}
else
{

	for(unsigned int i=0;i<ebclusters_Eta_.size();i++)
	{
		ebclusters_Eta_MC1_index.push_back(-1);
		ebclusters_Eta_MC2_index.push_back(-1);
	}
	for(unsigned int i=0;i<eeclusters_Eta_.size();i++)
	{
		eeclusters_Eta_MC1_index.push_back(-1);
		eeclusters_Eta_MC2_index.push_back(-1);
	}

//	cout<<"DEBUG   entering MCTruthAssoc for eta, gen size = "<<Gamma1MC_Eta_.size()<<endl;

	for(unsigned int i=0;i<Gamma1MC_Eta_.size();i++)
	{
//		cout<<"DEBUG  trying to match GEN etas with photon energies:  "<<Gamma1MC_Eta_[i].E()<<"  "<<Gamma2MC_Eta_[i].E()<<endl;
		double min_deltaR_G1 = 9999.9;
		double min_deltaR_G2 = 9999.9;
		bool G1_found_in_eb = false;
		bool G1_found_in_ee = false;
		bool G2_found_in_eb = false;
		bool G2_found_in_ee = false;
		unsigned int G1_cluster_ind = -1;
		unsigned int G2_cluster_ind = -1;

		//match G1
//		cout<<"DEBUG   eta, trying to match G1...."<<endl;
		for(unsigned int j=0;j<ebclusters_Eta_.size();j++)
		{
			if(ebclusters_Eta_MC1_index[j]>=0 || ebclusters_Eta_MC2_index[j]>=0) continue;
			double deltaR_this = GetDeltaR(Gamma1MC_Eta_[i].Eta(), ebclusters_Eta_[j].eta(), Gamma1MC_Eta_[i].Phi(), ebclusters_Eta_[j].phi());
			if(deltaR_this < min_deltaR_G1 && deltaR_this < deltaR) 
			{
				G1_found_in_eb	= true;
				G1_cluster_ind = j;
				min_deltaR_G1 = deltaR_this;
			}
		}	
		if(!G1_found_in_eb)
		{
			for(unsigned int j=0;j<eeclusters_Eta_.size();j++)
                	{
				if(eeclusters_Eta_MC1_index[j]>=0 || eeclusters_Eta_MC2_index[j]>=0) continue;
                        	double deltaR_this = GetDeltaR(Gamma1MC_Eta_[i].Eta(), eeclusters_Eta_[j].eta(), Gamma1MC_Eta_[i].Phi(), eeclusters_Eta_[j].phi());
                        	if(deltaR_this < min_deltaR_G1 && deltaR_this < deltaR)
                        	{       
                                	G1_found_in_ee  = true;
                                	G1_cluster_ind = j;
					min_deltaR_G1 = deltaR_this;
                        	}
                	} 
		}
		//match G2, if G1 not matched to any reco photon, then don't need to match G2
		if((!G1_found_in_eb) && (!G1_found_in_ee)) continue;
		else
		{
//			cout<<"DEBUG   eta, G1 matched, now trying to match G2..."<<endl;
			for(unsigned int j=0;j<ebclusters_Eta_.size() && j!=G1_cluster_ind ;j++)
                        {
				if(ebclusters_Eta_MC1_index[j]>=0 || ebclusters_Eta_MC2_index[j]>=0) continue;
                                double deltaR_this = GetDeltaR(Gamma2MC_Eta_[i].Eta(), ebclusters_Eta_[j].eta(), Gamma2MC_Eta_[i].Phi(), ebclusters_Eta_[j].phi());
                                if(deltaR_this < min_deltaR_G2 && deltaR_this < deltaR)
                                {
                                        G2_found_in_eb  = true;
                                        G2_cluster_ind = j;
					min_deltaR_G2 = deltaR_this;
                                }
                        }
			if(!G2_found_in_eb)
			{
				for(unsigned int j=0;j<eeclusters_Eta_.size() && j!=G1_cluster_ind;j++)
                        	{
					if(eeclusters_Eta_MC1_index[j]>=0 || eeclusters_Eta_MC2_index[j]>=0) continue;
                                	double deltaR_this = GetDeltaR(Gamma2MC_Eta_[i].Eta(), eeclusters_Eta_[j].eta(), Gamma2MC_Eta_[i].Phi(), eeclusters_Eta_[j].phi());
                                	if(deltaR_this < min_deltaR_G2 && deltaR_this < deltaR)
                                	{
                                        	G2_found_in_ee  = true;
                                        	G2_cluster_ind = j;
                                        	min_deltaR_G2 = deltaR_this;
                                	}
                        	}

			}	
		}
		if((G1_found_in_eb || G1_found_in_ee) && (G2_found_in_eb || G2_found_in_ee))
		{
			if(G1_found_in_eb) ebclusters_Eta_MC1_index[G1_cluster_ind] = i;
			if(G1_found_in_ee) eeclusters_Eta_MC1_index[G1_cluster_ind] = i;
			if(G2_found_in_eb) ebclusters_Eta_MC2_index[G2_cluster_ind] = i;
			if(G2_found_in_ee) eeclusters_Eta_MC2_index[G2_cluster_ind] = i;
			//cout<<"DEBUG   GEN eta matched to reco eta...."<<endl;
			N_Eta_match ++;
		}
	}

}	
}

//load all collection from cfg file
void Pi0Tuplizer::loadEvent(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	if(FillL1SeedFinalDecision_) iEvent.getByToken(uGtAlgToken_,uGtAlg);	
	edm::ESHandle<CaloGeometry> geoHandle;
  	iSetup.get<CaloGeometryRecord>().get(geoHandle);
  	geometry = geoHandle.product();
  	estopology_ = new EcalPreshowerTopology(geoHandle);
	esGeometry_ = (dynamic_cast<const EcalPreshowerGeometry*>( (CaloSubdetectorGeometry*) geometry->getSubdetectorGeometry (DetId::Ecal,EcalPreshower) ));

	if(MCAssoc_) 
	{
		iEvent.getByToken(g4_simTk_Token_, simTracks_h);
		iEvent.getByToken(g4_simVtx_Token_, simVert_h);
		iEvent.getByToken(genParticlesToken_,genParticles);
	}

}


void Pi0Tuplizer::loadEvent_Pi0(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	ebclusters_Pi0_.clear();
	eeclusters_Pi0_.clear();
	ebSeedTime_Pi0_.clear();
	eeSeedTime_Pi0_.clear();
	ebclusters_Pi0_MC1_index.clear();
	ebclusters_Pi0_MC2_index.clear();
	eeclusters_Pi0_MC1_index.clear();
	eeclusters_Pi0_MC2_index.clear();

	ebNxtal_Pi0_.clear();
	eeNxtal_Pi0_.clear();
	ebS4S9_Pi0_.clear();
	eeS4S9_Pi0_.clear();
	ebS2S9_Pi0_.clear();
	eeS2S9_Pi0_.clear();
	ebS1S9_Pi0_.clear();
	eeS1S9_Pi0_.clear();
	
	foundEB = false;
	foundEE = false;
	foundES = false;

	iEvent.getByToken ( EBRecHitCollectionToken_Pi0_, ebRecHit);
	foundEB = ebRecHit->size() > 0;
	iEvent.getByToken ( EERecHitCollectionToken_Pi0_, eeRecHit);
	foundEE = eeRecHit->size() > 0;
	iEvent.getByToken ( ESRecHitCollectionToken_Pi0_, esRecHit);
	foundES = esRecHit->size() > 0;
	N_ebRecHit_Pi0_ += ebRecHit->size();
	N_eeRecHit_Pi0_ += eeRecHit->size();
	N_esRecHit_Pi0_ += esRecHit->size();
	
	N_ebRecHit += ebRecHit->size();
	N_eeRecHit += eeRecHit->size();
	N_esRecHit += esRecHit->size();


//	if (N_esRecHit_Pi0_ > 0) cout<<"DEBUG esRecHit_Pi0.size  "<<N_esRecHit_Pi0_<<endl;
}


void Pi0Tuplizer::loadEvent_Eta(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	ebclusters_Eta_.clear();
	eeclusters_Eta_.clear();
	ebclusters_Eta_MC1_index.clear();
	ebclusters_Eta_MC2_index.clear();
	eeclusters_Eta_MC1_index.clear();
	eeclusters_Eta_MC2_index.clear();
	ebSeedTime_Eta_.clear();
	eeSeedTime_Eta_.clear();
	ebNxtal_Eta_.clear();
	eeNxtal_Eta_.clear();
	ebS4S9_Eta_.clear();
	eeS4S9_Eta_.clear();
	ebS2S9_Eta_.clear();
	eeS2S9_Eta_.clear();
	ebS1S9_Eta_.clear();
	eeS1S9_Eta_.clear();
	

	foundEB = false;
	foundEE = false;
	foundES = false;

	iEvent.getByToken ( EBRecHitCollectionToken_Eta_, ebRecHit);
	foundEB = ebRecHit->size() > 0;
	iEvent.getByToken ( EERecHitCollectionToken_Eta_, eeRecHit);
	foundEE = eeRecHit->size() > 0;
	iEvent.getByToken ( ESRecHitCollectionToken_Eta_, esRecHit);
	foundES = esRecHit->size() > 0;
	N_ebRecHit_Eta_ += ebRecHit->size();
	N_eeRecHit_Eta_ += eeRecHit->size();
	N_esRecHit_Eta_ += esRecHit->size();
	
	N_ebRecHit += ebRecHit->size();
	N_eeRecHit += eeRecHit->size();
	N_esRecHit += esRecHit->size();

//	if(N_esRecHit_Eta_ > 0) cout<<"DEBUG esRecHit_Eta.size  "<<N_esRecHit_Eta_<<endl;
  	

}

//load cut: pi0
void Pi0Tuplizer::loadCut_Pi0(const edm::ParameterSet& iConfig)
{
	EB_Seed_E_Pi0_ 			= iConfig.getUntrackedParameter<double>("EB_Seed_E_Pi0_",0.5);
	EE_Seed_E_Pi0_ 			= iConfig.getUntrackedParameter<double>("EE_Seed_E_Pi0_",0.5);

	pairPtCut_barrel1_Pi0_ 		= iConfig.getUntrackedParameter<double>("pairPtCut_barrel1_Pi0_",2.6);
	pairPtCut_barrel2_Pi0_ 		= iConfig.getUntrackedParameter<double>("pairPtCut_barrel2_Pi0_",2.6);
	pairPtCut_endcap1_Pi0_ 		= iConfig.getUntrackedParameter<double>("pairPtCut_endcap1_Pi0_",3.0);
	pairPtCut_endcap2_Pi0_ 		= iConfig.getUntrackedParameter<double>("pairPtCut_endcap2_Pi0_",1.5);
	
 	gPtCut_barrel1_Pi0_		= iConfig.getUntrackedParameter<double>("gPtCut_barrel1_Pi0_",1.3);
	gPtCut_barrel2_Pi0_ 		= iConfig.getUntrackedParameter<double>("gPtCut_barrel2_Pi0_",1.3);
	gPtCut_endcap1_Pi0_ 		= iConfig.getUntrackedParameter<double>("gPtCut_endcap1_Pi0_",0.95);
	gPtCut_endcap2_Pi0_ 		= iConfig.getUntrackedParameter<double>("gPtCut_endcap2_Pi0_",0.65);
	
 	s4s9Cut_barrel1_Pi0_ 		= iConfig.getUntrackedParameter<double>("s4s9Cut_barrel1_Pi0_",0.83);
	s4s9Cut_barrel2_Pi0_ 		= iConfig.getUntrackedParameter<double>("s4s9Cut_barrel2_Pi0_",0.83);
	s4s9Cut_endcap1_Pi0_ 		= iConfig.getUntrackedParameter<double>("s4s9Cut_endcap1_Pi0_",0.95);
	s4s9Cut_endcap2_Pi0_ 		= iConfig.getUntrackedParameter<double>("s4s9Cut_endcap2_Pi0_",0.95);
	
 	nxtal1Cut_barrel1_Pi0_ 		= iConfig.getUntrackedParameter<double>("nxtal1Cut_barrel1_Pi0_",0.);
	nxtal1Cut_barrel2_Pi0_ 		= iConfig.getUntrackedParameter<double>("nxtal1Cut_barrel2_Pi0_",0.);
	nxtal1Cut_endcap1_Pi0_ 		= iConfig.getUntrackedParameter<double>("nxtal1Cut_endcap1_Pi0_",0.);
	nxtal1Cut_endcap2_Pi0_ 		= iConfig.getUntrackedParameter<double>("nxtal1Cut_endcap2_Pi0_",0.);

	nxtal2Cut_barrel1_Pi0_ 		= iConfig.getUntrackedParameter<double>("nxtal2Cut_barrel1_Pi0_",0.);
	nxtal2Cut_barrel2_Pi0_ 		= iConfig.getUntrackedParameter<double>("nxtal2Cut_barrel2_Pi0_",0.);
	nxtal2Cut_endcap1_Pi0_ 		= iConfig.getUntrackedParameter<double>("nxtal2Cut_endcap1_Pi0_",0.);
	nxtal2Cut_endcap2_Pi0_ 		= iConfig.getUntrackedParameter<double>("nxtal2Cut_endcap2_Pi0_",0.);

	isoGammaBeltdR_Zone_Pi0_	= iConfig.getUntrackedParameter<double>("isoGammaBeltdR_Zone_Pi0_",0.2);
	isoPairBeltdR_Zone_Pi0_		= iConfig.getUntrackedParameter<double>("isoPairBeltdR_Zone_Pi0_",0.2);
	isoGammaBeltdEta_Zone_Pi0_	= iConfig.getUntrackedParameter<double>("isoGammaBeltdEta_Zone_Pi0_",0.05);
	isoPairBeltdEta_Zone_Pi0_	= iConfig.getUntrackedParameter<double>("isoPairBeltdEta_Zone_Pi0_",0.05);
	isoPairCut_	= iConfig.getUntrackedParameter<double>("isoPairCut_",0.05);
	isoGammaCut_	= iConfig.getUntrackedParameter<double>("isoGammaCut_",0.05);

}

void Pi0Tuplizer::loadCut_Eta(const edm::ParameterSet& iConfig)
{
	EB_Seed_E_Eta_ 			= iConfig.getUntrackedParameter<double>("EB_Seed_E_Eta_",0.5);
	EE_Seed_E_Eta_ 			= iConfig.getUntrackedParameter<double>("EE_Seed_E_Eta_",0.5);

	pairPtCut_barrel1_Eta_ 		= iConfig.getUntrackedParameter<double>("pairPtCut_barrel1_Eta_",2.6);
	pairPtCut_barrel2_Eta_ 		= iConfig.getUntrackedParameter<double>("pairPtCut_barrel2_Eta_",2.6);
	pairPtCut_endcap1_Eta_ 		= iConfig.getUntrackedParameter<double>("pairPtCut_endcap1_Eta_",3.0);
	pairPtCut_endcap2_Eta_ 		= iConfig.getUntrackedParameter<double>("pairPtCut_endcap2_Eta_",1.5);
	
 	gPtCut_barrel1_Eta_		= iConfig.getUntrackedParameter<double>("gPtCut_barrel1_Eta_",1.3);
	gPtCut_barrel2_Eta_ 		= iConfig.getUntrackedParameter<double>("gPtCut_barrel2_Eta_",1.3);
	gPtCut_endcap1_Eta_ 		= iConfig.getUntrackedParameter<double>("gPtCut_endcap1_Eta_",0.95);
	gPtCut_endcap2_Eta_ 		= iConfig.getUntrackedParameter<double>("gPtCut_endcap2_Eta_",0.65);
	
 	s4s9Cut_barrel1_Eta_ 		= iConfig.getUntrackedParameter<double>("s4s9Cut_barrel1_Eta_",0.83);
	s4s9Cut_barrel2_Eta_ 		= iConfig.getUntrackedParameter<double>("s4s9Cut_barrel2_Eta_",0.83);
	s4s9Cut_endcap1_Eta_ 		= iConfig.getUntrackedParameter<double>("s4s9Cut_endcap1_Eta_",0.95);
	s4s9Cut_endcap2_Eta_ 		= iConfig.getUntrackedParameter<double>("s4s9Cut_endcap2_Eta_",0.95);
	
 	nxtal1Cut_barrel1_Eta_ 		= iConfig.getUntrackedParameter<double>("nxtal1Cut_barrel1_Eta_",0.);
	nxtal1Cut_barrel2_Eta_ 		= iConfig.getUntrackedParameter<double>("nxtal1Cut_barrel2_Eta_",0.);
	nxtal1Cut_endcap1_Eta_ 		= iConfig.getUntrackedParameter<double>("nxtal1Cut_endcap1_Eta_",0.);
	nxtal1Cut_endcap2_Eta_ 		= iConfig.getUntrackedParameter<double>("nxtal1Cut_endcap2_Eta_",0.);

	nxtal2Cut_barrel1_Eta_ 		= iConfig.getUntrackedParameter<double>("nxtal2Cut_barrel1_Eta_",0.);
	nxtal2Cut_barrel2_Eta_ 		= iConfig.getUntrackedParameter<double>("nxtal2Cut_barrel2_Eta_",0.);
	nxtal2Cut_endcap1_Eta_ 		= iConfig.getUntrackedParameter<double>("nxtal2Cut_endcap1_Eta_",0.);
	nxtal2Cut_endcap2_Eta_ 		= iConfig.getUntrackedParameter<double>("nxtal2Cut_endcap2_Eta_",0.);

	isoGammaBeltdR_Zone_Eta_	= iConfig.getUntrackedParameter<double>("isoGammaBeltdR_Zone_Eta_",0.3);
	isoPairBeltdR_Zone_Eta_		= iConfig.getUntrackedParameter<double>("isoPairBeltdR_Zone_Eta_",0.3);
	isoGammaBeltdEta_Zone_Eta_	= iConfig.getUntrackedParameter<double>("isoGammaBeltdEta_Zone_Eta_",0.1);
	isoPairBeltdEta_Zone_Eta_	= iConfig.getUntrackedParameter<double>("isoPairBeltdEta_Zone_Eta_",0.1);

}

// ------------ method called once each job just before starting event loop  ------------
void Pi0Tuplizer::beginJob()
{

setBranches();

}

//------ Method called once each job just after ending the event loop ------//
void Pi0Tuplizer::endJob()
{

}

// ------------ method called when starting to processes a run  ------------
void Pi0Tuplizer::beginRun(edm::Run const&, edm::EventSetup const& iSetup)
{

}

// ------------ method called when ending the processing of a run  ------------
void Pi0Tuplizer::endRun(edm::Run const&, edm::EventSetup const&)
{

}


// set branch address
void Pi0Tuplizer::setBranches()
{
	//initiate vectors
	L1SeedBitFinalDecision = new std::vector<int>;
	L1SeedBitFinalDecision->clear();
	
//pi0 ntuple
  	Pi0Events->Branch("runNum", &runNum, "runNum/i");
  	Pi0Events->Branch("lumiNum", &lumiNum, "lumiNum/i");
  	Pi0Events->Branch("eventNum", &eventNum, "eventNum/i");
  	Pi0Events->Branch("eventTime", &eventTime, "eventTime/i");
  	Pi0Events->Branch( "allL1SeedFinalDecision", allL1SeedFinalDecision, "allL1SeedFinalDecision[300]/O");
	Pi0Events->Branch("L1SeedBitFinalDecision", "vector<int>", &L1SeedBitFinalDecision);

	Pi0Events->Branch( "N_ebRecHit", &N_ebRecHit, "N_ebRecHit/I");
	Pi0Events->Branch( "N_eeRecHit", &N_eeRecHit, "N_eeRecHit/I");
	Pi0Events->Branch( "N_esRecHit", &N_esRecHit, "N_esRecHit/I");
	Pi0Events->Branch( "N_ebRecHit_Pi0_", &N_ebRecHit_Pi0_, "N_ebRecHit_Pi0_/I");
	Pi0Events->Branch( "N_eeRecHit_Pi0_", &N_eeRecHit_Pi0_, "N_eeRecHit_Pi0_/I");
	Pi0Events->Branch( "N_esRecHit_Pi0_", &N_esRecHit_Pi0_, "N_esRecHit_Pi0_/I");
	Pi0Events->Branch( "N_ebRecHit_Eta_", &N_ebRecHit_Eta_, "N_ebRecHit_Eta_/I");
	Pi0Events->Branch( "N_eeRecHit_Eta_", &N_eeRecHit_Eta_, "N_eeRecHit_Eta_/I");
	Pi0Events->Branch( "N_esRecHit_Eta_", &N_esRecHit_Eta_, "N_esRecHit_Eta_/I");
	Pi0Events->Branch( "N_Pho_rec", &N_Pho_rec, "N_Pho_rec/I");
	Pi0Events->Branch( "N_ebPho_rec", &N_ebPho_rec, "N_ebPho_rec/I");
	Pi0Events->Branch( "N_eePho_rec", &N_eePho_rec, "N_eePho_rec/I");
	
	Pi0Events->Branch( "N_Pho_rec_Pi0_", &N_Pho_rec_Pi0_, "N_Pho_rec_Pi0_/I");
	Pi0Events->Branch( "N_ebPho_rec_Pi0_", &N_ebPho_rec_Pi0_, "N_ebPho_rec_Pi0_/I");
	Pi0Events->Branch( "N_eePho_rec_Pi0_", &N_eePho_rec_Pi0_, "N_eePho_rec_Pi0_/I");
	Pi0Events->Branch( "N_Pho_rec_Eta_", &N_Pho_rec_Eta_, "N_Pho_rec_Eta_/I");
	Pi0Events->Branch( "N_ebPho_rec_Eta_", &N_ebPho_rec_Eta_, "N_ebPho_rec_Eta_/I");
	Pi0Events->Branch( "N_eePho_rec_Eta_", &N_eePho_rec_Eta_, "N_eePho_rec_Eta_/I");

	Pi0Events->Branch( "N_Pair_rec", &N_Pair_rec, "N_Pair_rec/I");
	Pi0Events->Branch( "N_Pi0_rec", &N_Pi0_rec, "N_Pi0_rec/I");
	Pi0Events->Branch( "N_Pi0_gen", &N_Pi0_gen, "N_Pi0_gen/I");
	Pi0Events->Branch( "N_Pi0_genall", &N_Pi0_genall, "N_Pi0_genall/I");
	Pi0Events->Branch( "N_Pi0_match", &N_Pi0_match, "N_Pi0_match/I");
	Pi0Events->Branch( "N_Eta_rec", &N_Eta_rec, "N_Eta_rec/I");
	Pi0Events->Branch( "N_Eta_gen", &N_Eta_gen, "N_Eta_gen/I");
	Pi0Events->Branch( "N_Eta_genall", &N_Eta_genall, "N_Eta_genall/I");
	Pi0Events->Branch( "N_Eta_match", &N_Eta_match, "N_Eta_match/I");
	
	Pi0Events->Branch( "N_ebPair_rec", &N_ebPair_rec, "N_ebPair_rec/I");
	Pi0Events->Branch( "N_ebPi0_rec", &N_ebPi0_rec, "N_ebPi0_rec/I");
	Pi0Events->Branch( "N_ebEta_rec", &N_ebEta_rec, "N_ebEta_rec/I");
	
	Pi0Events->Branch( "N_eePair_rec", &N_eePair_rec, "N_eePair_rec/I");
	Pi0Events->Branch( "N_eePi0_rec", &N_eePi0_rec, "N_eePi0_rec/I");
	Pi0Events->Branch( "N_eeEta_rec", &N_eeEta_rec, "N_eeEta_rec/I");

	Pi0Events->Branch( "nIsoGamma0p3Pi0_genall", nIsoGamma0p3Pi0_genall, "nIsoGamma0p3Pi0_genall[N_Pi0_genall]/I");		
	Pi0Events->Branch( "nIsoGamma0p2Pi0_genall", nIsoGamma0p2Pi0_genall, "nIsoGamma0p2Pi0_genall[N_Pi0_genall]/I");		
	Pi0Events->Branch( "nIsoGamma0p1Pi0_genall", nIsoGamma0p1Pi0_genall, "nIsoGamma0p1Pi0_genall[N_Pi0_genall]/I");		
	Pi0Events->Branch( "ptPi0_genall", ptPi0_genall, "ptPi0_genall[N_Pi0_genall]/F");		
	Pi0Events->Branch( "etaPi0_genall", etaPi0_genall, "etaPi0_genall[N_Pi0_genall]/F");		
	Pi0Events->Branch( "phiPi0_genall", phiPi0_genall, "phiPi0_genall[N_Pi0_genall]/F");		

	Pi0Events->Branch( "nIsoGamma0p3Pi0_gen", nIsoGamma0p3Pi0_gen, "nIsoGamma0p3Pi0_gen[N_Pi0_gen]/I");		
	Pi0Events->Branch( "nIsoGamma0p2Pi0_gen", nIsoGamma0p2Pi0_gen, "nIsoGamma0p2Pi0_gen[N_Pi0_gen]/I");		
	Pi0Events->Branch( "nIsoGamma0p1Pi0_gen", nIsoGamma0p1Pi0_gen, "nIsoGamma0p1Pi0_gen[N_Pi0_gen]/I");		
	Pi0Events->Branch( "ptPi0_gen", ptPi0_gen, "ptPi0_gen[N_Pi0_gen]/F");		
	Pi0Events->Branch( "etaPi0_gen", etaPi0_gen, "etaPi0_gen[N_Pi0_gen]/F");		
	Pi0Events->Branch( "phiPi0_gen", phiPi0_gen, "phiPi0_gen[N_Pi0_gen]/F");		
	Pi0Events->Branch( "deltaRG1G2Pi0_gen", deltaRG1G2Pi0_gen, "deltaRG1G2Pi0_gen[N_Pi0_gen]/F");		
	Pi0Events->Branch( "enG1_Pi0_gen", enG1_Pi0_gen, "enG1_Pi0_gen[N_Pi0_gen]/F");		
	Pi0Events->Branch( "ptG1_Pi0_gen", ptG1_Pi0_gen, "ptG1_Pi0_gen[N_Pi0_gen]/F");		
	Pi0Events->Branch( "etaG1_Pi0_gen", etaG1_Pi0_gen, "etaG1_Pi0_gen[N_Pi0_gen]/F");		
	Pi0Events->Branch( "phiG1_Pi0_gen", phiG1_Pi0_gen, "phiG1_Pi0_gen[N_Pi0_gen]/F");		
	Pi0Events->Branch( "enG2_Pi0_gen", enG2_Pi0_gen, "enG2_Pi0_gen[N_Pi0_gen]/F");		
	Pi0Events->Branch( "ptG2_Pi0_gen", ptG2_Pi0_gen, "ptG2_Pi0_gen[N_Pi0_gen]/F");		
	Pi0Events->Branch( "etaG2_Pi0_gen", etaG2_Pi0_gen, "etaG2_Pi0_gen[N_Pi0_gen]/F");		
	Pi0Events->Branch( "phiG2_Pi0_gen", phiG2_Pi0_gen, "phiG2_Pi0_gen[N_Pi0_gen]/F");		

	Pi0Events->Branch( "nIsoGamma0p3Eta_genall", nIsoGamma0p3Eta_genall, "nIsoGamma0p3Eta_genall[N_Eta_genall]/I");		
	Pi0Events->Branch( "nIsoGamma0p2Eta_genall", nIsoGamma0p2Eta_genall, "nIsoGamma0p2Eta_genall[N_Eta_genall]/I");		
	Pi0Events->Branch( "nIsoGamma0p1Eta_genall", nIsoGamma0p1Eta_genall, "nIsoGamma0p1Eta_genall[N_Eta_genall]/I");		
	Pi0Events->Branch( "ptEta_genall", ptEta_genall, "ptEta_genall[N_Eta_genall]/F");		
	Pi0Events->Branch( "etaEta_genall", etaEta_genall, "etaEta_genall[N_Eta_genall]/F");		
	Pi0Events->Branch( "phiEta_genall", phiEta_genall, "phiEta_genall[N_Eta_genall]/F");		

	Pi0Events->Branch( "nIsoGamma0p3Eta_gen", nIsoGamma0p3Eta_gen, "nIsoGamma0p3Eta_gen[N_Eta_gen]/I");		
	Pi0Events->Branch( "nIsoGamma0p2Eta_gen", nIsoGamma0p2Eta_gen, "nIsoGamma0p2Eta_gen[N_Eta_gen]/I");		
	Pi0Events->Branch( "nIsoGamma0p1Eta_gen", nIsoGamma0p1Eta_gen, "nIsoGamma0p1Eta_gen[N_Eta_gen]/I");		
	Pi0Events->Branch( "ptEta_gen", ptEta_gen, "ptEta_gen[N_Eta_gen]/F");		
	Pi0Events->Branch( "etaEta_gen", etaEta_gen, "etaEta_gen[N_Eta_gen]/F");		
	Pi0Events->Branch( "phiEta_gen", phiEta_gen, "phiEta_gen[N_Eta_gen]/F");		
	Pi0Events->Branch( "deltaRG1G2Eta_gen", deltaRG1G2Eta_gen, "deltaRG1G2Eta_gen[N_Eta_gen]/F");		
	Pi0Events->Branch( "enG1_Eta_gen", enG1_Eta_gen, "enG1_Eta_gen[N_Eta_gen]/F");		
	Pi0Events->Branch( "ptG1_Eta_gen", ptG1_Eta_gen, "ptG1_Eta_gen[N_Eta_gen]/F");		
	Pi0Events->Branch( "etaG1_Eta_gen", etaG1_Eta_gen, "etaG1_Eta_gen[N_Eta_gen]/F");		
	Pi0Events->Branch( "phiG1_Eta_gen", phiG1_Eta_gen, "phiG1_Eta_gen[N_Eta_gen]/F");		
	Pi0Events->Branch( "enG2_Eta_gen", enG2_Eta_gen, "enG2_Eta_gen[N_Eta_gen]/F");		
	Pi0Events->Branch( "ptG2_Eta_gen", ptG2_Eta_gen, "ptG2_Eta_gen[N_Eta_gen]/F");		
	Pi0Events->Branch( "etaG2_Eta_gen", etaG2_Eta_gen, "etaG2_Eta_gen[N_Eta_gen]/F");		
	Pi0Events->Branch( "phiG2_Eta_gen", phiG2_Eta_gen, "phiG2_Eta_gen[N_Eta_gen]/F");		

	Pi0Events->Branch( "fromPi0", fromPi0, "fromPi0[N_Pair_rec]/O");		
	Pi0Events->Branch( "mPi0_rec", mPi0_rec, "mPi0_rec[N_Pair_rec]/F");		
	Pi0Events->Branch( "ptPi0_rec", ptPi0_rec, "ptPi0_rec[N_Pair_rec]/F");		
	Pi0Events->Branch( "isoPi0_rec", isoPi0_rec, "isoPi0_rec[N_Pair_rec]/F");		
	Pi0Events->Branch( "etaPi0_rec", etaPi0_rec, "etaPi0_rec[N_Pair_rec]/F");		
	Pi0Events->Branch( "phiPi0_rec", phiPi0_rec, "phiPi0_rec[N_Pair_rec]/F");		
	Pi0Events->Branch( "isoG1_rec", isoG1_rec, "isoG1_rec[N_Pair_rec]/F");		
	Pi0Events->Branch( "isoG2_rec", isoG2_rec, "isoG2_rec[N_Pair_rec]/F");		
	Pi0Events->Branch( "enG1_rec", enG1_rec, "enG1_rec[N_Pair_rec]/F");		
	Pi0Events->Branch( "enG2_rec", enG2_rec, "enG2_rec[N_Pair_rec]/F");		
	Pi0Events->Branch( "enG1_true", enG1_true, "enG1_true[N_Pair_rec]/F");		
	Pi0Events->Branch( "enG2_true", enG2_true, "enG2_true[N_Pair_rec]/F");		
	Pi0Events->Branch( "etaG1_rec", etaG1_rec, "etaG1_rec[N_Pair_rec]/F");		
	Pi0Events->Branch( "etaG2_rec", etaG2_rec, "etaG2_rec[N_Pair_rec]/F");		
	Pi0Events->Branch( "phiG1_rec", phiG1_rec, "phiG1_rec[N_Pair_rec]/F");		
	Pi0Events->Branch( "phiG2_rec", phiG2_rec, "phiG2_rec[N_Pair_rec]/F");		
	Pi0Events->Branch( "ptG1_rec", ptG1_rec, "ptG1_rec[N_Pair_rec]/F");		
	Pi0Events->Branch( "ptG2_rec", ptG2_rec, "ptG2_rec[N_Pair_rec]/F");		
	Pi0Events->Branch( "iEtaG1_rec", iEtaG1_rec, "iEtaG1_rec[N_Pair_rec]/I");		
	Pi0Events->Branch( "iXG1_rec", iXG1_rec, "iXG1_rec[N_Pair_rec]/I");		
	Pi0Events->Branch( "iEtaG2_rec", iEtaG2_rec, "iEtaG2_rec[N_Pair_rec]/I");		
	Pi0Events->Branch( "iXG2_rec", iXG2_rec, "iXG2_rec[N_Pair_rec]/I");		
	Pi0Events->Branch( "iPhiG1_rec", iPhiG1_rec, "iPhiG1_rec[N_Pair_rec]/I");		
	Pi0Events->Branch( "iYG1_rec", iYG1_rec, "iYG1_rec[N_Pair_rec]/I");		
	Pi0Events->Branch( "iPhiG2_rec", iPhiG2_rec, "iPhiG2_rec[N_Pair_rec]/I");		
	Pi0Events->Branch( "iYG2_rec", iYG2_rec, "iYG2_rec[N_Pair_rec]/I");		
	Pi0Events->Branch( "deltaRG1G2_rec", deltaRG1G2_rec, "deltaRG1G2_rec[N_Pair_rec]/F");		
	Pi0Events->Branch( "nxtalG1_rec", nxtalG1_rec, "nxtalG1_rec[N_Pair_rec]/I");		
	Pi0Events->Branch( "nxtalG2_rec", nxtalG2_rec, "nxtalG2_rec[N_Pair_rec]/I");		
	Pi0Events->Branch( "seedTimeG1_rec", seedTimeG1_rec, "seedTimeG1_rec[N_Pair_rec]/F");		
	Pi0Events->Branch( "seedTimeG2_rec", seedTimeG2_rec, "seedTimeG2_rec[N_Pair_rec]/F");		
	Pi0Events->Branch( "s4s9G1_rec", s4s9G1_rec, "s4s9G1_rec[N_Pair_rec]/F");		
	Pi0Events->Branch( "s4s9G2_rec", s4s9G2_rec, "s4s9G2_rec[N_Pair_rec]/F");		
	Pi0Events->Branch( "s2s9G1_rec", s2s9G1_rec, "s2s9G1_rec[N_Pair_rec]/F");		
	Pi0Events->Branch( "s2s9G2_rec", s2s9G2_rec, "s2s9G2_rec[N_Pair_rec]/F");		
	Pi0Events->Branch( "s1s9G1_rec", s1s9G1_rec, "s1s9G1_rec[N_Pair_rec]/F");		
	Pi0Events->Branch( "s1s9G2_rec", s1s9G2_rec, "s1s9G2_rec[N_Pair_rec]/F");		
//photon ntuple
	PhoEvents->Branch("runNum", &runNum, "runNum/i");
  	PhoEvents->Branch("lumiNum", &lumiNum, "lumiNum/i");
  	PhoEvents->Branch("eventNum", &eventNum, "eventNum/i");
  	PhoEvents->Branch("eventTime", &eventTime, "eventTime/i");
  	PhoEvents->Branch("pho_E", &pho_E, "pho_E/F");
  	PhoEvents->Branch("pho_seedE", &pho_seedE, "pho_seedE/F");
  	PhoEvents->Branch("pho_Eta", &pho_Eta, "pho_Eta/F");
  	PhoEvents->Branch("pho_iEta", &pho_iEta, "pho_iEta/I");
  	PhoEvents->Branch("pho_iX", &pho_iX, "pho_iX/I");
  	PhoEvents->Branch("pho_Phi", &pho_Phi, "pho_Phi/F");
  	PhoEvents->Branch("pho_iPhi", &pho_iPhi, "pho_iPhi/I");
  	PhoEvents->Branch("pho_iY", &pho_iY, "pho_iY/I");
  	PhoEvents->Branch("pho_Pt", &pho_Pt, "pho_Pt/F");
  	PhoEvents->Branch("pho_SeedTime", &pho_SeedTime, "pho_SeedTime/F");
  	PhoEvents->Branch("pho_ClusterTime", &pho_ClusterTime, "pho_ClusterTime/F");
  	PhoEvents->Branch("pho_S4S9", &pho_S4S9, "pho_S4S9/F");
  	PhoEvents->Branch("pho_S2S9", &pho_S2S9, "pho_S2S9/F");
  	PhoEvents->Branch("pho_S1S9", &pho_S1S9, "pho_S1S9/F");
  	PhoEvents->Branch("pho_Nxtal", &pho_Nxtal, "pho_Nxtal/I");
  	PhoEvents->Branch("pho_x", &pho_x, "pho_x/F");
  	PhoEvents->Branch("pho_y", &pho_y, "pho_y/F");
  	PhoEvents->Branch("pho_z", &pho_z, "pho_z/F");
	
}

// clear the content of all variables in the output ntuple
void Pi0Tuplizer::resetBranches()
{
	runNum = -1;
	lumiNum = -1;
	eventNum = -1;
	eventTime = -1;
	for(int i=0;i<NL1SEED;i++)
	{
	allL1SeedFinalDecision[i] = false;
	}	
	L1SeedBitFinalDecision->clear();
	N_Pair_rec = 0;
	N_ebPair_rec = 0;
	N_eePair_rec = 0;
	N_Pi0_rec = 0;
	N_Pi0_gen = 0;
	N_Pi0_genall = 0;
	N_Pi0_match = 0;
	N_ebPi0_rec = 0;
	N_eePi0_rec = 0;
	N_Eta_rec = 0;
	N_Eta_gen = 0;
	N_Eta_genall = 0;
	N_Eta_match = 0;
	N_ebEta_rec = 0;
	N_eeEta_rec = 0;
	N_Pho_rec = 0;
	N_ebPho_rec = 0;
	N_eePho_rec = 0;
	N_Pho_rec_Pi0_ = 0;
	N_ebPho_rec_Pi0_ = 0;
	N_eePho_rec_Pi0_ = 0;
	N_Pho_rec_Eta_ = 0;
	N_ebPho_rec_Eta_ = 0;
	N_eePho_rec_Eta_ = 0;

	N_ebRecHit = 0;
	N_eeRecHit = 0;
	N_esRecHit = 0;
	N_ebRecHit_Pi0_ = 0;
	N_eeRecHit_Pi0_ = 0;
	N_esRecHit_Pi0_ = 0;
	N_ebRecHit_Eta_ = 0;
	N_eeRecHit_Eta_ = 0;
	N_esRecHit_Eta_ = 0;

	for(int i=0;i<NPI0MAX;i++)
	{

		nIsoGamma0p3Pi0_genall[i] = 0;
		nIsoGamma0p2Pi0_genall[i] = 0;
		nIsoGamma0p1Pi0_genall[i] = 0;
		ptPi0_genall[i] = 0;
		etaPi0_genall[i] = 0;
		phiPi0_genall[i] = 0;

		nIsoGamma0p3Pi0_gen[i] = 0;
		nIsoGamma0p2Pi0_gen[i] = 0;
		nIsoGamma0p1Pi0_gen[i] = 0;
		ptPi0_gen[i] = 0;
		etaPi0_gen[i] = 0;
		phiPi0_gen[i] = 0;
		deltaRG1G2Pi0_gen[i] = 0;
		enG1_Pi0_gen[i] = 0;
		ptG1_Pi0_gen[i] = 0;
		etaG1_Pi0_gen[i] = 0;
		phiG1_Pi0_gen[i] = 0;
		enG2_Pi0_gen[i] = 0;
		ptG2_Pi0_gen[i] = 0;
		etaG2_Pi0_gen[i] = 0;
		phiG2_Pi0_gen[i] = 0;
		
		nIsoGamma0p3Eta_genall[i] = 0;
		nIsoGamma0p2Eta_genall[i] = 0;
		nIsoGamma0p1Eta_genall[i] = 0;
		ptEta_genall[i] = 0;
		etaEta_genall[i] = 0;
		phiEta_genall[i] = 0;

		nIsoGamma0p3Eta_gen[i] = 0;
		nIsoGamma0p2Eta_gen[i] = 0;
		nIsoGamma0p1Eta_gen[i] = 0;
		ptEta_gen[i] = 0;
		etaEta_gen[i] = 0;
		phiEta_gen[i] = 0;
		deltaRG1G2Eta_gen[i] = 0;
		enG1_Eta_gen[i] = 0;
		ptG1_Eta_gen[i] = 0;
		etaG1_Eta_gen[i] = 0;
		phiG1_Eta_gen[i] = 0;
		enG2_Eta_gen[i] = 0;
		ptG2_Eta_gen[i] = 0;
		etaG2_Eta_gen[i] = 0;
		phiG2_Eta_gen[i] = 0;


		fromPi0[i] = false;
		mPi0_rec[i] = 0;
		ptPi0_rec[i] = 0;
		isoPi0_rec[i] = 0;
		isoG1_rec[i] = 0;
		isoG2_rec[i] = 0;
		etaPi0_rec[i] = 0;
		phiPi0_rec[i] = 0;
		enG1_rec[i] = 0;
		enG2_rec[i] = 0;
		etaG1_rec[i] = 0;
		etaG2_rec[i] = 0;
		phiG1_rec[i] = 0;
		phiG2_rec[i] = 0;
		ptG1_rec[i] = 0;
		ptG2_rec[i] = 0;
		iEtaG1_rec[i] = -999;
		iXG1_rec[i] = -999;
		iEtaG2_rec[i] = -999;
		iXG2_rec[i] = -999;
		iPhiG1_rec[i] = -999;
		iYG1_rec[i] = -999;
		iPhiG2_rec[i] = -999;
		iYG2_rec[i] = -999;
		deltaRG1G2_rec[i] = 0;
		nxtalG1_rec[i] = 0;
		nxtalG2_rec[i] = 0;
		seedTimeG1_rec[i] = 0;
		seedTimeG2_rec[i] = 0;
		s4s9G1_rec[i] = 0;
		s4s9G2_rec[i] = 0;
		s2s9G1_rec[i] = 0;
		s2s9G2_rec[i] = 0;
		s1s9G1_rec[i] = 0;
		s1s9G2_rec[i] = 0;
	}

	pho_E=0;
	pho_seedE=0;
	pho_Eta=0;
	pho_iEta=-999;
	pho_iX=-999;
	pho_Phi=0;
	pho_iPhi=-999;
	pho_iY=-999;
	pho_Pt=0;
	pho_SeedTime=0;
	pho_ClusterTime=0;
	pho_S4S9=0;
	pho_S2S9=0;
	pho_S1S9=0;
	pho_Nxtal=0;
	pho_x=0;
	pho_y=0;
	pho_z=0;

	GammaIndex_Pi0_.clear();
	Gamma1MC_Pi0_.clear();
	Gamma2MC_Pi0_.clear();

	Gamma1MC_Eta_.clear();
	GammaIndex_Eta_.clear();
	Gamma2MC_Eta_.clear();
	
}

void Pi0Tuplizer::resetPhoBranches()
{
	pho_E=0;
	pho_seedE=0;
        pho_Eta=0;
        pho_iEta=-999;
        pho_iX=-999;
        pho_Phi=0;
        pho_iPhi=-999;
        pho_iY=-999;
        pho_Pt=0;
        pho_SeedTime=0;
        pho_ClusterTime=0;
        pho_S4S9=0;
        pho_S2S9=0;
        pho_S1S9=0;
        pho_Nxtal=0;
        pho_x=0;
        pho_y=0;
        pho_z=0;
}


void Pi0Tuplizer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
      edm::ParameterSetDescription desc;
      desc.setUnknown();
      descriptions.addDefault(desc);
}

int Pi0Tuplizer::diff_neta_s(int neta1, int neta2)
{
	int diff = neta1 - neta2;
	return diff;
}

int Pi0Tuplizer::diff_nphi_s(int nphi1, int nphi2)
{
	int diff;
	if(abs(nphi1-nphi2) < (360-abs(nphi1-nphi2))) diff = nphi1 - nphi2;
	
	else
	{
		diff=360-abs(nphi1-nphi2);
		if(nphi1>nphi2) diff=-diff;
	}	
	return diff;
}


float Pi0Tuplizer::GetDeltaR(float eta1, float eta2, float phi1, float phi2){

  return sqrt( (eta1-eta2)*(eta1-eta2)
        + DeltaPhi(phi1, phi2)*DeltaPhi(phi1, phi2) );

}


float Pi0Tuplizer::DeltaPhi(float phi1, float phi2){
  float diff = fabs(phi2 - phi1);
  while (diff >acos(-1)) diff -= 2*acos(-1);
  while (diff <= -acos(-1)) diff += 2*acos(-1);
  
  return diff;

}

PreshowerTools::PreshowerTools(const CaloGeometry* extGeom, CaloSubdetectorTopology* topology_p,  edm::Handle< ESRecHitCollection > & esHandle) : geom_(extGeom)
{
	estopology_ = topology_p;

    	for (ESRecHitCollection::const_iterator it = esHandle->begin(); it != esHandle->end(); it++) {
	rechits_map.insert(std::make_pair(it->id(), *it));
	}
}

PreshowerCluster PreshowerTools::makeOnePreshowerCluster(int stripwindow, ESDetId *strip)
{
	PreshowerCluster finalcluster;
	esroad_2d.clear();
	used_strips.clear();
	int plane = strip->plane();
	EcalRecHitCollection clusterRecHits;
	RecHitsMap recHits_pos;
	EcalPreshowerNavigator navigator(*strip, estopology_);
   	navigator.setHome(*strip);
	findESRoad(stripwindow,*strip,navigator,plane);

   	if ( plane == 1 ) {
      		ESDetId strip_north = navigator.north();
      		findESRoad(stripwindow,strip_north,navigator,plane);
      		navigator.home();
      		ESDetId strip_south = navigator.south();
      		findESRoad(stripwindow,strip_south,navigator,plane);
      		navigator.home();
   	}
   	if ( plane == 2 ) {
      		ESDetId strip_east = navigator.east();
      		findESRoad(stripwindow,strip_east,navigator,plane);
      		navigator.home();
      		ESDetId strip_west = navigator.west();
      		findESRoad(stripwindow,strip_west,navigator,plane);
      		navigator.home();
   	}

	float E_max = 0.;
   	bool found = false;
  	RecHitsMap::iterator max_it;	
	std::vector<ESDetId>::iterator itID;
//	cout<<"DEBUG preshower reco esroad_2d.size "<<esroad_2d.size()<<endl;
//	cout<<"DEBUG preshower reco rechits_map.size "<<rechits_map.size()<<endl;
   	for (itID = esroad_2d.begin(); itID != esroad_2d.end(); itID++) {
     		RecHitsMap::iterator strip_it = rechits_map.find(*itID);
     		if(!goodStrip(strip_it)) continue;

     		DetId nonblindstripid (itID->rawId());
		float E = strip_it->second.energy();
		
//		cout<<"DEBUG preshower reco E pId  "<< E <<endl;
     		if ( E > E_max) {
        		E_max = E;
        		found = true;
        		max_it = strip_it;
     		}
   	}
	
	if ( !found ) {
		for (itID = esroad_2d.begin(); itID != esroad_2d.end(); itID++) {
               	DetId blindstripid (itID->rawId());
		}

//		cout<<"DEBUG preshower reco return 00001"<<endl;
           	return finalcluster;
	}
	
	clusterRecHits.push_back(max_it->second);
   	recHits_pos.insert(std::make_pair(max_it->first, max_it->second));
   	used_strips.insert(max_it->first);
	ESDetId next, strip_1, strip_2;
   	navigator.setHome(max_it->first);
   	ESDetId startES = max_it->first;
	if (plane == 1) {
		int nadjacents_east = 0;
     		while ( (next=navigator.east()) != ESDetId(0) && next != startES && nadjacents_east < 2 ) {
       		++nadjacents_east;
       		RecHitsMap::iterator strip_it = rechits_map.find(next);

                if(!goodStrip(strip_it)) continue;
		clusterRecHits.push_back(strip_it->second);
		if ( nadjacents_east==1 ) strip_1 = next;
        	used_strips.insert(strip_it->first);
     		}
			
		navigator.home();
     		int nadjacents_west = 0;
     		while ( (next=navigator.west()) != ESDetId(0) && next != startES && nadjacents_west < 2 ) {
        	++nadjacents_west;
        	RecHitsMap::iterator strip_it = rechits_map.find(next);
        	if(!goodStrip(strip_it)) continue;
        	clusterRecHits.push_back(strip_it->second);
        	if ( nadjacents_west==1 ) strip_2 = next;
        	used_strips.insert(strip_it->first);
     		}
	}
	else if (plane == 2) {
		int nadjacents_north = 0;
    	 	while ( (next=navigator.north()) != ESDetId(0) && next != startES && nadjacents_north < 2 ) {
        		++nadjacents_north;  
        		RecHitsMap::iterator strip_it = rechits_map.find(next);
        		if(!goodStrip(strip_it)) continue;      
        		clusterRecHits.push_back(strip_it->second);
        		if ( nadjacents_north==1 ) strip_1 = next;
        		used_strips.insert(strip_it->first);
     		}
		navigator.home();
     		int nadjacents_south = 0;
     		while ( (next=navigator.south()) != ESDetId(0) && next != startES && nadjacents_south < 2 ) {
        		++nadjacents_south;   
        		RecHitsMap::iterator strip_it = rechits_map.find(next);
        		if(!goodStrip(strip_it)) continue;      
        		clusterRecHits.push_back(strip_it->second);
        		if ( nadjacents_south==1 ) strip_2 = next;
        		used_strips.insert(strip_it->first);
     		}
	}
	else {
		std::cout << " Wrong plane number" << plane <<", null cluster will be returned! " << std::endl;
     		return finalcluster;
	}
	RecHitsMap::iterator strip_it1, strip_it2;
   	if ( strip_1 != ESDetId(0)) {
     		strip_it1 = rechits_map.find(strip_1);
     		recHits_pos.insert(std::make_pair(strip_it1->first, strip_it1->second));
   	}
   	if ( strip_2 != ESDetId(0) ) {
     		strip_it2 = rechits_map.find(strip_2);
     		recHits_pos.insert(std::make_pair(strip_it2->first, strip_it2->second));
   	}

	RecHitsMap::iterator cp;
	double energy_pos = 0;
	double x_pos = 0;
	double y_pos = 0;
	double z_pos = 0;
	for (cp = recHits_pos.begin(); cp!=recHits_pos.end(); cp++ ) {
		double E = cp->second.energy();
		energy_pos += E;
		GlobalPoint position = geom_->getPosition(cp->first);
		x_pos += E * position.x();
		y_pos += E * position.y();
		z_pos += E * position.z();
	}

	if(energy_pos>0.) {
     		x_pos /= energy_pos;
     		y_pos /= energy_pos;
     		z_pos /= energy_pos;
  	}

	EcalRecHitCollection::iterator it;
	double Eclust = 0;
	int stripscounter = 0;

	for (it=clusterRecHits.begin(); it != clusterRecHits.end(); it++) {
		Eclust += it->energy();
		stripscounter++;
	}
	std::vector< std::pair<DetId, float> > usedHits;
 	PreshowerCluster output(Eclust,   math::XYZPoint(x_pos,y_pos, z_pos) , usedHits , plane);
	
//	cout<<"DEBUG makePreshower ... Eclust = "<<Eclust<<endl;
	return output;
}


void PreshowerTools::findESRoad(int stripwindow, ESDetId strip, EcalPreshowerNavigator theESNav, int plane) {

   	if ( strip == ESDetId(0) ) return;

    	ESDetId next;
    	theESNav.setHome(strip);
	esroad_2d.push_back(strip);

    	if (plane == 1) {
		int n_east= 0;
      		while ( ((next=theESNav.east()) != ESDetId(0) && next != strip) ) {
         	esroad_2d.push_back(next);
         	++n_east;
         	if (n_east == stripwindow) break;
      		}
		int n_west= 0;
      		theESNav.home();
      		while ( ((next=theESNav.west()) != ESDetId(0) && next != strip )) {
         	esroad_2d.push_back(next);
         	++n_west;
         	if (n_west == stripwindow) break;
      		}
   	}
   	else if (plane == 2) {
		int n_north= 0;
     		while ( ((next=theESNav.north()) != ESDetId(0) && next != strip) ) {
        	esroad_2d.push_back(next);
        	++n_north;
        	if (n_north == stripwindow) break;
     		}

		int n_south= 0;
     		theESNav.home();
     		while ( ((next=theESNav.south()) != ESDetId(0) && next != strip) ) {
        	esroad_2d.push_back(next);
        	++n_south;
        	if (n_south == stripwindow) break;
     		}
   	}

   	theESNav.home();
}
bool PreshowerTools::goodStrip(RecHitsMap::iterator candidate_it){
	if ( (used_strips.find(candidate_it->first) != used_strips.end())  ||        //...if it already belongs to a cluster
        (candidate_it == rechits_map.end() )                    ||        //...if it corresponds to a hit
        (candidate_it->second.energy() <= 0. ) )   // ...if it has a negative or zero energy
     	{
     		return false;
     	}

	return true;

}

//define this as a plug-in
DEFINE_FWK_MODULE(Pi0Tuplizer);

