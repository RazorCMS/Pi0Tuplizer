// -*- C++ -*-
/*
 *   Description: reconstruction of pi0/eta and save the ntuple
*/
// Author: Zhicai Zhang
// Created: Fri Mar 24 11:01:09 CET 2017


#include "Pi0Tuplizer.h"

using namespace std;

#define DEBUG


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
	//isPi0_	= iConfig.getUntrackedParameter<bool>("isPi0",false);
	FillL1SeedFinalDecision_	= iConfig.getUntrackedParameter<bool>("FillL1SeedFinalDecision",false);
	FillDiPhotonNtuple_	= iConfig.getUntrackedParameter<bool>("FillDiPhotonNtuple",false);
	FillPhotonNtuple_	= iConfig.getUntrackedParameter<bool>("FillPhotonNtuple",false);

        if(FillL1SeedFinalDecision_)   uGtAlgToken_ = consumes<BXVector<GlobalAlgBlk>>(iConfig.getUntrackedParameter<edm::InputTag>("uGtAlgInputTag"));	
	
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
	if(FillDiPhotonNtuple_) Pi0Events = fs->make<TTree>("Pi0Events", "reconstructed pi0/eta ntuple");	
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

	loadEvent_Pi0(iEvent, iSetup); 
	resetPhoBranches();
	if(foundEB) recoPhoCluster_EB(true);//reconstruct photon clusters in EB
	resetPhoBranches();
	if(foundEE) recoPhoCluster_EE(true);//reconstruct photon clusters in EE
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
	if(FillDiPhotonNtuple_) recoDiPhoEvents_EB(false);//reconstruct pi0/eta events from the photon clusters in EB
	N_ebPair_rec = N_Pair_rec - N_eePair_rec;
	if(FillDiPhotonNtuple_) recoDiPhoEvents_EE(false);//reconstruct pi0/eta events from the photon clusters in EE
	N_eePair_rec = N_Pair_rec - N_ebPair_rec;
	
	N_Eta_rec = N_Pair_rec - N_Pi0_rec;
	N_ebEta_rec = N_ebPair_rec - N_ebPi0_rec;
	N_eeEta_rec = N_eePair_rec - N_eePi0_rec;
//fill ntuple
#ifdef DEBUG
//	cout<<"N_Pho: "<<N_Pho_rec<<"  N_ebRecHit: "<<N_ebRecHit<<"   N_eeRecHit:  "<<N_eeRecHit<<endl;
#endif
	//if(FillDiPhotonNtuple_ && N_Pair_rec > 0) Pi0Events->Fill();
	Pi0Events->Fill();

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
	
		if(isPi0_)
		{	
		eeclusters_Pi0_.push_back( CaloCluster( e3x3, clusPos, CaloID(CaloID::DET_ECAL_ENDCAP), enFracs, CaloCluster::undefined, seed_id ) );
		eeSeedTime_Pi0_.push_back( itseed->time() );	
		eeNxtal_Pi0_.push_back(RecHitsInWindow.size());
		eeS4S9_Pi0_.push_back(s4s9);		
		eeS2S9_Pi0_.push_back(s2s9);		
		eeS1S9_Pi0_.push_back((maxEne)/e3x3);		
		}
		else
		{	
		eeclusters_Eta_.push_back( CaloCluster( e3x3, clusPos, CaloID(CaloID::DET_ECAL_ENDCAP), enFracs, CaloCluster::undefined, seed_id ) );
		eeSeedTime_Eta_.push_back( itseed->time() );	
		eeNxtal_Eta_.push_back(RecHitsInWindow.size());
		eeS4S9_Eta_.push_back(s4s9);		
		eeS2S9_Eta_.push_back(s2s9);		
		eeS1S9_Eta_.push_back((maxEne)/e3x3);		
		}

		//fill photon cluster
		pho_E = e3x3;
		pho_seedE = maxEne;
		pho_Eta = clusPos.eta();
		pho_iX = seed_id.ix();
		pho_Phi = clusPos.phi();
		pho_iY = seed_id.iy();
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
	
	for(std::vector<CaloCluster>::const_iterator g1  = (isPi0_? ebclusters_Pi0_.begin() : ebclusters_Eta_.begin() ); g1 != (isPi0_? ebclusters_Pi0_.end() : ebclusters_Eta_.end()); ++g1, ++i)
  	{
		if(g1->seed().subdetId()!=1) continue;

		int j=i+1;
		for(std::vector<CaloCluster>::const_iterator g2 = g1+1; g2 != (isPi0_? ebclusters_Pi0_.end() : ebclusters_Eta_.end()); ++g2, ++j ) 
		{
			if(g2->seed().subdetId()!=1) continue;
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

			if( pi0P4_nocor.mass()<0.03 && pi0P4.mass() < 0.03 ) continue;
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


			//fill pi0/eta ntuple
			if(N_Pair_rec >= NPI0MAX-1) break; // too many pi0s
			if( FillDiPhotonNtuple_ && pi0P4.mass() > ((isPi0_)?0.03:0.2) && pi0P4.mass() < ((isPi0_)?0.25:1.) )
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
	
	for(std::vector<CaloCluster>::const_iterator g1  = (isPi0_? eeclusters_Pi0_.begin() : eeclusters_Eta_.begin() ); g1 != (isPi0_? eeclusters_Pi0_.end() : eeclusters_Eta_.end()); ++g1, ++i)
  	{
		if(g1->seed().subdetId()!=2) continue;

		int j=i+1;
		for(std::vector<CaloCluster>::const_iterator g2 = g1+1; g2 != (isPi0_? eeclusters_Pi0_.end() : eeclusters_Eta_.end()); ++g2, ++j ) 
		{
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

			if( pi0P4_nocor.mass()<0.03 && pi0P4.mass() < 0.03 ) continue;
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


			//fill pi0/eta ntuple
			if(N_Pair_rec >= NPI0MAX-1) break; // too many pi0s
			if( FillDiPhotonNtuple_ && pi0P4.mass() > ((isPi0_)?0.03:0.2) && pi0P4.mass() < ((isPi0_)?0.25:1.) )
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


//load all collection from cfg file
void Pi0Tuplizer::loadEvent(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	if(FillL1SeedFinalDecision_) iEvent.getByToken(uGtAlgToken_,uGtAlg);	
	edm::ESHandle<CaloGeometry> geoHandle;
  	iSetup.get<CaloGeometryRecord>().get(geoHandle);
  	geometry = geoHandle.product();
  	estopology_ = new EcalPreshowerTopology(geoHandle);
}


void Pi0Tuplizer::loadEvent_Pi0(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	ebclusters_Pi0_.clear();
	eeclusters_Pi0_.clear();
	ebSeedTime_Pi0_.clear();
	eeSeedTime_Pi0_.clear();
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

}


void Pi0Tuplizer::loadEvent_Eta(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	ebclusters_Eta_.clear();
	eeclusters_Eta_.clear();
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
	Pi0Events->Branch( "N_Eta_rec", &N_Eta_rec, "N_Eta_rec/I");
	
	Pi0Events->Branch( "N_ebPair_rec", &N_ebPair_rec, "N_ebPair_rec/I");
	Pi0Events->Branch( "N_ebPi0_rec", &N_ebPi0_rec, "N_ebPi0_rec/I");
	Pi0Events->Branch( "N_ebEta_rec", &N_ebEta_rec, "N_ebEta_rec/I");
	
	Pi0Events->Branch( "N_eePair_rec", &N_eePair_rec, "N_eePair_rec/I");
	Pi0Events->Branch( "N_eePi0_rec", &N_eePi0_rec, "N_eePi0_rec/I");
	Pi0Events->Branch( "N_eeEta_rec", &N_eeEta_rec, "N_eeEta_rec/I");

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
	N_ebPi0_rec = 0;
	N_eePi0_rec = 0;
	N_Eta_rec = 0;
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


//define this as a plug-in
DEFINE_FWK_MODULE(Pi0Tuplizer);

