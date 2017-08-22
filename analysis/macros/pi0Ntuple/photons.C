#include <iostream>

void photons()
{

	string file_in_Name = "../../../pi0Ntuple_Eta.root";
	string file_out_Name = "photons.root";
	
	TFile *f_in = new TFile(file_in_Name.c_str(),"READ");
        TTree *tree_in = (TTree*)f_in->Get("ntuples/Pi0Events");


	int N_Pair_rec;
	
	bool fromPi0[100];
	float deltaRG1G2_rec[100];

	float enG1_rec[100];
	float enG2_rec[100];
	
	float enG1_true[100];
	float enG2_true[100];

	float etaG1_rec[100];
	float etaG2_rec[100];
	float phiG1_rec[100];
	float phiG2_rec[100];
	int iEtaG1_rec[100];
	int iEtaG2_rec[100];
	int iPhiG1_rec[100];
	int iPhiG2_rec[100];
	int iXG1_rec[100];
	int iXG2_rec[100];
	int iYG1_rec[100];
	int iYG2_rec[100];
	int nxtalG1_rec[100];
	int nxtalG2_rec[100];
	float s4s9G1_rec[100];
	float s4s9G2_rec[100];
	float s2s9G1_rec[100];
	float s2s9G2_rec[100];
	float s1s9G1_rec[100];
	float s1s9G2_rec[100];
	

	bool STr2_fromPi0;	
	float STr2_DeltaRG1G2;
	float STr2_eta;
	float STr2_phi;
	int STr2_iEtaiX;
	int STr2_iPhiiY;
	int   STr2_Nxtal;
	float STr2_S4S9;	
	float STr2_S2S9;	
	float STr2_S1S9;
	int STr2_SM_dist;
	int STr2_M_dist;
	float STr2_dist;
	float STr2_EOverEOther;
	float STr2_enG_rec;
	float STr2_enG_true;


	tree_in->SetBranchAddress( "N_Pair_rec",&N_Pair_rec);	
	tree_in->SetBranchAddress( "fromPi0", 	fromPi0);	
	tree_in->SetBranchAddress( "deltaRG1G2_rec", 	deltaRG1G2_rec);	
	tree_in->SetBranchAddress( "enG1_rec", 	enG1_rec);	
	tree_in->SetBranchAddress( "enG2_rec", 	enG2_rec);	
	tree_in->SetBranchAddress( "enG1_true", enG1_true);	
	tree_in->SetBranchAddress( "enG2_true", enG2_true);	
	tree_in->SetBranchAddress( "etaG1_rec", 	etaG1_rec);	
	tree_in->SetBranchAddress( "etaG2_rec", 	etaG2_rec);	
	tree_in->SetBranchAddress( "phiG1_rec", 	phiG1_rec);	
	tree_in->SetBranchAddress( "phiG2_rec", 	phiG2_rec);	
	tree_in->SetBranchAddress( "iEtaG1_rec", 	iEtaG1_rec);	
	tree_in->SetBranchAddress( "iEtaG2_rec", 	iEtaG2_rec);	
	tree_in->SetBranchAddress( "iXG1_rec", 	iXG1_rec);	
	tree_in->SetBranchAddress( "iXG2_rec", 	iXG2_rec);	
	tree_in->SetBranchAddress( "iYG1_rec", 	iYG1_rec);	
	tree_in->SetBranchAddress( "iYG2_rec", 	iYG2_rec);	
	tree_in->SetBranchAddress( "iPhiG1_rec", 	iPhiG1_rec);	
	tree_in->SetBranchAddress( "iPhiG2_rec", 	iPhiG2_rec);	
	tree_in->SetBranchAddress( "nxtalG1_rec", 	nxtalG1_rec);	
	tree_in->SetBranchAddress( "nxtalG2_rec", 	nxtalG2_rec);	
	tree_in->SetBranchAddress( "s4s9G1_rec", 	s4s9G1_rec);	
	tree_in->SetBranchAddress( "s4s9G2_rec", 	s4s9G2_rec);	
	tree_in->SetBranchAddress( "s2s9G1_rec", 	s2s9G1_rec);	
	tree_in->SetBranchAddress( "s2s9G2_rec", 	s2s9G2_rec);	
	tree_in->SetBranchAddress( "s1s9G1_rec", 	s1s9G1_rec);	
	tree_in->SetBranchAddress( "s1s9G2_rec", 	s1s9G2_rec);	
	
	
		
	int N_Entries_Pi0 = tree_in->GetEntries();
	
	TFile *f_out = new TFile(file_out_Name.c_str(),"RECREATE");
        TTree *Tree_Optim_gamma = new TTree("Tree_Optim_gamma","Output TTree gamma");
	
	TTree *Tree_Optim_gamma1 = new TTree("Tree_Optim_gamma1","Output TTree gamma1");
	TTree *Tree_Optim_gamma2 = new TTree("Tree_Optim_gamma2","Output TTree gamma2");

	Tree_Optim_gamma->Branch("STr2_fromPi0",	&STr2_fromPi0,	"STr2_fromPi0/O");
	Tree_Optim_gamma->Branch("STr2_enG_rec",	&STr2_enG_rec,	"STr2_enG_rec/F");
	Tree_Optim_gamma->Branch("STr2_enG_true",	&STr2_enG_true,	"STr2_enG_true/F");
	Tree_Optim_gamma->Branch("STr2_DeltaRG1G2",	&STr2_DeltaRG1G2,	"STr2_DeltaRG1G2/F");
	Tree_Optim_gamma->Branch("STr2_eta",	&STr2_eta,	"STr2_eta/F");
	Tree_Optim_gamma->Branch("STr2_phi",	&STr2_phi,	"STr2_phi/F");
	Tree_Optim_gamma->Branch("STr2_Nxtal",	&STr2_Nxtal,	"STr2_Nxtal/I");
	Tree_Optim_gamma->Branch("STr2_iEtaiX",	&STr2_iEtaiX,	"STr2_iEtaiX/I");
	Tree_Optim_gamma->Branch("STr2_iPhiiY",	&STr2_iPhiiY,	"STr2_iPhiiY/I");
	Tree_Optim_gamma->Branch("STr2_S4S9",	&STr2_S4S9,	"STr2_S4S9/F");
	Tree_Optim_gamma->Branch("STr2_S2S9",	&STr2_S2S9,	"STr2_S2S9/F");
	Tree_Optim_gamma->Branch("STr2_S1S9",	&STr2_S1S9,	"STr2_S1S9/F");
	Tree_Optim_gamma->Branch("STr2_SM_dist",	&STr2_SM_dist,	"STr2_SM_dist/I");
	Tree_Optim_gamma->Branch("STr2_M_dist",	&STr2_M_dist,	"STr2_M_dist/I");
	Tree_Optim_gamma->Branch("STr2_dist",	&STr2_dist,	"STr2_dist/F");
	Tree_Optim_gamma->Branch("STr2_EOverEOther",	&STr2_EOverEOther,	"STr2_EOverEOther/F");

	Tree_Optim_gamma1->Branch("STr2_fromPi0",	&STr2_fromPi0,	"STr2_fromPi0/O");
	Tree_Optim_gamma1->Branch("STr2_enG_rec",	&STr2_enG_rec,	"STr2_enG_rec/F");
	Tree_Optim_gamma1->Branch("STr2_enG_true",	&STr2_enG_true,	"STr2_enG_true/F");
	Tree_Optim_gamma1->Branch("STr2_DeltaRG1G2",	&STr2_DeltaRG1G2,	"STr2_DeltaRG1G2/F");
	Tree_Optim_gamma1->Branch("STr2_eta",	&STr2_eta,	"STr2_eta/F");
	Tree_Optim_gamma1->Branch("STr2_phi",	&STr2_phi,	"STr2_phi/F");
	Tree_Optim_gamma1->Branch("STr2_Nxtal",	&STr2_Nxtal,	"STr2_Nxtal/I");
	Tree_Optim_gamma1->Branch("STr2_iEtaiX",	&STr2_iEtaiX,	"STr2_iEtaiX/I");
	Tree_Optim_gamma1->Branch("STr2_iPhiiY",	&STr2_iPhiiY,	"STr2_iPhiiY/I");
	Tree_Optim_gamma1->Branch("STr2_S4S9",	&STr2_S4S9,	"STr2_S4S9/F");
	Tree_Optim_gamma1->Branch("STr2_S2S9",	&STr2_S2S9,	"STr2_S2S9/F");
	Tree_Optim_gamma1->Branch("STr2_S1S9",	&STr2_S1S9,	"STr2_S1S9/F");
	Tree_Optim_gamma1->Branch("STr2_SM_dist",	&STr2_SM_dist,	"STr2_SM_dist/I");
	Tree_Optim_gamma1->Branch("STr2_M_dist",	&STr2_M_dist,	"STr2_M_dist/I");
	Tree_Optim_gamma1->Branch("STr2_dist",	&STr2_dist,	"STr2_dist/F");
	Tree_Optim_gamma1->Branch("STr2_EOverEOther",	&STr2_EOverEOther,	"STr2_EOverEOther/F");

	Tree_Optim_gamma2->Branch("STr2_fromPi0",	&STr2_fromPi0,	"STr2_fromPi0/O");
	Tree_Optim_gamma2->Branch("STr2_enG_rec",	&STr2_enG_rec,	"STr2_enG_rec/F");
	Tree_Optim_gamma2->Branch("STr2_enG_true",	&STr2_enG_true,	"STr2_enG_true/F");
	Tree_Optim_gamma2->Branch("STr2_DeltaRG1G2",	&STr2_DeltaRG1G2,	"STr2_DeltaRG1G2/F");
	Tree_Optim_gamma2->Branch("STr2_eta",	&STr2_eta,	"STr2_eta/F");
	Tree_Optim_gamma2->Branch("STr2_phi",	&STr2_phi,	"STr2_phi/F");
	Tree_Optim_gamma2->Branch("STr2_Nxtal",	&STr2_Nxtal,	"STr2_Nxtal/I");
	Tree_Optim_gamma2->Branch("STr2_iEtaiX",	&STr2_iEtaiX,	"STr2_iEtaiX/I");
	Tree_Optim_gamma2->Branch("STr2_iPhiiY",	&STr2_iPhiiY,	"STr2_iPhiiY/I");
	Tree_Optim_gamma2->Branch("STr2_S4S9",	&STr2_S4S9,	"STr2_S4S9/F");
	Tree_Optim_gamma2->Branch("STr2_S2S9",	&STr2_S2S9,	"STr2_S2S9/F");
	Tree_Optim_gamma2->Branch("STr2_S1S9",	&STr2_S1S9,	"STr2_S1S9/F");
	Tree_Optim_gamma2->Branch("STr2_SM_dist",	&STr2_SM_dist,	"STr2_SM_dist/I");
	Tree_Optim_gamma2->Branch("STr2_M_dist",	&STr2_M_dist,	"STr2_M_dist/I");
	Tree_Optim_gamma2->Branch("STr2_dist",	&STr2_dist,	"STr2_dist/F");
	Tree_Optim_gamma2->Branch("STr2_EOverEOther",	&STr2_EOverEOther,	"STr2_EOverEOther/F");

	for(int i=0;i<N_Entries_Pi0;i++)
        {
		tree_in->GetEntry(i);
		if(i%100 == 0) cout<<"processing "<<i<<" / "<<N_Entries_Pi0<<endl;
		
		for(int j=0;j<N_Pair_rec;j++)
		{
			STr2_fromPi0 = fromPi0[j];	
			STr2_DeltaRG1G2 = deltaRG1G2_rec[j];

			STr2_enG_rec = enG1_rec[j];
			STr2_enG_true = enG1_true[j];
			STr2_eta = etaG1_rec[j];
			STr2_phi = phiG1_rec[j];
			STr2_Nxtal = nxtalG1_rec[j];
			STr2_iEtaiX = (abs(STr2_eta)<1.479) ? (iEtaG1_rec[j]) : (iXG1_rec[j]);
			STr2_iPhiiY = (abs(STr2_eta)<1.479) ? (iPhiG1_rec[j]) : (iYG1_rec[j]);
			STr2_S4S9 = s4s9G1_rec[j];	
			STr2_S2S9 = s2s9G1_rec[j];	
			STr2_S1S9 = s1s9G1_rec[j];	
			STr2_SM_dist = ((STr2_iPhiiY-1)%20<10)*((STr2_iPhiiY-1)%20) + (((STr2_iPhiiY-1)%20)>=10)*(19-(STr2_iPhiiY-1)%20);
                	STr2_M_dist  = (abs(STr2_iEtaiX)<=25)*(((abs(STr2_iEtaiX)-1)%25<12)*((abs(STr2_iEtaiX)-1)%25) + (((abs(STr2_iEtaiX)-1)%25)>=12)*(24-(abs(STr2_iEtaiX)-1)%25))
                            +(abs(STr2_iEtaiX)>25) * (((abs(STr2_iEtaiX)-26)%20<10)*((abs(STr2_iEtaiX)-26)%20) + (((abs(STr2_iEtaiX)-26)%20)>=10)*(19-(abs(STr2_iEtaiX)-26)%20));
			STr2_dist = sqrt(1.0*STr2_SM_dist*STr2_SM_dist + 1.0*STr2_M_dist*STr2_M_dist);
			STr2_EOverEOther = enG1_rec[j]/enG2_rec[j];
			
			Tree_Optim_gamma->Fill();	
			Tree_Optim_gamma1->Fill();	
			
			STr2_enG_rec = enG2_rec[j];
			STr2_enG_true = enG2_true[j];
			STr2_eta = etaG2_rec[j];
			STr2_phi = phiG2_rec[j];
			STr2_Nxtal = nxtalG2_rec[j];
			STr2_iEtaiX = (abs(STr2_eta)<1.479) ? (iEtaG2_rec[j]) : (iXG2_rec[j]);
			STr2_iPhiiY = (abs(STr2_eta)<1.479) ? (iPhiG2_rec[j]) : (iYG2_rec[j]);
			STr2_S4S9 = s4s9G2_rec[j];	
			STr2_S2S9 = s2s9G2_rec[j];	
			STr2_S1S9 = s1s9G2_rec[j];	
			STr2_SM_dist = ((STr2_iPhiiY-1)%20<10)*((STr2_iPhiiY-1)%20) + (((STr2_iPhiiY-1)%20)>=10)*(19-(STr2_iPhiiY-1)%20);
                	STr2_M_dist  = (abs(STr2_iEtaiX)<=25)*(((abs(STr2_iEtaiX)-1)%25<12)*((abs(STr2_iEtaiX)-1)%25) + (((abs(STr2_iEtaiX)-1)%25)>=12)*(24-(abs(STr2_iEtaiX)-1)%25))
                            +(abs(STr2_iEtaiX)>25) * (((abs(STr2_iEtaiX)-26)%20<10)*((abs(STr2_iEtaiX)-26)%20) + (((abs(STr2_iEtaiX)-26)%20)>=10)*(19-(abs(STr2_iEtaiX)-26)%20));
			STr2_dist = sqrt(1.0*STr2_SM_dist*STr2_SM_dist + 1.0*STr2_M_dist*STr2_M_dist);
			STr2_EOverEOther = enG2_rec[j]/enG1_rec[j];
	
			Tree_Optim_gamma->Fill();	
			Tree_Optim_gamma2->Fill();	
		}
		

	}

  	Tree_Optim_gamma->Write();
        Tree_Optim_gamma1->Write();
        Tree_Optim_gamma2->Write();
}
