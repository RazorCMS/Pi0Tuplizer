#ifndef PHOEVENTS_HH
#define PHOEVENTS_HH

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <vector>

#define NL1SEED 300
#define NPI0MAX 300


class PhoEvents {
public:
	TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   	Int_t           fCurrent;

	uint    runNum;
        uint    lumiNum;
        uint    eventNum;
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
        float pho_x;
        float pho_y;
        float pho_z;

 	PhoEvents(TTree *tree=0);
        virtual ~PhoEvents();
        virtual Int_t    Cut(Long64_t entry);
        virtual Int_t    GetEntry(Long64_t entry);
        virtual Long64_t LoadTree(Long64_t entry);
        virtual void     Init(TTree *tree);
        virtual void     Loop();
        virtual Bool_t   Notify();
        virtual void     Show(Long64_t entry = -1);
 
};

#endif
