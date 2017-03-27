#ifndef PI0EVENTS_HH
#define PI0EVENTS_HH

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <vector>

#define NL1SEED 300
#define NPI0MAX 300


class Pi0Events {
public:
	TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   	Int_t           fCurrent;
	
	uint    runNum;
        uint    lumiNum;
        uint    eventNum;
        bool    allL1SeedFinalDecision[NL1SEED];
        std::vector<int> * L1SeedBitFinalDecision;

        int     N_Pi0_rec;
        float   mPi0_rec[NPI0MAX];
        float   ptPi0_rec[NPI0MAX];
        float   etaPi0_rec[NPI0MAX];
        float   phiPi0_rec[NPI0MAX];
        float   enG1_rec[NPI0MAX];
        float   enG2_rec[NPI0MAX];
        float   etaG1_rec[NPI0MAX];
        float   etaG2_rec[NPI0MAX];
        float   phiG1_rec[NPI0MAX];
        float   phiG2_rec[NPI0MAX];
        float   ptG1_rec[NPI0MAX];
        float   ptG2_rec[NPI0MAX];
        int     iEtaG1_rec[NPI0MAX];
        int     iEtaG2_rec[NPI0MAX];
        int     iPhiG1_rec[NPI0MAX];
        int     iPhiG2_rec[NPI0MAX];
        float   deltaRG1G2_rec[NPI0MAX];
        int     nxtalG1_rec[NPI0MAX];
        int     nxtalG2_rec[NPI0MAX];
        float   seedTimeG1_rec[NPI0MAX];
        float   seedTimeG2_rec[NPI0MAX];
        float   s4s9G1_rec[NPI0MAX];
        float   s4s9G2_rec[NPI0MAX];
        float   s2s9G1_rec[NPI0MAX];
        float   s2s9G2_rec[NPI0MAX];
        float   s1s9G1_rec[NPI0MAX];
        float   s1s9G2_rec[NPI0MAX];

	Pi0Events(TTree *tree=0);
   	virtual ~Pi0Events();
   	virtual Int_t    Cut(Long64_t entry);
   	virtual Int_t    GetEntry(Long64_t entry);
   	virtual Long64_t LoadTree(Long64_t entry);
   	virtual void     Init(TTree *tree);
   	virtual void     Loop();
   	virtual Bool_t   Notify();
   	virtual void     Show(Long64_t entry = -1);
};

#endif
