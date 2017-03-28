#include "Pi0Events.hh"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>


void Pi0Events::Loop()
{

	if (fChain == 0) return;

   	Long64_t nentries = fChain->GetEntriesFast();

   	Long64_t nbytes = 0, nb = 0;
   	for (Long64_t jentry=0; jentry<nentries;jentry++) 
	{
      		Long64_t ientry = LoadTree(jentry);
      		if (ientry < 0) break;
      		nb = fChain->GetEntry(jentry);   nbytes += nb;
	}
}

Pi0Events::Pi0Events(TTree *tree) : fChain(0)
{
        if (tree == 0)
        {
                TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("test.root");
                if (!f || !f->IsOpen())
                {
                        f = new TFile("test.root");
                }
        TDirectory * dir = (TDirectory*)f->Get("test.root:/ntuples");
        dir->GetObject("Pi0Events",tree);

        }
   Init(tree);
}

Pi0Events::~Pi0Events()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t Pi0Events::GetEntry(Long64_t entry)
{
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}

Long64_t Pi0Events::LoadTree(Long64_t entry)
{
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void Pi0Events::Init(TTree *tree)
{
   if (!tree) return;
   fChain = tree;

   L1SeedBitFinalDecision = 0;

   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("runNum", &runNum);
   fChain->SetBranchAddress("lumiNum", &lumiNum);
   fChain->SetBranchAddress("eventNum", &eventNum);
   fChain->SetBranchAddress("allL1SeedFinalDecision", allL1SeedFinalDecision);
   fChain->SetBranchAddress("L1SeedBitFinalDecision", &L1SeedBitFinalDecision);
   fChain->SetBranchAddress("N_Pi0_rec", &N_Pi0_rec);
   fChain->SetBranchAddress("mPi0_rec", mPi0_rec);
   fChain->SetBranchAddress("ptPi0_rec", ptPi0_rec);
   fChain->SetBranchAddress("etaPi0_rec", etaPi0_rec);
   fChain->SetBranchAddress("phiPi0_rec", phiPi0_rec);
   fChain->SetBranchAddress("enG1_rec", enG1_rec);
   fChain->SetBranchAddress("enG2_rec", enG2_rec);
   fChain->SetBranchAddress("etaG1_rec", etaG1_rec);
   fChain->SetBranchAddress("etaG2_rec", etaG2_rec);
   fChain->SetBranchAddress("phiG1_rec", phiG1_rec);
   fChain->SetBranchAddress("phiG2_rec", phiG2_rec);
   fChain->SetBranchAddress("ptG1_rec", ptG1_rec);
   fChain->SetBranchAddress("ptG2_rec", ptG2_rec);
   fChain->SetBranchAddress("iEtaG1_rec", iEtaG1_rec);
   fChain->SetBranchAddress("iEtaG2_rec", iEtaG2_rec);
   fChain->SetBranchAddress("iPhiG1_rec", iPhiG1_rec);
   fChain->SetBranchAddress("iPhiG2_rec", iPhiG2_rec);
   fChain->SetBranchAddress("deltaRG1G2_rec", deltaRG1G2_rec);
   fChain->SetBranchAddress("nxtalG1_rec", nxtalG1_rec);
   fChain->SetBranchAddress("nxtalG2_rec", nxtalG2_rec);
   fChain->SetBranchAddress("seedTimeG1_rec", seedTimeG1_rec);
   fChain->SetBranchAddress("seedTimeG2_rec", seedTimeG2_rec);
   fChain->SetBranchAddress("s4s9G1_rec", s4s9G1_rec);
   fChain->SetBranchAddress("s4s9G2_rec", s4s9G2_rec);
   fChain->SetBranchAddress("s2s9G1_rec", s2s9G1_rec);
   fChain->SetBranchAddress("s2s9G2_rec", s2s9G2_rec);
   fChain->SetBranchAddress("s1s9G1_rec", s1s9G1_rec);
   fChain->SetBranchAddress("s1s9G2_rec", s1s9G2_rec);
   
   Notify();
}

Bool_t Pi0Events::Notify()
{
return kTRUE;
}

void Pi0Events::Show(Long64_t entry)
{
   if (!fChain) return;
   fChain->Show(entry);
}

Int_t Pi0Events::Cut(Long64_t entry)
{
   return 1;
}

