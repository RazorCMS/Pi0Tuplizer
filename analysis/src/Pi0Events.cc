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

