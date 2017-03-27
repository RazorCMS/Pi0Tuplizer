#include "PhoEvents.hh"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>


void PhoEvents::Loop()
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

PhoEvents::PhoEvents(TTree *tree) : fChain(0)
{
        if (tree == 0)
        {
                TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("test.root");
                if (!f || !f->IsOpen())
                {
                        f = new TFile("test.root");
                }
        TDirectory * dir = (TDirectory*)f->Get("test.root:/ntuples");
        dir->GetObject("PhoEvents",tree);

        }
   Init(tree);
}

PhoEvents::~PhoEvents()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t PhoEvents::GetEntry(Long64_t entry)
{
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}

Long64_t PhoEvents::LoadTree(Long64_t entry)
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

void PhoEvents::Init(TTree *tree)
{

   Notify();
}

Bool_t PhoEvents::Notify()
{
return kTRUE;
}

void PhoEvents::Show(Long64_t entry)
{
   if (!fChain) return;
   fChain->Show(entry);
}

Int_t PhoEvents::Cut(Long64_t entry)
{
   return 1;
}

