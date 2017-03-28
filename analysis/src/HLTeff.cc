#include "HLTeff.hh"

#include <map>
#include <fstream>
#include <sstream>
#include <string>
#include <assert.h>

//root
#include <TH1F.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TStyle.h>

using namespace std;

//HLTeff::HLTeff(TTree *tree) : Pi0Events(tree)
HLTeff::HLTeff(TTree *onTree, TTree *offTree)
{
	if(onTree==0)	cout<<"empty onTree!"<<endl;
	if(offTree==0)	cout<<"empty offTree!"<<endl;
	
	onT = new Pi0Events(onTree);
	offT = new Pi0Events(offTree);
}

HLTeff::~HLTeff()
{

	delete onT;
	delete offT;

}


void HLTeff::Analyze()
{
	if ( onT->fChain == 0 ) return;
	if ( offT->fChain == 0 ) return;
	int onT_nentries = onT->fChain->GetEntries();
	int offT_nentries = offT->fChain->GetEntries();
  	std::cout << "[INFO]: Total Entries (onT) = " << onT_nentries << "\n";
  	std::cout << "[INFO]: Total Entries (offT) = " << offT_nentries << "\n";

	TH1F *h_L1Seed_Online = new TH1F("h_L1Seed_Online","h_L1Seed_Online",NL1SEED,0,NL1SEED);
	TH1F *h_L1Seed_Offline = new TH1F("h_L1Seed_Offline","h_L1Seed_Offline",NL1SEED,0,NL1SEED);
	
	for(int i=0;i<onT_nentries;i++)
	{
		onT->fChain->GetEntry(i);
		//cout<<" entry "<<i<<" seed bit size  "<<onT->L1SeedBitFinalDecision->size()<<endl;
		for(int ind_bit=0;ind_bit<onT->L1SeedBitFinalDecision->size();ind_bit++)
		{
			h_L1Seed_Online->Fill(onT->L1SeedBitFinalDecision->at(ind_bit));
		}
	}

	for(int i=0;i<offT_nentries;i++)
	{
		offT->fChain->GetEntry(i);
		//cout<<" entry "<<i<<" seed bit size  "<<onT->L1SeedBitFinalDecision->size()<<endl;
		for(int ind_bit=0;ind_bit<offT->L1SeedBitFinalDecision->size();ind_bit++)
		{
			h_L1Seed_Offline->Fill(offT->L1SeedBitFinalDecision->at(ind_bit));
		}
	}

	gStyle->SetOptStat(0);

	TCanvas *c = new TCanvas("c_L1Seed","c_L1Seed",100,100,1000,700);

	c->SetFillColor(0);
	c->SetBorderMode(0);
	c->SetBorderSize(2);
	c->SetLeftMargin( leftMargin );
	c->SetRightMargin( rightMargin );
	c->SetTopMargin( topMargin );
	c->SetBottomMargin( bottomMargin );
	c->SetFrameBorderMode(0);
	c->SetFrameBorderMode(0);

	h_L1Seed_Online->SetTitle("");
	h_L1Seed_Online->GetXaxis()->SetTitle("L1 Bit");
	h_L1Seed_Online->GetYaxis()->SetTitle("Events");
	
	h_L1Seed_Online->SetLineColor(4);
	h_L1Seed_Offline->SetLineColor(2);
	
	h_L1Seed_Online->Draw();
	h_L1Seed_Offline->Draw("same");

	c->SaveAs("/afs/cern.ch/user/z/zhicaiz/www/sharebox/ECAL/HLT/h_L1Seed_Bit.pdf");	
	c->SaveAs("/afs/cern.ch/user/z/zhicaiz/www/sharebox/ECAL/HLT/h_L1Seed_Bit.png");	
	c->SaveAs("/afs/cern.ch/user/z/zhicaiz/www/sharebox/ECAL/HLT/h_L1Seed_Bit.C");

	TH1F* ratio = new TH1F( * h_L1Seed_Offline);
 	ratio->Divide( h_L1Seed_Online );
	ratio->SetTitle("");
	ratio->GetXaxis()->SetTitle("L1 Bit");
	ratio->GetYaxis()->SetTitle("Offline/Online");
	ratio->Draw("E");

	c->SaveAs("/afs/cern.ch/user/z/zhicaiz/www/sharebox/ECAL/HLT/h_L1Seed_Bit_Ratio.pdf");	
	c->SaveAs("/afs/cern.ch/user/z/zhicaiz/www/sharebox/ECAL/HLT/h_L1Seed_Bit_Ratio.png");	
	c->SaveAs("/afs/cern.ch/user/z/zhicaiz/www/sharebox/ECAL/HLT/h_L1Seed_Bit_Ratio.C");



}

