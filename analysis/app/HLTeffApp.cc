//C++ INCLUDES
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <map>

//ROOT INCLUDES
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TStyle.h>
#include <TGraphAsymmErrors.h>
#include <TEfficiency.h>
#include <TLegend.h>


//LOCAL INCLUDES
#include "Pi0Events.hh"

using namespace std;

int main( int argc, char* argv[])
{
	std::string inputfile_online = argv[1];
	std::string inputfile_offline = argv[2];

	if(inputfile_online == "")
	{
	  std::cerr << "[ERROR]: please provide an input file " << std::endl;
          return -1;
	}
	std::cout<<"using input file 1: "<<inputfile_online<<std::endl;


	TChain* onTree = new TChain();
	onTree->SetName("ntuples/Pi0Events");
	onTree->Add( inputfile_online.c_str() );
	if ( onTree == NULL ) return -1;	

	if(inputfile_offline == "")
	{
	  std::cerr << "[ERROR]: please provide another input file " << std::endl;
          return -1;
	}
	std::cout<<"using input file 2: "<<inputfile_offline<<std::endl;

	TChain* offTree = new TChain();
	offTree->SetName("ntuples/Pi0Events");
	offTree->Add( inputfile_offline.c_str() );
	if ( offTree == NULL ) return -1;	


	
	
	Pi0Events * onT = new Pi0Events(onTree);
        Pi0Events * offT = new Pi0Events(offTree);

	if ( onT->fChain == 0 ) return -1;
        if ( offT->fChain == 0 ) return -1;
        int onT_nentries = onT->fChain->GetEntries();
        int offT_nentries = offT->fChain->GetEntries();
        std::cout << "[INFO]: Total Entries (onT) = " << onT_nentries << "\n";
        std::cout << "[INFO]: Total Entries (offT) = " << offT_nentries << "\n";

        TH1F *h_L1Seed_Online = new TH1F("h_L1Seed_Online","h_L1Seed_Online",NL1SEED,0,NL1SEED);
        TH1F *h_L1Seed_Offline = new TH1F("h_L1Seed_Offline","h_L1Seed_Offline",NL1SEED,0,NL1SEED);


	int N_Online[NL1SEED] = {0};
	int N_Offline[NL1SEED] = {0};

	float frac_Online[NL1SEED] = {0.0};
	float frac_Offline[NL1SEED] = {0.0};

	float Off_vs_On[NL1SEED] = {0.0};
	float x_bit[NL1SEED] = {0.0};
	float exLow[NL1SEED] = {0.0};
	float exHigh[NL1SEED] = {0.0};
	float eyLow[NL1SEED] = {0.0};
	float eyHigh[NL1SEED] = {0.0};
       
	for(int i=0;i<onT_nentries;i++)
        {
                onT->fChain->GetEntry(i);
                for(int ind_bit=0;ind_bit<onT->L1SeedBitFinalDecision->size();ind_bit++)
                {
                        h_L1Seed_Online->Fill(onT->L1SeedBitFinalDecision->at(ind_bit));
			if( onT->L1SeedBitFinalDecision->at(ind_bit) < NL1SEED) N_Online[onT->L1SeedBitFinalDecision->at(ind_bit)] ++;		
                }
        }

        for(int i=0;i<offT_nentries;i++)
        {
                offT->fChain->GetEntry(i);
                for(int ind_bit=0;ind_bit<offT->L1SeedBitFinalDecision->size();ind_bit++)
                {
                        h_L1Seed_Offline->Fill(offT->L1SeedBitFinalDecision->at(ind_bit));
			if( offT->L1SeedBitFinalDecision->at(ind_bit) < NL1SEED) N_Offline[offT->L1SeedBitFinalDecision->at(ind_bit)] ++;		
                }
        }
		
	for(int i=0;i<NL1SEED;i++)
	{
		x_bit[i] = 1.0*i;

		if(N_Online[i] > 100)
		{
			Off_vs_On[i] = (N_Offline[i]*1.0)/(N_Online[i]*1.0);
			eyLow[i] = Off_vs_On[i] - TEfficiency::ClopperPearson((UInt_t)N_Online[i], (UInt_t)N_Offline[i], 0.68269, false);
			eyHigh[i] = - Off_vs_On[i] + TEfficiency::ClopperPearson((UInt_t)N_Online[i], (UInt_t)N_Offline[i], 0.68269, true);

		}
		frac_Online[i] = (N_Online[i]*1.0)/(onT_nentries);
		frac_Offline[i] = (N_Offline[i]*1.0)/(offT_nentries);
		Off_vs_On[i] *= 100.0;
		eyLow[i] *= 100.0;
		eyHigh[i] *= 100.0;
	}

        gStyle->SetOptStat(0);

        TCanvas *c = new TCanvas("c_L1Seed","c_L1Seed",100,100,1000,700);

	c->SetFillColor(0);
        c->SetBorderMode(0);
        c->SetBorderSize(2);
        c->SetFrameBorderMode(0);
        c->SetFrameBorderMode(0);

        h_L1Seed_Online->SetTitle("");
        h_L1Seed_Online->GetXaxis()->SetRangeUser(0,250);
        h_L1Seed_Online->GetXaxis()->SetTitle("L1 Bit");
        h_L1Seed_Online->GetYaxis()->SetTitle("Events");

        h_L1Seed_Online->SetLineColor(4);
        h_L1Seed_Online->SetLineWidth(2);
        h_L1Seed_Offline->SetLineColor(2);
        h_L1Seed_Offline->SetLineWidth(2);

        h_L1Seed_Online->Draw();
        h_L1Seed_Offline->Draw("same");


	TLegend *legend_h = new TLegend(0.50,0.6,0.85,0.88, NULL,"brNDC");
	legend_h->SetTextSize(0.05);
  	legend_h->SetBorderSize(0);
  	legend_h->SetLineColor(1);
  	legend_h->SetLineStyle(1);
  	legend_h->SetLineWidth(1);
  	legend_h->SetFillColor(0);
  	legend_h->SetFillStyle(1001);
	
	legend_h->AddEntry(h_L1Seed_Online, "online selection");
	legend_h->AddEntry(h_L1Seed_Offline, "offline selection");

	legend_h->Draw();

        c->SaveAs("/afs/cern.ch/user/z/zhicaiz/www/sharebox/ECAL/HLT/h_L1Seed_Bit_2016B.pdf");
        c->SaveAs("/afs/cern.ch/user/z/zhicaiz/www/sharebox/ECAL/HLT/h_L1Seed_Bit_2016B.png");
        c->SaveAs("/afs/cern.ch/user/z/zhicaiz/www/sharebox/ECAL/HLT/h_L1Seed_Bit_2016B.C");

	TGraphAsymmErrors *ratio = new TGraphAsymmErrors(NL1SEED, x_bit, Off_vs_On, exLow, exHigh, eyLow, eyHigh); 
	//TH1F* ratio = new TH1F( * h_L1Seed_Offline);
        //ratio->Divide( h_L1Seed_Online );
        ratio->SetTitle("");
        ratio->SetLineColor(4);
        ratio->SetLineWidth(2);
        ratio->GetXaxis()->SetTitle("L1 Bit");
        ratio->GetYaxis()->SetTitle("Offline/Online (%)");
        ratio->GetXaxis()->SetRangeUser(0,250);
        ratio->GetYaxis()->SetRangeUser(1.0,30.0);
        ratio->Draw("AP");

	//ratio->SetMarkerSize(2);
        ratio->SetMarkerStyle(7);

        c->SaveAs("/afs/cern.ch/user/z/zhicaiz/www/sharebox/ECAL/HLT/h_L1Seed_Bit_Ratio_2016B.pdf");
        c->SaveAs("/afs/cern.ch/user/z/zhicaiz/www/sharebox/ECAL/HLT/h_L1Seed_Bit_Ratio_2016B.png");
        c->SaveAs("/afs/cern.ch/user/z/zhicaiz/www/sharebox/ECAL/HLT/h_L1Seed_Bit_Ratio_2016B.C");


	//read the L1 table from txt file

	ifstream L1file("data/L1Table_2016v10.txt");

	std::string s_line;
		
	std::map <std::string, std::string> L1Table;//bit, name

	while(!L1file.eof())
	{
		getline(L1file,s_line);
		if(!s_line.empty())
		{
		std::stringstream ss(s_line);
		
		vector<std::string> tokens;
		std::string buf;
		while (ss >> buf)
       		{
			tokens.push_back(buf);
		}

		L1Table.insert(std::pair<std::string, std::string>(tokens[3], tokens[0]));
		}
	}


	//print out the efficienty table
	FILE* m_outfile = fopen("/afs/cern.ch/user/z/zhicaiz/www/sharebox/ECAL/HLT/HLT_eff_table.txt", "w+");	
	for(int i=0;i<NL1SEED;i++)
	{

		stringstream ss;
		ss << i;
		string str_bit = ss.str();
	
		if(L1Table.find(str_bit.c_str()) == L1Table.end()) continue;

		if(frac_Online[i]>0.05)
		{
			fprintf(m_outfile,"%2d & %s & %6.2f & %6.2f \\\\ \n", i, L1Table[str_bit].c_str(), 100.*frac_Online[i], Off_vs_On[i]);
		}
	}
	

	return 0;
}
	
