//C++ INCLUDES
#include <iostream>
//ROOT INCLUDES
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
//LOCAL INCLUDES
#include "HLTeff.hh"

using namespace std;

int main( int argc, char* argv[])
{
	std::string inputfile = argv[1];

	if(inputfile == "")
	{
	  std::cerr << "[ERROR]: please provide an input file " << std::endl;
          return -1;
	}
	std::cout<<"using input file: "<<inputfile<<std::endl;


	TChain* theChain = new TChain();
	theChain->SetName("ntuples/Pi0Events");
	theChain->Add( inputfile.c_str() );
	if ( theChain == NULL ) return -1;	

/********************add your analyzer here******************/
	
	HLTeff analyzer_HLTeff(theChain);
	analyzer_HLTeff.Analyze();	

/************************************************************/
	return 0;
}
	
