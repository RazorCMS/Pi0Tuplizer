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
	std::string inputfile1 = argv[1];
	std::string inputfile2 = argv[2];

	if(inputfile1 == "")
	{
	  std::cerr << "[ERROR]: please provide an input file " << std::endl;
          return -1;
	}
	std::cout<<"using input file 1: "<<inputfile1<<std::endl;


	TChain* theChain1 = new TChain();
	theChain1->SetName("ntuples/Pi0Events");
	theChain1->Add( inputfile1.c_str() );
	if ( theChain1 == NULL ) return -1;	

/********************add your analyzer here******************/

	if(inputfile2 == "")
	{
	  std::cerr << "[ERROR]: please provide another input file " << std::endl;
          return -1;
	}
	std::cout<<"using input file 2: "<<inputfile2<<std::endl;

	TChain* theChain2 = new TChain();
	theChain2->SetName("ntuples/Pi0Events");
	theChain2->Add( inputfile2.c_str() );
	if ( theChain2 == NULL ) return -1;	
	
	HLTeff analyzer_HLTeff(theChain1, theChain2);
	analyzer_HLTeff.Analyze();	

/************************************************************/
	return 0;
}
	
