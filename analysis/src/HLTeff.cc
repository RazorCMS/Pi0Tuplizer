#include "HLTeff.hh"

#include <map>
#include <fstream>
#include <sstream>
#include <string>
#include <assert.h>

using namespace std;

HLTeff::HLTeff(TTree *tree) : Pi0Events(tree)
{

}

HLTeff::~HLTeff()
{

}


void HLTeff::Analyze()
{
	if ( fChain == 0 ) return;
  	std::cout << "[INFO]: Total Entries = " << fChain->GetEntries() << "\n";

}

