#ifndef HLTEFF_HH
#define HLTEFF_HH

#include "Pi0Events.hh" 

#include <TROOT.h>
#include <TChain.h>
#include <TTree.h>
#include <TFile.h>
#include "TLorentzVector.h"

#include <map>
#include <string>
#include <vector>
#include <iostream>
using namespace std;

//class HLTeff: public Pi0Events {
class HLTeff{
    public :
	Pi0Events * onT;
	Pi0Events * offT;
        HLTeff(TTree *onTree=0, TTree *offTree=0);
	virtual ~HLTeff();
	virtual void Analyze();

	const float leftMargin   = 0.12;
	const float rightMargin  = 0.05;
	const float topMargin    = 0.07;
	const float bottomMargin = 0.12;

};

#endif

