// system include files

#include <TROOT.h>


#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"

#include <vector>
#include <string>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <cstdlib> //as stdlib.h    
#include <cstdio>
#include <cmath>
#include <sstream> //to use ostringstream to convert numbers to string in c++    


using namespace std;
//using namespace RooFit;

 
void  detID_iEtaiPhi_cmssw()
{


	// EB test
	cout<<"EB test..."<<endl;
        for(int ieta=-85;ieta<=85 && ieta!=0;ieta++)
        {
                for(int iphi=1;iphi<=360;iphi++)
                {
			EBDetId *myID = new EBDetId(ieta, iphi);	
                        cout<<ieta<<"  "<<iphi<<"  "<<myID->rawId()<<"  "<<myID->ieta()<<"  "<<myID->iphi()<<endl;

                }
        }


	//EE+ test
	cout<<"EE+ test..."<<endl;      
        for(int iX=1;iX<=100;iX++)
        {       
                for(int iY=1;iY<=100;iY++)
                {
			EEDetId *myID = new EEDetId(iX,iY, 1);
                        cout<<iX<<"  "<<iY<<"  "<<myID->rawId()<<"  "<<myID->ix()<<"  "<<myID->iy()<<endl;
                        
                }
        } 

	//EE- test
	cout<<"EE- test..."<<endl;      
        for(int iX=1;iX<=100;iX++)
        {       
                for(int iY=1;iY<=100;iY++)
                {
			EEDetId *myID = new EEDetId(iX,iY, -1);
                        cout<<iX<<"  "<<iY<<"  "<<myID->rawId()<<"  "<<myID->ix()<<"  "<<myID->iy()<<endl;
                        
                }
        } 


 } 
