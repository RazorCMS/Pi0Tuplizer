
int detID_from_iEtaiPhi(int iEta_or_iX=1, int iPhi_or_iY=1, bool isEB = true, bool isEEMinus = false)
{
	uint32_t detID = 0;
	int Ecal = 3;
	int EcalBarrel=1;
	int EcalEndcap=2;
	int iz = isEEMinus?-1:1;

	if(isEB) 
	{
		detID = ((Ecal&0xF)<<28)|((EcalBarrel&0x7)<<25);
		detID |= ((iEta_or_iX>0)?(0x10000|(iEta_or_iX<<9)):((-iEta_or_iX)<<9))|(iPhi_or_iY&0x1FF);
	}
	else
	{
		detID = ((Ecal&0xF)<<28)|((EcalEndcap&0x7)<<25);
		
		detID |=(iPhi_or_iY&0x7f)|((iEta_or_iX&0x7f)<<7)|((iz>0)?(0x4000):(0));
	}

	return int(detID);	
};


int iEta_or_iX_from_detID(int detID=1, bool isEB = true)
{
	int iEta_or_iX = 0;
	uint32_t id_ = uint32_t(detID);
	if(isEB)
	{
		int zside = (id_&0x10000)?(1):(-1);
		int ietaAbs = (id_>>9)&0x7F;
		iEta_or_iX =  zside*ietaAbs;
	}
	else
	{
		iEta_or_iX = (id_>>7)&0x7F;
	}
	return iEta_or_iX;
	
};

int iPhi_or_iY_from_detID(int detID=1, bool isEB = true)
{
	int iPhi_or_iY = 0;
	uint32_t id_ = uint32_t(detID);
	if(isEB)
	{
		iPhi_or_iY =  id_&0x1FF;
	}
	else
	{
		iPhi_or_iY = id_&0x7F;
	}
	return iPhi_or_iY;
};

 
void  detID_iEtaiPhi()
{
	//EB test
	cout<<"EB test..."<<endl;	
	for(int ieta=-85;ieta<=85 && ieta!=0;ieta++)
	{	
		for(int iphi=1;iphi<=360;iphi++)
		{
			int detID = detID_from_iEtaiPhi(ieta, iphi, true, false);
			cout<<ieta<<"  "<<iphi<<"  "<<detID<<"  "<<iEta_or_iX_from_detID(detID,true)<<"  "<<iPhi_or_iY_from_detID(detID,true)<<endl;
			
		}
	} 
	//EE+ test
	cout<<"EE+ test..."<<endl;	
	for(int iX=1;iX<=100;iX++)
	{	
		for(int iY=1;iY<=100;iY++)
		{
			int detID = detID_from_iEtaiPhi(iX, iY, false, false);
			cout<<iX<<"  "<<iY<<"  "<<detID<<"  "<<iEta_or_iX_from_detID(detID,false)<<"  "<<iPhi_or_iY_from_detID(detID,false)<<endl;
			
		}
	} 

	//EE- test
	cout<<"EE- test..."<<endl;	
	for(int iX=1;iX<=100;iX++)
	{	
		for(int iY=1;iY<=100;iY++)
		{
			int detID = detID_from_iEtaiPhi(iX, iY, false, true);
			cout<<iX<<"  "<<iY<<"  "<<detID<<"  "<<iEta_or_iX_from_detID(detID,false)<<"  "<<iPhi_or_iY_from_detID(detID,false)<<endl;
			
		}
	} 


} 
