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


//use of global variables associated to the tree can aviod copy of tree each time you call the function

uint start_run_tmp;
uint end_run_tmp;	
vector <float> *IC_time_all;
vector <int> *detID_all;

float getTimeCalibConstant(TTree *tree, vector <uint> & start_run, vector <uint> & end_run, uint run, int detID)
{
	float timeCalib = 0.0;

	int N_entries = tree->GetEntries();
	int i_entry=0;
	for(int i=0;i<start_run.size();i++)
	{
		if(run>= start_run[i] && run<= end_run[i])
		{
			i_entry = i;
			break;
		}
	}

	if(i_entry>= N_entries) return timeCalib;
	tree->GetEntry(i_entry);
	std::vector<int>::iterator p_id;
	p_id = std::find(detID_all->begin(), detID_all->end(), detID);
	if (p_id == detID_all->end()) return timeCalib;
	int idx = std::distance(detID_all->begin(), p_id);

	if(idx<IC_time_all->size()) timeCalib = IC_time_all->at(idx);	

	return timeCalib;
};


void DB_read_example()
{

/////////////////read time calibration constant tree
	std::string timeCalib_fileName = "/eos/cms/store/group/phys_susy/razor/EcalTiming/EcalTimeCalibConstants_Legacy2016_v1/EcalTimeCalibConstants_Legacy2016_v1.root";

	vector <uint> start_run;//start run of all IOV 
	vector <uint> end_run;//end run of all IOV
	 
	TFile f_timeCalib(timeCalib_fileName.c_str());
 	TTree *tree_timeCalib = (TTree*)f_timeCalib.Get("timeCalib");
	
	tree_timeCalib->SetBranchAddress("start_run", &start_run_tmp);
	tree_timeCalib->SetBranchAddress("end_run", &end_run_tmp);
	tree_timeCalib->SetBranchAddress("IC_time", &IC_time_all);
	tree_timeCalib->SetBranchAddress("detID", &detID_all);

	int N_entries_timeCalib = tree_timeCalib->GetEntries();

	for(int i=0;i<N_entries_timeCalib;i++)
	{
		tree_timeCalib->GetEntry(i);
		start_run.push_back(start_run_tmp);
		end_run.push_back(end_run_tmp);
	}
	
	//test 
	uint test_run = 273158;
	cout<<"EB test..."<<endl;
        for(int ieta=-85;ieta<=85 && ieta!=0;ieta++)
        {
                for(int iphi=1;iphi<=360;iphi++)
                {
                        int detID = detID_from_iEtaiPhi(ieta, iphi, true, false);
			cout<<test_run<<"  "<<ieta<<"  "<<iphi<<"  "<<detID;			
			float time_calib = getTimeCalibConstant(tree_timeCalib, start_run,end_run,test_run, detID);
			cout<<"   "<<time_calib<<endl;			
                }
        }

}
