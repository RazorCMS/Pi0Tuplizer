//root
#include <TTree.h>
#include <TFile.h>
#include <iostream>

//c++
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>

void timeCalibTag_perFile(std::string input_file_name="input.dat", std::string outfileName="output.root")
{

	TFile *outfile = new TFile(outfileName.c_str(),"RECREATE");
	TTree *outTree = new TTree("timeCalib", "timeCalibs for different IOVs");

	uint start_run;
	uint end_run;
	std:: vector <int> iEta;
	std:: vector <int> iPhi;
	std:: vector <float> IC_time;
	std:: vector <int> detID;
	
	outTree->Branch("start_run", &start_run);//, "start_run/i");
	outTree->Branch("end_run", &end_run);//, "end_run/i");
	outTree->Branch("iEta", &iEta);
	outTree->Branch("iPhi", &iPhi);
	outTree->Branch("IC_time", &IC_time);
	outTree->Branch("detID", &detID);

	if(input_file_name.empty()) return;	

	std::size_t since_pos = input_file_name.find("_since_");
	std::size_t till_pos = input_file_name.find("_till_");

	std::string substr_since = input_file_name.substr(since_pos+9, 6);	
	std::string substr_till = input_file_name.substr(till_pos+8, 6);	

	stringstream ss_since;
	ss_since<<substr_since;
	ss_since>>start_run;

	stringstream ss_till;
	ss_till<<substr_till;
	ss_till>>end_run;

	cout<<"processing IOV from "<<start_run<<"   to  "<<end_run<<endl;

	ifstream inputFile(input_file_name.c_str());
        std::string ped_line;

        while(!inputFile.eof())
        {
                getline(inputFile,ped_line);
                if(!ped_line.empty())
                {
                std::stringstream ss_line(ped_line);

                vector<std::string> tokens;
                std::string buf;
                while (ss_line >> buf)
                {
                        tokens.push_back(buf);
                }
		int iEta_tmp, iPhi_tmp, detID_tmp;
		float IC_time_tmp;

		stringstream ss_0, ss_1, ss_3, ss_4;

		ss_0<<tokens[0];
		ss_0>>iEta_tmp;
		
		ss_1<<tokens[1];
		ss_1>>iPhi_tmp;

		ss_3<<tokens[3];
		ss_3>>IC_time_tmp;
		
		ss_4<<tokens[4];
		ss_4>>detID_tmp;
		
		//cout<<iEta_tmp<<"   "<<iPhi_tmp<<"   "<<IC_time_tmp<<"   "<<rms_G12_tmp<<"   "<<mean_G6_tmp<<"   "<<rms_G6_tmp<<"   "<<mean_G1_tmp<<"   "<<rms_G1_tmp<<"   "<<detID_tmp<<endl;	
		iEta.push_back(iEta_tmp);
		iPhi.push_back(iPhi_tmp);
		IC_time.push_back(IC_time_tmp);
		detID.push_back(detID_tmp);
                }
        }	
	outTree->Fill();
	outTree->Write();
	outfile->Close();

}
