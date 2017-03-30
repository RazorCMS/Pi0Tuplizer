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

void pedestal_perFile(std::string input_file_name="input.dat", std::string outfileName = "output.root")
{

	TFile *outfile = new TFile(outfileName.c_str(),"RECREATE");
	TTree *outTree = new TTree("pedestal", "pedestals for different IOVs");

	uint start_time_second;
	uint end_time_second;
	std:: vector <int> iEta;
	std:: vector <int> iPhi;
	std:: vector <float> mean_G12;
	std:: vector <float> rms_G12;
	std:: vector <float> mean_G6;
	std:: vector <float> rms_G6;
	std:: vector <float> mean_G1;
	std:: vector <float> rms_G1;
	std:: vector <int> detID;
	
	outTree->Branch("start_time_second", &start_time_second);//, "start_time_second/i");
	outTree->Branch("end_time_second", &end_time_second);//, "end_time_second/i");
	outTree->Branch("iEta", &iEta);
	outTree->Branch("iPhi", &iPhi);
	outTree->Branch("mean_G12", &mean_G12);
	outTree->Branch("rms_G12", &rms_G12);
	outTree->Branch("mean_G6", &mean_G6);
	outTree->Branch("rms_G6", &rms_G6);
	outTree->Branch("mean_G1", &mean_G1);
	outTree->Branch("rms_G1", &rms_G1);
	outTree->Branch("detID", &detID);


	if(input_file_name.empty()) return;	

	std::size_t since_pos = input_file_name.find("_since_");
	std::size_t till_pos = input_file_name.find("_till_");

	std::string substr_since = input_file_name.substr(since_pos+7, 19);	
	std::string substr_till = input_file_name.substr(till_pos+6, 19);	

	long long int s_time_tmp;
	long long int e_time_tmp;
	
	stringstream ss_since;
	ss_since<<substr_since;
	ss_since>>s_time_tmp;

	stringstream ss_till;
	ss_till<<substr_till;
	ss_till>>e_time_tmp;
	
	start_time_second = s_time_tmp>>32;
	end_time_second = e_time_tmp>>32;

	cout<<"processing IOV from "<<s_time_tmp<<" (Unix "<< start_time_second<<" ) to  "<<e_time_tmp<<" (Unix "<<end_time_second<<"  )"<<endl;


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
		float mean_G12_tmp, rms_G12_tmp, mean_G6_tmp, rms_G6_tmp, mean_G1_tmp, rms_G1_tmp;

		stringstream ss_0, ss_1, ss_3, ss_4, ss_5, ss_6, ss_7, ss_8, ss_9;

		ss_0<<tokens[0];
		ss_0>>iEta_tmp;
		
		ss_1<<tokens[1];
		ss_1>>iPhi_tmp;

		ss_3<<tokens[3];
		ss_3>>mean_G12_tmp;
		
		ss_4<<tokens[4];
		ss_4>>rms_G12_tmp;
	
		ss_5<<tokens[5];
		ss_5>>mean_G6_tmp;
	
		ss_6<<tokens[6];
		ss_6>>rms_G6_tmp;
		
		ss_7<<tokens[7];
		ss_7>>mean_G1_tmp;
	
		ss_8<<tokens[8];
		ss_8>>rms_G1_tmp;
	
		ss_9<<tokens[9];
		ss_9>>detID_tmp;
		
		//cout<<iEta_tmp<<"   "<<iPhi_tmp<<"   "<<mean_G12_tmp<<"   "<<rms_G12_tmp<<"   "<<mean_G6_tmp<<"   "<<rms_G6_tmp<<"   "<<mean_G1_tmp<<"   "<<rms_G1_tmp<<"   "<<detID_tmp<<endl;	
		iEta.push_back(iEta_tmp);
		iPhi.push_back(iPhi_tmp);
		mean_G12.push_back(mean_G12_tmp);
		rms_G12.push_back(rms_G12_tmp);
		mean_G6.push_back(mean_G6_tmp);
		rms_G6.push_back(rms_G6_tmp);
		mean_G1.push_back(mean_G1_tmp);
		rms_G1.push_back(rms_G1_tmp);
		detID.push_back(detID_tmp);
                }
        }	
	outTree->Fill();
	
	outTree->Write();

}
