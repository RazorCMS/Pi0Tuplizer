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

void pedestal_cmd()
{
	std::string input_list="list_EcalPedestals_Legacy2016_time_v1.list";


	uint start_time_second;
	uint end_time_second;

		
	ifstream inputListFile(input_list.c_str());	


	std::string input_file_name;

	while(!inputListFile.eof())
	{
	getline(inputListFile, input_file_name);
	if(input_file_name.empty()) break;	

 	std::size_t since_pos = input_file_name.find("_since_");
        std::size_t till_pos = input_file_name.find("_till_");

        std::string substr_since = input_file_name.substr(since_pos+7, 19);
        std::string substr_till = input_file_name.substr(till_pos+6, 19);

        long long int s_time_tmp;
        long long int e_time_tmp;

        stringstream ss_since, ss_since_second;
        ss_since<<substr_since;
        ss_since>>s_time_tmp;

        stringstream ss_till, ss_till_second;
        ss_till<<substr_till;
        ss_till>>e_time_tmp;

        start_time_second = s_time_tmp>>32;
        end_time_second = e_time_tmp>>32;


        ss_since_second<<start_time_second;
        ss_till_second<<end_time_second;

        std::string substr_since_second = ss_since_second.str();
        std::string substr_till_second = ss_till_second.str();

	std::string outfileName = "/eos/cms/store/group/phys_susy/razor/EcalTiming/EcalPedestals_Legacy2016_time_v1/root_files/tree_since_"+substr_since_second+"_till_"+substr_till_second+".root";


	cout<<"root -l -b -q pedestal.C+(\\\""<<input_file_name<<"\\\",\\\""<<outfileName+"\\\")"<<endl;
	
	}


}
