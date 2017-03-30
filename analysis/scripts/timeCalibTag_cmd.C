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

void timeCalibTag_cmd()
{
	std::string input_list="list_EcalTimeCalibConstants_Legacy2016_v1.list";

		
	ifstream inputListFile(input_list.c_str());	


	std::string input_file_name;

	while(!inputListFile.eof())
	{
	getline(inputListFile, input_file_name);
	if(input_file_name.empty()) break;	

 	std::size_t since_pos = input_file_name.find("_since_");
        std::size_t till_pos = input_file_name.find("_till_");



	std::string substr_since = input_file_name.substr(since_pos+9, 6);
        std::string substr_till = input_file_name.substr(till_pos+8, 6);

	std::string outfileName = "/eos/cms/store/group/phys_susy/razor/EcalTiming/EcalTimeCalibConstants_Legacy2016_v1/root_files/tree_since_"+substr_since+"_till_"+substr_till+".root";


	cout<<"root -l -b -q timeCalibTag_perFile.C+(\\\""<<input_file_name<<"\\\",\\\""<<outfileName+"\\\")"<<endl;
	
	}


}
