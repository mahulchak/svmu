#include<iostream>
#include<fstream>
#include<cstdlib>
#include "indel.h"

using namespace std;


int main(int argc, char * argv[])
{
        if(argc==1)
        {cerr<<"Usage: "<<argv[0]<<" -d delta_file.out"<<endl;
        exit(EXIT_FAILURE);
        }

	ifstream fin;
	ofstream fout;

	asmMerge merge; 
	string header,ref_name,qu_name,tempname;
	
	int qu_st = 0;
	int qu_end = 0;
	int r_st = 0;
	int r_end = 0;
	
	fin.open(argv[2]);

	while(getline(fin,header))
	{
		if(header[0] =='>')
             	   {
		
            		ref_name = xtractcol(header,' ',1);	
			ref_name = ref_name.substr(1);
			merge.r_name.push_back(ref_name); 
			qu_name = xtractcol(header,' ',2); 
			merge.q_name.push_back(qu_name);
			merge.ref_len[ref_name] = atoi(xtractcol(header,' ',3).c_str()); //length is stored by chromosome name not by alignment pair
			tempname = ref_name.append(qu_name); // tempname is the index for the map. they describe alignment pairs( e.g. BackboneX ctgY). should be unique
			merge.q_len[tempname] = atoi(xtractcol(header,' ',4).c_str());			
			//storeIndex.push_back(tempname);
		}
		
		if(header[0] != '>' && header.size()>10)
		{
			r_st = atoi(xtractcol(header,' ',1).c_str());
			merge.ref_st[tempname].push_back(r_st); // storing the coordinates for each alignment hit
			r_end = atoi(xtractcol(header,' ',2).c_str());
			merge.ref_end[tempname].push_back(r_end);
			qu_st = atoi(xtractcol(header,' ',3).c_str());
			merge.q_st[tempname].push_back(qu_st);
			qu_end = atoi(xtractcol(header,' ',4).c_str());
			merge.q_end[tempname].push_back(qu_end);
		}
	}
	fin.close();

	writeToFile(merge);
	findIndel(merge);
	//fillChromPos(merge);
	//buildCoverage(merge);
	//filterInsCall(merge);
	for(map<string,vector<int> >::iterator it = merge.storeInsStart.begin();it != merge.storeInsStart.end(); it++)
	{
		for(unsigned int i =0; i<merge.storeInsStart[it->first].size();i++)
		{
			cout<<it->first<<"\t"<<merge.storeInsStart[it->first][i]<<"\t"<<merge.storeInsEnd[it->first][i]<<endl;
		}
	}

	return 0;
}
