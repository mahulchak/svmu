

#include<iostream>
#include<fstream>
#include<cstdlib>
#include "qmerge.h"

using namespace std;


int main(int argc, char * argv[])
{
        if(argc==1)
        {cerr<<"Usage: "<<argv[0]<<"-d1 delta_test.out -d2 delta_ref.out -q query.list"<<endl;
        exit(EXIT_FAILURE);
        }

	ifstream delta1,delta2,fin;
	

	asmMerge merge,merge1; 
	
	string header,ref_name,qu_name,tempname,name;
	
	int qu_st = 0;
	int qu_end = 0;
	int r_st = 0;
	int r_end = 0;
	
	delta1.open(argv[2]);

	while(getline(delta1,header))
	{
		if(header[0] =='>')
             	   {
		
            		ref_name = xtractcol(header,' ',1);	
			ref_name = ref_name.substr(4); //blast adds |c| so remove >|c|
			merge.r_name.push_back(ref_name); 
			qu_name = xtractcol(header,' ',2); 
			merge.q_name.push_back(qu_name);
			tempname = ref_name.append(qu_name); // tempname is the index for the map. they describe alignment pairs( e.g. BackboneX ctgY). should be unique. 
			merge.ref_len[tempname] = atoi(xtractcol(header,' ',3).c_str()); 
			merge.q_len[tempname] = atoi(xtractcol(header,' ',4).c_str()); 
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
	delta1.close();
	delta2.open(argv[4]);
	while(getline(delta2,header))
        {
                if(header[0] =='>')
                   {

                        ref_name = xtractcol(header,' ',1);
                        ref_name = ref_name.substr(4);
                        merge1.r_name.push_back(ref_name);
                        qu_name = xtractcol(header,' ',2);
//cout<<qu_name<<endl;
                        merge1.q_name.push_back(qu_name);
                        tempname = ref_name.append(qu_name); // tempname is the index for the map. they describe alignment pairs( e.g. BackboneX ctgY). should be unique. 
                        merge1.ref_len[tempname] = atoi(xtractcol(header,' ',3).c_str());
                        merge1.q_len[tempname] = atoi(xtractcol(header,' ',4).c_str());
                }

                if(header[0] != '>' && header.size()>10)
                {
                        r_st = atoi(xtractcol(header,' ',1).c_str());
                        merge1.ref_st[tempname].push_back(r_st); // storing the coordinates for each alignment hit
                        r_end = atoi(xtractcol(header,' ',2).c_str());
                        merge1.ref_end[tempname].push_back(r_end);
                        qu_st = atoi(xtractcol(header,' ',3).c_str());
                        merge1.q_st[tempname].push_back(qu_st);
                        qu_end = atoi(xtractcol(header,' ',4).c_str());
                        merge1.q_end[tempname].push_back(qu_end);
                }
	}
	delta2.close();
	fin.open(argv[6]);
	
	
	while(getline(fin,header))
	{
//cout<<header<<endl;
		countCopy(header,merge);
	}
	fin.close();
	fin.open(argv[6]);
	while(getline(fin,header))
	{
		countCopy(header,merge1);
	}

	for(map<string,int>::iterator it=merge.storeCount.begin();it!= merge.storeCount.end();it++)
	{
	 if((it->second) != merge1.storeCount[it->first])
		{
		cout<<it->first<<"\t"<<it->second<<"\t"<<merge1.storeCount[it->first]<<endl;
		}
	}
	
	fin.close();


	return 0;
}
