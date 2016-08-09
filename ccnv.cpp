#include<iostream>
#include<fstream>
#include<cstdlib>
#include "val.h"

using namespace std;


int main(int argc, char * argv[])
{
        if(argc==1)
        {
		cerr<<"Usage: "<<argv[0]<<" -d1 delta_candidate.out -d2 delta_ref.out -q query_name.list -c cutoff_for_repeats -qco cutoff_query_merging -rco cutoff_ref_merging -D master_delta.out"<<endl;
        	exit(EXIT_FAILURE);
        }

	ifstream Delta,delta1,delta2,fin;
	ofstream fout,chromMap,orphan;	
	
	asmMerge master,merge,merge1; 
	
	string header,ref_name,qu_name,tempname,name;
	
	int qu_st = 0;
	int qu_end = 0;
	int r_st = 0;
	int r_end = 0;
	int copyCutoff = 0;
	double rCO = 0;
	int qCO = 0;
	copyCutoff = strtod(argv[8],NULL);
	rCO = strtod(argv[12], NULL);
	qCO = stoi(argv[10],NULL);
//cout<<rCO <<"\t"<< qCO<<endl;
//cout<<string(argv[8])<<endl;
	if(string(argv[14]) != "")
	{
		Delta.open(argv[14]);
		while(getline(Delta,header))
       		{
        	        if(header[0] =='>')
                	{

                        	ref_name = xtractcol(header,' ',1);
	                        ref_name = ref_name.substr(1); //blast adds |c| so remove >|c|
        	                master.r_name.push_back(ref_name);
                	        qu_name = xtractcol(header,' ',2);
                	        master.q_name.push_back(qu_name);
                	        tempname = ref_name.append(qu_name); // tempname is the index for the map. they describe alignment pairs( e.g. BackboneX ctgY). should be unique.
                       		master.ref_len[tempname] = atoi(xtractcol(header,' ',3).c_str());
                	        master.q_len[tempname] = atoi(xtractcol(header,' ',4).c_str());
               		}

       		        if(header[0] != '>' && header.size()>10)
               		{
                       		r_st = atoi(xtractcol(header,' ',1).c_str());
                     		master.ref_st[tempname].push_back(r_st); // storing the coordinates for each alignment hit
               		      	r_end = atoi(xtractcol(header,' ',2).c_str());
                  		master.ref_end[tempname].push_back(r_end);
                     		qu_st = atoi(xtractcol(header,' ',3).c_str());
                  		master.q_st[tempname].push_back(qu_st);
                  		qu_end = atoi(xtractcol(header,' ',4).c_str());
                     		master.q_end[tempname].push_back(qu_end);
              		 }
       		}
	

		Delta.close();
	
		ovlStoreCalculator(master);//calculate the alignment length for each alignment
		findChromPartner(master);	
		chromMap.open("ctgmap.txt");
		for(map<string,string>::iterator it=master.storeHomolog.begin();it!= master.storeHomolog.end();it++)
		{
			chromMap<< it->first << "\t"<< it->second <<endl;
		}
		chromMap.close();
	}
	delta1.open(argv[2]);//open the second delta file

	while(getline(delta1,header))
	{
		if(header[0] =='>')
             	   {
		
            		ref_name = xtractcol(header,' ',1);	
			ref_name = ref_name.substr(5); //blast adds |c| so remove >|c|
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
                        ref_name = ref_name.substr(5);
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
collapseRange(merge,rCO,qCO);
collapseRange(merge1,rCO,qCO);	
//masterQlist(merge);
//masterQlist(merge1);
//reducList(merge);	
//reducList(merge1);
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
	fout.open("cnv_report.tsv");
	orphan.open("orphan.seq");
	for(map<string,int>::iterator it=merge.storeCount.begin();it!= merge.storeCount.end();it++)
	{
		 if(((it->second) > merge1.storeCount[it->first]) && (it->second < copyCutoff)) //if former has more copies.copyCutoff is for filtering TEs and repeats
		{
		//cout<<it->first<<"\t"<<it->second<<"\t"<<merge1.storeCount[it->first]<<endl;
			for(unsigned i=0;i<merge.storeName[it->first].size();i++) //access the query names 
			{
				tempname = merge.storeName[it->first][i]; //create the index
				for(unsigned j =0;j<merge.q_st[tempname].size();j++)
				{
					if((merge.q_end[tempname][j] != 0) && (merge.q_st[tempname][j]<merge.q_end[tempname][j]) && (merge1.storeCount[it->first] != 0))
					{
					
						fout<<tempname<<"\t"<<merge.q_st[tempname][j]<<"\t"<<merge.q_end[tempname][j]<<"\t"<<it->second<<"\t"<<merge1.storeCount[it->first]<<"\t"<<"+"<<endl;
					}
					if((merge.q_end[tempname][j] != 0) && (merge.q_st[tempname][j]>merge.q_end[tempname][j]) && (merge1.storeCount[it->first] != 0) )
					{
						fout<<tempname<<"\t"<<merge.q_end[tempname][j]<<"\t"<<merge.q_st[tempname][j]<<"\t"<<it->second<<"\t"<<merge1.storeCount[it->first]<<"\t"<<"-"<<endl;
					}
				}
			}
		}
		if((it->second ==0)||(merge1.storeCount[it->first]==0))
		{
			orphan<<it->first<<"\t"<<it->second<<"\t"<<merge1.storeCount[it->first]<<endl;
		}
	
	}
	
	
	fin.close();
	fout.close();
	orphan.close();
	return 0;
}
