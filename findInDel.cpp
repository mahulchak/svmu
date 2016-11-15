#include<iostream>
#include<fstream>
#include<cstdlib>
#include "indel.h"

using namespace std;


int main(int argc, char * argv[])
{
        if(argc==1)
        {cerr<<"Usage: "<<argv[0]<<" subcommand -d delta_file.out -m mutation_type(I/D) -p (qins/rins)"<<endl;
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
	
	fin.open(argv[3]);

	while(getline(fin,header))
	{
		if(header[0] =='>')
             	   {
		
            		ref_name = xtractcol(header,' ',1);	
			ref_name = ref_name.substr(1);
			merge.r_name.push_back(ref_name); 
			qu_name = xtractcol(header,' ',2); 
			merge.q_name.push_back(qu_name);
			//merge.ref_len[ref_name] = atoi(xtractcol(header,' ',3).c_str()); //length is stored by chromosome name not by alignment pair
			tempname = ref_name;
			tempname.append(qu_name); // tempname is the index for the map. they describe alignment pairs( e.g. BackboneX ctgY). should be unique
			merge.ref_len[tempname] = atoi(xtractcol(header,' ',3).c_str());
			merge.q_len[tempname] = atoi(xtractcol(header,' ',4).c_str());			
			fillChromPos(ref_name,merge,merge.ref_len[tempname]);
			fillChromPos(qu_name,merge,merge.q_len[tempname]);
		}
		
		if(header[0] != '>' && header.size()>10)
		{
			r_st = atoi(xtractcol(header,' ',1).c_str());
			merge.ref_st[tempname].push_back(r_st); // storing the coordinates for each alignment hit
			r_end = atoi(xtractcol(header,' ',2).c_str());
			merge.ref_end[tempname].push_back(r_end);
			addCoverage(merge,ref_name, r_st,r_end);
			qu_st = atoi(xtractcol(header,' ',3).c_str());
			merge.q_st[tempname].push_back(qu_st);
			qu_end = atoi(xtractcol(header,' ',4).c_str());
			merge.q_end[tempname].push_back(qu_end);
			if(qu_st<qu_end)
			{
				addCoverage(merge,qu_name,qu_st,qu_end);
			}
			if(qu_st>qu_end)
			{
				addCoverage(merge,qu_name,qu_end,qu_st);
			}
		}
	}
	fin.close();

	writeToFile(merge);
	//cout<<"summary file written"<<endl;
	
	if(string(argv[1]) == "indel")
	{
	double c = 0;
	float prop = stof(string(argv[7]),NULL);
	string name;
		 if(*argv[5] == 'D')
	        {
	                findIndel(merge,'I',prop);
     		}
	        if(*argv[5] == 'I')
	        {
       		        findIndel(merge,'D',prop);
                	collapseRange(merge);
		}
		cout<<"Ins Chr"<<"\t"<<"Ins Start"<<"\t"<<"Del End"<<"\t"<<"Del Chr"<<"\t"<<"Del Start"<<"\t"<<"Del End"<<"\t"<<"Mutation Type"<<"\t"<<"Coverage"<<endl;
		for(map<string,vector<int> >::iterator it = merge.storeInsStart.begin();it != merge.storeInsStart.end(); it++)
	        {
	                for(unsigned int i =0; i<merge.storeInsStart[it->first].size();i++)
      		        {
				name = it->first;
				c = cov(merge,name,merge.storeInsStart[it->first][i],merge.storeInsEnd[it->first][i]);
				if(c >0.1)
				{
                	        	cout<<it->first<<"\t"<<merge.storeInsStart[it->first][i]<<"\t"<<merge.storeInsEnd[it->first][i]<<"\t"<<merge.storeDelName[it->first][i]<<"\t"<<merge.storeDelStart[it->first][i]<<"\t"<<merge.storeDelEnd[it->first][i]<<"\t"<<"C"<<"\t"<<c<<endl;
				}
				else if (!(c > 0.1))
				{
					cout<<it->first<<"\t"<<merge.storeInsStart[it->first][i]<<"\t"<<merge.storeInsEnd[it->first][i]<<"\t"<<merge.storeDelName[it->first][i]<<"\t"<<merge.storeDelStart[it->first][i]<<"\t"<<merge.storeDelEnd[it->first][i]<<"\t"<<"S"<<"\t"<<c<<endl;
				}
               		}
		}
        }
	
	if(string(argv[1]) == "duplication")
	{
		findDup(merge);
	}

	//ovlStoreCalculator(merge);
	//fillChromPos(merge);
	//buildCoverage(merge);
//	for(unsigned int i = 0; i< merge.r_name.size();i++)
//	{
//		tempname = merge.r_name[i] + merge.q_name[i];
		 
//		for(unsigned int j = 0;j< merge.newrefSt[tempname].size();j++)
//		for(unsigned int j = 0; j<merge.anmlrefSt[tempname].size();j++)
//		 {
//if(merge.newrefSt[tempname].size()!=0){
//			cout<<tempname<<"\t"<< merge.newrefSt[tempname][j]<<"\t"<<merge.newrefEnd[tempname][j]<<"\t"<<merge.newq_St[tempname][j]<<"\t"<<merge.newq_end[tempname][j]<<endl;
//					}
//			cout<<tempname<<"\t"<<merge.anmlrefSt[tempname][j]<<"\t"<<merge.anmlrefEnd[tempname][j]<<"\t"<<merge.anmlq_St[tempname][j]<<"\t"<<merge.anmlq_end[tempname][j]<<endl;
//		 }
		//if(merge.refChromPos[tempname].size() >0)
		
//		for(unsigned int j=0;j<merge.refChromPos[merge.r_name[i]].size();j++)
//		{
//			cout<<merge.r_name[i]<<"\t"<<j<<"\t"<<merge.refChromPos[merge.r_name[i]][j]<<endl;
//		}
		
//	}

	if(string(argv[1]) == "coverage")
	{
		fin.open("range.txt");
		string str,job_name,chr_name,numeric;
		size_t pos;
		int x = 0,y=0;
		while(getline(fin,str))
		{
	        pos = str.find('\t');
	        chr_name = str.substr(0,pos);
	        job_name = str;
	        job_name[pos] = ':';
	        numeric = str.substr(pos);
        	x = stoi(numeric,&pos);
   		pos = str.find('\t',pos+1);
		job_name[pos] = '-';
		numeric = numeric.substr(pos-1);
		y = stoi(numeric,nullptr);
        //cout<<"mummerplot -fat -large -postscript -prefix "<<job_name<<" -r "<<chr_name<<" -x ["<<x<<":"<<y<<"] a2i.mdelta"<<endl;
//        cout<<"hello"<<endl;
		cout<<chr_name<<"\t"<<x<<"\t"<<y<<"\t"<<cov(merge,chr_name,x,y)<<endl;
	}

	
		
}
return 0;
}	
