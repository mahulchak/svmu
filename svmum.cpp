#include<iostream>
#include "sv.h"

int main(int argc, char * argv[])
{

	if(argc<3)
	{
	cerr<<"Usage: "<<argv[0]<<" foo.mgaps mumTobed.Chr.bed chrom"<<endl;
	exit(EXIT_FAILURE);
	}
	ifstream fin;
	ofstream fout;

	fin.open(argv[1]);
	fout.open(argv[2]);
	
	string str,name,temp,outFileName;
	int ref_st1,ref_st2,ref_end,aln_len,q_st1,q_st2,q_end,ctMum,totAln;
	
	mgapC cluster;
	
	map<int,vector<int>> Mref; // the first element would be the start coordinate of a master cluster
	map<int,vector<int>> Qref; //the first element would be the start and the 2nd elements would be master query
	map<int,string> ref_name;
	map<int,string> q_name;

	while(getline(fin,str))
	{
		if(str[0] == '>')
		{
			name=str.substr(1); //remove > from the string
		}
		ctMum =0; //this is the count of mums within a cluster.reset it here.

		while(str[0] != '#' && str[0] != '>' && !fin.eof())
		{
//cout<<q_st2<<endl;
			if(ctMum ==0)  //the first match.the only match when a cluster has a single match
			{
				temp = xtractcol(str,'\t',1);
				ref_st1 = stoi(temp,nullptr);
				temp = xtractcol(str,'\t',3);
				aln_len = stoi(temp,nullptr);
				temp = xtractcol(str,'\t',2);
				q_st1 = stoi(temp,nullptr);		
				totAln = aln_len;
			}
			else
			{
				temp = xtractcol(str,'\t',1);
                       		ref_st2 = stoi(temp,nullptr);
				temp = xtractcol(str,'\t',3);
				aln_len = stoi(temp,nullptr);
				temp = xtractcol(str,'\t',2);
				q_st2 = stoi(temp,nullptr);
				totAln = aln_len+totAln;
			}	
		ctMum++; //counts total mums inside a cluster		
		if(!fin.eof())
		{
			getline(fin,str); //this will also read the name of the next seq
		}
		}
		if(ctMum == 1) //if there is only one mum in the cluster
		{
			ref_end = ref_st1 + aln_len;
			q_end = q_st1 + aln_len;
		}
		if(ctMum >1) //if there are more than one mum in the cluster
		{
			ref_end = ref_st2 + aln_len;
			q_end = q_st2 + aln_len;
		}
		fout<<string(argv[3])<<"\t"<<ref_st1<<"\t"<<ref_end<<"\t"<<q_st1<<"\t"<<q_end<<"\t"<<name<<"\t"<<totAln<<"\t"<<(ref_end-ref_st1)<<endl;
		
		if(str[0] == '>' && !fin.eof()) //this is to extract the name when getline reads a seq names within the second while loop
		{
			name=str.substr(1); //remove > from the string
	        }

		ref_st1= 0; //resetting them so that they print 0 when an alignment is absent
		ref_end = 0;
		q_st1 = 0;
		q_end = 0;
		totAln = 0;
	}
	fout.close();
	fin.close();

	fin.open(argv[2]);
	ctMum = 0; //reset ctMum
		while(getline(fin,str))
		{
			if(ctMum>0)
			{
			cluster.refName[ctMum] = xtractcol(str,'\t',1);
			temp = xtractcol(str,'\t',6);
			cluster.qName[ctMum] = temp.substr(1);
			temp = xtractcol(str,'\t',2);
			cluster.refClust[ctMum].push_back(stoi(temp,nullptr));
			temp = xtractcol(str,'\t',3);
			cluster.refClust[ctMum].push_back(stoi(temp,nullptr));
			temp = xtractcol(str,'\t',4);
			cluster.qClust[ctMum].push_back(stoi(temp,nullptr));
			temp = xtractcol(str,'\t',5);
			cluster.qClust[ctMum].push_back(stoi(temp,nullptr));
			}
		ctMum++;
		}
	name = string(argv[3]);
	comparClust(cluster);
	outFileName = "SV_report."+ name +".tsv";
	fout.open(outFileName.c_str());
	filterDup(cluster);
	for(unsigned int k=0;k<cluster.dupCord.size();k++)
	{
		if(cluster.dupCord[k][0] != 0)
		{
			fout<<cluster.dupName[k][0]<<"\t"<<cluster.dupCord[k][0]<<"\t"<<cluster.dupCord[k][1]<<"\t"<<cluster.dupName[k][1]<<"\t"<<cluster.dupCord[k][2]<<"\t"<<cluster.dupCord[k][3]<<"\t"<<cluster.dupName[k][2]<<"\t"<<cluster.dupCord[k][4]<<"\t"<<cluster.dupCord[k][5]<<endl;
		}
	}
	fin.close();
	fout.close();



return 0;
}
