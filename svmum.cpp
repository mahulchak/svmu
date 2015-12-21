#include<iostream>
#include "sv.h"

int main(int argc, char * argv[])
{

	if(argc<3)
	{
	cerr<<"Usage: "<<argv[0]<<"foo.mgaps mumTobed.Chr.bed chrom"<<endl;
	exit(EXIT_FAILURE);
	}
	ifstream fin;
	ofstream fout;

	fin.open(argv[1]);
	fout.open(argv[2]);
	
	string str,name,temp;
	int ref_st1,ref_st2,ref_end,aln_len,q_st1,q_st2,q_end,ctMum;

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
//cout<<temp<<endl;
				ref_st1 = stoi(temp,nullptr);
//cout<<ref_st1<<endl;
				temp = xtractcol(str,'\t',3);
				aln_len = stoi(temp,nullptr);
//cout<<aln_len<<endl;
				temp = xtractcol(str,'\t',2);
				q_st1 = stoi(temp,nullptr);
//cout<<q_st1<<endl;
			}
			else
			{
//cout<<str<<endl;
				temp = xtractcol(str,'\t',1);
                       		ref_st2 = stoi(temp,nullptr);
				temp = xtractcol(str,'\t',3);
				aln_len = stoi(temp,nullptr);
				temp = xtractcol(str,'\t',2);
				q_st2 = stoi(temp,nullptr);
//cout<<q_st2<<endl;
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
		fout<<string(argv[3])<<"\t"<<ref_st1<<"\t"<<ref_end<<"\t"<<q_st1<<"\t"<<q_end<<"\t"<<name<<endl;

		if(str[0] == '>' && !fin.eof()) //this is to extract the name when getline reads a seq names within the second while loop
		{
			name=str.substr(1); //remove > from the string
	        }

		ref_st1= 0; //resetting them so that they print 0 when an alignment is absent
		ref_end = 0;
		q_st1 = 0;
		q_end = 0;
	}
	fout.close();
	fin.close();
//cout<<"crossed first hurdle"<<endl;
//may be add a system command here for bedtools?
	fin.open(argv[2]);
	ctMum = 0; //reset ctMum
		while(getline(fin,str))
		{
			if(ctMum>0)
			{
			ref_name[ctMum] = xtractcol(str,'\t',1);
			temp = xtractcol(str,'\t',6);
			q_name[ctMum] = temp.substr(1);
			temp = xtractcol(str,'\t',2);
			Mref[ctMum].push_back(stoi(temp,nullptr));
//cout<<Mref[ctMum][0]<<endl;
			temp = xtractcol(str,'\t',3);
			Mref[ctMum].push_back(stoi(temp,nullptr));
			temp = xtractcol(str,'\t',4);
			Qref[ctMum].push_back(stoi(temp,nullptr));
			temp = xtractcol(str,'\t',5);
			Qref[ctMum].push_back(stoi(temp,nullptr));
			}
		ctMum++;
		}
	comparClust(ref_name,q_name,Mref,Qref);
//2. store the cluster data in a map
//3. generate a merged cluster
//4.compare 2 with 3 and call duplicates. if the query clusters don't overlap but reference clusters do. May be use a gap of 100bp between the query clusters
//write the results
//generate a formatted delta output that can be used by bedtools

return 0;
}
