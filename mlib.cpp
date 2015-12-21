#include<iostream>
#include "sv.h"

string xtractcol(string str, char c, int n)
{
int count =0;
int j =0; // j is the index for the string that'll be returned
string elem;
for(unsigned int i=0;i<str.size() && count<n;i++)
        {
        if(str[i] == c)
                {
                count++; //keep a count of the field separator
                }
        if(count == n-1)
                {
		elem.push_back(str[i]);
                j++;
                }
        }
return elem;
}


void comparClust(map<int,string>& ref_name,map<int,string>& q_name, map<int,vector<int> > & mRef, map<int,vector<int> > & mQ) //this function compares the terminal coordinates of a cluster (v1 = ref cords, v2 = query cords)  with exisiting list of clusters (the first two arguments)
{
vector<int> rv;
vector<int> qv;
string name;
ofstream fout;
fout.open("SV_report.ref.tsv");
int refD,qD;
	for(unsigned int i =1; i<ref_name.size();i++)
	{
		rv = mRef[i];
		qv = mQ[i];
		name = q_name[i];
//cout<<mRef[i][0]<<endl;
		qD = rv[1] -rv[0];
		for(unsigned int j=1;j<mRef.size();j++)
		{
			refD = mRef[j][1] - mRef[j][0];//the length of the cluster
//cout << int(refD*0.5)<<endl;
			if((!(rv[1]<mRef[j][0]) && !(rv[0]>mRef[j][1])) && (q_name[j] == name) && (refD>500)) //if the rv ref cluster does not fall completely outside a master ref cluster range. same query. ref cluster length should be at least 500bp
			{
//if(qv[0] == 7104120){cout<<refD<<overlapD(rv,mRef[j])<<endl;
				if(((qv[1] < (mQ[j][0]+5)) ||(qv[0] >(mQ[j][1]-5))) && (overlapD(rv,mRef[j])>50)) //if the qv cluster falls completely outside the master query cluster range;5 is added for random matches
				{
					fout<<ref_name[j]<<"\t"<<rv[0]<<"\t"<<rv[1]<<"\t"<<q_name[j]<<"\t"<<qv[0]<<"\t"<<qv[1]<<"\t"<<mQ[j][0]<<"\t"<<mQ[j][1]<<endl;
				}
			}			
 		}
	}
fout.close();
}
	 
	
int overlapD(vector<int>& rv, vector<int>& mRef) //computes overlap between two reference sequence ranges
{
int D;
	if((rv[0] > mRef[0]) && (rv[0] <mRef[1])) // if the query start is within the reference range
	{
		D = (mRef[1] - rv[0]);
	}
	if((rv[1] > mRef[0]) && (rv[1] < mRef[1]))//if query end is within the reference range
	{
		D = (rv[1] - mRef[0]);
	}
	if((rv[0] <mRef[0]) && (rv[1] > mRef[1]))
	{
		D = (mRef[1] - rv[0]);
	} 
return D;
}

