#include<iostream>
#include "val.h"

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
//////////////////////////////////////////////////////////////////////////////////
void countCopy(string & str, asmMerge & merge)
{
	string tempname;
	int tot_count = 0;

	for(unsigned int i = 0;i<merge.r_name.size();i++)
	{
		tempname = merge.r_name[i] + merge.q_name[i];

		if(tempname.find(str) != string::npos)
		{
		merge.seqCount[str]++;
			for(unsigned int j=0;j<merge.ref_st[tempname].size();j++)
			{
				if((merge.ref_end[tempname][j] - merge.ref_st[tempname][j]) > int(merge.ref_len[tempname]*0.9)) //the length of the duplicates have to be at least half
				{
					//tot_count = merge.ref_st[tempname].size() + tot_count;
					tot_count++;
				}
				else
				{
					merge.ref_end[tempname][j] = 0; // experimental
					merge.ref_st[tempname][j] = 0;
					merge.q_st[tempname][j] = 0;
					merge.q_end[tempname][j] = 0;
				}
			}
		merge.storeName[str].push_back(tempname);
		}
	}

	merge.storeCount[str] = tot_count;
	
}
////////////////////////////////////////////////////////////////////////////////
void collapseRange(asmMerge & merge,double & rCO,int & qCO) // this collapses two consecutive query ranges if they are close to each other
{
	string tempname;
	int len =0;
	
	for(unsigned int i =0; i<merge.r_name.size();i++)
	{	
		tempname = merge.r_name[i] + merge.q_name[i];
		len = merge.ref_len[tempname];
		for(unsigned j = 1; j<merge.ref_st[tempname].size();j++)
		{
			if((chkOvl(merge,tempname,j,rCO,qCO) == 0) && ((merge.ref_end[tempname][j]-merge.ref_st[tempname][j-1])<len)) // if the ranges don't overlap
			{
				merge.ref_end[tempname][j-1] = merge.ref_end[tempname][j];
				merge.q_end[tempname][j-1] = merge.q_end[tempname][j];
				merge.ref_st[tempname].erase(merge.ref_st[tempname].begin()+j);
				merge.ref_end[tempname].erase(merge.ref_end[tempname].begin()+j);
				merge.q_st[tempname].erase(merge.q_st[tempname].begin()+j);
				merge.q_end[tempname].erase(merge.q_end[tempname].begin()+j);
				j--;
			}
		}
	}
}
///////////////////////////////////////////////////////////////////////////////
bool chkOvl(asmMerge & merge, string & str,unsigned int & j, double & rCO, int & qCO) //to provide the condition for merging two query ranges
{
	int end1 = merge.ref_end[str][j-1];
	int st2 = merge.ref_st[str][j];
	int qend1 = merge.q_end[str][j-1];
	int qst2 = merge.q_st[str][j];
	int r_len = merge.ref_len[str];

	//if((!(st2 < end1)) && (st2 - end1 <int(merge.ref_len[str]*0.12)) && ( abs(qst2 - qend1) <int(merge.ref_len[str]*0.20))) // the gap between the ranges has to be less than 12% of the reference length and 20% for the query
	if((abs(st2-end1)< int (rCO*r_len)) && ( abs(qst2 - qend1) < qCO)) //cutoff distance for merging reference and query
	{
		return 0;
	}
	else
	{
		return 1;
	}
	 	
}
/////////////////////////////////////////////////////////////////////////////
int ovlCalculator(vector<int>& q_st, vector<int>& q_end)
{
int ovl = 0;

        for(unsigned int j =0;j<q_st.size();j++)
        {
        ovl = ovl + abs(q_st[j] - q_end[j]);
        }


return ovl;
}
/////////////////////////////////////////////////////////////////////////////////////
void ovlStoreCalculator(asmMerge & merge)
{
	string tempname;
		for(unsigned int i =0;i<merge.r_name.size();i++)
		{
			tempname = merge.r_name[i]+merge.q_name[i];
			merge.ovlStore[tempname] = ovlCalculator(merge.q_st[tempname],merge.q_end[tempname]);
		}
}
////////////////////////////////////////////////////////////////////////////////////

