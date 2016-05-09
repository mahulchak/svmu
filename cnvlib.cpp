#include<iostream>
#include "qmerge.h"

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
			for(unsigned int j=0;j<merge.ref_st[tempname].size();j++)
			{
				if((merge.ref_end[tempname][j] - merge.ref_st[tempname][j]) > int(merge.ref_len[tempname]*0.5)) //the length of the duplicates have to be at least half
				{
					//tot_count = merge.ref_st[tempname].size() + tot_count;
					tot_count++;
			//cout<< merge.r_name[i]<<"\t"<<merge.q_name[i]<<"\t"<<merge.ref_st[tempname].size()<<endl;
				}
			}
		}
	}

	merge.storeCount[str] = tot_count;
}
	
