#include<iostream>
#include "indel.h"

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

void writeToFile(asmMerge & merge)
	{
	ofstream fout;
	fout.open("aln_summary.tsv");
	fout<<"REF"<<'\t'<<"QUERY"<<'\t'<<"REF-LEN"<<'\t'<<"Q-LEN"<<'\t'<<"REF-ST"<<'\t'<<"REF-END"<<'\t'<<"Q-ST"<<'\t'<<"Q-END"<<endl;
		string tempname;

		for(unsigned int i =0;i<merge.r_name.size();i++)
		{
			tempname = merge.r_name[i] + merge.q_name[i];
			for(unsigned int j =0; j<merge.ref_st[tempname].size();j++)
			{
				fout<<merge.r_name[i]<<'\t'<<merge.q_name[i]<<'\t'<<merge.ref_len[merge.r_name[i]]<<'\t'<<merge.q_len[tempname]<<'\t'<<merge.ref_st[tempname][j]<<'\t'<<merge.ref_end[tempname][j]<<'\t'<<merge.q_st[tempname][j]<<'\t'<<merge.q_end[tempname][j]<<'\t'<<endl;
			}
		}
	fout.close();
	} 
		
			
/////////////////////////////////////////////////////////////////////////this function checks whether a query is completely contained within a reference

void findIndel(asmMerge & merge)
{
	string tempname;
	int gapRef =0;
	int gapQ = 0;
	char indel = 'n';//deletion is 'd' and insertion is 'i'. 'n' is not found
	for(unsigned int i =0;i<merge.r_name.size();i++)
	{
		tempname = merge.r_name[i] + merge.q_name[i];
		for(unsigned int j = 0; j<merge.ref_st[tempname].size();j++)
		{
			for(unsigned int k=0; k<merge.ref_st[tempname].size();k++)
			{
				//jth cluster is downstream of kth cluster;
				//condition 1: there should be gap between start and end of the two clusters. umbrella condition.

				if(((merge.ref_st[tempname][j] - merge.ref_end[tempname][k] >0) || (merge.q_st[tempname][j] - merge.q_end[tempname][k]) >0))
				{
					//calculate the distance between the reference and query clusters

					gapRef = merge.ref_st[tempname][j] - merge.ref_end[tempname][k];
					gapQ = merge.q_st[tempname][j] - merge.q_end[tempname][k];
		// don't use a translocated cluster for indel calls. a bunch of translocated clusters is probably okay though
					if (((abs(gapRef) <50000) && (abs(gapQ) <50000)) && ((gapRef>0) || (gapQ>0)))
					{
						indel =	checkIndel(tempname,merge,k,j); 
					}
					if(indel == 'i') //this means deletion in the reference or insertion in query
					{
					//	cout<<"Ins"<<"\t"<<merge.r_name[i]<<"\t"<<merge.ref_end[tempname][k]<<"\t"<<merge.ref_st[tempname][j]<<"\t"<<gapRef<<"\t"<<merge.q_end[tempname][k]<<"\t"<<merge.q_st[tempname][j]<<endl;
						merge.storeInsStart[merge.r_name[i]].push_back(merge.ref_end[tempname][k]);
						merge.storeInsEnd[merge.r_name[i]].push_back(merge.ref_st[tempname][j]);
					}
					if(indel == 'd') //this means deletion in the query or insertion in the reference
					{
						//cout<<"Del"<<"\t"<<merge.q_name[i]<<"\t"<<merge.q_end[tempname][k]<<"\t"<<gapQ<<"\t"<<merge.r_name[i]<<"\t"<<merge.ref_end[tempname][k]<<endl;
				
		
					}
					indel = 'n';
				}	
			}		
		}
	}
		
}
//////////////////////////////////////////////////////////////////////////////////////////
char checkIndel(string tempname, asmMerge & merge, int k, int j)
{
//don't want to use clusters where gap between two clusters lacking the deletion is more than 5 bp
//if the deletion is in the reference then gap between the two clusters will be bigger

	int gapRef = merge.ref_st[tempname][j] - merge.ref_end[tempname][k];
        int gapQ = merge.q_st[tempname][j] - merge.q_end[tempname][k];

//either the queries should overlap or the queries should not be more than 5bp apart
	if (((chkOvlQ(tempname,merge,k,j) == 1) || (gapQ >0 && gapQ<5)) && (gapRef >0 && gapRef>gapQ))
	{
		 return 'i';
	}
//either the reference clusters should overlap or they should not be more than 5 bp apart
	if(((chkOvlR(tempname,merge,k,j) == 1) || (gapRef >0 && gapRef<5)) && (gapQ >0 && gapQ>gapRef))
	{
		return 'd';
	}
	else
	{
		return 'n';
	}
}
//////////////////////////////////////////////////////////////////////////////////////
void fillChromPos(asmMerge & merge)
{
	for (map<string,int>::iterator it = merge.ref_len.begin();it!= merge.ref_len.end();it++)
	{
		for(int i = 0;i<it->second ;i++)
		{
			merge.refChromPos[it->first].push_back(0); //because chromosome positions begin at 1
		}
//cout<<"Done with chromosome "<<it->first<<endl;
	}
	
	//for(map<string,int>::iterator it1 = merge.q_len.begin();it1 != merge.q_len.end();it1++)
	//{
	//	for(int i =0;i<it1->second; i++)
	//	{
	//		merge.qChromPos[it1->first].push_back(0);
	//	}
      //	}
}
//////////////////////////////////////////////////////////////////////////////////////

void buildCoverage(asmMerge & merge)
{
	string tempname;
	for(unsigned int i = 0;i<merge.r_name.size();i++)
	{
			
		tempname = merge.r_name[i] + merge.q_name[i];
		for(unsigned int j =0; j<merge.ref_st[tempname].size();j++)
		{
			addCoverage(merge,merge.r_name[i],merge.ref_st[tempname][j],merge.ref_end[tempname][j]);
		}
//	cout<<"finished processing "<<i<<" th alignment"<<endl;
	}
}
////////////////////////////////////////////////////////////////////////////////////////////

void addCoverage(asmMerge & merge,string & str, int ref_st, int ref_end)
{

	for(int i = ref_st;i<ref_end;i++)
	{
		if(i>0)//i should always be that but lets be careful
		{
		merge.refChromPos[str][i-1] = merge.refChromPos[str][i-1]+1; //added 1 to coverage
		}
	}
}


///////////////////////////////////////////////////////////////////////////////////////////
bool chkOvlQ(string tempname,asmMerge & merge,int k,int j)
{
	int D = abs(merge.q_st[tempname][k] - merge.q_end[tempname][k]) + abs (merge.q_st[tempname][j]-merge.q_end[tempname][j]);
	int d = maxD(merge.q_st[tempname][k],merge.q_end[tempname][k],merge.q_st[tempname][j],merge.q_end[tempname][j]);

	if(D>d) //when they overlap, sum of the two intervals would be bigger than the interval between the endpoints
	{
		return 1;
	}
	else
	{
		return 0;
	}
}
/////////////////////////////////////////////////////////////////////////////////////////
bool chkOvlR(string tempname,asmMerge & merge,int k,int j)
{
	int D = abs(merge.ref_st[tempname][k] - merge.ref_end[tempname][k]) + abs(merge.ref_st[tempname][j] - merge.ref_end[tempname][j]);
  	int d = maxD(merge.ref_st[tempname][k],merge.ref_end[tempname][k],merge.ref_st[tempname][j],merge.ref_end[tempname][j]);
        
	if(D>d)
        {
                return 1;
        }
        else
        {
                return 0;
        }
}

////////////////////////////////////////////////////////////////////////////////////////
int maxD(int & qf1,int & qe1, int & qf2, int & qe2)
{
        int d1,d2,d3,d4;
        vector<int> d;

        d1 = abs(qf1 - qe2);
        d.push_back(d1);
        
        d2 = abs(qf1 - qf2);
        d.push_back(d2);

        d3 = abs(qe1 - qf2);
        d.push_back(d3);

        d4 = abs(qe1 - qe2);
        d.push_back(d4);
        sort(d.begin(),d.end());
        return d[3];
}
///////////////////////////////////////////////////////////////////////////////////////////

void filterInsCall(asmMerge & merge)
{
	string refName;
	int start =0;
	int end=0;
	int tempStart = 0;
	int tempEnd =0;
	for(map<string,vector<int> >::iterator it=merge.storeInsStart.begin();it != merge.storeInsEnd.end();it++)
	{
		refName = it->first;
		for(unsigned int i =0;i<merge.storeInsStart[refName].size();i++)
		{
			start = merge.storeInsStart[refName][i];
			end = merge.storeInsEnd[refName][i];
			tempStart = start;
			tempEnd = end;
cout<<refName<<"\t"<<start<<"\t"<<end<<" became ";
			for(int j = 0; j< int(abs(end-start)*0.5);j++)
			{
//if the reference chromosome has coverage at the insertion range then subtract it from the range
				if(merge.refChromPos[refName][start+j] != 0)
				{
					tempStart++;
				}
				if(merge.refChromPos[refName][end-j] != 0)
				{
					tempEnd--;
				}
			}
			merge.storeInsStart[refName][i] = tempStart;
			merge.storeInsEnd[refName][i] = tempEnd;
cout<<tempStart<<"\t"<<tempEnd<<endl;
		}
	}
}








	
	
