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
bool rsort(int i,int j){return (i>j);}
/////////////////////////////////////////////////////////////////////////////////
double cov(asmMerge & merge,string & str,int st, int end)
{
	double sum =0.0,cov=0.0;

	if(find(merge.q_name.begin(),merge.q_name.end(),str) != merge.q_name.end())
	{
		if(st>end)
		{
			swap(st,end);
		}
		for(unsigned int i=st;i<end;i++)
		{
			sum = sum + merge.refChromPos[str][i];
		}
	}
	if(find(merge.r_name.begin(),merge.r_name.end(),str) != merge.r_name.end())
	{
		for(unsigned int i=st;i<end;i++)
		{
			sum = sum + merge.refChromPos[str][i];
		}
	}
	cov = sum/(end-st);
	return cov;
}	
				
//////////////////////////////////////////////////////////////////////////////////

void writeToFile(asmMerge & merge)
	{
	ofstream fout;
	fout.open("aln_summary.tsv");
	fout<<"REF"<<'\t'<<"QUERY"<<'\t'<<"REF-LEN"<<'\t'<<"Q-LEN"<<'\t'<<"REF-ST"<<'\t'<<"REF-END"<<'\t'<<"Q-ST"<<'\t'<<"Q-END"<<endl;
		string tempname;
		double covRef = 0, covQ =0;
		for(unsigned int i =0;i<merge.r_name.size();i++)
		{
			tempname = merge.r_name[i] + merge.q_name[i];
			for(unsigned int j =0; j<merge.ref_st[tempname].size();j++)
			{
				covRef = cov(merge,merge.r_name[i],merge.ref_st[tempname][j],merge.ref_end[tempname][j]);	
				covQ = cov(merge,merge.q_name[i],merge.q_st[tempname][j],merge.q_end[tempname][j]);
				fout<<merge.r_name[i]<<'\t'<<merge.q_name[i]<<'\t'<<merge.ref_len[tempname]<<'\t'<<merge.q_len[tempname]<<'\t'<<merge.ref_st[tempname][j]<<'\t'<<merge.ref_end[tempname][j]<<"\t"<<covRef<<'\t'<<merge.q_st[tempname][j]<<'\t'<<merge.q_end[tempname][j]<<'\t'<<covQ<<endl;
			}
		}
	fout.close();
	} 
		
			
/////////////////////////////////////////////////////////////////////////this function checks whether a query is completely contained within a reference

void findIndel(asmMerge & merge, char mutType, float & prop)
{
	string tempname;
	int gapRef =0;
	int gapQ = 0;
	char indel = 'n';//deletion is 'd' and insertion is 'i'. 'n' is not found
	for(unsigned int i =0;i<merge.r_name.size();i++)
	{
if(merge.r_name[i] == merge.q_name[i].substr(1))
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
		// don't use a translocated cluster for indel calls. a big chunk of translocated clusters is probably okay though
					if (((abs(gapRef) <50000) && (abs(gapQ) <50000)) && ((gapRef>0) || (gapQ>0)))
					{
						indel =	checkIndel(tempname,merge,k,j, prop); 
					}
					if((indel == 'i') && (mutType == 'D')) //this means insertion in the reference or deletion in query
					{
						//cout<<"Ins"<<"\t"<<merge.r_name[i]<<"\t"<<merge.ref_end[tempname][k]<<"\t"<<merge.ref_st[tempname][j]<<"\t"<<gapRef<<"\t"<<merge.q_end[tempname][k]<<"\t"<<merge.q_st[tempname][j]<<endl;
						merge.storeInsStart[merge.r_name[i]].push_back(merge.ref_end[tempname][k]);
						merge.storeInsEnd[merge.r_name[i]].push_back(merge.ref_st[tempname][j]);
						merge.storeDelStart[merge.r_name[i]].push_back(merge.q_end[tempname][k]);
						merge.storeDelEnd[merge.r_name[i]].push_back(merge.q_st[tempname][j]);
						merge.storeDelName[merge.r_name[i]].push_back(merge.q_name[i]);
					}
					if((indel == 'i') && (mutType == 'I')) // this means deletion in the reference and insertion in the query
					{
							merge.storeInsStart[merge.q_name[i]].push_back(merge.q_end[tempname][k]);
							merge.storeInsEnd[merge.q_name[i]].push_back(merge.q_st[tempname][j]);
												
					}
					indel = 'n';
				}	
			}		
		}
}
	}
		
}
//////////////////////////////////////////////////////////////////////////////////////////
char checkIndel(string tempname, asmMerge & merge, int k, int j, float & prop)
{
//don't want to use clusters where gap between two clusters lacking the deletion is more than 5 bp
//if the deletion is in the reference then gap between the two clusters will be bigger
	int gapQ = 0;
        int gapRef = merge.ref_st[tempname][j] - merge.ref_end[tempname][k];
        int st1 =0,st2 = 0,end1 =0,end2 =0, ovlD =0;
        if(merge.q_st[tempname][k] < merge.q_end[tempname][k])
        {
                st1 = merge.q_st[tempname][k];
                end1 = merge.q_end[tempname][k];
        }
        if(merge.q_st[tempname][k] > merge.q_end[tempname][k])  //inverting the ends 
        {
                st1 = merge.q_end[tempname][k];
                 end1 = merge.q_st[tempname][k];
        }
        if(merge.q_st[tempname][j] < merge.q_end[tempname][j])
        {
                st2 = merge.q_st[tempname][j];
                end2 = merge.q_end[tempname][j];
        }
        if(merge.q_st[tempname][j] > merge.q_end[tempname][j])
        {
                st2 = merge.q_end[tempname][j];
                end2 = merge.q_st[tempname][j];
        }
         if(chkOvlQ(tempname,merge,k,j) == 1)
         {
                ovlD = abs(max(st2,st1) - min(end1,end2));
	 }


	if((st1<end1) && (st2<end2)) // both are forward oriented
	{
        	gapQ = merge.q_st[tempname][j] - merge.q_end[tempname][k];//gapQ will be negative when they overlap

//either the queries should overlap or the queries should not be more than 5bp apart
		//if ((((chkOvlQ(tempname,merge,k,j) == 1) && (ovlD<100))|| (gapQ>0 && gapQ<5)) && (gapRef >100))
		if ((((chkOvlQ(tempname,merge,k,j) == 1) && (ovlD<100))|| (gapQ>0 && gapQ<int(prop*gapRef))) && (gapRef >100))
		{
		 return 'i';
		}
	//either the reference clusters should overlap or they should not be more than 5 bp apart
		//if(((chkOvlR(tempname,merge,k,j) == 1) || (gapRef >0 && gapRef<5)) && (gapQ >0 && gapQ>gapRef))
		if(((chkOvlR(tempname,merge,k,j) == 1) || (gapRef >0 && gapRef<int(prop*gapRef))) && (gapQ >0 && gapQ>gapRef))
		{
		return 'd';
		}
		else
		{
		return 'n';
		}
	}
	if((st1>end1) && (st2>end2)) // if both are reverse oriented
	{
		gapQ = merge.q_end[tempname][k] - merge.q_st[tempname][j]; // gapQ is negative when they overlap
		if (((chkOvlQ(tempname,merge,k,j) == 1) || (gapQ<5)) && (gapRef >100))
                {
                 return 'i';
                }
		if(((chkOvlR(tempname,merge,k,j) == 1) || (gapRef >0 && gapRef<5)) && (gapQ >0 && gapQ>gapRef))
                {
                return 'd';
                }
                else
                {
                return 'n';
                }
	}
	//if(((st1>end1) && (st2<end2)) || ((st1>end1) && (st2<end2)))
	else
	{
		return 'n';
	}
}
///////////////////////this function creates a chromosome///////////////////////////////
void fillChromPos(string & tempname,asmMerge & merge,int & length)
{
	for(int i = 0;i<length;i++)
		{
			merge.refChromPos[tempname].push_back(0); //because chromosome positions begin at 1
		}
//cout<<"Done with chromosome "<<it->first<<endl;
}
/////this function creates coverage for each position on each reference chromosome/////

void buildCoverage(asmMerge & merge)
{
	string tempname;
	for(unsigned int i = 0;i<merge.r_name.size();i++)
	{
			
		tempname = merge.r_name[i] + merge.q_name[i];
		//for(unsigned int j =0; j<merge.newrefSt[tempname].size();j++)
		for(unsigned int j = 0;j<merge.anmlrefSt[tempname].size();j++)
		{
			addCoverage(merge,tempname,merge.anmlrefSt[tempname][j],merge.anmlrefEnd[tempname][j]);
		}
//	cout<<"finished processing "<<i<<" th alignment"<<endl;
	}
}
/////////////calculate coverage for each reference position//////////////////////////////

void addCoverage(asmMerge & merge,string & str, int ref_st, int ref_end)
{

	for(int i = ref_st;i<ref_end;i++)
	{
		if(i>0)//i should always be that but lets be careful
		{
		merge.refChromPos[str][i-1] = merge.refChromPos[str][i-1]+1; //added 1 to coverage
//cout<<str<<"\t"<<ref_st<<"\t"<<ref_end<<endl;
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

	d1 = abs(qf1 - qe1);
	d.push_back(d1);

        d1 = abs(qf1 - qe2);
        d.push_back(d1);
        
        d2 = abs(qf1 - qf2);
        d.push_back(d2);

        d3 = abs(qe1 - qf2);
        d.push_back(d3);

        d4 = abs(qe1 - qe2);
        d.push_back(d4);
        
	d4 = abs(qf2 - qe2);
	d.push_back(d4);
	sort(d.begin(),d.end());
        return d[5];
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
//cout<<tempStart<<"\t"<<tempEnd<<endl;
		}
	}
}

////////////////////////////////////////////////////////
void collapseRange(asmMerge & merge)
{

int start = 0, end =0;
for(map<string,vector<int> >::iterator it = merge.storeInsStart.begin();it != merge.storeInsStart.end();it++)
{
	for(unsigned int j = 0; j<merge.storeInsStart[it->first].size();j++)
	{
		start = merge.storeInsStart[it->first][j];
		end = merge.storeInsEnd[it->first][j];
		for(unsigned int k = 0;k<merge.storeInsStart[it->first].size();k++)
		{
			if(!(merge.storeInsStart[it->first][k] > end) && !(merge.storeInsEnd[it->first][k] < start))//if they overlap
			{
				if(!((merge.storeInsStart[it->first][k] == start) && (merge.storeInsEnd[it->first][k] == end)))//if they are not identical
				{
					merge.storeInsStart[it->first][j] = max(start,merge.storeInsStart[it->first][k]);//the bigger of the two start
					merge.storeInsStart[it->first][k] = max(start,merge.storeInsStart[it->first][k]);//the bigger of the two start
					merge.storeInsEnd[it->first][j] = min(end,merge.storeInsEnd[it->first][k]);
					merge.storeInsEnd[it->first][k] = min(end,merge.storeInsEnd[it->first][k]);
				}
			}
			
		}
	}
}
}
///////////////////////////////////////////////////////
void lisCalculator(asmMerge & merge,string & tempname,vector<int>& ref_st,vector<int>& ref_end,vector<int>& q_st, vector<int>& q_end,int dist)
{
int ovl = 0;
long int storedQend = 0,storedRefEnd = 0;
vector<int> v;
map<int,vector<int> > ovl2refSt,ovl2anmlrefSt;
map<int,vector<int> > ovl2qSt,ovl2anmlqSt;
map<int,vector<int> > ovl2refEnd,ovl2anmlrefEnd;
map<int,vector<int> > ovl2qEnd,ovl2anmlqEnd;
vector<int> temprefSt,anmlrefSt;
vector<int> temprefEnd,anmlrefEnd;
vector<int> tempqSt,anmlqSt;
vector<int> tempqEnd,anmlqEnd;
vector<int> temp;
vector<int> tempR;
        if(q_st.size() >1)
        {
                for(unsigned int j =0;j<q_st.size();j++)
                {
                        if(q_st[j]<q_end[j]) //if this MUM is forward oriented
                        {
                                storedQend = 0;
                                storedRefEnd = 0;
				temp = q_st;
				tempR = ref_st;
				sort(temp.begin()+j,temp.end());
				sort(tempR.begin()+j,tempR.end());
                                for(unsigned k = j;k<q_st.size();k++)
                                {

                                       if(!(q_st[k] < storedQend) && !(ref_st[k] < storedRefEnd) && !(pos(q_st[k],temp)>k) && !(pos(ref_st[k],tempR)>k)) 
                                        {
                                                ovl = ovl + abs(q_st[k] - q_end[k]);
                                                storedQend = q_end[k];
                                                storedRefEnd = ref_end[k];
						temprefSt.push_back(ref_st[k]);
						temprefEnd.push_back(ref_end[k]);
						tempqSt.push_back(q_st[k]);
						tempqEnd.push_back(q_end[k]);
						if(k==j)
						{
							j = q_st.size();
						}
						
                                        }
					 else
                                        {
                                                anmlrefSt.push_back(ref_st[k]);
                                                anmlrefEnd.push_back(ref_end[k]);
                                                anmlqSt.push_back(q_st[k]);
                                                anmlqEnd.push_back(q_end[k]);
                                         }
                                }

                        }
                        if(q_st[j]>q_end[j]) //if this MUM is reverse oriented
                        {
				temp = q_st;
				tempR = ref_st;
				sort(temp.begin(),temp.end(),rsort);
				sort(tempR.begin(),tempR.end(),rsort);
                                storedQend = 100000000000;
                                storedRefEnd = 0;
                                for(unsigned k = j;k<q_st.size();k++)
                                {
                                        if(j == k)
                                        {
                                                if((q_st[k] < storedQend) && !(pos(q_st[k],temp)>k) && (ref_st[k] > storedRefEnd) && !(pos(ref_st[k],tempR)>k))
                                                {
                                                        ovl = ovl + abs(q_st[k] - q_end[k]);
                                                        storedQend = q_end[k]; 
                                                        storedRefEnd = ref_end[k];
							temprefSt.push_back(ref_st[k]);
							temprefEnd.push_back(ref_end[k]);
							tempqSt.push_back(q_st[k]);
							tempqEnd.push_back(q_end[k]);
                                                }
						else
						{
							anmlrefSt.push_back(ref_st[k]);
							anmlrefEnd.push_back(ref_end[k]);
							anmlqSt.push_back(q_st[k]);
							anmlqEnd.push_back(q_end[k]);
						}
                                        }
                                        if(k>j)
                                        {
                                                if((q_st[k] < storedQend) && (abs(q_st[k]-storedQend)<dist) && (ref_st[k] > storedRefEnd))
                                                {
                                                        ovl = ovl + abs(q_st[k] - q_end[k]);
                                                        storedQend = q_end[k];
                                                        storedRefEnd = ref_end[k];
							temprefSt.push_back(ref_st[k]);
							temprefEnd.push_back(ref_end[k]);
							tempqSt.push_back(q_st[k]);
							tempqEnd.push_back(q_end[k]);
                                                }
						else
						{
							anmlrefSt.push_back(ref_st[k]);
                                                        anmlrefEnd.push_back(ref_end[k]);
                                                        anmlqSt.push_back(q_st[k]);
                                                        anmlqEnd.push_back(q_end[k]);
						}

                                        }
                                }
                        }
//cout<<tempname<<"\t"<<q_st[j]<<"\t"<<ovl<<endl;
			ovl2refSt[ovl] = temprefSt;
			ovl2anmlrefSt[ovl] = anmlrefSt;
			anmlrefSt.clear();
			temprefSt.clear();
			
			ovl2qSt[ovl] = tempqSt;
			ovl2anmlqSt[ovl] = anmlqSt;
			anmlqSt.clear();
			tempqSt.clear();

			ovl2refEnd[ovl] = temprefEnd;
			ovl2anmlrefEnd[ovl] = anmlrefEnd;
			anmlrefEnd.clear();
			temprefEnd.clear();

			ovl2qEnd[ovl] = tempqEnd;
			ovl2anmlqEnd[ovl] = anmlqEnd;
			anmlqEnd.clear();	
			tempqEnd.clear();
			
			
                        v.push_back(ovl);
                        ovl = 0;

                }
        }
        if(q_st.size() == 1)
        {
                ovl = abs(q_st[0] - q_end[0]);
                ovl2refSt[ovl].push_back(ref_st[0]);
		ovl2qSt[ovl].push_back(q_st[0]);
		ovl2refEnd[ovl].push_back(ref_end[0]);
		ovl2qEnd[ovl].push_back(q_end[0]);

                v.push_back(abs(q_st[0] - q_end[0]));
        }
                sort(v.begin(),v.end());
                ovl = v[v.size()-1];
	
	merge.newrefSt[tempname] = ovl2refSt[ovl];
	merge.anmlrefSt[tempname] = ovl2anmlrefSt[ovl];
        ovl2refSt.clear();
	ovl2anmlrefSt.clear();
        merge.newrefEnd[tempname] = ovl2refEnd[ovl];
   	merge.anmlrefEnd[tempname] = ovl2anmlrefEnd[ovl];
	ovl2refEnd.clear();
	ovl2anmlrefEnd.clear();
        merge.newq_St[tempname] = ovl2qSt[ovl];
	merge.anmlq_St[tempname] = ovl2anmlqSt[ovl];
	ovl2anmlqSt.clear();
   	ovl2qSt.clear();
        merge.newq_end[tempname] = ovl2qEnd[ovl];
	merge.anmlq_end[tempname] = ovl2anmlqEnd[ovl];
	ovl2anmlqEnd.clear();
	ovl2qEnd.clear();
}

//////////////////////////////////////////////////
void ovlStoreCalculator(asmMerge & merge)
{
        string tempname;
	int dist = 0;
        for(unsigned int i =0;i<merge.r_name.size();i++)
        {
                tempname = merge.r_name[i]+merge.q_name[i];
		dist = midDist(merge.ref_st[tempname],merge.ref_end[tempname],merge.q_st[tempname],merge.q_end[tempname]);
                //merge.ovlStore[tempname] = ovlCalculator(merge.q_st[tempname],merge.q_end[tempname]);
                lisCalculator(merge,tempname,merge.ref_st[tempname],merge.ref_end[tempname],merge.q_st[tempname],merge.q_end[tempname],dist);
        }
}
/////////////////////////////////////////////
int midDist(vector<int>& ref_st,vector<int>& ref_end,vector<int>& q_st, vector<int>& q_end)
{
	vector<int> dist;
	int length = 0;
	for(unsigned int i =0; i<ref_st.size()-1;i++)
	{
		dist.push_back(abs(q_end[i] - q_st[i+1]));		
	}
	
sort(dist.begin(),dist.end());
length = dist.size()-1;
	if(length == 1)
	{
		length = 0;
	}
	if(length > 1)
	{
		length = int((dist.size()-1) *0.80) +1;
	}
//cout<<q_st[0]<<"\t"<<dist[length]<<endl;
return dist[length];

}	
////////////////////////////////////////////////
void addvTol(vector<int> & q_st,vector <int> & temp,int k)
{
	for(unsigned int i = k;i<q_st.size();i++)
	{
		temp.push_back(q_st[i]);
	}
	sort(temp.begin(),temp.end());
}
/////////////////////////////////////////////////
unsigned int pos(int & elem, vector<int> & v)
{
	unsigned int i = 0;
	while(elem != v[i])
	{
		i++;
	}
return i;
}
///////////////////////////////////////////////
void findDup(asmMerge & merge)
{
	string tempname,refName,qName;
	int refStart =0,refEnd =0;
	for(map<string,vector<int> >::iterator it = merge.anmlrefSt.begin();it!= merge.anmlrefSt.end();it++)
	{
		tempname = it->first;
		refName = tempname.substr(0,tempname.find(" ")-1);
		qName = tempname.substr(tempname.find(" "));
		for(unsigned int i =0; i < merge.anmlrefSt[tempname].size();i++)
		{
			refStart = merge.anmlrefSt[tempname][i];
			refEnd = merge.anmlrefEnd[tempname][i];
			//refName = tempname.substr(0,(tempname.find(" ")-1));
cout<<	refName<<"\t"<<refStart<<"\t"<<refEnd<<endl;
		}
	}
}
					
