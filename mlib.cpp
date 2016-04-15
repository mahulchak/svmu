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


//this function compares the terminal coordinates of a cluster (v1 = ref cords, v2 = query cords)  with exisiting list of clusters (the first two arguments)
void comparClust(mgapC & cluster)
{
vector<int> rv;
vector<int> qv;
vector<int> dup_term; // stores vector returned by dup_ends
vector<int> filter_term;
string name;
int refD; //qD
int dupCt = 0;
int filterCt = 0;
int dup_length_Ref,dup_length_Q1,dup_length_Q2,countCopy,falseCount;

	for(unsigned int i =1; i<cluster.refName.size();i++)
	{
		rv = cluster.refClust[i];
		qv = cluster.qClust[i];
		name = cluster.qName[i];
		countCopy = 0; //reset countCopy for each range
		falseCount = 0;
		for(unsigned int j=1;j<cluster.refClust.size();j++)
		{
			refD = cluster.refClust[j][1] - cluster.refClust[j][0];//the length of the cluster
//cout << int(refD*0.5)<<endl;
			if((!(rv[1]<cluster.refClust[j][0]) && !(rv[0]>cluster.refClust[j][1])) && (refD>100))
			{
				//if(((qv[1] < (cluster.qClust[j][0]+5)) ||(qv[0] >(cluster.qClust[j][1]-5))) && (overlapD(rv,cluster.refClust[j])>100)) //if the qv cluster falls completely outside the master query cluster range;5 is added for random matches
				if(overlapD(rv,cluster.refClust[j])>100)
				{	
					dup_term = findDupEnds(rv[0],rv[1],cluster.refClust[j][0],cluster.refClust[j][1],qv[0],qv[1],cluster.qClust[j][0],cluster.qClust[j][1]);
					if(!dup_term.empty())
					{
						dup_length_Ref = (dup_term[1] - dup_term[0]);
						dup_length_Q1 = (dup_term[3]-dup_term[2]);
						dup_length_Q2 = (dup_term[5] - dup_term[4]);
if((rv[0] == 25372) || (cluster.refClust[j][0] == 25372))
{
cout<<"Duplicate "<<"\t"<<rv[0]<<"\t"<<rv[1]<<"\t"<<cluster.refClust[j][0]<<"\t"<<cluster.refClust[j][1]<<"\t"<<qv[0]<<"\t"<<qv[1]<<"\t"<<cluster.qClust[j][0]<<"\t"<<cluster.qClust[j][1]<<"\t"<<dup_term[4]<<"\t"<<dup_term[5]<<"\t"<<dup_length_Ref<<"\t"<<dup_length_Q1<<"\t"<<dup_length_Q2<<endl;
}
						if((abs(dup_length_Ref - dup_length_Q1) < int(dup_length_Ref*0.5)) && (abs(dup_length_Ref - dup_length_Q2) < int(dup_length_Ref*0.5))) //inferred lengths from ref and query are similar, they may differ by up to 0.5* ref length
						{

							cluster.dupName[dupCt].push_back(cluster.refName[j]);
							cluster.dupName[dupCt].push_back(name); // the query being compared to
							cluster.dupName[dupCt].push_back(cluster.qName[j]); // the comparing query
							cluster.dupCord[dupCt].push_back(dup_term[0]);
							cluster.dupCord[dupCt].push_back(dup_term[1]);
							cluster.dupCord[dupCt].push_back(dup_term[2]);
							cluster.dupCord[dupCt].push_back(dup_term[3]);
							cluster.dupCord[dupCt].push_back(dup_term[4]);
							cluster.dupCord[dupCt].push_back(dup_term[5]);
							dupCt++;
							countCopy++; //countCopy is incremented eveery time a new hit is found
						}
					}
					
				}
			}
			
			if((!(qv[1]<cluster.qClust[j][0]) && !(qv[0] > cluster.qClust[j][1])) && (cluster.qName[j] == name)) //if the query clusters are overlapping
			{
				if(((rv[1]<cluster.refClust[j][0]) || (rv[0] > cluster.refClust[j][1])) && (overlapD(qv,cluster.qClust[j]) > 100)) //if query clusters overlap but references don't
				//if(abs(abs(qv[0] - cluster.qClust[j][0])-abs(rv[0] - cluster.refClust[j][0]))>65)
				{
					filter_term = findDupEnds(qv[0],qv[1],cluster.qClust[j][0],cluster.qClust[j][1],rv[0],rv[1],cluster.refClust[j][0],cluster.refClust[j][1]);
//if((rv[0] == 19312) || (cluster.refClust[j][0] == 19312)) {cout<<"Added as filter"<<endl;
//cout<<"Filter "<<"\t"<<filter_term[2]<<"\t"<<filter_term[3]<<"\t"<<filter_term[4]<<"\t"<<filter_term[5]<<"\t"<<(filter_term[1] - filter_term[0])<<"\t"<<(filter_term[3]-filter_term[2])<<"\t"<<(filter_term[5] - filter_term[4])<<endl;}
					if(!filter_term.empty() && (abs((filter_term[3]-filter_term[2])-dup_length_Ref)<int(0.05*dup_length_Ref)))//count it as filter only if the size condition is met
					{
						cluster.filterName[filterCt] = cluster.qName[j]; // there is only one reference so storing query name
						cluster.filterList[filterCt].push_back(filter_term[0]);
						cluster.filterList[filterCt].push_back(filter_term[1]);
						cluster.filterList[filterCt].push_back(filter_term[2]);
						cluster.filterList[filterCt].push_back(filter_term[3]);
						cluster.filterList[filterCt].push_back(filter_term[4]);
						cluster.filterList[filterCt].push_back(filter_term[5]);
						filterCt++;
						falseCount++;
//cout<<dup_length_Ref<<"\t"<<falseCount<<"\t"<<filter_term[3]-filter_term[2]<<"\t"<<filter_term[5] - filter_term[4]<<endl;
					}
				}
			}
//do the same thing for clusters that are duplicated in reference but single copy in query (to remove FPs due to duplicated gene segments)			
 		}
			if(countCopy > 3)
			{
//cout<<"dupCt is "<<dupCt<<" countCopy is "<<countCopy<<endl;
				for (int k = (dupCt-countCopy);k<dupCt;k++) //if the query hits more than 3 times then mark all the hits as false positives
				{
					cluster.dupCord[k][0] = 0;
				}
			}
			if(countCopy <4)//the querying coords which found these dups has hit these many times in the genome
			{
				for(int k = (dupCt-countCopy);k<dupCt;k++)
				{
					cluster.dupCord[k].push_back(countCopy);
					cluster.dupCord[k].push_back(falseCount);
				}
			}
	}
}
	 
	
int overlapD(vector<int>& rv, vector<int>& mRef) //computes overlap between two sequence ranges but does not check if there is an overlap
{
int D;
	if(!(rv[0]==mRef[0]) && !(rv[1]==mRef[1])) //if the ranges are not identical
	{
		D = min(mRef[1],rv[1]) - max(mRef[0],rv[0]);
	}
return D;
}

vector<int> findDupEnds(int & ref_st1, int & ref_end1,int & ref_st2, int & ref_end2, int & q_st1,int & q_end1, int & q_st2, int & q_end2) //reports the conservative estimate of the duplicate span 
{
	vector<int> dup_ends;
	int start,end,startQ1,endQ1,startQ2,endQ2;
	start = max(ref_st1,ref_st2);
	end = min(ref_end1,ref_end2);
	dup_ends.push_back(start);
	dup_ends.push_back(end);
	startQ1 = q_st1 + abs(start - ref_st1);//mapping ref start on query start
	endQ1 = q_end1 - abs(ref_end1 - end);
	dup_ends.push_back(startQ1);
	dup_ends.push_back(endQ1);
	startQ2 = q_st2 + abs(start - ref_st2);
	endQ2 = q_end2 - abs(ref_end2 - end);
	dup_ends.push_back(startQ2);
	dup_ends.push_back(endQ2);
	return dup_ends;
}

bool ovlChk(vector<int> &v1, vector<int> & v2)
{
	if((!(v2[1]<v1[0])) && (!(v2[0]>v1[1])))
	{
		return 1;
	}
	else
	{
		return 0;
	}
}

void filterDup(mgapC & cluster) // this will compare the reciprocal calls to remove false calls from the dup calls list
{
	vector<int> tempFilter1(2);
	vector<int> tempFilter2(2);
	vector<int> masterDup(2);
	vector<int> qFilterChk1(2);
	
	vector<int> qMasterFilter1(2);
	vector<int> qMasterFilter2(2);
	string qName;
	int masterQDist1,masterQDist2,qDist;
	for(unsigned int i = 0;i<cluster.filterList.size();i++)
	{
		tempFilter1[0] = cluster.filterList[i][2];//start of a duplicate in reference (FP)
		tempFilter1[1] = cluster.filterList[i][3]; // end of a duplicate in reference
		tempFilter2[0] = cluster.filterList[i][4];//start of a duplicate in reference (FP)
		tempFilter2[1] = cluster.filterList[i][5];//end of a duplicate in reference
		qFilterChk1[0] = cluster.filterList[i][0]; //start of a duplicate in query coordinates (FP)
		qFilterChk1[1] = cluster.filterList[i][1]; //end of a duplicate in query coordinates (FP)
		//refDist1 = tempFilter1[1] - tempFilter1[0];
		//refDist2 = tempFilter2[1] - tempFilter2[0];
		qDist = qFilterChk1[1] - qFilterChk1[0];
		qName = cluster.filterName[i];

		for(unsigned int j = 0;j<cluster.dupCord.size();j++)
		{
			if(cluster.dupCord[j][0] != 0) // if it is not already identified as a FP
			{
				masterDup[0] = cluster.dupCord[j][0]; //denotes the reference start cord for the dup
				masterDup[1] = cluster.dupCord[j][1]; //reference end for the dup
				qMasterFilter1[0] = cluster.dupCord[j][2]; //dup start in the query
				qMasterFilter1[1] = cluster.dupCord[j][3]; // dup end in the query
				qMasterFilter2[0] = cluster.dupCord[j][4];
				qMasterFilter2[1] = cluster.dupCord[j][5];
				//masterRefDist = cluster.dupCord[j][1] - cluster.dupCord[j][0];
				masterQDist1 = cluster.dupCord[j][3] - cluster.dupCord[j][2];
				masterQDist2 = qMasterFilter2[1] - qMasterFilter2[0];
				//if(((((ovlChk(tempFilter1,masterDup) == 1) && (abs(masterRefDist-refDist1)<int(masterRefDist*.05))) || ((ovlChk(tempFilter2,masterDup) == 1) && (abs(masterRefDist-refDist2)<int(masterRefDist*.05))))))
				//if(((((ovlChk(tempFilter1,masterDup) == 1) && (refDist1 > int(masterRefDist*.95))) || ((ovlChk(tempFilter2,masterDup) == 1) && (refDist2 >int(masterRefDist*.95))))))
				//if((ovlChk(tempFilter1,masterDup) == 1) || (ovlChk(tempFilter2,masterDup) == 1))
				//{
					if(((ovlChk(qFilterChk1,qMasterFilter1) == 1) && (qName ==cluster.dupName[j][0])) || ((ovlChk(qFilterChk1,qMasterFilter2) == 1) && (qName == cluster.dupName[j][1])))
					{
						if((abs(masterQDist1 - qDist) < int(masterQDist1*.05)) ||(abs(masterQDist2 - qDist) < int(masterQDist2*.05))) 
						{ 
							cluster.dupCord[j][7] = cluster.dupCord[j][7] + 1; // add to show the count of duplication in the reference corresponding to the dup call
cout<<cluster.dupCord[j][0]<<"\t"<<cluster.dupCord[j][1]<<"\t"<<cluster.dupCord[j][2]<<"\t"<<cluster.dupCord[j][3]<<"\t"<<cluster.dupCord[j][4]<<"\t"<<cluster.dupCord[j][5]<<endl;
						}
					}
				//}
			}
		}
			
	}	
}

/////////////////////////////////////////////////////////////////////////////////////////

void removeExactDups(mgapC & cluster)
{

	for(unsigned int i = 0; i<cluster.dupCord.size();i++)
	{
		for(unsigned int j =0; j< cluster.dupCord.size();j++)
		{
			if((cluster.dupCord[i][0] == cluster.dupCord[j][0]) && (cluster.dupCord[i][1] == cluster.dupCord[j][1]) && (cluster.dupCord[i][2] == cluster.dupCord[j][2]) && (cluster.dupCord[i][3]==cluster.dupCord[j][3]) && (cluster.dupCord[i][4] == cluster.dupCord[j][4]) && (cluster.dupCord[i][5] == cluster.dupCord[j][5]) &&(i!=j)) //if two entries are identican then its unlikely that they are real duplicates. instead of doing hard filtering, copy number can be reset to 0
			{
				cluster.dupCord[i][0] = 0;
				cluster.dupCord[j][0] = 0;
			}
		}
	}
}
				
				 
//void checkTE(mgapC & cluster) //pass the result from the comparclust function
//{
 //map<int,vector<int> > tempRef;
//tempRef = cluster.dupRef1;
//string key1,key2;
//	for(unsigned int i=0; i<cluster.dupRef1.size();i++)
//	{
//		for(unsigned int j =0; j<tempRef.size();j++)
//		{
//			if((cluster.dupRef1[i][0] == tempRef[j][0]) && (cluster.dupRef1[i][1] == tempRef[j][1])) //if the clusters are same
//			{
//				key1 = to_string(cluster.dupRef1[i][0]) ;
//				key1.append(" ");
//				key2 = to_string(cluster.dupRef1[i][1]);
//				key1.append(key2);
//				if(cluster.dupCount[key1] == 0) //if encountering for the first time
//				{
//					cluster.dupCount[key1] =1;
//				}
//				if(cluster.dupCount[key1] > 0)
//				{
//					cluster.dupCount[key1] = cluster.dupCount[key1] +1;
//				}
//			}
 
//		}
//	}
//}		
		

