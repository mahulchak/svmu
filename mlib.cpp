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
string name;
int refD; //qD
int dupCt = 0;
int filterCt = 0;
	for(unsigned int i =1; i<cluster.refName.size();i++)
	{
		rv = cluster.refClust[i];
		qv = cluster.qClust[i];
		name = cluster.qName[i];
		//qD = rv[1] -rv[0];
		for(unsigned int j=1;j<cluster.refClust.size();j++)
		{
			refD = cluster.refClust[j][1] - cluster.refClust[j][0];//the length of the cluster
//cout << int(refD*0.5)<<endl;
			//if((!(rv[1]<cluster.refClust[j][0]) && !(rv[0]>cluster.refClust[j][1])) && (cluster.qName[j] == name) && (refD>500)) //if the rv ref cluster does not fall completely outside a master ref cluster range. same query. ref cluster length should be at least 500bp
			if((!(rv[1]<cluster.refClust[j][0]) && !(rv[0]>cluster.refClust[j][1])) && (refD>500))
			{
//if(qv[0] == 7104120){cout<<refD<<overlapD(rv,mRef[j])<<endl;
				if(((qv[1] < (cluster.qClust[j][0]+5)) ||(qv[0] >(cluster.qClust[j][1]-5))) && (overlapD(rv,cluster.refClust[j])>50)) //if the qv cluster falls completely outside the master query cluster range;5 is added for random matches
				{	
					cluster.dupName[dupCt].push_back(cluster.refName[j]);
					cluster.dupName[dupCt].push_back(name); // the query being compared to
					cluster.dupName[dupCt].push_back(cluster.qName[j]); // the comparing query
					dup_term = findDupEnds(rv[0],rv[1],cluster.refClust[j][0],cluster.refClust[j][1],qv[0],qv[1],cluster.qClust[j][0],cluster.qClust[j][1]);
					cluster.dupCord[dupCt].push_back(dup_term[0]);
					cluster.dupCord[dupCt].push_back(dup_term[1]);
					cluster.dupCord[dupCt].push_back(dup_term[2]);
					cluster.dupCord[dupCt].push_back(dup_term[3]);
					cluster.dupCord[dupCt].push_back(dup_term[4]);
					cluster.dupCord[dupCt].push_back(dup_term[5]);
					dupCt++;
					
				}
			}
			
			if((!(qv[1]<cluster.qClust[j][0]) && !(qv[0] > cluster.qClust[j][1])) && (cluster.qName[j] == name))
			{
				if(((rv[1]<cluster.refClust[j][0]) || (rv[0] > cluster.refClust[j][1])) && (overlapD(qv,cluster.qClust[j]) > 50)) //if query clusters overlap but references don't
				{
					cluster.filterList[filterCt].push_back(rv[0]);
					cluster.filterList[filterCt].push_back(rv[1]);
					cluster.filterList[filterCt].push_back(cluster.refClust[j][0]);
					cluster.filterList[filterCt].push_back(cluster.refClust[j][1]);
					filterCt++;
				}
			}
//do the same thing for clusters that are duplicated in reference but single copy in query (to remove FPs due to duplicated gene segments)			
 		}
	}
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

vector<int> findDupEnds(int & ref_st1, int & ref_end1,int & ref_st2, int & ref_end2, int & q_st1,int & q_end1, int & q_st2, int & q_end2) //reports the conservative estimate of the duplicate span 
{
	vector<int> dup_ends;
	int dup_in_query1;
	int dup_in_query2;
	if(ref_st2<ref_end1)
	{
		dup_ends.push_back(ref_st2); // predicted start point of the duplicate
		dup_ends.push_back(ref_end1); //conservative estimate of the end point of the duplicate
		dup_in_query1 = q_end1 - (ref_end1-ref_st2); //subtract dup length from end of query cluster
		dup_ends.push_back(dup_in_query1); //dup_start
		dup_ends.push_back(q_end1);
		dup_in_query2 = q_st2 + (ref_end1-ref_st2);
		dup_ends.push_back(q_st2);
		dup_ends.push_back(dup_in_query2);
	}
	if(ref_end2 > ref_st1) // if the querying cluser is upstream of the resident cluster
	{
		dup_ends.push_back(ref_st1);
		dup_ends.push_back(ref_end2);
		dup_in_query1 = q_st1 + (ref_end2-ref_st1);
		dup_ends.push_back(q_st1);
		dup_ends.push_back(dup_in_query1);
		dup_in_query2 = q_end2 - (ref_end2-ref_st1);
		dup_ends.push_back(dup_in_query2);
		dup_ends.push_back(q_end2);
	}
	
	return dup_ends;
}

void filterDup(mgapC & cluster) // this will compare the reciprocal calls to remove false calls from the dup calls list
{
	vector<int> tempFilter1;
	vector<int> tempFilter2;
	vector<int> masterDup;
	
	for(unsigned int i = 0;i<cluster.filterList.size();i++)
	{
		tempFilter1.push_back(cluster.filterList[i][0]);
		tempFilter1.push_back(cluster.filterList[i][1]);
		tempFilter2.push_back(cluster.filterList[i][2]);
		tempFilter2.push_back(cluster.filterList[i][3]);
		for(unsigned int j = 0;j<cluster.dupCord.size();j++)
		{
			masterDup.push_back(cluster.dupCord[j][0]);
			masterDup.push_back(cluster.dupCord[j][1]);
			if((overlapD(tempFilter1,masterDup) > 0)||(overlapD(tempFilter2,masterDup) > 0))
			{
				cluster.dupCord[j][0] = 0; // reset the first number to 0, to identify FP
			}
			masterDup.clear();
		}
			
		tempFilter1.clear();
		tempFilter2.clear();
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
		

