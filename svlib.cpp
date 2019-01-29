#include<iostream>
#include "sv.h"
using namespace std;
using chroms = map<string,chromPair>;
using ccov = vector<int>;
using vq = vector<qord>;

/////////////////////////////////////////////////////////
bool qusort(mI mi1, mI mi2)
{
	//return (min(mi1.y1,mi1.y2) < min(mi2.y1,mi2.y2)) ||((min(mi1.y1,mi1.y2) == min(mi2.y1,mi2.y2)) && (max(mi1.y1,mi1.y2)<max(mi2.y1,mi2.y2)));
	return(max(mi1.y1,mi1.y2) < max(mi2.y1,mi2.y2)) || ((max(mi1.y1,mi1.y2) == max(mi2.y1,mi2.y2)) && (min(mi1.y1,mi1.y2) > min(mi1.y1,mi1.y2)));
}
//////////////////////////////////////////////////////
ccov makeChromBucket(int refLen)
{
	ccov v;
	for(int i=0;i<refLen;i++)
	{
		v.push_back(0);
	}
return v;
}
/////////////////////////////////////////////////////////
bool msort(mI mi1, mI mi2)
{
	return	(mi1.x2 < mi2.x2) || ((mi1.x2 == mi2.x2) && (mi1.x1 > mi2.x1));
}
/////////////////////////////////////////////////////////
bool isort(mI m1, mI m2)
{
	return (max(m1.mv[0],m1.mv[1])<max(m2.mv[0],m2.mv[1]));
}
//////////////////////////////////////////////////////////
bool iqsort(mI m1, mI m2)
{
	return (min(m1.mv[0],m1.mv[1]) < min(m2.mv[0],m2.mv[1]));
}
///////////////////////////////////////////////////////////
bool dsort(mI m1,mI m2)
{
	if((m1.x2 != m2.x2)||(m1.x1 != m2.x1))
		return m1.x2 < m2.x2;
	return m1.y2 < m2.y2;
}	
////////////////////////////////////////////////////////////
bool findInnie(vector<mI> & mums,mI mi)
{
	int i = 0;
	//int i = int(mums.size()) -1;
	bool found = false;
	//while(!(mi.x2 >mums[i].x2))//until they become just on more than equal
	while((mi.x2 > mums[i].x1) && (i<mums.size()))
	//while((mi.x2 < mums[i].x1) && (i>0))
	{
//cout<<"debug\t"<<mi.rn<<'\t'<<mi.x1<<'\t'<<mi.x2<<'\t'<<mums[i].rn<<'\t'<<mums[i].x1<<'\t'<<mums[i].x2<<endl;
		if((mi.x1 > (mums[i].x1-1)) && (mi.x2 < (mums[i].x2+1)))
		//if((mi.x1 > (mums[i].x1)) && (mi.x2 < (mums[i].x2)))
		{
//cout<<"debug\t"<<mi.rn<<'\t'<<mi.x1<<'\t'<<mi.x2<<'\t'<<mums[i].rn<<'\t'<<mums[i].x1<<'\t'<<mums[i].x2<<endl;
			if(!(mi == mums[i]))
			{
//cout<<"debug\t"<<mi.rn<<'\t'<<mi.x1<<'\t'<<mi.x2<<'\t'<<mums[i].rn<<'\t'<<mums[i].x1<<'\t'<<mums[i].x2<<endl;
				found = true;
				break;
			}
		}
		++i;
		//--i;
	}
	return found;
}
/////////////////////////////////////////////////////////////
bool findInnieQ(vector<mI> & mums,mI mi) // checks if mi query embeds into another query
{
	int i =0;
	bool found = false;
	while((max(mi.y2,mi.y1)) > (min(mums[i].y1,mums[i].y2)) && (i<mums.size()))
	{
		if((min(mi.y1,mi.y2) > (min(mums[i].y1,mums[i].y2)-1)) && (max(mi.y2,mi.y1) <(max(mums[i].y2,mums[i].y1)+1)))
		{
			if(!(mi == mums[i]))
			{
				found = true;
				break;
			}
		}
		++i;
	}
	return found;
}
////////////////////////////////////////////////////////////////////
void storeCords(ccov & masterRef,ccov & masterQ, mI & mi)
{

	int ty1 = 0, ty2 =0;
	
	if(mi.y1 > mi.y2)//if reverse oriented
	{
		ty1 = mi.y2;
		ty2 = mi.y1;
	}
	if(mi.y1 < mi.y2)//forward oriented
	{
		ty1 = mi.y1;
		ty2 = mi.y2;
	}
	for(int i = mi.x1-1; i<mi.x2;i++)
	{
		masterRef[i]++;
	}
	
	for(int j = ty1-1; j<ty2;j++)
	{
		masterQ[j]++;
	}
}
///////////////////////////////////////////////////////
vector<double> getCoverage(mI & mi, ccov & masterRef, ccov & masterQ, float p)
{
	int d = 0, cov = 0,covCount=0, medCov =0;
	double c;
	vector<double> cc;
	map<int,int> covFreq;//holds coverage frequency for the genomic interval
//cout<<"0\t"<<mi.rn<<'\t'<<mi.x1<<"\t"<<mi.x2<<"\t"<<mi.y1<<"\t"<<mi.y2<<"\n";
	d = mi.x2 - mi.x1;
	for(int i = mi.x1-1;i<mi.x2;i++)
	{
//cout<<"1\t"<<mi.x1<<"\t"<<mi.x2<<"\t"<<mi.y1<<"\t"<<mi.y2<<'\t'<<i<<"\n";
		cov = cov + masterRef[i];	
		covFreq[masterRef[i]]++;
	}
	for(map<int,int>::iterator it = covFreq.begin();it!= covFreq.end();it++)
	{
		if((covCount <int(d*p)) && (d>5))
		{
			covCount = covCount + it->second;
			medCov = it->first;
		}
	}
//cout<<"2\t"<<mi.x1<<"\t"<<mi.x2<<"\t"<<mi.y1<<"\t"<<mi.y2<<"\t"<<covCount<<"\t"<<medCov<<"\n";
	c = cov/double(d);
	if(d>5)
	{
		cc.push_back(double(medCov));
	}
	else
	{
		cc.push_back(c);
	}
	
	cov = 0;
	medCov = 0;
	covFreq.erase(covFreq.begin(),covFreq.end());//destroying the previous map
	covCount = 0;
	d = abs(mi.y1-mi.y2);
	for(int i = min(mi.y1,mi.y2)-1;i<max(mi.y1,mi.y2);i++)
	{
		cov = cov + masterQ[i];
		covFreq[masterQ[i]]++;
	}
	for(map<int,int>::iterator it = covFreq.begin();it!= covFreq.end();it++)
	{
		if((covCount <int(d*p)+1) && (d>5))
		{
			covCount = covCount + it->second;
			medCov = it->first;
		}
	}
//cout<<"3\t"<<medCov<<endl;
	c = cov/double(d);
	if(d>5)
	{
		cc.push_back(double(medCov));
	}
	else
	{
		cc.push_back(c);	
	}
return cc;
//cout<<c<<endl;
}
		
//////////////////////////////////////////////////////
vector<double> getCoverage(mI & mi, ccov & masterRef, ccov & masterQ)
{
	int d = 0, cov = 0;
	double c;
	vector<double> cc;
	d = mi.x2 - mi.x1;
	for(int i = mi.x1-1;i<mi.x2;i++)
	{
		cov = cov + masterRef[i];
	}
	c = cov/double(d);
	cc.push_back(c);
	cov = 0;
	d = abs(mi.y1-mi.y2);
	for(int i = min(mi.y1,mi.y2)-1;i<max(mi.y1,mi.y2);i++)
	{
		cov = cov + masterQ[i];
	}
	c = cov/double(d);
	cc.push_back(c);
	return cc;
}		
////////////////////////////////////////////////////
mI findClosest(mI & mi, vector<mI> & mums)
{
	map<double,mI> storeDist;
	map<double,mI>::iterator it;
	double d1 =0, d2 =0, d = 0,Dist= 0,rd1 =0,rd2=0,rd=0;
	mI invmi;//store the swapped gapmi
	Dist = sqrt(pow(abs(mi.x2-mi.x1),2)+pow(abs(mi.y1-mi.y2),2));
	invmi = mi;
	invmi.y1 = mi.y2;//invmi is reverse complement of gapmi
	invmi.y2 = mi.y1;
	for(unsigned int j=0; j <mums.size();j++)
	{
		d1 = pow(abs(mi.x1 - mums[j].x1),2);//if start of mum preceeds the gap start then effective mum start is the gap start. will do it for query too.
		d2 = pow(abs(mi.y1 - mums[j].y1),2);
		d = abs(sqrt(d1+d2));
		rd1 = pow(abs(invmi.x1 - mums[j].x1),2);
		rd2 = pow(abs(invmi.y1 - mums[j].y1),2);
		rd = abs(sqrt(rd1 + rd2));
		if(d < rd) // if forward orientation is closer
		{
			storeDist[d] = mums[j];
		}
		if(rd < d)
		{
			storeDist[rd] = mums[j];
		}
	}
	it = storeDist.begin();
	//if(it->first >(2*Dist))
//cout << mi.rn<<"\t"<<mi.x1<<"\t"<<mi.x2<<"\t"<<mi.qn<<"\t"<<mi.y1<<"\t"<<mi.y2<<"\t"<<it->second.rn<<"\t"<<it->second.x1<<"\t"<<it->second.x2<<"\t"<<it->second.y1<<"\t"<<it->second.y2<<"\t"<<it->first<<"\t"<<Dist<<"\t"<<d1<<"\t"<<d2<<endl;
//	if(it->first > Dist)
//	{
//		it->second.y1 = 0;
//	}
return it->second;
}
////////////////////////////////////////////////////
void gapCloser(mI & mi, vector<mI> ncm, vector<mI> & cm)
{
	mI tempmi;
	vector<mI> smum; //selected mums that overlap with the gap.
	if(mi.x2 - mi.x1>0)
	{
		for(unsigned int i = 0;i<ncm.size();i++)
		{
			if((!(ncm[i].x2 < mi.x1)) && (!(max(ncm[i].y1,ncm[i].y2)<min(mi.y1,mi.y2))) && (!(ncm[i].x1>mi.x2)) && (!(min(ncm[i].y1,ncm[i].y2)>max(mi.y2,mi.y1)))) //ncm mum does not fall outside
			{
				smum.push_back(ncm[i]);
//if(mi.rn == "2L")
//{
//cout<<ncm.size()<<"\t"<<i<<"\t"<<mi.rn<<"\t"<<mi.qn<<mi.x1<<"\t"<<mi.x2<<"\t"<<mi.y1<<"\t"<<mi.y2<<"\t"<<ncm[i].x1<<"\t"<<ncm[i].x2<<"\t"<<ncm[i].y1<<"\t"<<ncm[i].y2<<endl;
//cout<<ncm.size()<<endl;
//}
			}
		}
		if(smum.size()>0)
		{	
			tempmi = findClosest(mi,smum); //find the closest mum from the ncm pool
			if(tempmi.y1 != 0) // add the condition that when it is forward this happens
			{
				cm.push_back(tempmi);
//if(tempmi.rn == "2R")
//{
//	cout<<"cm\t"<<tempmi.rn<<" "<<tempmi.x1<<" "<<tempmi.x2<<" "<<tempmi.qn<<" "<<tempmi.y1<<" "<<tempmi.y2<<endl;
//}
				mi.x1 = tempmi.x2+1; //adjust the gap coordinates
				if(tempmi.y1 < tempmi.y2)//forward oriented
				{
					mi.y1 = tempmi.y2+1; //adjust the gap coordinates
				}
				if(tempmi.y1 > tempmi.y2)//reverse oriented
				{
					if(abs(mi.y1-tempmi.y1) < abs(mi.y2 - tempmi.y2))//if tempmi is closer to the leftgap end
					{
						mi.y1 = tempmi.y1 +1;
					}
					if(abs(mi.y1-tempmi.y1) > abs(mi.y2 - tempmi.y2))
					{
						mi.y2 = tempmi.y2 -1;
					}
				}
				gapCloser(mi,smum,cm);//need to make it dependent on the size of ncm if ncm size does not change, that mean no more solution is there
			}
					
		}
	}
	else
	{
		//return;
	}
	
}
/////////////////////////////////////////////////////
void gapCloserRev(mI & mi, vector<mI> ncm, vector<mI> & cm)
{
	mI tempmi;
	vector<mI> smum; //selected mums that overlap with the gap.
	if((mi.x2 - mi.x1>0) && (mi.y1 - mi.y2 >0)) //checking if both of them has gaps. this is for inverted sequences
	{
		for(unsigned int i = 0;i<ncm.size();i++)
		{
			//if ncm[i] is also on reverse strand
			if((!(ncm[i].x2 < mi.x1) && !(min(ncm[i].y1,ncm[i].y2)>mi.y1)) && (!(ncm[i].x1>mi.x2)) && (!(max(ncm[i].y1,ncm[i].y2)<mi.y2))) //ncm mum does not fall outside
			{
				smum.push_back(ncm[i]);
			}
		}
		if(smum.size()>0)
		{
			tempmi = findClosest(mi,smum); //find the closest mum from the ncm pool
			cm.push_back(tempmi);
			mi.x1 = tempmi.x2+1; //adjust the gap coordinates
			mi.y1 = min(tempmi.y1,tempmi.y2)-1; //adjust the gap coordinates
			gapCloserRev(mi,smum,cm);//need to make it dependent on the size of ncm if ncm size does not change, that mean no more solution is there
		}
	}
	else
{
		//return;
	}	
}
/////////////////////////////////////////////////	
int nearestInt(double d)//returns the nearest integer
{	
	int in;
	in = int(d);
	if(in +0.5 >d) //if d was less than in.5
	{
		return in;
	}
	else
	{
		return in+1;
	}
}
/////////////////////////////////////////////	
void xtracTrans(map<int,vq> & mRef,vector<mI> & cm, ofstream & ftest)
{
	//vector<mI> qcm = cm;
	int k = cm.size()-1;//to use for traversing the vector in reverse direction
	//sort(qcm.begin(),qcm.end(),qusort); //sorted cm based on query coordinates
	mI gapmi;
	for(unsigned int i= 1;i<cm.size()-1;i++)
	{
		gapmi.rn = cm[i-1].rn;
		gapmi.qn = cm[i-1].qn;
		gapmi.x1 = cm[i-1].x2;
		gapmi.x2 = cm[i].x1;
                gapmi.y1 = min(cm[i-1].y2,cm[i].y1);
                gapmi.y2 = max(cm[i-1].y2,cm[i].y1);

		if((cm[i].y2 < cm[i-1].y1) && (cm[i].y1 <cm[i].y2) && (cm[i-1].y1 < cm[i-1].y2)) //if out of order and neither are reverse oriented
		{
			if(!(!(cm[i].x1 < cm[i-1].x1) && !(cm[i].x2 > cm[i-1].x2)))
			{
				ftest<<cm[i].rn<<"\t"<<cm[i].x1<<"\t"<<cm[i].x2<<"\tTRANSLOC\t"<<cm[i].qn<<"\t"<<cm[i].y1<<"\t"<<cm[i].y2<<endl;
			}
			if(!(cm[i].x1 < cm[i-1].x1) && !(cm[i].x2 > cm[i-1].x2))
			{
				cm[i] = cm[i-1];
			}	
		
		}
		if((cm[k-i].y1 > cm[k-i+1].y2) && (cm[k-i+1].y1 < cm[k-i+1].y2) && (cm[k-i].y1 < cm[k-i].y2)) //if out of order and neither are reverse oriented
		{
			if(!(!(cm[k-i+1].x1 < cm[k-i].x1) && !(cm[k-i+1].x2 > cm[k-i].x2)))
			{
				ftest<<cm[k-i].rn<<"\t"<<cm[k-i].x1<<"\t"<<cm[k-i].x2<<"\tTRANSLOC\t"<<cm[k-i].qn<<"\t"<<cm[k-i].y1<<"\t"<<cm[k-i].y2<<endl;
			}
		}
	}
}
