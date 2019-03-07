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
	return (max(mi1.y1,mi1.y2) < max(mi2.y1,mi2.y2)) ||((max(mi1.y1,mi1.y2) == max(mi2.y1,mi2.y2)) && (min(mi1.y1,mi1.y2)>min(mi2.y1,mi2.y2)));

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
bool lsort(mI m1,mI m2)
{
	return (m1.l>m2.l) || ((m1.l == m2.l) && (m1.x1<m2.x1));
}	
////////////////////////////////////////////////////////////
void findInnie(vector<mI> & mums,mI & mi)
{
	int i = 0;
	while((mi.x2 > mums[i].x1) && (i<mums.size()))
	{
//cout<<"debug\t"<<mi.rn<<'\t'<<mi.x1<<'\t'<<mi.x2<<'\t'<<mums[i].rn<<'\t'<<mums[i].x1<<'\t'<<mums[i].x2<<endl;
		if((mi.x1 > (mums[i].x1-1)) && (mi.x2 < (mums[i].x2+1)))
		{
//cout<<"debug\t"<<mi.rn<<'\t'<<mi.x1<<'\t'<<mi.x2<<'\t'<<mums[i].rn<<'\t'<<mums[i].x1<<'\t'<<mums[i].x2<<endl;
			if(!(mi == mums[i]))
			{
//cout<<"debug\t"<<mi.rn<<'\t'<<mi.x1<<'\t'<<mi.x2<<'\t'<<mums[i].rn<<'\t'<<mums[i].x1<<'\t'<<mums[i].x2<<endl;
				if(mi.c == 'q')
				{
					mi.c = 'd';
					break;
				}
				else
				{
					mi.c = 'r';
				}
			}
		}
		if((min(mi.y1,mi.y2) > (min(mums[i].y1,mums[i].y2)-1)) && (max(mi.y2,mi.y1) <(max(mums[i].y2,mums[i].y1)+1)))
		{
			if(!(mi == mums[i]))
			{
				if(mi.c == 'r')
				{
					mi.c = 'd';
					break;
				}
				else
				{
					mi.c = 'q';
				}
			}
		}
		++i;
	}
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
void storeNameCount(ccov & chromDensityRef,ccov & chromDensityQ,map<string,int> & lookUpRef,map<string,int> & lookUpQ,mI & mi)
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
		if((chromDensityRef[i] >0) && (chromDensityRef[i] != lookUpQ[mi.qn]))//if the existing position does not have the same chrom and has not already been mapped to >2 chroms
		{
			chromDensityRef[i] = -2;
		}
		if(chromDensityRef[i] == 0)//no chrom name has been recorded yet
		{
			chromDensityRef[i] = lookUpQ[mi.qn];
		}
			
	}

	for(int i = ty1-1; i<ty2;i++)
	{
		if((chromDensityQ[i] >0) && (chromDensityQ[i] != lookUpRef[mi.rn]))//if the existing position does not have the same chrom and has not already been mapped to >2 chroms
		{
			chromDensityQ[i] = -2;

		}
		if(chromDensityQ[i] == 0)//no chrom name has been recorded yet
		{
			chromDensityQ[i] = lookUpRef[mi.rn];
		}
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
	d = max(abs(mi.x2 - mi.x1),1);//to avoid using 0
	for(int i = mi.x1-1;i<mi.x2;i++)
	{
		cov = cov + masterRef[i];
	}
	c = cov/double(d);
	cc.push_back(c);
	cov = 0;
	d = 0;
	d = max(abs(mi.y1-mi.y2),1);//to avoid using 0
	for(int i = (min(mi.y1,mi.y2)-1);i<max(mi.y1,mi.y2);i++)
	{
		cov = cov + masterQ[i];
	}
	c = cov/double(d);
	cc.push_back(c);
	return cc;
}
////////////////////////////////////////////////////////////
vector<double> getChromCount(mI & mi, ccov & chromDensityRef, ccov & chromDensityQ)
{
	vector<double> cc;
	double c;
	int d = 0, cov = 0;
	d = max(abs(mi.x1 - mi.x2),1);
	for(int i = mi.x1-1;i<mi.x2;i++)
	{
		if(chromDensityRef[i] == -2)
		{
                	++cov;
		}
	}
	c = double(cov)/double(d);
	cc.push_back(c);
	cov = 0;
	d = max(abs(mi.y1-mi.y2),1);
	for(int i = min(mi.y1,mi.y2)-1;i<max(mi.y1,mi.y2);i++)
	{
		if(chromDensityQ[i] == -2)
		{
			++cov;
		}
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
	double d1 =0, d2 =0, d = 0,Dist= 0,rd2=0,rd=0;
	mI invmi;//store the swapped gapmi
	sort(mums.begin(),mums.end(),lsort);//sort the mums by length
	//Dist = sqrt(pow(abs(mi.x2-mi.x1),2)+pow(abs(mi.y1-mi.y2),2));
	invmi = mi;
	invmi.y1 = mi.y2;//invmi is reverse complement of gapmi
	invmi.y2 = mi.y1;
	for(unsigned int j=0; j <mums.size();j++)
	{
		if((mums[j].rn == mi.rn) && (mums[j].qn == mi.qn))
		{
			d1 = pow(abs(mi.x1 - max(mums[j].x1,mi.x1)),2);//if start of mum preceeds the gap start then effective mum start is the gap start. will do it for query too.
			d2 = pow(abs(mi.y1 - mums[j].y1),2);
			d = abs(sqrt(d1+d2));
			rd2 = pow(abs(invmi.y1 - mums[j].y1),2);
			rd = abs(sqrt(d1 + rd2));
			if(d < rd) // if forward orientation is closer
			{
				storeDist[d] = mums[j];
			}
			else
			{
				storeDist[rd] = mums[j];
			}
		}
//cout<<"dist\t"<<d<<'\t'<<rd<<'\t'<<mums[j].rn<<'\t'<<mums[j].x1<<'\t'<<mums[j].x2<<'\t'<<mums[j].qn<<'\t'<<mums[j].y1<<'\t'<<mums[j].y2<<endl;
	}
	it = storeDist.begin();
	if((mums.size()>1)&& (mums[0].l == mums[1].l))//if both lengths are same, get the closest
	{
		return it->second;
	}
	else
	{
		return mums[0];//return the longest mum
	}
}
////////////////////////////////////////////////////
void gapCloser(mI & mi, vector<mI> ncm, vector<mI> & cm)
{
	mI tempmi;
	mI gapRight,gapLeft;//two sides of the new split gap
	vector<mI> smum; //selected mums that overlap with the gap.
	int refOvl =0, qOvl =0;
	double refProp, qProp;
	if((mi.x2 > mi.x1) && (mi.y2 > mi.y1)) //gaps in both ref and query exist, and in forward strand. inverted coordinates are converted to forward coordinates. so gaps are always in forward direction
	{
//cout<<"GAPS\t"<<mi.rn<<'\t'<<mi.x1<<'\t'<<mi.x2<<'\t'<<mi.qn<<'\t'<<mi.y1<<'\t'<<mi.y2<<endl;
		for(unsigned int i = 0;i<ncm.size();i++)
		{
			//if((!(ncm[i].x2 < mi.x1)) && (!(max(ncm[i].y1,ncm[i].y2)<min(mi.y1,mi.y2))) && (!(ncm[i].x1>mi.x2)) && (!(min(ncm[i].y1,ncm[i].y2)>max(mi.y2,mi.y1)))) //ncm mum does not fall outside
			if((!(ncm[i].x2 < mi.x1)) && (!(ncm[i].x1 > mi.x2)) && (mi.rn == ncm[i].rn) && (mi.qn == ncm[i].qn))
			{
				if((!(max(ncm[i].y2,ncm[i].y1) < mi.y1)) && (!(min(ncm[i].y1,ncm[i].y2) > mi.y2)))
				{
					refOvl = min(ncm[i].x2,mi.x2) - max(ncm[i].x1,mi.x1);
					refProp = double(refOvl)/double(ncm[i].x2-ncm[i].x1);
					qOvl = min(max(ncm[i].y1,ncm[i].y2),mi.y2) - max(min(ncm[i].y1,ncm[i].y2),mi.y1);
					qProp = double(qOvl)/double(abs(ncm[i].y2-ncm[i].y1));
					if((refProp>0.5) && (qProp>0.5))
					{				
						smum.push_back(ncm[i]);
//cout<<mi.rn<<'\t'<<mi.x1<<'\t'<<mi.x2<<'\t'<<mi.qn<<'\t'<<mi.y1<<'\t'<<mi.y2<<'\t'<<ncm[i].rn<<'\t'<<ncm[i].x1<<'\t'<<ncm[i].x2<<'\t'<<ncm[i].qn<<'\t'<<'\t'<<ncm[i].y1<<'\t'<<ncm[i].y2<<endl;
					}
				}
			}
		}
		if(smum.size()>0)
		{	
			tempmi = findClosest(mi,smum); //find the closest mum from the ncm pool
			if(tempmi.y1 != 0) // add the condition that when it is forward this happens
			{
				cm.push_back(tempmi);
				gapRight.rn = mi.rn;
				gapRight.qn = mi.qn;
				gapRight.x1 = min(tempmi.x2+1,mi.x2); //adjust the gap coordinates
				gapRight.x2 = mi.x2;
				gapLeft.x1 = mi.x1;
				gapLeft.x2 = max(tempmi.x1,mi.x1);
				gapLeft.rn = mi.rn;
				gapLeft.qn = mi.qn;

				if(tempmi.y1 < tempmi.y2)//forward oriented
				{
					gapRight.y1 = min(max(tempmi.y2+1,tempmi.y1),mi.y2); //adjust the gap coordinates
					gapRight.y2 = mi.y2;
					gapLeft.y1 = mi.y1;
					gapLeft.y2 = max(min(tempmi.y1,tempmi.y2+1),mi.y1);
				}
				if(tempmi.y1 > tempmi.y2)//reverse oriented
				{
					gapRight.y1 = min(max(tempmi.y1 +1,tempmi.y2),mi.y2);
					gapRight.y2 = mi.y2;
					gapLeft.y1 = mi.y1;
					gapLeft.y2 = max(min(tempmi.y2,tempmi.y1+1),mi.y1);	
				}
//cout<<"tempmi\t"<<tempmi.rn<<'\t'<<tempmi.x1<<'\t'<<tempmi.x2<<'\t'<<tempmi.qn<<'\t'<<tempmi.y1<<'\t'<<tempmi.y2<<endl;
//cout<<"Right\t"<<gapRight.rn<<'\t'<<gapRight.x1<<'\t'<<gapRight.x2<<'\t'<<gapRight.qn<<'\t'<<gapRight.y1<<'\t'<<gapRight.y2<<endl;
//cout<<"Left\t"<<gapLeft.rn<<'\t'<<gapLeft.x1<<'\t'<<gapLeft.x2<<'\t'<<gapLeft.qn<<'\t'<<gapLeft.y1<<'\t'<<gapLeft.y2<<endl;
				gapCloser(gapRight,smum,cm);//need to make it dependent on the size of ncm if ncm size does not change, that mean no more solution is there
				gapCloser(gapLeft,smum,cm);
			}
					
		}
	}
		
}
/////////////////////////////////////////////////////
mI returnMumByQ1(int & y1,vector<mI> & mums)//returns the mum from query sorted mums which has the same y1 as y1
{
	unsigned int i =0;
	unsigned int k = mums.size()-1;
	mI tempmi;
	while((mums[i].y1 != y1) && (mums[k-i].y1 != y1) && (i<mums.size()))
	{
		++i;
	}
	if(mums[mums.size()-1].y1 == y1)
	{
		tempmi = mums[mums.size()-1];//there is no number after this so send it
	}	
	else
	{
		if(mums[i].y1 == y1)
		{
			tempmi = mums[i+1];
		}
		if(mums[k-i].y1 == y1)
                {
			tempmi = mums[(k-i)+1];
                }

	}
	return tempmi;
}
///////////////////////////////////////////////////////////
mI returnMumByQ2(int & y1,vector<mI> & mums)//returns the mum from query sorted mums which has the same y1 as y1
{
	unsigned int i = 0;
	unsigned int k = mums.size()-1;
	mI tempmi;
	while((mums[i].y1 != y1) && (mums[k-i].y1 != y1) && (i<mums.size()))
	{
		++i;
	}
	if(mums[0].y1 == y1) //if the first element
	{
		tempmi = mums[0];//there is no element before this so send the same
	}

	else
	{
		if(mums[i].y1 == y1)
		{
			tempmi = mums[i-1];
		}
		if(mums[k-i].y1 == y1)
		{
			tempmi = mums[(k-i)-1];
		}
	}
	return tempmi;
}
///////////////////////////////////////////////////////
void gapCloserRev(mI & mi, vector<mI> ncm, vector<mI> & cm)
{
	mI tempmi;
	vector<mI> smum; //selected mums that overlap with the gap.
	if((mi.x2 > mi.x1) && (mi.y1 > mi.y2)) //checking if both of them has gaps. this is for inverted sequences
	{
		for(unsigned int i = 0;i<ncm.size();i++)
		{
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
