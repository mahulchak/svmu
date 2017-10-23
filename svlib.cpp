#include<iostream>
#include "sv.h"
using namespace std;
using chroms = map<string,chromPair>;
using ccov = vector<int>;
using vq = vector<qord>;

/////////////////////////////////////////////////////////
bool qusort(mI mi1, mI mi2)
{
	return (min(mi1.y1,mi1.y2) < min(mi2.y1,mi2.y2)) ||((min(mi1.y1,mi1.y2) == min(mi2.y1,mi2.y2)) && (max(mi1.y1,mi1.y2)<max(mi2.y1,mi2.y2)));
}
ccov makeChromBucket(int refLen)
{
	ccov v;
	for(int i=0;i<refLen;i++)
	{
		v.push_back(0);
	}
return v;
}	
///////////////////////////////////////////////////////////
int findDist(int & x1,int & y1, int & c) //c is the intercept of the absolute diagonal
{
	int c1 = 0, dist =0;
	
	c1 = abs(x1 - y1);
	
	dist = abs(c1 -c);
	dist = int(dist/sqrt(2));

	return dist;
}
///////////////////////////////////////////////////////////
//bool chkOverlap(int & x1, int & x2,int & y1, int &y2) 
bool detectShadow(mI & mum, vector<mI> & mums, unsigned int i) //return whether mum is shadow mum or not
{
	unsigned int count = 0;
	bool found = false;
	while((found == false) && (count<mums.size())) // end of mum reference cannot be beyond mums[i] start if mum is a shadow
	{
		if((!(mum.x1<mums[count].x1)) && (!(mum.x2>mums[count].x2)) && (count != i))
		{
			found =true;
		}
		if((!(mum.y1 < mums[count].y1)) && (!(mum.y2 > mums[count].y2)) && (count != i))
		{
			found = true;
		}
	count++;
	}

return found;
}
/////////////////////////////////////////////////////////
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
void storeCords(map<int,vq> & mRef, mI & mi)
{	
	int refC = mi.x1;
	int ci = refC * (-1); //ci is minus i
	qord temp;
	if(mi.y1 < mi.y2 ) //both are on the same strand
	{
		int qC = mi.y1;
		while( refC<mi.x2+1)
		{
			if((find(mi.mv.begin(),mi.mv.end(),refC) == mi.mv.end()) && (find(mi.mv.begin(),mi.mv.end(),ci) == mi.mv.end())) //if this position does not have a indel
			{
				refC++;
				qC++;
				temp.name = mi.qn;
				temp.cord = qC-1;
				mRef[refC-1].push_back(temp);
			}
			if((find(mi.mv.begin(),mi.mv.end(),refC) != mi.mv.end()) && (find(mi.mv.begin(),mi.mv.end(),ci) == mi.mv.end())) //position has insertion
			{
				refC++;
				temp.name = mi.qn;
				temp.cord = qC-1;
				mRef[refC-1].push_back(temp);
			}
			if((find(mi.mv.begin(),mi.mv.end(),refC) == mi.mv.end()) && (find(mi.mv.begin(),mi.mv.end(),ci) != mi.mv.end())) //position has deletion
			{
				qC++;
			}
		}
	}
	if(mi.y1 > mi.y2 )//if two are on different strands
	{
		int qC = mi.y1; //y1 is bigger than y2
		while(refC<mi.x2+1)
		{
			if((find(mi.mv.begin(),mi.mv.end(),refC) == mi.mv.end()) && (find(mi.mv.begin(),mi.mv.end(),ci) == mi.mv.end())) //if this position does not have a indel
			{
				refC++;
				qC--;
				temp.name = mi.qn;
				temp.cord = qC+1;	
				mRef[refC-1].push_back(temp);
			}
			if((find(mi.mv.begin(),mi.mv.end(),refC) != mi.mv.end()) && (find(mi.mv.begin(),mi.mv.end(),ci) == mi.mv.end())) //position has insertion
			{
				refC++;
				temp.name = mi.qn;
				temp.cord = qC+1;
				mRef[refC-1].push_back(temp);				
			}
			if((find(mi.mv.begin(),mi.mv.end(),refC) == mi.mv.end()) && (find(mi.mv.begin(),mi.mv.end(),ci) != mi.mv.end())) //position has deletion
			{
				qC--;
			}
		}
	}
}
		
//////////////////////////////////////////////////////////
vector<double> getCoverage(mI & mi, ccov & masterRef, ccov & masterQ)
{
	int d = 0, cov = 0;
	double c;
	vector<double> cc;
//cout<<mi.x1<<"\t"<<mi.x2<<"\t"<<mi.y1<<"\t"<<mi.y2<<"\t";
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
//cout<<c<<endl;
}
		
		
		
//////////////////////////////////////////////////////
mI findClosest(mI & mi, vector<mI> & mums, unsigned int i,ccov & masterRef, ccov & masterQ)
{
	map<double,mI> storeDist;
	map<double,mI>::iterator it;
	double d1 =0, d2 =0, d = 0;
	int ty1 =0, ty2 =0; //to switch the coordinates of the inverted MUMs
	vector<double> vd;
	ty1 = mi.y2;
	
	if(mi.y1 > mi.y2) // if reverse oriented
	{
		ty1 = mi.y1;
	}

	for(unsigned int j=i+1; j <mums.size();j++)
	{
		vd = getCoverage(mums[j],masterRef,masterQ);
		if(mums[j].y1 > mums[j].y2) //on the other strand
		{
			ty2 = mums[j].y2; //swap the values
			d1 = pow(abs(mi.x2 - mums[j].x1),2);
			d2 = pow(abs(ty1 - ty2),2);
              		d = sqrt(d1+d2);
		}
		if(mums[j].y1 < mums[j].y2)
		{
			d1 = pow(abs(mi.x2 -mums[j].x1),2);
			d2 = pow(abs(ty1 - mums[j].y1),2);
			d = sqrt(d1+d2);
		}
//cout << mi.x2<<"\t"<<mi.y2<<"\t"<<mums[j].x1<<"\t"<< mums[j].y1<<endl;
		storeDist[d] = mums[j];
		
	}
	
	it = storeDist.begin();
	 		
return it->second;
}
////////////////////////////////////////////////////
mI findClosest(mI & mi, vector<mI> & mums)
{
	map<double,mI> storeDist;
	map<double,mI>::iterator it;
	double d1 =0, d2 =0, d = 0;
	int ty1 =0, ty2 =0; //to switch the coordinates of the inverted MUMs
	ty1 = mi.y1;
	if(mi.y1 > mi.y2) // if reverse oriented
	{
		ty1 = mi.y1;
	}
	for(unsigned int j=0; j <mums.size();j++)
	{
		if(mums[j].y1 > mums[j].y2) //on the other strand
		{
			ty2 = mums[j].y2; //swap the values
			d1 = pow(abs(mi.x1 - mums[j].x1),2);
			d2 = pow(abs(ty1 - ty2),2);
			 d = sqrt(d1+d2);
		}
		if(mums[j].y1 < mums[j].y2)
		{
			d1 = pow(abs(mi.x1 -mums[j].x1),2);
			d2 = pow(abs(ty1 - mums[j].y1),2);
			d = sqrt(d1+d2);
		}
		storeDist[d] = mums[j];
	}
	it = storeDist.begin();
return it->second;
}
////////////////////////////////////////////////////
//void splitByCoverage(chromPair & cp, ccov & chrom,vector<mI> & mums, ccov & masterRef, ccov & masterQ) // returns the percentage of gap filled by the mi in mums
void splitByCoverage(chromPair & cp, ccov & chrom, ccov & masterQ,ofstream & findel)
{
	int cov=0, lastcov=0, nextcov =0;
	//vector<mI> mum;
	mI mi,tempmi,gapmi;
	vector<double> vd;
	mi.x1 =1;
	mi.x1 = cp.cm[0].x1;
	for(unsigned int j=0;j<cp.cm.size();j++)
	{
		//for(unsigned int i =1; i<chrom.size()-1;i++)
		mi.x1 = cp.cm[j].x1;
		if(j>0)
		{
			gapmi.x1 = min(cp.cm[j-1].x2,cp.cm[j].x1);
			gapmi.x2 = max(cp.cm[j-1].x2,cp.cm[j].x1);
			gapmi.y1 = min(cp.cm[j-1].y2, cp.cm[j].y1);
			gapmi.y2 = max(cp.cm[j-1].y2, cp.cm[j].y1);
			vd = getCoverage(gapmi,chrom,masterQ);
			if(vd[0] < 0.1)
			{
				findel<<cp.cm[j].rn<<"\t"<<gapmi.x1<<"\t"<<gapmi.x2<<"\tREF_INS\t"<<cp.cm[j].qn<<"\t"<<gapmi.y1<<"\t"<<gapmi.y2<<endl;
			}
			if(vd[1]<0.1)
			{
				findel<<cp.cm[j].qn<<"\t"<<gapmi.y1<<"\t"<<gapmi.y2<<"\tQUERY_INS\t"<<cp.cm[j].rn<<"\t"<<gapmi.x1<<"\t"<<gapmi.x2<<endl;
			}
		for(int i = cp.cm[j].x1-1; i<cp.cm[j].x2;i++)
		{
			cov = chrom[i];
			lastcov = chrom[i-1];
			nextcov = chrom[i+1];
			if((cov != lastcov) && (cov == nextcov))
			{
				mi.x1 = i+1;			
			}
			if((cov == lastcov) && (cov != nextcov))
			{
				mi.x2 = i+1;
			mi.rn = cp.cm[0].rn;
			mi.qn = cp.cm[0].qn;
				if(chrom[mi.x1-1] >1)
				//if(chrom[mi.x1-1] > 0) //this is to count for those which are fewer in query than reference
				{
					if((cp.cc.size() == 0) && (mi.x2 -mi.x1 >20)) //at least 20 bp or more should show cnv
					{
						cp.cc.push_back(mi);
//cout<<mi.rn<<"\t"<<mi.x1<<"\t"<<mi.x2<<"\t"<<chrom[mi.x1-1]<<"\t"<<chrom[mi.x2-1]<<endl;
					}
					if((cp.cc.size() >0) && !(mi == cp.cc[cp.cc.size()-1]) && (mi.x2-mi.x1>20))
					{
						cp.cc.push_back(mi);
					}
//cout<<mi.rn<<"\t"<<mi.x1<<"\t"<<mi.x2<<"\t"<<chrom[mi.x1-1]<<"\t"<<chrom[mi.x2-1]<<endl;
				}
				//if(chrom[mi.x1-1] ==0) 
				//{
				//	cp.in.push_back(mi);
//cout<<mi.rn<<"\t"<<mi.x1<<"\t"<<mi.x2<<"\t"<<chrom[mi.x1-1]<<endl;
				}
						
//cout<<mi.x1<<"\t"<<mi.x2<<"\t"<<mi.y1<<"\t"<<mi.y2<<"\t"<<vd[0]<<"\t"<<vd[1]<<endl;
			}
//cout<<i<<"\t"<<mi.x1<<"\t"<<mi.x2<<"\t"<<chrom[mi.x1-1]<<"\t"<<chrom[mi.x2-1]<<endl;//because coverage is 0 based but coordinate is 1 based
		}
	}
	
//return mum;	
}						
/////////////////////////////////////////////////
void gapCloser(mI & mi, vector<mI> ncm, vector<mI> & cm)
{
	mI tempmi;\
	vector<mI> smum; //selected mums that overlap with the gap.
	if((mi.x2 - mi.x1>0) && (mi.y2 - mi.y1 >0)) //checking if both of them has gaps. needs to do this for inverted sequences
	{
		for(unsigned int i = 0;i<ncm.size();i++)
		{
			if((!(ncm[i].x2 < mi.x1) && !(max(ncm[i].y1,ncm[i].y2)<mi.y1)) && (!(ncm[i].x1>mi.x2) && !(min(ncm[i].y1,ncm[i].y2)>mi.y2))) //ncm mum does not fall outside
			{
				smum.push_back(ncm[i]);
//cout<<ncm.size()<<"\t"<<i<<"\t"<<mi.rn<<"\t"<<mi.qn<<mi.x1<<"\t"<<mi.x2<<"\t"<<mi.y1<<"\t"<<mi.y2<<"\t"<<ncm[i].x1<<"\t"<<ncm[i].x2<<"\t"<<ncm[i].y1<<"\t"<<ncm[i].y2<<endl;
			}
		}
		if(smum.size()>0)
		{	
			tempmi = findClosest(mi,smum); //find the closest mum from the ncm pool
			cm.push_back(tempmi);
//cout<<tempmi.rn<<" "<<tempmi.x1<<" "<<tempmi.x2<<" "<<tempmi.y1<<" "<<tempmi.y2<<endl;
			mi.x1 = tempmi.x2+1; //adjust the gap coordinates
			mi.y1 = max(tempmi.y1,tempmi.y2)+1; //adjust the gap coordinates
			gapCloser(mi,smum,cm);//need to make it dependent on the size of ncm if ncm size does not change, that mean no more solution is there
			
		}
		//return;
	}
	else
	{
		//return;
	}
	
}
/////////////////////////////////////////////////	
vector<mI> findQuery(map<int,vq> & mRef, mI & mi,ccov & masterRef, ccov & masterQ)
{
	
	//vector<mI> mums(mRef[mi.x1].size());//creating the vector of the coverage size
	vector<string> qnames;//will be used to screen TEs based on number of contigs they are present on
	qord temp;
	vector<double> vd;
	//vd = getCoverage(mi,masterRef,masterQ);
	vector<mI> mums(masterRef[mi.x1-1]);
	
	int rcov =0,cov1 =0;
	vector<int> vi;
	bool found =false;
	sort(mRef[mi.x1].begin(),mRef[mi.x1].end());
	sort(mRef[mi.x2].begin(),mRef[mi.x2].end());
	//for(unsigned int j=0; j<mRef[mi.x1].size();j++)
	for(unsigned int j=0; j<masterRef[mi.x1-1];j++)
	{
		mums[j].x1 = mi.x1;
		mums[j].x2 = mi.x2;
		mums[j].y1 = mRef[mi.x1][j].cord;	
		mums[j].y2 = mRef[mi.x2][j].cord;
		mums[j].qn = mRef[mi.x1][j].name;
		mums[j].rn = mi.rn;
		if(j ==0)
		{
			qnames.push_back(mums[j].qn);
		}
		if(find(qnames.begin(),qnames.end(),mums[j].qn) == qnames.end() && (j>0)) //if the qname hasn't already been entered
		{
			qnames.push_back(mums[j].qn);
		}
	}
	for(unsigned int i = 0; i< mums.size();i++)
	{
		vd = getCoverage(mums[i],masterRef,masterQ); //add these in the function arguments
		rcov = nearestInt(vd[0]);
		cov1 = nearestInt(vd[1]);	
		//if(rcov != cov1) //if they are unequal. counts both less and more copies
		if(rcov > cov1) //counts only copies which are more
		{
			found = true;
		}
		if((qnames.size() > 1)) //present in more than 1 chromosomes/contigs
		{
			mums[i].qn = mums[i].qn + " trans";
		}	
	}
	if(found == false)
	{
		mums.clear();
	}	
	qnames.clear();
return mums;
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
void xtracTrans(vector<mI> & cm, ofstream & ftest)
{
	//vector<mI> qcm = cm;
	int k = cm.size()-1;//to use for traversing the vector in reverse direction
	//sort(qcm.begin(),qcm.end(),qusort); //sorted cm based on query coordinates
	for(unsigned int i= 1;i<cm.size()-1;i++)
	{
		if((cm[i].y2 < cm[i-1].y1) && (cm[i].y1 <cm[i].y2) && (cm[i-1].y1 < cm[i-1].y2)) //if out of order and neither are reverse oriented
		{
			if(!(!(cm[i].x1 < cm[i-1].x1) && !(cm[i].x2 > cm[i-1].x2)))
			{
				ftest<<"1 "<<cm[i].rn<<"\t"<<cm[i].x1<<"\t"<<cm[i].x2<<"\t"<<cm[i].qn<<"\t"<<cm[i].y1<<"\t"<<cm[i].y2<<endl;
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
				ftest<<"2 "<<cm[k-i].rn<<"\t"<<cm[k-i].x1<<"\t"<<cm[k-i].x2<<"\t"<<cm[k-i].qn<<"\t"<<cm[k-i].y1<<"\t"<<cm[k-i].y2<<endl;
			}
		}
	}
}
