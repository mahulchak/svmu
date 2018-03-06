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
void storeCords(map<int,vq> & mRef, mI & mi,ofstream & fout)
{	
	int refC = mi.x1;
	int ci = 0; //assuming that the first base in a mum will not have an insertion in query
	qord temp;
	vector<int>:: iterator it;
	if(mi.y1 < mi.y2 ) //both are on the same strand
	{
		int qC = mi.y1;
		while( refC<mi.x2+1)
		{
			if((find(mi.mv.begin(),mi.mv.end(),refC) == mi.mv.end()) && (find(mi.mv.begin(),mi.mv.end(),ci) == mi.mv.end())) //if this position does not have a indel
			{
				temp.name = mi.qn;
				temp.cord = qC;
				mRef[refC-1].push_back(temp);
				refC++;
				qC++;
				ci = refC * (-1);
				if(find(mi.mv.begin(),mi.mv.end(),ci) != mi.mv.end())//if the next is a del
				{
					temp.name = mi.qn;
					temp.cord = qC;
					mRef[refC-1].push_back(temp);
				}
//fout<<mi.rn<<"\t"<<refC<<"\t"<<mi.qn<<"\t"<<qC<<endl;
			}
			if((find(mi.mv.begin(),mi.mv.end(),refC) != mi.mv.end())) //position has insertion
			{
				temp.name = mi.qn;
				temp.cord = qC-1;//because when insertion begins, query cord does not increase
				mRef[refC-1].push_back(temp);
				refC++;
				ci = refC * (-1);
//fout<<mi.rn<<"\t"<<refC<<"\t"<<mi.qn<<"\t"<<qC<<endl;
			}
			if((find(mi.mv.begin(),mi.mv.end(),ci) != mi.mv.end())) //position has deletion
			{
				qC++;
				it = find(mi.mv.begin(),mi.mv.end(),ci);
				it++;
				--ci;
				while((*it == -1) && (it != mi.mv.end()))
				{
					qC++;
					it++;
				}
				if((*it != -1) || (it == mi.mv.end()))
				{
					refC++;
					qC++;
					ci = refC * (-1);
				}
			}
//fout<<mi.rn<<"\t"<<refC<<"\t"<<mi.qn<<"\t"<<qC<<endl;
		}
	}
	if(mi.y1 > mi.y2 )//if two are on different strands
	{
		int qC = mi.y1; //y1 is bigger than y2
		while(refC<mi.x2+1)
		{
			if((find(mi.mv.begin(),mi.mv.end(),refC) == mi.mv.end()) && (find(mi.mv.begin(),mi.mv.end(),ci) == mi.mv.end())) //if this position does not have a indel
			{
				temp.name = mi.qn;
				temp.cord = qC;	
				mRef[refC-1].push_back(temp);
				refC++;
				qC--;
				ci = refC * (-1);
				if(find(mi.mv.begin(),mi.mv.end(),ci) != mi.mv.end())
				{
					temp.name = mi.qn;
					temp.cord = qC;
					mRef[refC-1].push_back(temp);
				}
//fout<<mi.rn<<"\t"<<refC<<"\t"<<mi.qn<<"\t"<<qC<<endl;
			}
			if((find(mi.mv.begin(),mi.mv.end(),refC) != mi.mv.end())) //position has insertion
			{
				temp.name = mi.qn;
				temp.cord = qC+1;
				mRef[refC-1].push_back(temp);				
				refC++;
				ci = refC * (-1);
//fout<<mi.rn<<"\t"<<refC<<"\t"<<mi.qn<<"\t"<<qC<<endl;
			}
			if((find(mi.mv.begin(),mi.mv.end(),ci) != mi.mv.end())) //position has deletion
			{
				--qC;
				it = find(mi.mv.begin(),mi.mv.end(),ci);
				it++;
				--ci;
				while((*it == -1) && (it != mi.mv.end()))
				{
					--qC;
					it++;
				}
				if((*it != -1) || (it == mi.mv.end()))//stretch of -1 has ended
				{
					refC++;
					--qC;
					ci = refC * (-1);
				}
//fout<<mi.rn<<"\t"<<refC<<"\t"<<mi.qn<<"\t"<<qC<<endl;
			}
		}
	}
}
		
//////////////////////////////////////////////////////////
vector<double> getCoverage(mI & mi, ccov & masterRef, ccov & masterQ, float p)
{
	int d = 0, cov = 0,covCount=0, medCov =0;
	double c;
	vector<double> cc;
	map<int,int> covFreq;//holds coverage frequency for the genomic interval
//cout<<mi.x1<<"\t"<<mi.x2<<"\t"<<mi.y1<<"\t"<<mi.y2<<"\t";
	d = mi.x2 - mi.x1;
	for(int i = mi.x1-1;i<mi.x2;i++)
	{
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
//cout<<mi.x1<<"\t"<<mi.x2<<"\t"<<mi.y1<<"\t"<<mi.y2<<"\t"<<covCount<<"\t"<<medCov<<"\t";
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
//cout<<medCov<<endl;
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
	double d1 =0, d2 =0, d = 0,Dist= 0;
	int ty1 =0, ty2 =0; //to switch the coordinates of the inverted MUMs
	ty1 = mi.y1;
	Dist = sqrt(pow(abs(mi.x2-mi.x1),2)+pow(abs(mi.y1-mi.y2),2));
	//Dist = abs(mi.x2 - mi.x1);
	if(mi.y1 > mi.y2) // if reverse oriented
	{
		ty1 = mi.y1;
	}
	for(unsigned int j=0; j <mums.size();j++)
	{
		if(mums[j].y1 > mums[j].y2) //on the other strand
		{
			//if(mi.c != 'r')//if both ends are not inverted
			//{
				//ty2 = mums[j].y2; //swap the values
			//	ty2 = mums[j].y1;
			//}
			//if(mi.c == 'r')
			//{
			//	ty2 = mums[j].y1;
			//}
			if(abs(mums[j].y2-mi.y1) > abs(mums[j].y1 - mi.y1))
			{
				ty2 = mums[j].y1;
			}
			if(abs(mums[j].y2-mi.y1) < abs(mums[j].y1 - mi.y1))
			{
				ty2 = mums[j].y2;
			}
			d1 = pow(abs(mi.x1 - max(mums[j].x1,mi.x1)),2);//if start of mum preceeds the gap start then effective mum start is the gap start. will do it for query too.
			d2 = pow(abs(ty1 - ty2),2);
			 d = abs(sqrt(d1+d2));
		}
		if(mums[j].y1 < mums[j].y2)
		{
			d1 = pow(abs(mi.x1 -max(mums[j].x1,mi.x1)),2);
			d2 = pow(abs(ty1 - mums[j].y1),2);
			d = abs(sqrt(d1+d2));
		}
		//d = abs(mums[j].x1 - mi.x1); 
		storeDist[d] = mums[j];
	}
	it = storeDist.begin();
	//if(it->first >(2*Dist))
//	cout << mi.rn<<"\t"<<mi.x1<<"\t"<<mi.x2<<"\t"<<mi.qn<<"\t"<<mi.y1<<"\t"<<mi.y2<<"\t"<<it->second.rn<<"\t"<<it->second.x1<<"\t"<<it->second.x2<<"\t"<<it->second.y1<<"\t"<<it->second.y2<<"\t"<<it->first<<"\t"<<Dist<<"\t"<<d1<<"\t"<<d2<<endl;
	if(it->first > Dist)
	{
		it->second.y1 = 0;
	}
return it->second;
}
////////////////////////////////////////////////////
void splitByCoverageSen(chromPair & cp, map<int,vq> & mRef,ccov & chrom, ccov & masterQ)
{
	int cov=0, lastcov=0, nextcov =0;
	//vector<mI> mum;
	mI mi,tempmi,gapmi;
	vector<double> vd;
	mi.x1 =1;
		for(unsigned int i= 1;i<chrom.size()-1;i++)
		{
			cov = chrom[i];
			lastcov = chrom[i-1];
			nextcov = chrom[i+1];
			if((cov != lastcov) && (cov == nextcov))
			{
				mi.x1 = i+1;			
				while(mRef[mi.x1-1].size() != chrom[mi.x1-1])
				{
					mi.x1++;
				}
			}
			if((cov == lastcov) && (cov != nextcov))
			{
				mi.x2 = i+1;
				while(mRef[mi.x2-1].size() != chrom[mi.x2-1])
				{
					--mi.x2;
				}
				mi.rn = cp.cm[0].rn;
				mi.qn = cp.cm[0].qn;
				//if(chrom[mi.x1-1] >1)
				if(chrom[mi.x1-1] > 0)
				{
					mi.y1 = 0;
					mi.y2 = 0;
					if((cp.cc.size() == 0) && (mi.x2 -mi.x1 >20)) //at least 20 bp or more should show cnv
					{
						cp.cc.push_back(mi);
//cout<<mi.rn<<"\t"<<mi.x1<<"\t"<<mi.x2<<"\t"<<chrom[mi.x1-1]<<"\t"<<chrom[mi.x2-1]<<endl;
					}
					if((cp.cc.size() >0) && !(mi == cp.cc[cp.cc.size()-1]) && (mi.x2-mi.x1>20) && (find(cp.cc.begin(),cp.cc.end(),mi) == cp.cc.end()))
					{
						cp.cc.push_back(mi);
					}
//cout<<mi.rn<<"\t"<<mi.x1<<"\t"<<mi.x2<<"\t"<<chrom[mi.x1-1]<<"\t"<<chrom[mi.x2-1]<<endl;
				}
			}
						
//cout<<mi.x1<<"\t"<<mi.x2<<"\t"<<mi.y1<<"\t"<<mi.y2<<"\t"<<vd[0]<<"\t"<<vd[1]<<endl;
		}
		
}
///////////////////////////////////////////////////
void splitByCoverage(chromPair & cp,map<int,vq> & mRef, ccov & chrom, ccov & masterQ)
//void splitByCoverage(chromPair & cp, ccov & chrom, ccov & masterQ)
{
	int cov=0, lastcov=0, nextcov =0;
	mI mi,tempmi,gapmi;
	vector<double> vd;
	mi.x1 =1;
	mi.x1 = cp.cm[0].x1;
	for(unsigned int j=0;j<cp.cm.size();j++)
	{	mi.x1 = cp.cm[j].x1;
		for(int i = cp.cm[j].x1-1; i<cp.cm[j].x2;i++)
		{
			cov = chrom[i];
			//cov = mRef[i].size();
			lastcov = chrom[i-1];
			//lastcov = mRef[i-1].size();
			nextcov = chrom[i+1];
			//nextcov = mRef[i+1].size();
			if((cov != lastcov) && (cov == nextcov))
			{
				mi.x1 = i+1;
				while(mRef[mi.x1-1].size() != chrom[mi.x1-1])
				{
					mi.x1++;
				}
			}
			if((cov == lastcov) && (cov != nextcov))
			{
				mi.x2 = i+1;
				while(mRef[mi.x2-1].size() != chrom[mi.x2-1])
				{
					--mi.x2;
				}
				mi.rn = cp.cm[0].rn;
				mi.qn = cp.cm[0].qn;
				if(chrom[mi.x1-1] >1)
				//if(mRef[mi.x1-1].size()>1)
				{
					mi.y1 = 0;
					mi.y2 = 0;
					if((cp.cc.size() == 0) && (mi.x2 -mi.x1 >20)) //at least 20 bp or more should show cnv
					{
						cp.cc.push_back(mi);
					}
					if((cp.cc.size() >0) && !(mi == cp.cc[cp.cc.size()-1]) && (mi.x2-mi.x1>20) && (find(cp.cc.begin(),cp.cc.end(),mi) == cp.cc.end()))
					{
						cp.cc.push_back(mi);
					}
				}
			}
		}
	}
}	
						
/////////////////////////////////////////////////
void gapCloser(mI & mi, vector<mI> ncm, vector<mI> & cm)
{
	mI tempmi;\
	vector<mI> smum; //selected mums that overlap with the gap.
	//if((mi.x2 - mi.x1>0) && (mi.y2 - mi.y1 >0)) //checking if both of them has gaps. needs to do this for inverted sequences
	if(mi.x2 - mi.x1>0)
	{
		for(unsigned int i = 0;i<ncm.size();i++)
		{
			if((!(ncm[i].x2 < mi.x1)) && (!(max(ncm[i].y1,ncm[i].y2)<min(mi.y1,mi.y2))) && (!(ncm[i].x1>mi.x2)) && (!(min(ncm[i].y1,ncm[i].y2)>max(mi.y2,mi.y1)))) //ncm mum does not fall outside
			//if((!(ncm[i].x2 < mi.x1)) && (!(ncm[i].x1>mi.x2)))
			{
				smum.push_back(ncm[i]);
//cout<<ncm.size()<<"\t"<<i<<"\t"<<mi.rn<<"\t"<<mi.qn<<mi.x1<<"\t"<<mi.x2<<"\t"<<mi.y1<<"\t"<<mi.y2<<"\t"<<ncm[i].x1<<"\t"<<ncm[i].x2<<"\t"<<ncm[i].y1<<"\t"<<ncm[i].y2<<endl;
			}
		}
		if(smum.size()>0)
		{	
			tempmi = findClosest(mi,smum); //find the closest mum from the ncm pool
			if(tempmi.y1 != 0)
			{
				cm.push_back(tempmi);
//cout<<tempmi.rn<<" "<<tempmi.x1<<" "<<tempmi.x2<<" "<<tempmi.y1<<" "<<tempmi.y2<<endl;
				mi.x1 = tempmi.x2+1; //adjust the gap coordinates
				mi.y1 = max(tempmi.y1,tempmi.y2)+1; //adjust the gap coordinates
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
vector<mI> findQuery(map<int,vq> & mRef, mI & mi,ccov & masterRef, ccov & masterQ, ccov & masterHQ)
{
	
	//vector<mI> mums(mRef[mi.x1].size());//creating the vector of the coverage size
	vector<string> qnames;//will be used to screen TEs based on number of contigs they are present on
	qord temp;
	vector<double> vd;
	//vd = getCoverage(mi,masterRef,masterQ);
	vector<mI> mums(masterRef[mi.x1-1]);//the last element is for coverage
	mI tmi;//this is to add the coverage info as the last element of mums
	int rcov =0,cov1 =0, tdcov = 0;
	vector<int> vi;
	bool found =true;
//	sort(mRef[mi.x1].begin(),mRef[mi.x1].end());
//	sort(mRef[mi.x2].begin(),mRef[mi.x2].end());
	//for(unsigned int j=0; j<mRef[mi.x1].size();j++)
	for(int j=0; j<masterRef[mi.x1-1];j++)
	{
cout<<mi.x1<<"\t"<<mi.x2<<"\ttotal requested coverage\t"<<masterRef[mi.x1-1]<<"\tand\t"<<masterRef[mi.x2-1]<<"\texisting coverage\t"<<mRef[mi.x1-1].size()<<"\tand\t"<<mRef[mi.x2-1].size()<<endl;
		mums[j].x1 = mi.x1;
		mums[j].x2 = mi.x2;
		mums[j].y1 = mRef[mi.x1-1][j].cord;	
		mums[j].y2 = mRef[mi.x2-1][j].cord;
		mums[j].qn = mRef[mi.x1-1][j].name;
		mums[j].rn = mi.rn;
//cout<<mums[j].rn<<"\t"<<mums[j].x1<<"\t"<<mums[j].x2<<"\t"<<mums[j].qn<<"\t"<<mums[j].y1<<"\t"<<mums[j].y2<<endl;
		if(j ==0)
		{
			qnames.push_back(mums[j].qn);
		}
		if(find(qnames.begin(),qnames.end(),mums[j].qn) == qnames.end() && (j>0)) //if the qname hasn't already been entered
		{
			qnames.push_back(mums[j].qn);
		}
	}
	for(unsigned int i=0;i<mums.size();i++)
	{
		vd = getCoverage(mums[i],masterRef,masterQ);
		rcov = nearestInt(vd[0]);
		cov1 = nearestInt(vd[1]);
		vd = getCoverage(mums[i],masterRef,masterHQ);//being done for masterHQ and NOT masterQ
		tdcov = nearestInt(vd[1]);
		if(rcov != cov1) //if they are unequal. counts both less and more copies
	//	if(rcov > cov1) //counts only copies which are more
		{
			found = true;
			tmi.x1 = rcov;
			tmi.x2 = cov1;
			tmi.y1 = tdcov;
		}
		else
		{
			found = false;
			break;
		}
		if((qnames.size() > 1)) //present in more than 1 chromosomes/contigs
		{
			mums[i].qn = mums[i].qn + " trans";
		}	
	}
//	tmi.x1 = rcov;
//	tmi.x2 = cov1;
//	tmi.y1 = tdcov;
	mums.push_back(tmi);//add the coverage info as the last element. 
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
