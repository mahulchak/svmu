#include "sv.h"
#include<iostream>

using namespace std;
using chroms = map<string,chromPair>;
using ccov = vector<int>;
using vq = vector<qord>;

void annotGaps(vector<mI> & cm, map<int,vq> & mRef, ccov & masterRef, ccov & masterQ, vector<mI> & cnv, map<int,vector<qord> > & umRef, string & refseq, string & qseq, vector<int> & seqLen,ofstream & fout, ofstream & fsmall)
{
	int refOvl =0; //overlap between reference intervals
	int qOvl = 0; //overlap between query intervals
	vector<double> cov(2);
	vector<int> vi;
	mI gapmi,cnvmi,tempmi;
	vector<mI> storedCNV;
	tempmi = cm[0]; // because the loop does not go into this
	callSmall(tempmi,umRef,refseq,qseq,seqLen, fsmall);
	bool fs = false;//fs = found SV
	for(unsigned int i=1; i< cm.size();i++)
	{
//cout<<"cm\t"<<cm[i].rn<<"\t"<<cm[i].x1<<"\t"<<cm[i].x2<<"\t"<<cm[i].qn<<"\t"<<cm[i].y1<<"\t"<<cm[i].y2<<endl;
		tempmi = cm[i];
		callSmall(tempmi,umRef,refseq,qseq,seqLen,fsmall);
		refOvl = cm[i].x1 - cm[i-1].x2;
		qOvl = cm[i].y1 - cm[i-1].y2;
		
		if((cm[i-1].y1 > cm[i-1].y2) && (cm[i].y1 >cm[i].y2)) //two subsequent mums are inverted
		{
			qOvl = cm[i-1].y2 - cm[i].y1;
		}
		if((cm[i-1].y1 > cm[i-1].y2) && (cm[i].y1 < cm[i].y2)) //only previous mum is inverted
		{
			if((i>1) && (cm[i-2].y1 < cm[i-2].y2))
			{
				qOvl = cm[i].y1 - cm[i-1].y1;
			}
			else
			{
				refOvl = 0;
				qOvl = 0;
			}
		}
		if((cm[i-1].y1 < cm[i-1].y2) && (cm[i].y1 >cm[i].y2)) //only second mum is inverted
		{
			if((i < (cm.size()-1)) && (cm[i+1].y1 < cm[i+1].y2))//third mum isn't inverted
			{
				qOvl = cm[i].y2 - cm[i-1].y2; //only if this is is a single inverted mum
			}
			else
			{
				qOvl = 0;
				refOvl = 0;
			}
			vi = findInvertSpan(cm,i);
			//fout<<"INV\t"<<cm[i].rn<<"\t"<<cm[i].x1<<"\t"<<cm[i].x2<<"\t"<<cm[i].qn<<"\t"<<cm[i].y2<<"\t"<<cm[i].y1<<endl;//inverted the query start and end
			fout<<cm[i].rn<<"\t"<<cm[i].x1<<"\t"<<cm[i].x2<<"\tINV\t"<<cm[i].qn<<"\t"<<vi[0]<<"\t"<<vi[1]<<endl;
		} 	
		cnvmi = findDup(cm[i-1],cm[i]);
		if(cnvmi.x1 != 0)
		{
			fout<<cnvmi.rn<<"\t"<<cnvmi.x1<<" "<<cnvmi.x2<<"\tCNV\t"<<cnvmi.qn<<"\t"<<cnvmi.y1<<"\t"<<cnvmi.y2<<endl;
			cm[i].y1 = cnvmi.y1; // changing the interval so that it is not used for calling CNV again			
		}

		if(refOvl < 0)
		{
			if(!(qOvl>0))//if qOvl is 0 or less
			{
				//cout<<"BD "<<cm[i-1].x2<<" "<<cm[i].x1<<" "<<max(cm[i-1].y1,cm[i-1].y2)<<" "<<min(cm[i].y1,cm[i].y2)<<endl;
			}
			else
			{
				gapmi.rn = cm[i-1].rn;
				gapmi.x1 = min(cm[i-1].x2,cm[i].x1);
				gapmi.x2 = max(cm[i-1].x2,cm[i].x1);
				gapmi.y1 = min(cm[i-1].y2,cm[i].y1);
				gapmi.y2 = max(cm[i].y1,cm[i-1].y2);
				findCnvOverlap(cnv,gapmi,storedCNV,fout); //find if there is any cnv in this query interval
				cov = getCoverage(gapmi,masterRef,masterQ);
				//fout<<"INS\t"<<cm[i-1].rn<<"\t"<<cm[i-1].x2<<"\t"<<cm[i].x1<<"\t"<<cm[i-1].qn<<"\t"<<max(cm[i-1].y1,cm[i-1].y2) +1<<"\t"<<min(cm[i].y1,cm[i].y2) -1<<"\t"<<cov[0]<<"\t"<<cov[1]<<endl;
				if(max(gapmi.y1,gapmi.y2) > min(gapmi.y1,gapmi.y2)) //they still differ after cnv has been subtracted
				{
					fout<<cm[i-1].rn<<"\t"<<cm[i-1].x2<<"\t"<<cm[i].x1<<"\tMOR INS\t"<<cm[i-1].qn<<"\t"<<min(gapmi.y1,gapmi.y2) + 1<<"\t"<<max(gapmi.y1,gapmi.y2) -1<<"\t"<<cov[0]<<"\t"<<cov[1]<<endl;
					fs = true; //it reported the insertion
				}
				if( fs == false)
				{
					if(!((cm[i-1].y1 > cm[i-1].y2) && (cm[i].y1 >cm[i].y2)))
					{
						fout<<cm[i-1].rn<<"\t"<<cm[i-1].x2<<"\t"<<cm[i].x1<<"\tINS\t"<<cm[i-1].qn<<"\t"<<max(cm[i-1].y1,cm[i-1].y2) +1<<"\t"<<min(cm[i].y1,cm[i].y2) -1<<"\t"<<cov[0]<<"\t"<<cov[1]<<"\t"<<qOvl<<endl;
					}
					if(((cm[i-1].y1 > cm[i-1].y2) && (cm[i].y1 >cm[i].y2))) //inverted
					{
						fout<<cm[i-1].rn<<"\t"<<cm[i-1].x2<<"\t"<<cm[i].x1<<"\tINS/INV\t"<<cm[i-1].qn<<"\t"<<min(cm[i-1].y1,cm[i-1].y2) +1<<"\t"<<max(cm[i].y1,cm[i].y2) -1<<"\t"<<cov[0]<<"\t"<<cov[1]<<endl;
					}
				}	
			}
		}
		if(refOvl >0)
		{
			if(refOvl > qOvl)
			{
				gapmi.rn = cm[i-1].rn;
				gapmi.x1 = min(cm[i-1].x2,cm[i].x1);
				gapmi.x2 = max(cm[i-1].x2,cm[i].x1);
				gapmi.y1 = min(cm[i].y1,cm[i].y2);
				gapmi.y2 = max(cm[i-1].y1,cm[i-1].y2);
				//perform a cnv search for reference and subtract it from refOvl
				cov = getCoverage(gapmi,masterRef,masterQ);
		//		findCnvOverlapInRef(cnv,gapmi,storedCNV,fout);
				if(!((cm[i-1].y1 > cm[i-1].y2) && (cm[i].y1 >cm[i].y2)))
				{
					fout<<cm[i-1].rn<<"\t"<<cm[i-1].x2 + 1 <<"\t"<<cm[i].x1 -1<<"\tDEL\t"<<cm[i-1].qn<<"\t"<<max(cm[i-1].y1,cm[i-1].y2)<<"\t"<<min(cm[i].y1,cm[i].y2)<<"\t"<<cov[0]<<"\t"<<cov[1]<<endl;
				}
				if(((cm[i-1].y1 > cm[i-1].y2) && (cm[i].y1 >cm[i].y2))) //inverted
				{
					fout<<cm[i-1].rn<<"\t"<<cm[i-1].x2 + 1 <<"\t"<<cm[i].x1 -1<<"\tDEL/INV\t"<<cm[i-1].qn<<"\t"<<min(cm[i-1].y1,cm[i-1].y2)<<"\t"<<max(cm[i].y1,cm[i].y2)<<"\t"<<cov[0]<<"\t"<<cov[1]<<endl;
				}
			}
			if(refOvl < qOvl)
			{
				gapmi.rn = cm[i-1].rn;
				gapmi.x1 = min(cm[i-1].x2,cm[i].x1);
                                gapmi.x2 = max(cm[i-1].x2,cm[i].x1);
				gapmi.y1 = min(cm[i-1].y2,cm[i].y1);
                                gapmi.y2 = max(cm[i].y1,cm[i-1].y2);
                                findCnvOverlap(cnv,gapmi,storedCNV,fout);
				if(max(gapmi.y1,gapmi.y2) > min(gapmi.y1,gapmi.y2)) //they still differ after cnv has been computed
				{
					cov = getCoverage(gapmi,masterRef,masterQ);
					fout<<cm[i-1].rn<<"\t"<<cm[i-1].x2 <<"\t"<<cm[i].x1<<"\tMOD INS\t"<<cm[i-1].qn<<"\t"<<min(gapmi.y1,gapmi.y2) +1 <<"\t"<<max(gapmi.y1,gapmi.y2) -1<<"\t"<<cov[0]<<"\t"<<cov[1]<<endl;
					fs = true;
				}
				cov = getCoverage(gapmi,masterRef,masterQ);
				if(fs == false)
				{
				//	if(!((cm[i-1].y1 > cm[i-1].y2) && (cm[i].y1 >cm[i].y2)))
				//	{
						fout<<cm[i-1].rn<<"\t"<<cm[i-1].x2 <<"\t"<<cm[i].x1<<"\tINS\t"<<cm[i-1].qn<<"\t"<<min(gapmi.y1,gapmi.y2) +1 <<"\t"<<max(gapmi.y1,gapmi.y2) -1<<"\t"<<cov[0]<<"\t"<<cov[1]<<endl;
				//	}
				//	if(((cm[i-1].y1 > cm[i-1].y2) && (cm[i].y1 >cm[i].y2))) //inverted
				//	{
				//		fout<<"INS\t"<<cm[i-1].rn<<"\t"<<cm[i-1].x2 <<"\t"<<cm[i].x1<<"\t"<<cm[i-1].qn<<"\t"<<max(gapmi.y1,gapmi.y2) +1 <<"\t"<<max(gapmi.y1,gapmi.y2) -1<<"\t"<<cov[0]<<"\t"<<cov[1]<<endl;
				  //      }

				}
			}
		}
		fs = false; //reset it

	}
} 
//////////////////////////////////////////////////////////////////////////////////
	
	
void findCnvOverlap(vector<mI> & cnv,mI & mi,vector<mI> & storedCNV,ofstream & fout) //returns a cnv if it overlaps a gap
{
	mI tempmi;
	int ovl =0;
	unsigned int count = 0;
	for(unsigned int i=0;i<cnv.size();i++)
	{
		tempmi = cnv[i];
		if((tempmi.y1 <tempmi.y2) && (cnv[i].qn.find("trans") == string::npos) && (mi.rn == tempmi.rn)) //forward strand. assumes that the gap is also forward oriented,does not have predicted te
		{
			if(!(tempmi.y2<mi.y1) && !(tempmi.y1>mi.y2))
			{
				ovl = min(tempmi.y2,mi.y2) - max(tempmi.y1,mi.y1);
			
				if((double(ovl)/double(tempmi.y2 -tempmi.y1)) >0.9)
				{
					if((find(storedCNV.begin(),storedCNV.end(),tempmi) == storedCNV.end()) && (mi.y1 != mi.y2))
					{
fout<<tempmi.rn<<"\t"<<tempmi.x1<<"\t"<<tempmi.x2<<"\tnCNV\t"<<tempmi.qn<<"\t"<<tempmi.y1<<"\t"<<tempmi.y2<<endl;
						storedCNV.push_back(tempmi);
						mi.y1 = min(tempmi.y2 + 1,mi.y2);//reduce the gap
						//mi.y2 = tempmi.y2 + 1;
						findCnvOverlap(cnv,mi,storedCNV,fout);			
						break;
					}
				}
				
			}
		}
		if((tempmi.y1 > tempmi.y2) && (cnv[i].qn.find("trans") == string::npos) && (mi.rn == tempmi.rn)) // reverse strand. assumes that the gap is also reverse oriented
		{
			if(!(tempmi.y2 > mi.y1) && !(tempmi.y1 < mi.y2))
			{
				ovl = min(tempmi.y1,mi.y1) - max(tempmi.y2,mi.y2);
			
				if((double(ovl)/double(tempmi.y1 -tempmi.y2))>0.9) 
				{
					if((find(storedCNV.begin(),storedCNV.end(),tempmi) == storedCNV.end()) && (mi.y1 != mi.y2))
					{
fout<<tempmi.rn<<"\t"<<tempmi.x1<<"\t"<<tempmi.x2<<"\tnCNV\t"<<tempmi.qn<<"\t"<<tempmi.y1<<"\t"<<tempmi.y2<<endl;

						storedCNV.push_back(tempmi);
						mi.y2 = max(tempmi.y2 - 1,mi.y1);
						//mi.y2 = tempmi.y2 -1;
						findCnvOverlap(cnv,mi,storedCNV,fout);
						break;
					}
				}
			}
		}
		count++;
	}

}			
/////////////////////////////////////////////////////////////////
mI findDup(mI & mi1, mI & mi2)
{
	mI tempmi;
	tempmi.x1= 0;//to use it as a filter for later
	if((mi2.x1 < mi1.x2) && !(mi1 == mi2)) //i.e. mi2.x1 falls within the previous interval. this works because mums are sorted by reference cord.
	{
//cout<<"candidate "<<mi1.x1<<" "<<mi1.x2<<" "<<mi2.x1<<" "<<mi2.x2<<endl;
		if(mi2.x2 > mi1.x2) //mi1 amd mi2 just overlaps but one is not contained within another
		{
			tempmi.rn = mi1.rn;
			tempmi.x1 = mi2.x1;
			tempmi.qn = mi1.qn;
			tempmi.x2 = min(mi1.x2,mi2.x2); //smallest of the two
			if(mi2.y1 <mi2.y2) //forward oriented
			{
				tempmi.y1 = mi2.y1;
				tempmi.y2 = mi2.y1 + (mi1.x2 - mi2.x1);
			}
			if(mi2.y1 > mi2.y2)//reverse oriented
			{
				tempmi.y1 = mi2.y1;
				tempmi.y2 = mi2.y1 - (mi1.x2 - mi2.x1);
			}
//cout<<"cnv "<<tempmi.rn<<" "<<tempmi.x1<<" "<<tempmi.x2<<" "<<tempmi.qn<<" "<<tempmi.y1<<" "<<tempmi.y2<<endl;				
		}
		if(!(mi2.x2>mi1.x2))//mi2 is contained within mi1
		{
			tempmi.rn = mi1.rn;
			tempmi.x1 = mi2.x1;
			tempmi.x2 = min(mi1.x2,mi2.x2);
			tempmi.qn = mi1.qn;
			tempmi.y1 = mi2.y1;
			tempmi.y2 = mi2.y2;
//cout<<"cnv "<<tempmi.rn<<" "<<tempmi.x1<<" "<<tempmi.x2<<" "<<tempmi.qn<<" "<<tempmi.y1<<" "<<tempmi.y2<<endl;
		}
	}
	return tempmi;
}
/////////////////////////////////////////////////////////////////////////////
void findCnvOverlapInRef(vector<mI> & cnv,mI & mi,vector<mI> & storedCNV,ofstream & fout) //returns a cnv if it overlaps a gap
{
	mI tempmi;
	int ovl =0;
//	unsigned int count = 0;
	for(unsigned int i=0;i<cnv.size();i++)
	{
		tempmi = cnv[i];
		if((cnv[i].qn.find("trans") == string::npos) && (mi.rn == tempmi.rn)) //forward strand. assumes that the gap is also forward oriented,does not have predicted te
		{
			if(!(tempmi.x2<mi.x1) && !(tempmi.x1>mi.x2))
			{
				ovl = min(tempmi.x2,mi.x2) - max(tempmi.x1,mi.x1);
				if((double(ovl)/double(tempmi.x2 -tempmi.x1)) >0.9)
				{
					if((find(storedCNV.begin(),storedCNV.end(),tempmi) == storedCNV.end()) && (mi.x1 != mi.x2))
					{
fout<<"rnCNV\t"<<tempmi.rn<<"\t"<<tempmi.x1<<"\t"<<tempmi.x2<<"\t"<<tempmi.qn<<"\t"<<tempmi.y1<<"\t"<<tempmi.y2<<endl;
						storedCNV.push_back(tempmi);
						findCnvOverlap(cnv,mi,storedCNV,fout);
						break;
					}
				}
			}
		}
	}
}
/////////////////////////////////////////////////////////////////////////////	
vector<int> findInvertSpan(vector<mI> & cm, int i)
{
	bool fi = true; //found invert abbreviated
	vector<int> vi(2);

	while(fi == true)
	{
		i++;
		if(cm[i].y1 > cm[i].y2)//ith mum is inverted
		{	
			fi = true;
			vi[0] = cm[i].y2;
			vi[1] = cm[i].y1;
		}
		if(cm[i].y1 < cm[i].y2)//ith mum is not inverted
		{
			fi = false;
		}
	
	}
	return vi;
}
