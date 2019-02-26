#ifndef SV_H_
#define SV_H_

#include<iostream>
#include<fstream>
#include<string>
#include<cstdlib>
#include<vector>
#include<map>
#include<algorithm>
#include<list>
#include<climits>
#include<cmath>
#include<utility>
#include<iomanip>

using namespace std;

//to store coordinates at alignment level
struct mI {
	string rn;
	string qn;
        int x1;//reference start
        int x2;//reference end
        int y1;//query start
        int y2;//query end
        char c;//qualifier. special comments or information can be added to this:i=inversion;
	vector<int> mv;        
	int l;//length of the MUM
        bool operator < (const mI& mum1) const
        {
//sort the mums by ref chrom name and then by start cords of refs 
                return (rn < mum1.rn)|| ((rn == mum1.rn) && (x1 < mum1.x1)) || ((rn == mum1.rn) && (x1 == mum1.x1) && (x2 < mum1.x2));
        }
        bool operator == (const mI& mum1) const
        {
                return x1 == mum1.x1 && x2 == mum1.x2 && y1 == mum1.y1 && y2 == mum1.y2;
        }
        };

struct qord {
	string name;
	int cord;
	bool operator < (const qord& q1) const
	{
		return(name<q1.name) || ((name == q1.name) && (cord<q1.cord));
	}
	};

class chromPair {
	public:
	vector<mI> mums;	
	vector<mI> cm; //conserved mems
	vector<mI> ncm; //conserved mems from reverse side
	vector<mI> gap;
};

bool qusort(mI mi1, mI mi2); //to sort the mI based on query coordinates
vector<int> makeChromBucket(int refLen);
bool msort(mI mi1, mI mi2);
bool isort(mI m1, mI m2);
bool lsort(mI m1,mI m2);
void storeCords(vector<int> & masterRef,vector<int> & masterQ, mI & mi);
void storeCords(vector<int> & masterQ, mI & mi);//overloaded
void storeNameCount(vector<int> & chromDensityRef,vector<int> & chromDensityQ,map<string,int> & lookUpRef,map<string,int> & lookUpQ, mI & mi);
void storeCords(map<int,vector<qord> > & mRef, mI & mi, ofstream & fout); //overloaded
void storeCordsCm(map<int,vector<qord>> & mRef, mI & mi);
mI findClosest(mI & mi, vector<mI> & mums);//overloaded function
vector<double> getCoverage(mI & mi, vector<int> & masterRef,vector<int> & masterQ);
vector<double> getCoverage(mI & mi, vector<int> & masterRef,vector<int> & masterQ,float p);
vector<double> getChromCount(mI & mi, vector<int> & chromDensityRef, vector<int> & chromDensityQ);
void gapCloser(mI & mi, vector<mI> ncm, vector<mI>& cm);
mI returnMumByQ1(int & y1,vector<mI> & mums);
mI returnMumByQ2(int & y1,vector<mI> & mums);
void gapCloserRev(mI & mi, vector<mI> ncm, vector<mI> & cm);
int nearestInt(double d);
void annotGaps(vector<mI> & cm,vector<int> & masterRef, vector<int> & masterQ,vector<int> & chromDensityRef, vector<int> & chromDensityQ,vector<mI> & cnv,map<int,vector<qord> > & umRef, string & refseq, string & qseq,vector<int> & seqlen,ofstream & fout, ofstream & fsmall,int & id);
void readUniq(ifstream & fin,vector<mI> & cm, map<int,vector<qord> > & umRef,vector<int> & masterHQ);
void callSmall(mI & mi,map<int,vector<qord> > & umRef, string & refseq, string & qseq,vector<int> & seqlen,ofstream & fsmall);
void findCnvOverlap(mI & gapmi,vector<mI> ncm, vector<mI> cnv, vector<int> & masterRef, vector<int> & masterQ,vector<int> & chromDensityRef,vector<int> & chromDensityQ,ofstream & fout, int & id);
mI findDupRef(mI & mi1, mI & mi2);
mI findDupQ(mI & m1, mI & m2);
char comp(char & N);
void findInnie(vector<mI> & mums,mI & mi);

#endif
