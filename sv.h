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

using namespace std;

//to store coordinates at alignment level
struct mI {
	string rn;
	string qn;
        int x1;//reference start
        int x2;//reference end
        int y1;//query start
        int y2;//query end
        //int m;//number of mutations
	vector<int> mv;        
        bool operator < (const mI& mum1) const
        {
                return(x1 < mum1.x1) || ((x1 == mum1.x1) && (x2 < mum1.x2));
        }
        bool operator == (const mI& mum1) const
        {
                return x1 == mum1.x1 && x2 == mum1.x2 && y1 == mum1.y1 && y2 == mum1.y2;
        }
        };
//to store coordinates at base pair level
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
	vector<mI> gap; //gaps are represented as mums
	vector<mI> cc; //cnv candidates
	vector<mI> in;//stores insertion mums in reference
	vector<mI> del; //stores deletion mums in query
};

bool qusort(mI mi1, mI mi2); //to sort the mI based on query coordinates
vector<int> makeChromBucket(int refLen);
void storeCords(vector<int> & masterRef,vector<int> & masterQ, mI & mi);
void storeCords(map<int,vector<qord> > & mRef, mI & mi); //overloaded
void storeCordsCm(map<int,vector<qord>> & mRef, mI & mi);
int findDist(int & x1, int & y1, int & c);//distance between the diagonal and the other MUMs
bool detectShadow(mI & mum, vector<mI> & mums, unsigned int n);
mI findClosest(mI & mi, vector<mI> & mums,unsigned int i, vector<int> & masterRef,vector<int> & masterQ);
mI findClosest(mI & mi, vector<mI> & mums);//overloaded function
void recordShadow(unsigned int i, unsigned int j, vector<mI> & mums, vector<mI> & sm);
vector<double> getCoverage(mI & mi, vector<int> & masterRef,vector<int> & masterQ);
void findPartnerCord(mI & mi, vector<mI> & mums,char c);
//void splitByCoverage(chromPair & cp,vector<int> & chrom, vector<mI> & mums,vector<int> & masterRef, vector<int> & masterQ);
void splitByCoverage(chromPair & cp,vector<int> & rchrom,vector<int> & qchrom);
void gapCloser(mI & mi, vector<mI> ncm, vector<mI>& cm);
vector<mI> findQuery(map<int,vector<qord> > & mRef, mI & mi,vector<int> & masterRef, vector<int> & masterQ);
int nearestInt(double d);
void annotGaps(vector<mI> & cm,map<int,vector<qord> > & mRef,vector<int> & masterRef, vector<int> & masterQ,vector<mI> & cnv,map<int,vector<qord> > & umRef, string & refseq, string & qseq,vector<int> & seqlen,ofstream & fout, ofstream & fsmall);
void readUniq(ifstream & fin,vector<mI> & cm, map<int,vector<qord> > & umRef);
void callSmall(mI & mi,map<int,vector<qord> > & umRef, string & refseq, string & qseq,vector<int> & seqlen,ofstream & fsmall);
void findCnvOverlap(vector<mI> & cnv,mI & mi,vector<mI> & storedCNV,ofstream & fout);
void findCnvOverlapInRef(vector<mI> & cnv,mI & mi,vector<mI> & storedCNV,ofstream & fout);
mI findDup(mI & mi1, mI & mi2);
char comp(char & N);
void xtracTrans(vector<mI> & cm,ofstream & ftest);
vector<int> findInvertSpan(vector<mI> & cm, int i);
#endif
