
#ifndef INDEL_H_
#define INDEL_H_

#include<vector>
#include<string>
#include<map>
#include<fstream>
#include<algorithm>
#include<list>
using namespace std;


class asmMerge {

public:

vector<string> r_name;
vector<string> q_name;
map<string,int> ref_len;
map<string,int> q_len;
map<string,vector<int> > ref_st;
map<string,vector<int> > ref_end;
map<string,vector<int> > q_st;
map<string,vector<int> > q_end;
map<string,vector<int> > storeDelStart;
map<string,vector<int> >storeDelEnd;
map<string,vector<string> >storeDelName;
map<string, vector<int> >storeInsStart;
map<string,vector<int> >storeInsEnd;
vector<string> refChrom;
vector<string> qChrom;
map<string,vector<int> >refChromPos;
map<string,vector<int> >qChromPos;
map<string,vector<int> > newrefSt;
map<string,vector<int> > newrefEnd;
map<string,vector<int> > newq_St;
map<string,vector<int> > newq_end;
map<string,vector<int> > anmlrefSt;
map<string,vector<int> > anmlrefEnd;
map<string,vector<int> > anmlq_St;
map<string,vector<int> > anmlq_end;
};

string xtractcol(string str,char c, int n);
void writeToFile(asmMerge & merge);
void findIndel(asmMerge & merge, char mutType, float & prop);
char checkIndel(string tempname,asmMerge & merge,int k,int j,float & prop);
void fillChromPos(string & tempname,asmMerge & merge,int & length);
void addCoverage(asmMerge & merge,string & str, int ref_st, int ref_end);
void buildCoverage(asmMerge & merge);
bool chkOvlQ(string tempname,asmMerge & merge,int k,int j);
bool chkOvlR(string tempname,asmMerge & merge,int k,int j);
int maxD(int & qf1,int & qe1, int & qf2, int & qe2);
void filterInsCall(asmMerge & merge);
void collapseRange(asmMerge & merge);
void lisCalculator(asmMerge & merge,string & tempname,vector<int>& ref_st,vector<int>& ref_end,vector<int>& q_st, vector<int>& q_end, int dist);
void ovlStoreCalculator(asmMerge & merge);
int midDist(vector<int>& ref_st,vector<int>& ref_end,vector<int>& q_st, vector<int>& q_end);
void addvTol(vector<int> & q_st,vector<int>& temp,int k);
unsigned int pos(int & elem, vector<int> & v);
double cov(asmMerge & merge,string & str,int st, int end);
void findDup(asmMerge & merge);
#endif
