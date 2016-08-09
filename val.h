
#ifndef VAL_H_
#define VAL_H_

#include<vector>
#include<string>
#include<map>
#include<fstream>
#include<algorithm>
#include<cstdlib>

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
map<string,int> storeCount;
map<string,int> seqCount;//stores the number of sequences to which a SV finds a hit
map<string,vector<string> > storeName;
map<string,int> ovlStore;
map<string,string> storeHomolog;
map<string,int> storeHomAln;
map<string,vector<int> > masterQst;
map<string,vector<int> > masterQend;
map<string,vector<double> >cov;
};



string xtractcol(string str,char c, int n);
void countCopy(string& fin, asmMerge & merge);
void collapseRange(asmMerge & merge, double & rCO, int & qCO);
bool chkOvl(asmMerge & merge, string & str,unsigned int & j, double & rCO, int & qCO);
int ovlCalculator(vector<int>& q_st, vector<int>& q_end);
void ovlStoreCalculator(asmMerge & merge);
void findChromPartner(asmMerge & merge);
void masterQlist(asmMerge & merge);
void reducList(asmMerge & merge);
bool innie(vector<int> & masterQst, vector<int> & masterQend,vector<double> & cov, int & qSt, int & qEnd);
#endif
