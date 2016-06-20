
#ifndef INDEL_H_
#define INDEL_H_

#include<vector>
#include<string>
#include<map>
#include<fstream>
#include<algorithm>

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
map<string, vector<int> >storeInsStart;
map<string,vector<int> >storeInsEnd;
vector<string> refChrom;
vector<string> qChrom;
map<string,vector<int> >refChromPos;
map<string,vector<int> >qChromPos;

};

string xtractcol(string str,char c, int n);
void writeToFile(asmMerge & merge);
void findIndel(asmMerge & merge);
char checkIndel(string tempname,asmMerge & merge,int k,int j);
void fillChromPos(asmMerge & merge);
void addCoverage(asmMerge & merge,string & str, int ref_st, int ref_end);
void buildCoverage(asmMerge & merge);
bool chkOvlQ(string tempname,asmMerge & merge,int k,int j);
bool chkOvlR(string tempname,asmMerge & merge,int k,int j);
int maxD(int & qf1,int & qe1, int & qf2, int & qe2);
void filterInsCall(asmMerge & merge);
#endif
