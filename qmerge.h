
#ifndef MASM_H_
#define MASM_H_

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
map<string,int> storeCount;


};



string xtractcol(string str,char c, int n);
void countCopy(string& fin, asmMerge & merge);
#endif
