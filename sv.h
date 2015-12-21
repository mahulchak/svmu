#ifndef SV_H_
#define SV_H_

#include<iostream>
#include<fstream>
#include<vector>
#include<string>
#include<map>
#include<cstdlib>
using namespace std;
string xtractcol(string str, char c, int n); //str: the string substring will be lifted, c: delimiting character, n = # of col

class repeats {

public:
vector<int> masterRepeat; //reference sequences that correspond to the repeats
vector<int> childRepeat; //query sequences that correspond to the repeats
};

void comparClust(map<int,string>& ref_name,map<int,string>& q_name, map<int,vector<int> >& mRef, map<int,vector<int> >& mQ);
int overlapD(vector<int>& rv, vector<int>& mRef);
//void chkOvl(map<int,string>& ref_name,map<int,string>& q_name,map<int,vector<int> >& mRef, map<int,vector<int> >& mQ, vector<int> & rv,vector<int> & qv, string & name, ofstream & fout);

#endif
