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

class mgapC {

public:
map<int,string> refName;
map<int,string> qName;
map<int,vector<int> > refClust; //stores all the clusters from the reference sequence
map<int,vector<int> > qClust; //stores all the clsuters from the query sequence
map<int,vector<int> > dupList;
map<int,vector<string> > dupName; //stores names of dup reference, query 1, query 2 
map<int,vector<int> > dupCord; //stores ends of dups in reference ,query1,query 2,coordinates
map<int,vector<int> > filterList; //first four coordinates are reference, last four are query
//map<string,int> dupCount;
};

void comparClust(mgapC & cluster);
int overlapD(vector<int>& rv, vector<int>& mRef);
//void checkTE(mgapC & cluster);
vector<int> findDupEnds(int & ref_st1, int & ref_end1,int & ref_st2, int & ref_end2, int & q_st1,int & q_end1, int & q_st2, int & q_end2);
void filterDup(mgapC & cluster);
bool ovlChk(vector<int> &v1, vector<int> & v2);
#endif
