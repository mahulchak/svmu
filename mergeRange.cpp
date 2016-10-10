#include<vector>
#include<cstdlib>
#include<iostream>
#include<fstream>
#include<algorithm>
using namespace std;

string xtractcol(string str,char c,int n);

int main(int argc, char * argv[])
{
	if(argc == 1)
	{
		cerr<<"Usage: "<<argv[0]<<" Sorted range file"<<endl;
		exit(EXIT_FAILURE);
	}

	string str;

	vector<string> name;
	vector<int> start;
	vector<int> end;
	ifstream fin;
	fin.open(argv[1]);
	while(getline(fin,str))
	{
		name.push_back(xtractcol(str,'\t',1));
	
		start.push_back(stoi(xtractcol(str,'\t',2),nullptr));
		
		end.push_back(stoi(xtractcol(str,'\t',3),nullptr));
	}
	fin.close();

	for(unsigned int i =0;i<name.size()-1;i++)
	{
			
		if((abs(start[i] - start[i+1]) <10) && (abs(end[i] - end[i+1])<10) && (name[i] == name[i+1]))
		{
			start[i] = min(start[i],start[i+1]);
			end[i] = max(end[i],end[i+1]);
			name.erase(name.begin()+i+1);
			start.erase(start.begin()+i+1);
			end.erase(end.begin()+i+1);
			i--;
		}

		//cout<<name[i]<<"\t"<<start[i]<<"\t"<<end[i]<<endl;
		
	}
	
	for(unsigned int i=0;i<name.size();i++)
	{
		cout<<name[i]<<"\t"<<start[i]<<"\t"<<end[i]<<endl;
	}
		
return 0;
}

string xtractcol(string str, char c, int n)
{
int count =0;
int j =0; // j is the index for the string that'll be returned
string elem;
for(unsigned int i=0;i<str.size() && count<n;i++)
        {
        if(str[i] == c)
                {
                count++; //keep a count of the field separator
                }
        if(count == n-1)
                {
                elem.push_back(str[i]);
                j++;
                }
        }
return elem;
}

