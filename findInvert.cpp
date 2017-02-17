#include<iostream>
#include<fstream>
#include<string>
#include<cstdlib>
#include<vector>
using namespace std;

int main(int argc, char *argv[])
{
	if(argc <2)
	{
		cerr<<"Usage: "<<argv[0]<<" foo.mgaps"<<endl;
		exit(EXIT_FAILURE);
	}
	
	
	string str = "Reverse";
	string foo = string(argv[1]);
	string line;
	int refStart = 0, refEnd = 0, qTemp = 0, range =0;
	size_t pos1,pos2;
	ifstream fin;
	bool record = 0;
	fin.open(argv[1]);
	vector<int> tempCluster;
	while(getline(fin,line))
	{
		if((line.find(str) != string::npos) && (line.find('>') != string::npos))//Reverse is present in the fasta header
		{
			record = true;
		}
		if((line.find(str) == string::npos) && (line.find('>') != string::npos))
		{
			record = false;
		}
		if(record == true)//if Reverse was in the fasta header
		{
			if((line.find('>') == string::npos) && (line.find('#') == string::npos))
			{
				refStart = stoi(line,&pos1);
				qTemp = stoi(line.substr(pos1),&pos2);
				range = stoi(line.substr(pos1+pos2));
				refEnd = refStart +range;
				tempCluster.push_back(refStart);
				tempCluster.push_back(refEnd);
			//cout<<line<<"\t"<<refStart<<"\t"<<refEnd<<"\t"<<range<<endl;	
			}
			else 
			{
				//cout<<line<<endl;
				if(tempCluster.size() != 0)
				{
					cout<<foo.substr(0,foo.find('.'))<<"\t"<<tempCluster[0]<<"\t"<<tempCluster[tempCluster.size()-1]<<endl;
				}
				tempCluster.clear();
			}
		}
	}
	fin.close();
	
	return 0;
}
		
