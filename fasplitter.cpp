//this takes each fasta sequence from a fasta file and generates new fasta file for each sequence. The name of the new fasta files are derived from the sequence names.
//useful for CNVnator and other stuff

#include<iostream>
#include<fstream>
#include<string>
#include<cstring>
#include<cstdlib>

using namespace std;

string clip_string(string str, char c);

int main(int argc, char *argv[])

{
ifstream fin;
ofstream fout;
	if(argc ==1)
	{
	cerr<<"Usage: "<<argv[0]<<" foo.fasta '.fa'(Y/N)?" <<endl;
	exit(EXIT_FAILURE);
	}
fin.open(argv[1]);
string str,str_hist,temp;
str_hist="empty";
while(getline(fin,str))
{
if(str[0] == '>')
	{
//	cout<<str<<endl;
	if(str_hist != str && str_hist != "empty")
		{fout.close();
		}
	if(*argv[2] == 'Y')
	{
		temp = clip_string(str,' ')+".fa"; //original was ' '.
//cout<<temp<<endl;
		fout.open(temp.c_str());
	}
	else
	{
		fout.open(clip_string(str,' ').c_str());
	}
	
	str_hist = str;

	}

fout<<str<<endl;
//fout.close();


}
return 0;
}


string clip_string(string str, char c)

{
string str1;
str1 = str;
for(unsigned int i =1;i<str.size();i++)
	{ str1[i-1] = str1[i];
		if (str1[i] == c)
		{ str1.resize(i-1); //changed this for quivered sequence names. original was i-1
		 
		 i = str.size();
		}
	}
//cout<<str1<<endl;
return str1;
}

