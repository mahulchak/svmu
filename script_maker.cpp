#include<iostream>
#include<fstream>
#include<string>
#include<cstdlib>

using namespace std;

int main(int argc, char * argv[])
{
	if ((argc == 1) || (argc <2))
	{
		cerr<<"Usage: "<<argv[0]<<" -q your_assembly.fasta -l list_of_query_fasta -s cluster_separation -ln minimum_length"<<endl;;
		exit(EXIT_FAILURE);
	} 
size_t pos;
string str,script_file,name;
ofstream fout;
ifstream fin;
fin.open(argv[4]);
int sep = stoi(argv[6],NULL);
int minOvl = stoi(argv[8],NULL);
	while(getline(fin,str)) //str is chromosome file name
	{
		pos = str.find(".fa");
		name =	str.substr(0,pos); // name is chromosome name
		script_file = "job_"+name;
		fout.open(script_file.c_str());	//script file is named after the chrom name
		fout<<"mummer -b -l 20 "<<str<<" "<<argv[2]<<" | mgaps -l "<<minOvl<<" -d 5 -f .12 -s "<<sep<<" > "<<name<<".mgaps"<<endl;
		fout<<"sed -i 's/^ //g' "<<name<<".mgaps"<<endl;
		fout<<"sed -i 's/^ //g' "<<name<<".mgaps"<<endl;
		fout<<"sed -i 's/^ //g' "<<name<<".mgaps"<<endl;
		fout<<"sed -i 's/^> />/g' "<<name<<".mgaps"<<endl;
		fout<<"sed -i 's/Reverse//g' "<<name<<".mgaps"<<endl;
		fout<<"gawk '{print $1\"\\t\"$2\"\\t\"$3\"\\t\"$4\"\\t\"$5\"\\t\"$6}' "<<name<<".mgaps > "<<name<<".awk.mgaps"<<endl;
		fout<<"./svmu "<<name<<".awk.mgaps mumToBed."<<name<<".bed "<<name<<endl; // check if the last name should be query contig or reference contig name
		fout.close();
	}
fin.close();
return 0;
}
