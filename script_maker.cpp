#include<iostream>
#include<fstream>
#include<string>
#include<cstdlib>

using namespace std;

int main(int argc, char * argv[])
{
	if ((argc == 1) || (argc <2))
	{
		cerr<<"Usage: "<<argv[0]<<" reference.fasta list_of_query_fasta"<<endl;;
		exit(EXIT_FAILURE);
	} 
size_t pos;
string str,script_file,name;
ofstream fout;
ifstream fin;
fin.open(argv[2]);
	while(getline(fin,str)) //str is chromosome file name
	{
		pos = str.find(".fa");
		name =	str.substr(0,pos); // name is chromosome name
		script_file = "job_"+name;
		fout.open(script_file);	//script file is named after the chrom name
		fout<<"mummer -b -l 20 "<<argv[1]<<" "<<str<<" | mgaps -l 100 -d 5 -f .12 -s 200 > "<<name<<".mgaps"<<endl;
		fout<<"sed -i 's/^ //g' "<<name<<".mgaps"<<endl;
		fout<<"sed -i 's/^ //g' "<<name<<".mgaps"<<endl;
		fout<<"sed -i 's/^ //g' "<<name<<".mgaps"<<endl;
		fout<<"sed -i 's/^> />/g' "<<name<<".mgaps"<<endl;
		fout<<"sed -i 's/Reverse//g' "<<name<<".mgaps"<<endl;
		fout<<"gawk '{print $1\"\\t\"$2\"\\t\"$3\"\\t\"$4\"\\t\"$5\"\\t\"$6}' "<<name<<".mgaps > "<<name<<".awk.mgaps"<<endl;
		fout<<"./svmum "<<name<<".awk.mgaps mumToBed."<<name<<".bed "<<name<<endl; // check if the last name should be query contig or reference contig name
		fout.close();
	}
fin.close();
return 0;
}
