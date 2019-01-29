#include<iostream>
#include "sv.h"
#include "seqIO.h"

using namespace std;

using chroms = map<string,chromPair>;
using ccov = vector<int>;
using vq = vector<qord>;
int main(int argc, char *argv[])
{
	if(argc <2)
	{
		cerr<<"Usage: "<<argv[0]<<" foo.delta ref.fasta query.fasta cutoff mode(h/l)"<<endl;
		exit(EXIT_FAILURE);
	}

	chroms allChrom;
	
	map<string,ccov> masterRef; //stores sequence coverage but it can also be used to find reference chromosome lengths
	map<string,ccov>masterQ; //stores sequence coverage but it can also be used to find query chromosome lengths
	map<string,ccov>masterHQ; //same as masterQ but records coverage only for homologous pairs
	map<string,vector<string> > cp; //cp is an alias for Chromosome partner. Each reference name index has a vector of unqiue alignments which are part of these
	map<string,vector<string> > hcp;//hcp stands for homologous cp
	
	map<string,map<int,vq> > mRef; //stores the coordinates of query on reference chromosomes
	map<string,map<int,vq> > umRef;//stores the coordinates of unique reference to query map; requires re-reading the file
	map<string,vector<mI> > mir;//stores all reference mums, but does not store query coordinates
	map<string,vector<mI> > miq; // stores all query mums, but does not store reference coordinates
	map<string,string> refseq;
	map<string,string> qseq;
	map<string,vector<int> > seqLen;//length of sequences.first element is ref and second is query
	map<string,bool> qStrand; //stores whether query strand is forward strand or reverse strand
	mI tempmi,prevmi,temprmi,tempmi2;

	string foo = string(argv[1]);
	string line, chromName,refName,qName,indexAln;
	int refStart = 0, refEnd = 0, qStart = 0, qEnd = 0, refLen =0, qLen =0, count = -1,indelPos =0, forCount =0, revCount = 0,index =0;
	unsigned int cutoff = 0;
	cutoff = stoi(argv[4]);
	double hcr;//homolog cutoff ratio
	vector<double> vd;
	vector<int> vi;
	vector<mI> vmi,tempVmi,vm,qvm,gapmi;//qvm is query sorted vm
	size_t pos1,pos2,namePos;
	
	ifstream fin, refFasta, qFasta;
	ofstream fout,fcnv,fsmall,ftrans,findel,fcords;
	fin.open(argv[1]);
	fcords.open("cords.txt");
	while(getline(fin,line))
	{
		
		if(line.find('>') != string::npos)//start of an aligning chromosome description
		{
						
			refName = line.substr(1,line.find(' ')-1);
			pos1 = line.find(' '); //position of the first space
			pos2 = line.find(' ',pos1+1);//position of the second space
			qName = line.substr(pos1+1, pos2-pos1-1); //up to the second space
//cout<<qName<<endl;
			pos1 = line.find(' ',pos2+1); //recycling pos1 to find pos3
			refLen = stoi(line.substr(pos2+1,pos1-pos2));//reference length
			qLen = stoi(line.substr(pos1));//from last space till end 
			indexAln = refName + qName;
			count = -1;
			seqLen[indexAln].push_back(refLen);
			seqLen[indexAln].push_back(qLen);
			cp[refName].push_back(indexAln); //adding the alignment to the list of refName alignments
			if(masterRef[refName].size() == 0)//if they have not been created
			{
				masterRef[refName] = makeChromBucket(refLen);
			}
			if(masterQ[qName].size() == 0)//if they have not been created
			{
				masterQ[qName] = makeChromBucket(qLen);
			}
		}
		if((line.size() <10) && (refName != "") && (count > -1))
		{
			
			indelPos = stoi(line);
			if(indelPos < -1)
			{
				refStart = refStart + abs(indelPos) -1;
				vi.push_back(refStart*-1);
			}
			if(indelPos == -1)
			{
				vi.push_back(-1);
			}
			if(indelPos > 0)
			{
				refStart = refStart + abs(indelPos);	
				vi.push_back(refStart);
			}
					
			if(indelPos ==0) //reached the end of the indel description
			{
				tempmi.mv = vi;
				storeCords(masterRef[refName],masterQ[qName],tempmi);
				//storeCords(mRef[refName],tempmi,fcords);
				tempmi.mv.clear();//delete this?
				allChrom[indexAln].mums.push_back(tempmi);
//cout<<indexAln<<tempmi.x1<<'\t'<<tempmi.x2<<'\t'<<tempmi.qn<<'\t'<<tempmi.y1<<'\t'<<tempmi.y2<<endl;
				vi.clear();//reset it once its values are used
			}
				
			count++;
			
		}
		if((line.find('>') == string::npos) && (line.size() >10) && (refName != "")) //when describing alignment segments
		{
//cout<<line<<endl;
				tempmi.rn = refName;
				tempmi.qn = qName;		
				refStart = stoi(line,&pos1);
				refEnd = stoi(line.substr(pos1),&pos2);
				qStart = stoi(line.substr(pos1+pos2), &namePos);
				qEnd = stoi(line.substr(pos1+pos2+namePos));
				tempmi.x1 = refStart;
				tempmi.x2 = refEnd;
				tempmi.y1 = qStart;
				tempmi.y2 = qEnd;
//cout<<refName<<"\t"<<refStart<<"\t"<<refEnd<<"\t"<<qName<<"\t"<<qStart<<"\t"<<qEnd<<"\t"<<allChrom[indexAln].mums.size()<<endl;
				count = 0;
	//			--refStart;//to count the mutation distance

		}
	}
	fin.close();
	for(chroms::iterator it = allChrom.begin();it!= allChrom.end();it++)
	{
		vm.clear();
		qvm.clear();
		tempVmi.clear();
		indexAln = it->first;
		count = 0;
		sort(allChrom[indexAln].mums.begin(),allChrom[indexAln].mums.end());
		if(allChrom[indexAln].mums.size() >2)
		{
			for(unsigned int i = 0; i<allChrom[indexAln].mums.size();i++)
			{
				tempmi = allChrom[indexAln].mums[i];
				fcords<<tempmi.rn<<'\t'<<tempmi.x1<<'\t'<<tempmi.x2<<'\t'<<tempmi.qn<<'\t'<<tempmi.y1<<'\t'<<tempmi.y2<<endl;
				if(findInnie(allChrom[indexAln].mums,tempmi) == false)
				{
					vd = getCoverage(tempmi,masterRef[tempmi.rn],masterQ[tempmi.qn],0.3);
					if((vd[0] <2) && (vd[1]<2))
					{
						qvm.push_back(tempmi);
						count = tempmi.x2 - tempmi.x1; //count the alignment length
//cout<<"cm\t"<<tempmi.rn<<'\t'<<tempmi.x1<<'\t'<<tempmi.x2<<'\t'<<tempmi.qn<<'\t'<<tempmi.y1<<'\t'<<tempmi.y2<<'\t'<<vd[0]<<'\t'<<vd[1]<<endl;
					}
					else
					{
						vm.push_back(tempmi);
						allChrom[indexAln].ncm.push_back(tempmi);
					}
				}
				else
				{
					tempVmi.push_back(tempmi);
				}
			}
			sort(allChrom[indexAln].mums.begin(),allChrom[indexAln].mums.end(),qusort);
			for(unsigned int i =0;i<tempVmi.size();i++)
			{
				tempmi = tempVmi[i];
				if(findInnieQ(allChrom[indexAln].mums,tempmi) == true)
				{
//cout<<"shadow\t"<<tempmi.rn<<'\t'<<tempmi.x1<<'\t'<<tempmi.x2<<'\t'<<tempmi.qn<<'\t'<<tempmi.y1<<'\t'<<tempmi.y2<<endl;
				}
				else
				{
cout<<"noshadow\t"<<tempmi.rn<<'\t'<<tempmi.x1<<'\t'<<tempmi.x2<<'\t'<<tempmi.qn<<'\t'<<tempmi.y1<<'\t'<<tempmi.y2<<endl;
				}
			}
			if(qvm.size()>1)
			{
				//vm = tempVmi;//vm has the ncm
				sort(vm.begin(),vm.end());
				sort(tempVmi.begin(),tempVmi.end(),qusort);
				sort(qvm.begin(),qvm.end());
	//		cout<<indexAln<<"\tLength is\t"<<count<<"\t"<<qvm[qvm.size()-1].x2-qvm[0].x1<<"\t"<<hcr<<"\t"<<qvm[qvm.size()-1].x2<<'\t'<<qvm[0].x1<<endl;
				for(unsigned int i = 0;i<qvm.size();i++)
				{
					prevmi = qvm[i];
					if(allChrom[indexAln].cm.size() > 0)
					{
						temprmi = allChrom[indexAln].cm[allChrom[indexAln].cm.size()-1];
					}
					allChrom[indexAln].cm.push_back(qvm[i]);
					if(prevmi.x1 - temprmi.x1 > 100)
					{
						if((prevmi.y1 > prevmi.y2) && (temprmi.y1 > temprmi.y2) && \
						 (temprmi.y2 - prevmi.y1 >100)) // inverted
						{
							tempmi2.rn = prevmi.rn;
							tempmi2.x1 = temprmi.x2;
							tempmi2.x2 = prevmi.x1;
							tempmi2.qn = prevmi.qn;
							tempmi2.y1 = temprmi.y2;
							tempmi2.y2 = prevmi.y1;
						}
						if((prevmi.y1 < prevmi.y2) && (temprmi.y1 < temprmi.y2) && \
						(prevmi.y1 - temprmi.y2 > 100)) //forward
						{
							tempmi2.rn = prevmi.rn;
							tempmi2.x1 = temprmi.x2;
							tempmi2.x2 = prevmi.x1;
							tempmi2.qn = prevmi.qn;
							tempmi2.y1 = temprmi.y1;
							tempmi2.y2 = prevmi.y1;
						}
						if(tempmi2.x1 != 0)
						{
//cout<<tempmi2.rn<<"\t"<<tempmi2.x1<<"\t"<<tempmi2.x2<<"\t"<<tempmi2.qn<<"\t"<<tempmi2.y1<<"\t"<<tempmi2.y2<<endl;
							gapmi.push_back(tempmi2);
							tempmi2.x1 = 0; // reset it
						}
					}//this is for x>100				
				}//this is for the for loop
				if(allChrom[indexAln].cm.size() > 0)
				{
					for(unsigned int i =0;i<gapmi.size();i++)
					{	
						tempmi = gapmi[i];
						gapCloser(tempmi,vm,allChrom[indexAln].cm);
					}
					hcp[allChrom[indexAln].cm[0].rn].push_back(indexAln);//homologous alignment
					sort(allChrom[indexAln].cm.begin(),allChrom[indexAln].cm.end());
				}					
				
			}//qvm.size()>0
			vm.clear();
			qvm.clear();
		}//vm.size()>2
		count = 0; //reset count for the next alignment
		allChrom[indexAln].mums.clear();
	}
cout<<"Done with gap filling "<<endl;	
	fcords.close();
	//ftrans.open("trans.txt");		
	if(argv[5][0] == 'h')
	{
	for(map<string,vector<string> >::iterator it = hcp.begin(); it != hcp.end();it++)
	{
		refName = it->first;
		for(unsigned int i = 0; i<hcp[refName].size();i++)
		{
			fin.open(argv[1]);
			indexAln = hcp[refName][i];
			qName = allChrom[indexAln].cm[0].qn;
			sort(allChrom[indexAln].cm.begin(),allChrom[indexAln].cm.end());
			//xtracTrans(mRef[refName],allChrom[indexAln].cm,ftrans);//sort it inside the function
			readUniq(fin,allChrom[indexAln].cm,umRef[refName],masterHQ[qName]);
			fin.close();
		}
	}
	}
	//ftrans.close();
	
	refFasta.open(argv[2]);//read in the reference fasta
	readfasta(refFasta,refseq);//load them into memory
	refFasta.close();
	qFasta.open(argv[3]);//read in the query fasta	
	readfasta(qFasta,qseq);
	qFasta.close();
	int id = 0;	
	fout.open("sv.txt");
	fcnv.open("cnv_all.txt");
	fsmall.open("small.txt");
	findel.open("indel.txt");
	fout<<"REF_CHROM\tREF_START\tREF_END\tSV_TYPE\tQ_CHROM\tQ_START\tQ_END\tID\tLEN\tCOV_REF\tCOV_Q"<<endl;
	for(map<string,vector<string> >::iterator it = hcp.begin(); it != hcp.end();it++)
	{
		refName = it->first;
		for(unsigned int i = 0; i<hcp[refName].size();i++)
		{		
			indexAln = hcp[refName][i];
			qName = allChrom[indexAln].mums[i].qn;
			sort(allChrom[indexAln].cm.begin(),allChrom[indexAln].cm.end());
			sort(allChrom[indexAln].ncm.begin(),allChrom[indexAln].ncm.end());
			allChrom[indexAln].gap.clear();//flushing the gaps vector
			vmi = allChrom[indexAln].ncm;
			annotGaps(allChrom[indexAln].cm,masterRef[refName],masterQ[qName],vmi,umRef[refName],refseq[refName],qseq[qName],seqLen[indexAln],fout,fsmall,id);			
			
			
		}
	}
	fout.close();
	fcnv.close();	
	fsmall.close();
	findel.close();
	
return 0;
}
			


