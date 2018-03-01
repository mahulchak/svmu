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
	map<string,string> refseq;
	map<string,string> qseq;
	map<string,vector<int> > seqLen;//length of sequences.first element is ref and second is query
	map<string,bool> qStrand; //stores whether query strand is forward strand or reverse strand
	mI tempmi,gapmi,prevmi;

	string foo = string(argv[1]);
	string line, chromName,refName,qName,indexAln;
	int refStart = 0, refEnd = 0, qStart = 0, qEnd = 0, refLen =0, qLen =0, count = -1,indelPos =0, forCount =0, revCount = 0;
	unsigned int cutoff = 0;
	cutoff = stoi(argv[4]);
	vector<double> vd;
	vector<int> vi;
	vector<mI> vmi,tempVmi;
	size_t pos1,pos2,namePos;
	
	ifstream fin, refFasta, qFasta,fcm;
	ofstream fout,fcnv,fsmall,ftrans,findel,fcords;
	fin.open(argv[1]);
	fcords.open("cords.txt");
	fcm.open("cm.txt");
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
			
			indelPos = abs(stoi(line));
			refStart = refStart + indelPos;
			if(indelPos <0)
			{	
				refStart = refStart * (-1);

			}
			vi.push_back(refStart);
		
//cout<<refName<<"\t"<<indelPos<<" " <<refStart<<"\t"<<refEnd<<"\t"<<qName<<"\t"<<qStart<<"\t"<<qEnd<<endl;
			if(indelPos ==0) //reached the end of the indel description
			{
				tempmi.mv = vi;
				//allChrom[indexAln].mums.push_back(tempmi);
				storeCords(masterRef[refName],masterQ[qName],tempmi);
				storeCords(mRef[refName],tempmi,fcords);
				tempmi.mv.clear();//delete this?
				allChrom[indexAln].mums.push_back(tempmi);
//cout<<refName<<"\t"<<refStart<<"\t"<<refEnd<<"\t"<<qName<<"\t"<<qStart<<"\t"<<qEnd<<"\t"<<allChrom[indexAln].mums.size()<<endl;
				vi.clear();//reset it once its values are used
			}
				
			count++;
			
		}
		if((line.find('>') == string::npos) && (line.size() >10) && (refName != "")) //when describing alignment segments
		{
		
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

				count = 0;
	//			--refStart;//to count the mutation distance

		}
	}
	fin.close();
	for(chroms::iterator it = allChrom.begin();it!= allChrom.end();it++)
	{
		indexAln = it->first;
		sort(allChrom[indexAln].mums.begin(),allChrom[indexAln].mums.end());
		for(unsigned int i = 0; i< allChrom[indexAln].mums.size();i++)
		{
			tempmi = allChrom[indexAln].mums[i];
			vd = getCoverage(tempmi,masterRef[tempmi.rn],masterQ[tempmi.qn]);
			if((vd[0] <1.3) && (vd[1]<1.3))
			{
				if(tempmi.y1>tempmi.y2)
				{
					revCount++;
				}
				else
				{
					forCount++;
				}
			}
		}
		if(revCount > forCount)
		{
			qStrand[indexAln] = true;//yes it is reverse strand
		}
		else
		{
			qStrand[indexAln] = false;
		}
		forCount = 0;
		revCount = 0;
		for(unsigned int i= 0; i<allChrom[indexAln].mums.size();i++)
		{
			tempmi = allChrom[indexAln].mums[i];
			vd = getCoverage(tempmi,masterRef[tempmi.rn],masterQ[tempmi.qn]);
			if((vd[0] <1.3) && (vd[1]<1.3))
			{
				if(allChrom[indexAln].cm.size() == 0)
				{
					gapmi.rn = tempmi.rn;
					gapmi.qn = tempmi.qn;
					gapmi.x1 = 1;
					gapmi.x2 = tempmi.x1;
					if(qStrand[indexAln] == false)
					{
						gapmi.y1 = 1;
						gapmi.y2 = min(tempmi.y1,tempmi.y2);
						gapmi.c = 'f';
					}
					else	
					{
						gapmi.y1 = seqLen[indexAln][1];//the largest length. correct this after the next element
						gapmi.y2 = max(tempmi.y1,tempmi.y2);
						gapmi.c = 'r';
					}
					allChrom[indexAln].gap.push_back(gapmi);
				}
					
				if(allChrom[indexAln].cm.size() >0)//once more than one element has been entered
				{	
					prevmi = allChrom[indexAln].cm[allChrom[indexAln].cm.size() -1];//the previous mi
					refStart = allChrom[indexAln].cm[allChrom[indexAln].cm.size() -1].x2;
					//if((prevmi.y1 <prevmi.y2) && (tempmi.y1 <tempmi.y2))
					//{
					//	qStart = max(allChrom[indexAln].cm[allChrom[indexAln].cm.size() -1].y2,allChrom[indexAln].cm[allChrom[indexAln].cm.size() -1].y1);
					//}
					if((prevmi.y1 > prevmi.y2) && (tempmi.y1 > tempmi.y2))
					{
						qStart = min(allChrom[indexAln].cm[allChrom[indexAln].cm.size() -1].y2,allChrom[indexAln].cm[allChrom[indexAln].cm.size() -1].y1);
					}
					else
					{
						qStart = max(allChrom[indexAln].cm[allChrom[indexAln].cm.size() -1].y2,allChrom[indexAln].cm[allChrom[indexAln].cm.size() -1].y1);
					}
					gapmi.rn = tempmi.rn;
					gapmi.qn = tempmi.qn;
					if(refStart < tempmi.x1)
					{
						if(!((prevmi.y1 >prevmi.y2) && (tempmi.y1 >tempmi.y2))) //not both of them are inverted
						{	
							gapmi.x1 = refStart;
							gapmi.x2 = tempmi.x1;
							gapmi.y1 = qStart;
							gapmi.y2 = min(tempmi.y1,tempmi.y2);
							gapmi.c = 'f';
							if(gapmi.y2 - gapmi.y1 >0)
							{
								allChrom[indexAln].gap.push_back(gapmi);
							}
						}
						if((prevmi.y1 > prevmi.y2) && (tempmi.y1 > tempmi.y2))
						{
							gapmi.x1 = refStart;
							gapmi.x2 = tempmi.x1;
							gapmi.y1 = qStart;
							gapmi.y2 = max(tempmi.y1,tempmi.y2);
							gapmi.c = 'r';
							if(gapmi.y1 - gapmi.y2 >0)
							{
								allChrom[indexAln].gap.push_back(gapmi);
							}
						}
					}
					
				}
				allChrom[indexAln].cm.push_back(tempmi);
				count = count + tempmi.x2 - tempmi.x1; //keeping a count of total unique alignment
//				cout<<"cm\t"<<indexAln<<"\t"<<tempmi.x1<<"\t"<<tempmi.x2<<"\t"<<tempmi.y1<<"\t"<<tempmi.y2<<"\t"<<vd[0]<<"\t"<<vd[1]<<endl;
			}
			
			else
			{
				allChrom[indexAln].ncm.push_back(tempmi);
				
			}
		}
		allChrom[indexAln].mums.clear(); //free up the memory
		if(allChrom[indexAln].cm.size()>cutoff)
		{
			hcp[allChrom[indexAln].cm[0].rn].push_back(indexAln);//homologous alignment
			
			for(unsigned int i=0;i<allChrom[indexAln].gap.size();i++)
			{
				if(allChrom[indexAln].gap[i].c != 'r')
				{
					gapCloser(allChrom[indexAln].gap[i], allChrom[indexAln].ncm, allChrom[indexAln].cm);
				}
				if(allChrom[indexAln].gap[i].c == 'r')
				{
					gapCloserRev(allChrom[indexAln].gap[i], allChrom[indexAln].ncm, allChrom[indexAln].cm);
				}
			}
			//sort(allChrom[indexAln].cm.begin(),allChrom[indexAln].cm.end(),qusort);
			allChrom[indexAln].gap.clear();//free up memory;will create gaps again later
				
			sort(allChrom[indexAln].cm.begin(),allChrom[indexAln].cm.end());
			for(unsigned int i=0;i<allChrom[indexAln].cm.size()-1;i++)
			{
				tempmi = allChrom[indexAln].cm[i];
				fcm<<tempmi.rn<<"\t"<<tempmi.x1<<"\t"<<tempmi.x2<<"\t"<<tempmi.qn<<"\t"<<tempmi.y1<<"\t"<<tempmi.y2<<endl;
					
			}
					
		}

		allChrom[indexAln].mums.clear(); //free up the memory
		allChrom[indexAln].gap.clear();//free up memory;will create gaps again later
		count = 0; //reset count for the next alignment
	}
	fcm.close();
//cout<<"Done with gap filling "<<endl;	
	ftrans.open("trans.txt");	
	for(map<string,vector<string> >::iterator it = hcp.begin(); it != hcp.end();it++)
	{
		refName = it->first;
		for(unsigned int i = 0; i<hcp[refName].size();i++)
		{
			fin.open(argv[1]);
			indexAln = hcp[refName][i];
			qName = allChrom[indexAln].mums[i].qn;
			sort(allChrom[indexAln].cm.begin(),allChrom[indexAln].cm.end());
			xtracTrans(mRef[refName],allChrom[indexAln].cm,ftrans);//sort it inside the function
			readUniq(fin,allChrom[indexAln].cm,umRef[refName],masterHQ[qName]);
			fin.close();
		}
	}
	ftrans.close();
	
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
			if(argv[5][0] =='h')
			{
				splitByCoverageSen(allChrom[indexAln],masterRef[refName],masterQ[qName]);
			}
			else
			{
				splitByCoverage(allChrom[indexAln],masterRef[refName],masterQ[qName]);	
			}
			allChrom[indexAln].gap.clear();//flushing the gaps vector
			for(unsigned int j=0; j<allChrom[indexAln].cc.size();j++)
			{
				tempVmi = findQuery(mRef[refName],allChrom[indexAln].cc[j],masterRef[refName],masterQ[qName], masterHQ[qName]);
				if(tempVmi.size()>0)
				{
					vmi.insert(vmi.end(),tempVmi.begin(),tempVmi.end());
					for(unsigned int i=0; i< tempVmi.size()-1;i++)//last element of tempVmi has the coverage info
					{
						fcnv<<tempVmi[i].rn<<"\t"<<tempVmi[i].x1<<"\t"<<tempVmi[i].x2<<"\tCNV\t"<<tempVmi[i].qn<<"\t"<<tempVmi[i].y1<<"\t"<<tempVmi[i].y2<<"\t"<<setfill('0')<<setw(10)<<tempVmi[i].x1<<tempVmi[i].rn<<"\t"<<tempVmi[i].x2-tempVmi[i].x1<<"\t"<<tempVmi[tempVmi.size()-1].x1<<"\t"<<tempVmi[tempVmi.size()-1].x2<<"\t"<<tempVmi[tempVmi.size()-1].y1<<endl;
					}
				}
			}
			annotGaps(allChrom[indexAln].cm,mRef[refName],masterRef[refName],masterQ[qName],vmi,umRef[refName],refseq[refName],qseq[qName],seqLen[indexAln],fout,fsmall,id);			
			
			
		}
	}
	fout.close();
	fcnv.close();	
	fsmall.close();
	findel.close();
	fcords.close();
return 0;
}
			


