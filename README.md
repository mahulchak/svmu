# svmu

The long molecule technologies have ushered in the era of extremely contiguous genome assemblies. These assemblies are empowering us to detect all structural variants present in a genome, a feat that can be accomplished via alignment of two  genome assemblies. Towards that goal, we present SVMU (Structural Variants from MUmmer): a pipeline to call duplicates,indels, and inversions from whole genome alignments using MUMmer. <b>It is still under development</b>. For questions and comments, please email me at mchakrab@uci.edu. 

Note: svmu is not a read mapping based copy number variation detection tool. It calls duplicates and TEs from alignments between two genome assemblies. A manuscript describing SVMU is under preparation but if you are publishing results obtained with this pipeline, please cite SVMU as "https://github.com/mahulchak/svmu".

Download and compile the programs -

 ```
 	g++ -Wall -std=c++0x fasplitter.cpp -o fasplitter
	g++ -Wall -std=c++0x mlib.cpp svmu.cpp -o svmu
	g++ -Wall -std=c++0x script_maker.cpp -o scriptmaker
	g++ -Wall -std=c++0x cnvlib.cpp ccnv.cpp -o checkCNV
 ```

Other programs needed for this pipeline:

  * <a href="http://mummer.sourceforge.net/">MUMmer</a>,  and <a href="https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download/"> BLAST</a>. Additionally, a program is needed to split the chromosomes into individual fasta files (supplied here). BLAST and MUMmer should be in your path.

  * Reference and a query genome assembly are needed in fasta file format. The program will report the sequences that are n copy in the reference genome but >n copy in the query genome.

Here is an example of how to use svmu pipeline to obtain a list of duplicates sequences from whole gnome alignment.

1. Add a prefix to the sequence names to avoid a bug in BLAST.
 ```
	sed 's/^>/>svmu/g' foo.fasta
 ```
 Run this on your assembly and the reference assembly.

1. Run fasplitter on the reference genome to split the contigs/chromosomes/scaffolds into component fasta files. This is necessary because mummer can be run on a single reference sequence at a time.

 ``` 
	./fasplitter reference_assembly.fasta Y
 ```
Using the 'Y' switch in fasplitter will ensure that the new fasta files have '.fa' in their names.

2. Create a list of the new fasta files

 ```
	ls *.fa > list_of_fa

 ```

  TIP: to avoid listing the reference assembly fasta file in the list, use ".fasta" as the file extension for the reference assembly. 
  This list can be edited to remove chromosomes/contigs you don't want to use for duplicate finding.
 
3. Run scriptmaker to generate duplicate calling scripts for all the component fasta files.

 ```   
	./scriptmaker -q your_assembly.fasta -l list_of_ref_fasta -s cluster_separation -ln minimum_length
 ```
  cluster separation is same as the -s option in mgaps (use 1000 for divergent species and 200-600 for individuals from the same strain). -ln denotes the minimum length of duplicates you want to detect.

4. Generate the list of all scripts and then add bash command to each of them.

 ```
	ls job_* > list_of_jobs

	sed -i 's/^/bash /g' list_of_jobs

 ```

5. Now run all the job scripts in parallel or serial mode depending on how many cores you have. If you are running in serial mode (i.e. want to use single processor) , do -

 ```
	bash list_of_jobs
 ```
 For using the parallel mode, you can use <a href="http://www.gnu.org/software/parallel/">GNU parallel</a> -

 ```
	cat list_of_jobs | parallel -j NPROC
 ```
 TIP: Replace 'NPROC' with the number of processors you want to use. If you have 4 processors, use 4.


6. Concatenate all the output files.

 ```
	cat *.tsv > all_chrom.tsv
 ```

7. The output tsv file has the following columns -
 
 ```
	REF_NAME REF_ST REF_END Q_NAME1 Q_ST1 Q_END1 Q_NAME2 Q_ST2 Q_END2 
 ```
  REF_NAME:reference chromosome where the parental gene is located.

  REF_ST: start coordinate of the reference genome segment that is duplicated.

  REF_END: end coordinate of the reference genome segment that is duplicated.

  Q_NAME1, Q_ST1, Q_END1: chromosome name, start and end coordinates of one of the copies in the query genome. 

  Q_NAME2, Q_ST2, Q_END2: chromosome name, start and end coordinates of the second copy in the query genome.
  
  

8. The 'all_chrom.tsv' file has both TE and duplicates in it. It also has false positives because ancestral duplicates may often be marked as polymorphic duplicate by svmu. These false positives are inevitable because aligning two genomes, especially at complex repeat regions, is not always perfect. To filter the false positives, we will use nucmer and a second program called checkCN. Next, extract reference names and reference coordinates from "all_chrom.tsv" to create a list of coordinates that can be mined using BLAST. You can use other tools or your own script.
 
 ```
	awk '{print $1" "$2"-"$3}' all_chrom.tsv|sort -k1,1 -k2,2n -u > blast.query.list
 ```

 TIP: If you wish to save on computation time, you can run bedtools merge (e.g. bedtools merge -d -100) on the "all_chrom.tsv" to merge overlapping coordinates before you create the blast query list.


9. Obtain the sequences corresponding to these reference coordinates using BLAST. Please see the BLAST manual to create a local BLAST database for your genome assemblies.

 ```
	blastdbcmd -db REF_BLAST_DB -dbtype nucl -entry_batch blast.query.list -out nucmer.query.fasta
 ```
 
10. Align all sequences in nucmer.query.fasta to the query genome and the reference genome using nucmer.

 ```
	nucmer -mumreference -prefix out.q nucmer.query.fasta your_assembly.fasta
	
	nucmer -mumreference -prefix out.r nucmer.query.fasta reference_assembly.fasta
 ```

11. Reformat the blast.query.list and then check the differences in copy number for each sequence using the program checkCN.

 ```
	sed -i 's/ /:/g' blast.query.list
	
	./checkCNV -d1 out.q.delta -d2 out.r.delta -q blast.query.list -c cutoff_for_repeats -qco cutoff_query_merging -rco cutoff_ref_merging
 ```
	qco determines how much divergence to allow within a duplicated query sequence. E.g. a 10Kb TE insertion can be accomodated by setting qco as 10000.
	rco determines how much divergence to allow within a duplicated sequence in the reference. it uses fractions and 0.05 typically works well.

 <u>Description of the "cnv_report.tsv" file</u>
 Column 1 = name of the CNV with the reference sequence name and the coordinates
 Column 2 = name of the query chromosome/contig
 Column 3 = start coordinate in the query chromosome
 Column 4 = end coordinate in the query chromosome
 Column 5 = copy number in the query chromosome
 Column 6 = copy number in the reference chromosome
 Column 7 = orientation of the copy in the query genome
 
Now you can use your own cutoff for copy number to separate  duplicates from TEs or you can use TE annotations in your reference genome to identify the TEs.
