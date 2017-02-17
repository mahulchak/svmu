# svmu

The long molecule technologies have ushered in the era of extremely contiguous genome assemblies. These assemblies are empowering us to detect all structural variants present in a genome, a feat that can be accomplished via alignment of two  genome assemblies. Towards that goal, we present SVMU (Structural Variants from MUmmer): a pipeline to call duplicates,indels, and inversions from whole genome alignments using MUMmer. <b>It is still under development</b>. If you encounter an issue, please email me at mchakrab@uci.edu. 

Note: svmu is not a read mapping based copy number variation detection tool. It calls duplicates and TEs from alignments between two genome assemblies. A manuscript describing SVMU is under preparation but if you are publishing results obtained with this pipeline, please cite SVMU as "https://github.com/mahulchak/svmu".

Download and compile the programs -

 ```
	make

 ```

Other programs needed for this pipeline:

  * <a href="http://mummer.sourceforge.net/">MUMmer</a>,  <a href="https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download/"> BLAST</a> (Note: current version of checkCNV is compatible only with BLAST versions prior to 2.3), and <a href="http://www.repeatmasker.org/">Repeatmasker</a> (optional). BLAST and MUMmer should be in your path.

  * Reference and a query genome assembly are needed in fasta file format. The program will report the sequences that are n copy in the reference genome but >n copy in the query genome.

Here is an example of how to use svmu pipeline to obtain a list of duplicates sequences from whole gnome alignment.

1. Add a prefix to the sequence names to avoid a bug in BLAST.
 ```
	sed -i 's/^>/>svmu/g' foo.fasta
 ```
 Run this on your assembly and the reference assembly.

2. Run fasplitter on the reference genome to split the contigs/chromosomes/scaffolds into component fasta files. This is necessary because mummer can be run on a single reference sequence at a time.

 ``` 
	./fasplitter reference_assembly.fasta Y
 ```
Using the 'Y' switch in fasplitter will ensure that the new fasta files have '.fa' in their names.

3. Create a list of the new fasta files

 ```
	ls *.fa > list_of_fa

 ```

  TIP: to avoid listing the reference assembly fasta file in the list, use ".fasta" as the file extension for the reference assembly. 
  This list can be edited to remove chromosomes/contigs you don't want to use for duplicate finding.
 
4. Run scriptmaker to generate duplicate calling scripts for all the component fasta files.

 ```   
	./scriptmaker -q your_assembly.fasta -l list_of_ref_fasta -s cluster_separation -ln minimum_length
 ```
  cluster separation is same as the -s option in mgaps (use 1000 for divergent species and 200-600 for individuals from the same strain). -ln denotes the minimum length of duplicates you want to detect.

5. Generate the list of all scripts and then add bash command to each of them.

 ```
	ls job_* > list_of_jobs

	sed -i 's/^/bash /g' list_of_jobs

 ```

6. Now run all the job scripts in parallel or serial mode depending on how many cores you have. If you are running in serial mode (i.e. want to use single processor) , do -

 ```
	bash list_of_jobs
 ```
 For the parallel mode, you can use <a href="http://www.gnu.org/software/parallel/">GNU parallel</a> -

 ```
	cat list_of_jobs | parallel -j NPROC
 ```
 TIP: Replace 'NPROC' with the number of processors you want to use. If you have 4 processors, use 4.


7. Concatenate all the output files.

 ```
	cat *.tsv > all_chrom.tsv
 ```

8. The output tsv file has the following columns -
 
 ```
	REF_NAME REF_ST REF_END Q_NAME1 Q_ST1 Q_END1 Q_NAME2 Q_ST2 Q_END2 
 ```
  REF_NAME:reference chromosome where the parental sequence is located.

  REF_ST: start coordinate of the reference genome segment that is duplicated.

  REF_END: end coordinate of the reference genome segment that is duplicated.

  Q_NAME1, Q_ST1, Q_END1: chromosome name, start and end coordinates of one of the copies in the query genome. 

  Q_NAME2, Q_ST2, Q_END2: chromosome name, start and end coordinates of the second copy in the query genome.
  
  

9. The 'all_chrom.tsv' file has both TE and duplicates in it. It also has false positives because ancestral duplicates may often be marked as polymorphic duplicate by svmu. These false positives are inevitable because aligning two genomes, especially at complex repeat regions, is not always perfect. To remove TEs, use your own TE annotation file (repeats.te.bed) for the reference sequence or create one with Repeatmasker.  Reference names and reference coordinates are taken from "all_chrom.tsv" to create a list of coordinates that can be fetched using BLAST. You can use other tools or your own script.
 
 ```	
	cat SV_report.l20.* | cut -f1-3 | bedtools subtract -a stdin -b repeats.te.bed | sort -u | awk '{if($3-$2>100) print $0}' | sort -k1,1 -k2,2n | awk '{print $1" "$2"-"$3}' > input_for_blast
	cat SV_report.l20.* | cut -f1-3 | bedtools subtract -a stdin -b repeats.te.bed | sort -u | awk '{if($3-$2>100) print $0}' | sort -k1,1 -k2,2n | awk '{print $1":"$2"-"$3}' > input_for_nucmer
	
 ```

10. Obtain the sequences corresponding to these reference coordinates using BLAST. Please see the BLAST manual to create a local BLAST database for your genome assemblies.

 ```
	blastdbcmd -db REF_BLAST_DB -dbtype nucl -entry_batch blast.query.list -out nucmer.query.fasta
 ```
 
11. Align all sequences in nucmer.query.fasta to the query genome and the reference genome using nucmer.

 ```
	nucmer -maxmatch -g cluster_sep -prefix out.q nucmer.query.fasta your_assembly.fasta
	
	nucmer -maxmatch -g cluster_sep -prefix out.r nucmer.query.fasta reference_assembly.fasta
 ```

12. Next check the differences in copy number for each sequence using the program checkCNV.

 ```
	./checkCNV -d1 out.q.delta -d2 out.r.delta -q input_for_nucmer -c cutoff_for_repeats -qco cutoff_query_merging -rco cutoff_ref_merging

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

12. Remove svmu from sequence names.
 ```
	sed -i 's/svmu//g' cnv_report.tsv
 ``` 
13. <b>Insertion-Deletion</b>: Align "your_assembly.fasta" to the "reference_assembly.fasta" using nucmer.

 ```
	nucmer -mumreference -prefix q2r ref.fasta your_asm.fasta
	delta-filter -m q2r.delta > q2r.mdelta
	findDel indel -d q2r.mdelta -m I -p p> ins.in.a
 ```
	This will give you the list of insertions (or deletions in your assembly) in the reference genome. Run nucmer in reverse orientation, i.e. set "your_asm.fasta" as the reference for nucmer instead of query, to obtain the list of insertions (deletions in the reference) in your assembly.

	m =  tells the program that we are looking for insertions.

	p = proportion between the gaps betweem query and reference MUMs. 0.1 is a good value (i.e. if the query has a 100bp sequence and the reference has a 1000bp sequence, surrounded by syntenic blocks, an insertion event will be inferred for the reference).


