# svMum

svMum is a pipeline to call duplicates from whole genome alignments using mummer. The pipeline is under development. For issues and bugs, please contact me at mchakrab@uci.edu.

Download and compile the programs -

 ```
	g++ -Wall -std=c++0x mlib.cpp svmum.cpp -o svmum
	g++ -Wall script_maker.cpp -o scriptmaker
 ```

Other programs needed for this pipeline:

  * You will need <a href="http://mummer.sourceforge.net/">MUMmer</a>,  <a href="https://github.com/arq5x/bedtools2/blob/master/README.md">bedtools</a> , and <a href="http://www.repeatmasker.org/"> Repeatmasker</a> to use this pipeline. Additionally, the program fasplitter from <a href = "https://github.com/mahulchak/Assembly-utils">Assembly-utils</a> is required to split the fasta file.

  * You need to have a reference and a query genome assembly (in fasta file format). If desired, the assemblies could be processed through repeatmasker before running the pipeline. The program will report the sequences that are single copy in the reference genome but >1 copy in the query genome.

Here is an example of how to use svMUM pipeline to obtain a list of duplicates sequences from whole gnome alignment.

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
  Edit the list if you don't want to use certain sequences for SV detection.E.g. you could ignore sequences that are < 50 kb long and hence remove sequence names which correspond to sequences <50kb.
 
3. Run scriptmaker to generate duplicate calling scripts for all the component fasta files.

 ```   
	./scriptmaker your_assembly.fasta list_of_fa
 ```

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
 Replace 'NPROC' with the number of processors you want to use. If you have 4 processors, use 4.

6. Concatenate all the output files.

 ```
	cat *.tsv > all_chrom.tsv
 ```

7. The output tsv file has the following columns -
 
 ```
	REF_NAME REF_ST REF_END Q_NAME1 Q_ST1 Q_END1 Q_ST2 Q_NAME2 Q_ST2 Q_END2
 ```
  REF_NAME:reference chromosome where the parental gene is located.

  REF_ST: start coordinate of the reference genome segment that is duplicated.

  REF_END: end coordinate of the reference genome segment that is duplicated.

  Q_NAME1, Q_ST1, Q_END1: chromosome name, start and end coordinates of one of the copies in the query genome. 

  Q_NAME2, Q_ST2, Q_END2: chromosome name, start and end coordinates of the second copy in the query genome.

8. The 'all_chrom.tsv' file has both TE (if repeatmasker was not used on the assembly) and duplicates in it. To separate TE from duplicates, you will need <a href="https://github.com/arq5x/bedtools2/blob/master/README.md">bedtools</a> and a file with TE annotations for the reference genome.
 
 TIP: You can use Repeatmasker to generate the TE annotation file if you already don't have a TE annotation file. The TE annotation file needs to be in this format:
 
 ```
	CHROM_NAME	TE_START	TE_END
 ```

 After you get the TE annotation file, you can use the following commands to obtain the list of sequences that are single copy in the reference genome but more than one copy in the other genome.

 ```
	cat all_chrom.tsv| awk '{if($10-$9 >100) print $0}' | sort -u -k9,9 |sort -k1,1 -k2,2n >name_all_singleton

	cat name_all_singleton | sort -k1,1 -k2,2n |bedtools subtract -A -a stdin -b TE.bed > list_of_all_dups

	awk '{ if ($3-$2>100 && $4 == $7) print $0}' name_all_singleton | sort -k1,1 -k2,2n |bedtools subtract -A -a stdin -b TE.bed |bedtools merge -i stdin >list_of_TD
 ```
 
  Replace 'name' with a unique identifier for the query genome.

  Explanations of the avove mentioned steps:

   * In the first step, only duplicates longer than 100 bp are kept. Only the unique duplicates are kept.
   * The second command filters all TEs from the duplicate calls and reports all duplicates.
   * The third command does the same thing as the second command but reports only duplicates within same chromosome or contig.

  These are provided only as examples. Other commands could also be used to partition the duplicate calls.



 

