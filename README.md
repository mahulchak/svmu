# svmu

SVMU (Structural Variants from MUmmer) 0.2beta is a byproduct of our (Emerson and Long labs at UCI) efforts to identify comprehensive sequence variants via alignment of two  genome assemblies. It calls SNPs, small indels, as well as duplicates, large indels, and inversions from whole genome alignments using MUMmer. 
<b>It is still under development</b>. If you encounter an issue, please email me at mchakrab@uci.edu. This version replaces the earlier version of svmu.

If you publish results obtained with this pipeline, please cite SVMU as described here https://www.nature.com/articles/s41588-017-0010-y.

1. Download and compile the programs -

 ```
	make

 ```

2. Align the reference and your sample genomes using nucmer: 

 ```
	nucmer -maxmatch --noextend -prefix sam2ref.mm ref.fasta sample.fasta
	
 ```
Unfortunately, svmu has a high memory footprint (we are working to reduce it) so if your svmu run crashes due to memory, run nucmer as follows -
 ```
	nucmer -mumreference --noextend -prefix sam2ref.mr ref.fasta sample.fasta

 ```

3. Run svmu on the delta file.

 ```
	svmu sam2ref.mm.delta ref.fasta sample.fasta n > sample.small.txt

 ```
  n represents the number of unique mum/syntenic blocks that should be present between two sequences to find the SVs between them. It can be 5, or 10, or 100 (a future update will likely get rid of this parameter). The program generates several files as output: 

	sv.txt: A bed file that summarizes structural mutations (indels, CNVs, inversions) in the sample genome with respect to the reference genome.  

	small.txt: A bed file containing SNPs and small indels that occur within syntenic blocks (or MUMs).

	cnv_all.txt: A bed file with all the reference genomic regions that are present in higher copy numbers (>1) in the sample genome. Those with "trans" in their names mean either it is a transposable element or non-TE copies of a gene in different chromosomes.

	indel.txt: A bed file with all the reference and quesry genomic regions that are missing in the other genome. 
	
	trans.txt: A bed file with the reference genomic regions that have been putatively translocated (may not include TEs). 

We are continuously working on this and other tools that facilitate variant detection from a population samples of high quality genomes. If you have an idea or suggestion (including collaboration ideas, write to me or J.J. Emerson or Anthony Long).
