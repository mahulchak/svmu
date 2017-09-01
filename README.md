# svmu

The long molecule technologies have ushered in the era of extremely contiguous genome assemblies. These assemblies are enabling us comprehensive detection of structural variants in a genome, a feat that can be accomplished via alignment of two  genome assemblies. Towards that goal, we present SVMU 0.2beta (Structural Variants from MUmmer): a program to call duplicates,indels, and inversions from whole genome alignments using MUMmer. <b>It is still under development</b>. If you encounter an issue, please email me at mchakrab@uci.edu. This version replaces the earlier version of svmu.

Note: svmu is not a read mapping based copy number variation detection tool. It calls duplicates and TEs from alignments between two genome assemblies. A manuscript describing SVMU is under preparation but if you are publishing results obtained with this pipeline, please cite SVMU as described here http://www.biorxiv.org/content/early/2017/03/08/114967.

1. Download and compile the programs -

 ```
	make

 ```

2. Align the reference and your sample genomes using nucmer: 

 ```
	nucmer -maxmatch -prefix sam2ref.mm ref.fasta sample.fasta
	
 ```

3. Run svmu on the delta file.

 ```
	svmu sam2ref.mm.delta ref.fasta sample.fasta n > sample.small.txt

 ```
  n represents the number of unique mum/syntenic blocks that should be present between two sequences to find the SVs between them. It can 5 or 10 10 or 100 (a new update will likely get rid of this parameter but bear with it for now). The program will generate "sv.txt" and "cnv_all.txt" files as output. The former contains all indels,inversions, and CNVs in the sample genome with respect to the reference genome. The later contains all CNV coordinates.



