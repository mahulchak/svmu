# svmu

SVMU (Structural Variants from MUmmer) 0.2beta is a byproduct of our (Emerson and Long labs at UCI) efforts to identify comprehensive sequence variants via alignment of two  genome assemblies. It calls SNPs, small indels, as well as duplicates, large indels, and inversions from whole genome alignments using MUMmer. 
<b>It is still under development</b>. If you encounter an issue, please email me at mchakrab@uci.edu. This version replaces the earlier version of svmu.

If you publish results obtained with this pipeline, please cite SVMU as described here http://www.biorxiv.org/content/early/2017/03/08/114967.

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
  n represents the number of unique mum/syntenic blocks that should be present between two sequences to find the SVs between them. It can be 5, or 10, or 100 (a new update will likely get rid of this parameter but bear with it for now). The program will generate "sv.txt" and "cnv_all.txt" files as output. The former contains all indels,inversions, and CNVs in the sample genome with respect to the reference genome. The latter contains all CNV coordinates.


We are continuously working on this and other tools that facilitate variant detection from a population samples of high quality genomes. If you have an idea or suggestion (including collaboration ideas, write to me or J.J. Emerson or Anthony Long).
