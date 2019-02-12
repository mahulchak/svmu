# svmu

SVMU (Structural Variants from MUmmer) 0.2beta is part of our (Emerson lab at UCI) ongoing efforts to identify comprehensive sequence variants via alignment of two  contiguous genome assemblies. It calls SNPs, small indels, as well as duplicates, large indels, and inversions from whole genome alignments using MUMmer. 
<b>It is still under development</b>. If you encounter an issue, please email me at mchakrab@uci.edu. This version replaces the earlier version of svmu. SVMU works with both MUMmer v3.23 and MUMmer v4.x.

NOTE: SVMU is currently in the middle of a big change, which will reduce memory usage drastically and increase accuracy of SV detection. Feel free to try the new version but be cautious with the results. However, if you are coming here looking for the svmu versions used in the A4 and DSPR papers, see below:

If you publish results obtained with this pipeline, please cite SVMU as described here https://www.nature.com/articles/s41588-017-0010-y. the version used in this paper can be found here : https://github.com/mahulchak/svmu/releases/tag/v0.1beta. If you looking for the version that was used in the DSPR paper, please download commits prior to March 6,2018.

1. Download and compile the programs -

 ```
	make

 ```

2. Align the reference and your sample genomes using nucmer: 

 ```
	nucmer -maxmatch -prefix sam2ref.mm ref.fasta sample.fasta
	
 ```
Older version of svmu (prior to Mar 6, 2018) had a high memory footprint so if your svmu run crashes due to memory, you can try to run nucmer as follows (or use the latest,buggy version) -
 ```
	nucmer -mumreference -noextend -prefix sam2ref.mr ref.fasta sample.fasta

 ```

3. Run svmu on the delta file.

 ```
	svmu sam2ref.mm.delta ref.fasta sample.fasta n snp_mode> sample.small.txt

 ```
 <b>n</b> represents the number of unique mum/syntenic blocks that should be present between two sequences to find the SVs between them. It can be 5, or 10, or 100. We will get rid of this parameter very soon. 
  
 <b>snp_mode</b> should be 'h' or 'l', depending on whether you want to get the SNPs and small indels or not. 'h' will lead to higher memory usage and longer run times. The program generates several files as output: 

	sv.txt: A tab delimited file that summarizes structural mutations (indels, CNVs, inversions) in the sample genome with respect to the reference genome.  

	small.txt: A tab delimited file containing SNPs and small indels that occur within syntenic blocks (or MUMs).

	cnv_all.txt: A tab delimited file with all the reference genomic regions that are present in higher copy numbers (>1) in the sample genome. Those with "trans" in their names mean either it is a transposable element or non-TE copies of a gene in different chromosomes.

	cm.txt: A bed file with the reference genomic regions that have been putatively translocated (may not include TEs). 

This is work in progress so do examine the output. Final goal is to facilitate variant detection from a population samples of high quality genomes. If you have an idea or suggestion (including collaboration ideas, write to me).


Finally, If you are using SVMU for your research, please keep in mind that SVMU has not been extensively tested on genomes bigger than Drosophila. So there is no gurantee that it will work well with other genomes. Currently it requires ~2.5G memory for the <i>D. melanogaster</i> genome.

KNOWN ISSUES:
1. If the reverse complementary strand of a chromosome was sequenced (relative to the reference), svmu will identify whole sequence as 'INV' or inverted. A simple workaround is to do grep -v 'INV' on your sv.txt file. A future fix will take care of this issue.
2. White space in fasta headers will cause segfault in svmu because nucmer strips all text following white space or tab present in the fasta headers.

