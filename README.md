# svmu

SVMU (Structural Variants from MUmmer) 0.3 is a revised version of SVMU 0.2. SVMU attempts to identify comprehensive sequence variants via alignment of two  contiguous genome assemblies. It calls duplicates, large indels, inversions, small indels, SNPs from whole genome alignments using MUMmer. 
<b>It is still under active development. We are incorporating new features and fixing bugs. One major feature we are planning to add is incorporation of LASTZ output in svmu. So please bear with us during this transition</b>. If you encounter an issue, feel free to email me at mchakrab@uci.edu. SVMU works with both MUMmer v3.23 and MUMmer v4.x. Support for LASTZ is coming soon.

NOTE: Feel free to try the newest version but be cautious with the results. If you have an idea or suggestion (including collaboration ideas), write to me. However, if you are coming here looking for the svmu versions used in the A4 and DSPR papers, see below:

If you publish results obtained with this pipeline, please cite SVMU as described here https://www.nature.com/articles/s41588-017-0010-y. the version used in this paper can be found here : https://github.com/mahulchak/svmu/releases/tag/v0.1beta. If you looking for the version that was used in the DSPR paper, please download the commits prior to March 6,2018.

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
	svmu sam2ref.mm.delta ref.fasta sample.fasta n snp_mode> sample.small.txt

 ```
 <b>n</b> represents the number of unique mum/syntenic blocks that should be present between two sequences to find the SVs between them. It can be 5, or 10, or 100. 
  
 <b>snp_mode</b> should be 'h' or 'l', depending on whether you want to get the SNPs and small indels or not. 'h' will lead to higher memory usage and longer run times. The program generates several files as output: 

	sv.txt: A tab delimited file that summarizes structural mutations (indels, CNVs, inversions) in the sample genome with respect to the reference genome.  

	small.txt: A tab delimited file containing SNPs and small indels that occur within syntenic blocks (or MUMs).

	cnv_all.txt: A tab delimited file with all the reference genomic regions that are present in higher copy numbers (>1) in the sample genome. Those with "trans" in their names mean either it is a transposable element or non-TE copies of a gene in different chromosomes.

	cm.txt: A bed file with the reference genomic regions that have been putatively translocated (may not include TEs). 

Finally, If you are using SVMU for your research, please keep in mind that SVMU has not been extensively tested on genomes bigger than Drosophila. So there is no gurantee that it will work well with other genomes. Currently it requires ~2.5G memory for the <i>D. melanogaster</i> genome.

KNOWN BUGS/BUGS WE ARE FIXING:
1. If the reverse complementary strand of a chromosome was sequenced (relative to the reference), svmu will identify whole sequence as 'INV' or inverted. A simple workaround is to do grep -v 'INV' on your sv.txt file. A future fix will take care of this issue.
2. White space in fasta headers will cause segfault in svmu because nucmer strips all text following white space or tab present in the fasta headers.
3. Translocated segments may show up as large indels. We are working on an experimental fix.
4. NUCMER is not very sensitive at divergent regions of the genome so we are bringing LASTZ into SVMU. MUMmer is not going away soon. LASTZ will complement the alignments NUCMER finds.

