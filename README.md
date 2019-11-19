# svmu

SVMU (Structural Variants from MUmmer) attempts to identify comprehensive sequence variants via alignment of two  contiguous genome assemblies. It combines the strengths of two powerful aligners MUMmer and LASTZ to annotate duplicates, large indels, inversions, small indels, SNPs from whole genome alignments.
 
<b>It is still under active development. We are incorporating new features and fixing bugs. One experimental feature we have added is processing of LASTZ output in svmu</b>. If you encounter an issue, feel free to email me at mchakrab@uci.edu. SVMU works with both MUMmer v3 and MUMmer v4. Support for LASTZ is experimental.

NOTE: Feel free to try the latest version but be cautious with the results. If you have an idea or suggestion (including collaboration ideas), write to me. However, if you are coming here looking for the svmu versions used in the A4 and DSPR papers, see below:

If you publish results obtained with this pipeline, please cite SVMU as described here https://www.nature.com/articles/s41467-019-12884-1. The version used in the paper are available through commits prior to March 6,2018.

1. Download and compile the programs -

 ```
	make

 ```

2. Obtain the mummer and lastz alignments: 

 ```
	nucmer --maxmatch --prefix sam2ref.mm ref.fasta sample.fasta

	lastz ref.fasta[multiple] sample.fasta[multiple] --chain --format=general:name1,strand1,start1,end1,name2,strand2,start2,end2 > sam_lastz.txt
	
 ```
	(LASTZ output that svmu reads should have only six columns as mentioned in the lastz command above)
	(--maxmatch can be omitted if two large genomes are being compared)

3. Run svmu:

 ```
	svmu sam2ref.mm.delta ref.fasta sample.fasta snp_mode sam_lastz.txt prefix 

 ```
 <b>snp_mode</b> should be 'h' or 'l'. h = report SNPs; l = no SNPs. currently, this option is turned off [will be activated soon].

 <b>prefix</b> provides a prefix that will be added to the output files. 

	sv.txt: A tab delimited file that summarizes structural mutations (indels, CNVs, inversions) in the sample genome with respect to the reference genome.  

	small.txt: A tab delimited file containing SNPs and small indels that occur within syntenic blocks (or MUMs).

	cnv_all.txt: A tab delimited file with all the reference genomic regions that are present in higher copy numbers (>1) in the sample genome. Those with "trans" in their names mean either it is a transposable element or non-TE copies of a gene in different chromosomes.

	cm.txt: A bed file with the reference genomic regions that are syntenic between the two genomes. 

Finally, If you are using SVMU for your research, please keep in mind that SVMU has not been extensively tested on genomes bigger than Drosophila. So there is no gurantee that it will work well with other genomes. Currently it requires ~2.5G memory for the <i>D. melanogaster</i> genome.

KNOWN BUGS/PLANNED FUTURE IMPROVEMENTS:
1. SVMU currently reports the inversion breakpoints and may not report the length of the inversion. Do inspect the reported inversions before you trust them fully. A future fix will take care of this issue.
2. White space in fasta headers will cause segfault in svmu because nucmer strips all text following white space or tab present in the fasta headers.
3. Translocated segments may show up as large indels. 
4. A lastz mode will be introduced that will use only lastz output to report SVs.
5. a sam to delta format converter is coming soon. That way, svmu will be compatible with minimap2.
