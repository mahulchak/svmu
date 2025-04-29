# svmu

SVMU (Structural Variants from MUmmer) attempts to identify comprehensive sequence variants via the alignment of two  contiguous genome assemblies. Its goal is to combine the strengths of different aligners to annotate duplicates, large indels, inversions, small indels, and SNPs from whole genome alignments.
 
<b>SVMU is being re-developed by the Chakraborty Lab at Texas A&M. We are incorporating new features (e.g., a VCF output). We'll make sure to announce here when the new version is available</b>. Feel free to continue using the existing version until then. If you encounter an issue, send an email to mahul@tamu.edu. SVMU currently works with MUMmer; both versions 3 and 4 are supported.

NOTE: If you have an idea or suggestion, including collaboration ideas, please email us at the address above. However, if you are interested in the svmu versions used in the 2018 Nature Genetics (A4) and 2019 Nature Communications (DSPR) papers, see below:

If you publish results obtained with this pipeline, please cite SVMU as described here https://www.nature.com/articles/s41467-019-12884-1. The version used in the paper is available through commits made before March 6, 2018.

1. Download and compile the programs -

 ```
	make

 ```

2. Obtain the mummer and lastz alignments: 

 ```
	nucmer --threads n --prefix sam2ref ref.fasta sample.fasta

	lastz ref.fasta[multiple] sample.fasta[multiple] --chain --format=general:name1,strand1,start1,end1,name2,strand2,start2,end2 > sam_lastz.txt
	
 ```
	(LASTZ output that svmu reads should have only six columns as mentioned in the lastz command above)
	(for relatively small genomes, --maxmatch can also be used for nucmer)

3. Run svmu:

 ```
	svmu sam2ref.mm.delta ref.fasta sample.fasta snp_mode sam_lastz.txt prefix 

 ```
 <b>snp_mode</b> should be 'h' or 'l'. h = report SNPs; l = no SNPs. Currently, this option is turned off [will be activated in the near future].

 <b>prefix</b> = a prefix that will be added to the output files. 

	sv.prefix.txt = A tab-delimited file that summarizes structural mutations (indels, CNVs, inversions) in the sample genome with respect to the reference genome.  

	small.prefix.txt: A tab-delimited file containing SNPs and small indels that occur within syntenic blocks (or MUMs).

	cnv_all.prefix.txt: A tab-delimited file with all the reference genomic regions that are present in higher copy numbers (>1) in the sample genome. Those with "trans" in their names mean either that it is a transposable element or non-TE copies of a gene in different chromosomes.

	cm.prefix.txt: A bed file with the reference genomic regions that are syntenic between the two genomes. 

Finally, if you are using SVMU for your research, please keep in mind that SVMU has not been extensively tested on genomes bigger than Drosophila. So, there is no guarantee that it will work well with other genomes. Currently, it requires ~2.5G memory for the <i>D. melanogaster</i> genome.

KNOWN BUGS/PLANNED FUTURE IMPROVEMENTS:
1. SVMU currently reports the inversion breakpoints but may not report the length of the inversion. Inspect the reported inversions before you trust them fully. A future update will address this issue.
2. White space in FASTA headers will cause a segfault in SVMU because NUCmer discards all text following white space or tabs present in the FASTA headers.
3. Translocated segments may appear as large insertions or deletions (indels). 
4. Using lastz is not currently supported, but it may be added in the new version of svmu.
5. The new version of svmu will be compatible with minimap2.
