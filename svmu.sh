./fasplitter reference_assembly.fasta Y

ls *.fa > list_of_fa

./scriptmaker your_assembly.fasta list_of_fa

ls job_* > list_of_jobs

sed -i 's/^/bash /g' list_of_jobs

bash list_of_jobs

cat list_of_jobs | parallel -j NPROC

cat *.tsv > all_chrom.tsv

awk '{print $1" "$2"-"$3}' all_chrom.tsv|sort -k1,1 -k2,2n -u > input_for_blast

[add a line to generate blast database for the genome]

blastdbcmd -db REF_BLAST_DB -dbtype nucl -entry_batch input_for_blast -out nucmer.query.fasta

sed 's/ /:/g' input_for_blast > input_for_nucmer

nucmer -maxmatch -prefix out.q nucmer.query.fasta [your_assembly.fasta]

nucmer -maxmatch -prefix out.r nucmer.query.fasta [reference_assembly.fasta]

./checkCN -d1 out.q.delta -d2 out.r.delta -q blast.query.list -c [cutoff] > my_cnvs.$CUTOFF

nucmer -mumreference -l 100 [reference_assembly.fasta] [your_assembly.fasta]
delta-filter -m out.delta > out.m.delta
./findDel -d out.m.delta > dels.in.[your_assembly]
