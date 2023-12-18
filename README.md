# Split_Chr_gff
##Script to split large chromosome sequence into sub-chromosomes (500Mb), then output .fasta sequence file and updated genes anno .gff3 file of the genome with sub-chromosomes.

# 1. ENV creates and activate
conda create -n genome_split -c bioconda biopython bcbio-gff
conda activate genome_split

# 2. Split chromosome and output new gff
python split.py <fasta_file> <gff_file> <output_prefix> <gff_output_file>
