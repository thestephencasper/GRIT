#!/usr/bin/env bash

# Before running this shell, the UNIX command line utilities for Entrez must be installed. See.
# https://www.ncbi.nlm.nih.gov/books/NBK179288/

# For info, see https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.39/

echo Please be patient. This may take a while. Several GB of data are downloading....

echo chr1...
esearch -db nucleotide -query "GenBank: CM000663.2" | efetch -format fasta > chr1.fasta  # chromosome 1
esearch -db nucleotide -query "NC_000001.11" | efetch -format fasta_cds_aa > chr1_cds.fasta  # chromosome 1 cds

echo chr2...
esearch -db nucleotide -query "GenBank: CM000664.2" | efetch -format fasta > chr2.fasta  # chromosome 2
esearch -db nucleotide -query "NC_000002.12" | efetch -format fasta_cds_aa > chr2_cds.fasta  # chromosome 2 cds

echo chr3...
esearch -db nucleotide -query "GenBank: CM000665.2" | efetch -format fasta > chr3.fasta  # chromosome 3
esearch -db nucleotide -query "NC_000003.12" | efetch -format fasta_cds_aa > chr3_cds.fasta  # chromosome 3 cds

echo chr4...
esearch -db nucleotide -query "GenBank: CM000666.2" | efetch -format fasta > chr4.fasta  # chromosome 4
esearch -db nucleotide -query "NC_000004.12" | efetch -format fasta_cds_aa > chr4_cds.fasta  # chromosome 4 cds

echo chr5...
esearch -db nucleotide -query "GenBank: CM000667.2" | efetch -format fasta > chr5.fasta  # chromosome 5
esearch -db nucleotide -query "NC_000005.10" | efetch -format fasta_cds_aa > chr5_cds.fasta  # chromosome 5 cds

echo chr6...
esearch -db nucleotide -query "GenBank: CM000668.2" | efetch -format fasta > chr6.fasta  # chromosome 6
esearch -db nucleotide -query "NC_000006.12" | efetch -format fasta_cds_aa > chr6_cds.fasta  # chromosome 6 cds

echo chr7...
esearch -db nucleotide -query "GenBank: CM000669.2" | efetch -format fasta > chr7.fasta  # chromosome 7
esearch -db nucleotide -query "NC_000007.14" | efetch -format fasta_cds_aa > chr7_cds.fasta  # chromosome 7 cds

echo chr8...
esearch -db nucleotide -query "GenBank: CM000670.2" | efetch -format fasta > chr8.fasta  # chromosome 8
esearch -db nucleotide -query "NC_000008.11" | efetch -format fasta_cds_aa > chr8_cds.fasta  # chromosome 8 cds

echo chr9...
esearch -db nucleotide -query "GenBank: CM000671.2" | efetch -format fasta > chr9.fasta  # chromosome 9
esearch -db nucleotide -query "NC_000009.12" | efetch -format fasta_cds_aa > chr9_cds.fasta  # chromosome 9 cds

echo chr10...
esearch -db nucleotide -query "GenBank: CM000672.2" | efetch -format fasta > chr10.fasta  # chromosome 10
esearch -db nucleotide -query "NC_000010.11" | efetch -format fasta_cds_aa > chr10_cds.fasta  # chromosome 10 cds

echo chr11...
esearch -db nucleotide -query "GenBank: CM000673.2" | efetch -format fasta > chr11.fasta  # chromosome 11
esearch -db nucleotide -query "NC_000011.10" | efetch -format fasta_cds_aa > chr11_cds.fasta  # chromosome 11 cds

echo chr12...
esearch -db nucleotide -query "GenBank: CM000674.2" | efetch -format fasta > chr12.fasta  # chromosome 12
esearch -db nucleotide -query "NC_000012.12" | efetch -format fasta_cds_aa > chr12_cds.fasta  # chromosome 12 cds

echo chr13...
esearch -db nucleotide -query "GenBank: CM000675.2" | efetch -format fasta > chr13.fasta  # chromosome 13
esearch -db nucleotide -query "NC_000013.11" | efetch -format fasta_cds_aa > chr13_cds.fasta  # chromosome 13 cds

echo chr14...
esearch -db nucleotide -query "GenBank: CM000676.2" | efetch -format fasta > chr14.fasta  # chromosome 14
esearch -db nucleotide -query "NC_000014.9" | efetch -format fasta_cds_aa > chr14_cds.fasta  # chromosome 14 cds

echo chr15...
esearch -db nucleotide -query "GenBank: CM000677.2" | efetch -format fasta > chr15.fasta  # chromosome 15
esearch -db nucleotide -query "NC_000015.10" | efetch -format fasta_cds_aa > chr15_cds.fasta  # chromosome 15 cds

echo chr16...
esearch -db nucleotide -query "GenBank: CM000678.2" | efetch -format fasta > chr16.fasta  # chromosome 16
esearch -db nucleotide -query "NC_000016.10" | efetch -format fasta_cds_aa > chr16_cds.fasta  # chromosome 16 cds

echo chr17...
esearch -db nucleotide -query "GenBank: CM000679.2" | efetch -format fasta > chr17.fasta  # chromosome 17
esearch -db nucleotide -query "NC_000017.11" | efetch -format fasta_cds_aa > chr17_cds.fasta  # chromosome 17 cds

echo chr18...
esearch -db nucleotide -query "GenBank: CM000680.2" | efetch -format fasta > chr18.fasta  # chromosome 18
esearch -db nucleotide -query "NC_000018.10" | efetch -format fasta_cds_aa > chr18_cds.fasta  # chromosome 18 cds

echo chr19...
esearch -db nucleotide -query "GenBank: CM000681.2" | efetch -format fasta > chr19.fasta  # chromosome 19
esearch -db nucleotide -query "NC_000019.10" | efetch -format fasta_cds_aa > chr19_cds.fasta  # chromosome 19 cds

echo chr20...
esearch -db nucleotide -query "GenBank: CM000682.2" | efetch -format fasta > chr20.fasta  # chromosome 20
esearch -db nucleotide -query "NC_000020.11" | efetch -format fasta_cds_aa > chr20_cds.fasta  # chromosome 20 cds

echo chr21...
esearch -db nucleotide -query "GenBank: CM000683.2" | efetch -format fasta > chr21.fasta  # chromosome 21
esearch -db nucleotide -query "NC_000021.9" | efetch -format fasta_cds_aa > chr21_cds.fasta  # chromosome 21 cds

echo chr22...
esearch -db nucleotide -query "GenBank: CM000684.2" | efetch -format fasta > chr22.fasta  # chromosome 22
esearch -db nucleotide -query "NC_000022.11" | efetch -format fasta_cds_aa > chr22_cds.fasta  # chromosome 22 cds

echo chrX...
esearch -db nucleotide -query "GenBank: CM000685.2" | efetch -format fasta > chrX.fasta  # chromosome X
esearch -db nucleotide -query "NC_000023.11" | efetch -format fasta_cds_aa > chrX_cds.fasta  # chromosome X cds

echo chrY...
esearch -db nucleotide -query "GenBank: CM000686.2" | efetch -format fasta > chrY.fasta  # chromosome Y
esearch -db nucleotide -query "NC_000024.10" | efetch -format fasta_cds_aa > chrY_cds.fasta  # chromosome Y cds
