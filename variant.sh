#!/bin/bash

# -------- readme
# Name of script: rnxtract.sh
# Purpose: a pipeline to extract a gene of interest from cDNA
# Dependencies: hisat2, samtools, bcftools, and EMBOSS
# Author: Dylan O. Burge
# Contact: dylan.o.burge@gmail.com
# Last modified: 20 April 2017

# -------- setup
# make a directory for your run: mkdir rnxtract
# go to this directory: cd rnxtract

# required software accessible from here: hisat2, samtools, bcftools, EMBOSS.
# if you are on a Mac, the package manager Homebrew will install these globally.
# find your reference genome (either cDNA or nuclear genome); must be .fasta.
# below, I call this "genome."
# alternatively, index using hisat2: hisat2-build genome.fasta root_name
# find your gene(s) of interest in this genome and define then in a BED file.
# give a FASTA nucleotide file with sequences that match the BED file.
# get a set of paired-end reads based on cDNA sequencing; must be .fastq.
# put the reference genome, sequencing reads, BAM, and FASTA into rnxtract.

# -------- definitions

# define the genes of interest (use your own file in place of genes.bed)
        bed=genes.bed
# define the fasta file for same genes (use your file name in place of genes.fasta):
        fasta=genes.fasta
# define sequencing reads (use your names in place of "forward" and "reverse"):
        reads1=forward.fastq
        reads2=reverse.fastq

# -------- usage
# script is built to run interactively, but in a single run for a set of reads.
# make script executable by running this: chmod +x rnxtract.sh
# run the script by typing the following into the terminal: ./rnxtract.sh

# -------- script (no need to edit this)
# use hisat2 to map the reads to your reference genome (takes some time)
        echo "mapping the genome"
        hisat2 -x genome -1 $reads1 -2 $reads2 --no-unal -p 4 -S organism.sam
        echo "finished mapping the genome"
# use samtools to create intermediate files based on region of interest:
        echo "begin gene-specific work"
        samtools view -b -L $bed organism.sam > organism.bam
        samtools sort -o sorted.organism.bam -O bam organism.bam
        samtools index sorted.organism.bam
# use samtools mpileup to get likelihoods
        samtools mpileup -uDl $bed -f genome.fasta sorted.organism.bam -o likelihood.bcf
# use bcftools to call variants against genes in the bam file:
        bcftools call likelihood.bcf -m -Ov -o organism.vcf
        bgzip organism.vcf
        tabix organism.vcf.gz
# use bcftools to get a consensus:
        bcftools consensus -f $fasta organism.vcf.gz > organism_nuc.fasta
# use transeq from EMBOSS to translate protein using first reading frame:
        transeq organism_nuc.fasta organism_protein.fasta

echo "done getting your genes of interest"
