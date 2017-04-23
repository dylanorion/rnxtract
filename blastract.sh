#!/bin/bash

# -------- readme

# Name of script: blastract.sh
# Purpose: blast for genes of interest in a set of contigs and obtain sequences.
# Dependencies: blast, strip.py (python script)
# Author: Dylan O. Burge
# Contact: dylan.o.burge@gmail.com
# Last modified: 20 April 2017

# -------- setup

# make a directory for your run: mkdir blast
# go to this directory: cd rnxtract
# required software and scripts accessible from here: blast, strip.py (python script).
# make sure the python script is executable: chmod +x strip.py
# if you are on a Mac, the package manager Homebrew will install blast globally.
# find your reference genome (either cDNA or nuclear genome); must be .fasta.
# below, I call this "genome."

# -------- definitions

# define the query sequences (use your file name in place of query.fasta):

fasta=query.fasta

# -------- usage

# script is built to run interactively, but in a single run for a set of query sequences.
# make script executable by running this: chmod +x blastoff.sh
# run the script by typing the following into the terminal: ./blastout.sh

# -------- script (no need to edit this)

echo "start blast script"

# use blast to make a database out of the subject (set of contigs):

makeblastdb -in genome.fasta -dbtype nucl

# use blast to query the reference genome:

tblastn -query $fasta -db genome.fasta -evalue 0.0000000000000000000001 -best_hit_overhang 0.25 -max_hsps 1 -outfmt 6 -out results.tbln.tab

# manipulate table and output a list of proteins:

awk '$9<$10 {print $2"\t"$9"\t"$10}' results.tbln.tab > tmp.fwd.bed
awk '$10<$9 {print $2"\t"$10"\t"$9}' results.tbln.tab > tmp.rev.bed
cat tmp.fwd.bed tmp.rev.bed > concat.bed
sortBed -i concat.bed > sorted.bed
bedtools merge -i sorted.bed -d 1000 > merged.bed
awk '{print $1"\t"$2"\t"$3"\t"$3-$2}' merged.bed > print_merged.bed
sort -k4nr print_merged.bed > genes.bed
awk '{print $1}' genes.bed > genes.list
echo "done getting your genes of interest"
