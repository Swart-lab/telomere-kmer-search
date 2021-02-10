#!/bin/bash
set -e

# filename prefix for output files
PREFIX=prefix
# number of processors for jellyfish
CPU=12
# k-mer length for counting; should be prime
KMERLEN=19
# min coverage for k-mers to use for repeat finding
MINCOV=100
# filename for generators file
GENERATORS=${PREFIX}.generators

if [ ${#@} == 0 ]; then
  # no positional arguments, help message
  echo "List names of fastq.gz files to be processed as arguments after script name"
  echo "Example:"
  echo "  bash $0 file1.fastq.gz file2.fastq.gz"
else
  # Generate "generators" file, listing gzip decompress commands for the input
  # fastq.gz files
  for FQGZ_FILE in $@; do
    echo "gzip -cd $FQGZ_FILE" >> $GENERATORS
  done
  # Count k-mers; memory and CPU requirements may need to be adjusted depending
  # on your data and available computing resources
  jellyfish count
    -g $GENERATORS -G 4 -m $KMERLEN -s 250M -t $CPU --bf-size 100G -C \
    -o ${PREFIX}.k${KMERLEN}.jf

  # Dump k-mer sequences and counts to a Fasta file
  jellyfish dump ${PREFIX}.k${KMERLEN}.jf > ${PREFIX}.k${KMERLEN}.dumps.fa

  # Discard singletons and low-coverage k-mers
  # Cutoff should be slightly over the expected average coverage of the genome,
  # to retain only repetitive elements
  ./discard_singletons.py --input ${PREFIX}.k${KMERLEN}.dumps.fa \
    --cutoff $MINCOV \
    --histo ${PREFIX}.k${KMERLEN}.hist.json \
    --output ${PREFIX}.k${KMERLEN}.cov${MINCOV}.json
fi
