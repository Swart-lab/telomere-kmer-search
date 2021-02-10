#!/usr/bin/env python3

import json
from Bio import Seq
from math import floor
from statistics import pstdev, mean
from collections import defaultdict
import argparse

parser = argparse.ArgumentParser("Identify repeat elements from k-mer counts")
parser.add_argument("--counts", "-c", help="Dict of counts per k-mer, in JSON format")
parser.add_argument("--output", "-o", help="Prefix for output files", default="repeat_info")
parser.add_argument("--maxzeroes", type=int, default=1, help="Maximum number zero-count expected k-mers to allow in reporting repeats")
args = parser.parse_args()

# functions -----------------------------------------------------------------------

def get_repeat_frame(seq: str, frame = 1):
    """Find tandem repeat length within sequence.
    
    Returns:
    tuple of repeat length (int) and repeat sequence (str)

    Edge cases:
    - If sequence is not repeating, repeat length equals length of input seq
    - If sequence is empty, repeat length is 1, repeat sequence is empty
    """
    # I love recursive functions
    if seq[0:-frame] == seq[frame:]:
        return(frame, seq[0:frame])
    else:
        return(get_repeat_frame(seq, frame+1))

def recurse_repeat_arr(seq:str, seqs = []):
    """Sequence of a repeat element in all reading frames.
    """
    if len(seqs) == len(seq):
        return(seqs)
    else:
        seqs.append(seq)
        return(recurse_repeat_arr(seq[-1]+seq[0:-1], seqs))

def get_canonical_repeat(seq: str):
    """For a tandem repeat element, report its canonical form.
    
    The canonical form of a repeat is the sequence in a reading frame and
    strand that minimizes the lexicographical sort order.
    """
    # Get all reading frames of the repeat, and their reverse complements
    seqs = recurse_repeat_arr(seq, [])
    seqs_rc = [Seq.reverse_complement(s) for s in seqs]
    # Report first after lexicographical sort
    seqs_all = seqs + seqs_rc
    seqs_all.sort()
    return(seqs_all[0])

def canonical_kmer(seq:str):
    """Report canonical version of a k-mer
    """
    ll = [seq, Seq.reverse_complement(seq)]
    ll.sort()
    return(ll[0])

def get_expected_kmers_from_repeat(seq: str, k: int, canonical=True):
    """For a tandem repeat element and a k-mer length, report all k-mers
    expected to be observed. k must be longer than repeat length.
    k is recommended to be prime
    """
    if k <= len(seq):
        return ValueError("k-mer length k must be longer than sequence length")
    else:
        seqs = recurse_repeat_arr(seq, [])
        quot = floor(k / len(seq)) # integer part of quotient
        rmdr = k%len(seq)          # remainder
        out = []
        for s in seqs:
            outseq = s*quot + s[0:rmdr]
            out.append(outseq)
        if canonical:
            out = [canonical_kmer(o) for o in out]
        return(out)

# main ----------------------------------------------------------------------

# read dict containing k-mer counts
with open(args.counts, "r") as fh:
    cnts = json.load(fh)

# reverse-sort k-mers by counts
cnts_sort = sorted(cnts, key=cnts.get, reverse=True)
print(f'Total of {len(cnts_sort)} k-mers in input')

# write output files
outjson = args.output + ".kmers.json"
outtbl = args.output + ".report.tsv"

fh = open(outtbl, "w")
fh.write("\t".join(['repeat','total','cv','zeroes']) + "\n") # header

# process k-mers 
out = defaultdict(lambda: defaultdict (int))
tracker = dict(zip(cnts_sort,([0]*len(cnts_sort))))
counter=0
for kmer in cnts_sort:
    counter += 1
    if counter % 100000 == 0:
        print(f'Processed {counter} k-mers from input')
    if tracker[kmer] == 0: # avoid kmers already processed
        (replen, rep) = get_repeat_frame(kmer)
        if replen < 19: # ignore sequences that are evidently not repeats
            rep = get_canonical_repeat(rep)
            outrec = {}
            expected_kmers = get_expected_kmers_from_repeat(rep, 19)
            for kk in expected_kmers:
                try:
                    # out[rep][kk] = cnts[kk]
                    outrec[kk] = cnts[kk]
                    tracker[kk] += 1
                except KeyError: 
                    # out[rep][kk] = 0
                    outrec[kk] = 0
            if list(outrec.values()).count(0) < int(args.maxzeroes):
                out[rep] = outrec
                fh.write("\t".join([rep,
                                    str(sum(outrec.values())),
                                    str(pstdev(outrec.values()) / mean(outrec.values())),
                                    str(list(outrec.values()).count(0))
                                    ])
                          + "\n")               

fh.close()

# dump k-mers for each repeat seq above cutoff into a json file
print(f'Writing output to file {outjson}')
with open(outjson, "w") as fh:
    fh.write(json.dumps(out,indent=2))

"""
# report counts and coverage
print(f'Writing output to file {outtbl}')
with open(outtbl, "w") as fh:
    fh.write("\t".join(['repeat','total','pstdev','zeroes']) + "\n") # header
    for rep in out:
        if list(out[rep].values()).count(0) < int(args.maxzeroes):
            fh.write("\t".join([rep,
                                str(sum(out[rep].values())),
                                str(pstdev(out[rep].values())),
                                str(list(out[rep].values()).count(0))
                                ])
                     + "\n")
# dump everything into a json file
print(f'Writing output to file {outjson}')
with open(outjson, "w") as fh:
    fh.write(json.dumps(out,indent=2))

"""
