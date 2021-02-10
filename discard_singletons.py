#!/usr/bin/env python3

"""Parse jellyfish dumps file and discard singletons/low coverage"""

import argparse
import json
import re
from collections import defaultdict

parser = argparse.ArgumentParser(
    description="Parse jellyfish dumps file and discard singletons/low coverage")
parser.add_argument("--input", "-i",
    help="Dumps file", action="store")
parser.add_argument("--histo",
    help="File to store histogram of k-mer counts", default="counts_histo.json")
parser.add_argument("--cutoff", type=int, default=2,
    help="Cutoff (inclusive) for k-mer coverage to include")
parser.add_argument("--output", "-o", default="counts_kmers.json",
    help="Output file containing kmers and counts above cutoff")

args = parser.parse_args()

histo = defaultdict(int) # Histogram of frequencies of k-mer counts
keeper = defaultdict(int) # dict of kmers above cutof

cur_count = 0

with open(args.input, "r") as fh:
    for line in fh:
        line = line.rstrip()
        matcher = re.match(r"^>(\d+)$", line)
        if matcher:
            cur_count = int(matcher.group(1))
            histo[cur_count] += 1
        else:
            if cur_count > args.cutoff:
                keeper[line] = cur_count

with open(args.histo, "w") as fh:
    fh.write(json.dumps(histo, indent=4))

with open(args.output, "w") as fh:
    fh.write(json.dumps(keeper, indent=4))
