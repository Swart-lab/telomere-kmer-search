Telomere repeat prediction from k-mer counts
============================================

Predict telomere repeats using k-mer counts from unassembled read libraries.

This assumes that telomeres are the most abundant direct tandem repeat
sequences in a given genome library. The output is a list of direct tandem
repeat unit sequences, sorted by frequency.

k-mers are first counted with [jellyfish
count](https://github.com/gmarcais/Jellyfish) before being k-mer sequences and
counts are dumped to Fasta format. Those above the minimum coverage are then
retained, and used to search for the most abundant tandem direct repeats. An
example bash script for this task is given in `do_jellyfish.sh`. These scripts
have been tested with jellyfish v2.2.10.

NB:

 * k-mer length used for k-mer counting should be prime and longer than the
   maximum expected telomere repeat unit length.
 * Minimum coverage cutoff supplied to `discard_singletons.py` should be
   slightly above the expected average coverage of the genome, to retain only
   repeat elements.

The JSON file containing counts of high-frequency k-mers is then supplied to
the script `find_repeats_from_kmers.py` along with the k-mer length used for
`jellyfish count`. Example:

```bash
./find_repeats_from_kmers.py \
  --maxzeroes 1 \
  -k 19 \
  --counts library.k19.cov100.json \
  --output library.k19.repeats
```
