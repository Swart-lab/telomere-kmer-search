Telomere repeat prediction from k-mer counts
============================================

Predict telomere repeats using k-mer counts from unassembled read libraries.

This assumes that telomeres are the most abundant direct tandem repeat
sequences in a given genome library. The output is a list of direct tandem
repeat unit sequences, sorted by frequency.

k-mers are first counted with [Jellyfish](https://github.com/gmarcais/Jellyfish)
before being k-mer sequences and counts are dumped to Fasta format. Those above
the minimum coverage are then retained, and used to search for the most
abundant tandem direct repeats. An example bash script for this task is given in
`do_jellyfish.sh`.

k-mer length used for k-mer counting should be prime and longer than the
maximum expected telomere repeat unit length.
