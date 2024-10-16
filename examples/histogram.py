from dna_hash import DNAHash
from dna_hash.tokenizers import Kmer, Canonical

from Bio import SeqIO
from matplotlib import pyplot as plt

hash_table = DNAHash(max_false_positive_rate=0.001)

tokenizer = Canonical(Kmer(6))

with open('covid-19-virus.fasta', 'r') as file:
    for record in SeqIO.parse(file, 'fasta'):
        for token in tokenizer.tokenize(str(record.seq)):
            hash_table.increment(token)

counts, bins = hash_table.histogram(20)

plt.stairs(counts, bins)
plt.title('Histogram of SARS-CoV-2 Genome')
plt.xlabel('Counts')
plt.ylabel('Frequency')
plt.show()
