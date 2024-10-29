from dna_hash import DNAHash
from dna_hash.tokenizers import Kmer, Canonical

from Bio import SeqIO
from matplotlib import pyplot as plt

hash_table = DNAHash(max_false_positive_rate=0.001)

tokenizer = Canonical(Kmer(6))

with open("covid-19-virus.fasta", "r") as file:
    for record in SeqIO.parse(file, "fasta"):
        for token in tokenizer.tokenize(str(record.seq)):
            hash_table.increment(token)

for sequence, count in hash_table.top(25):
    print(f"{sequence}: {count}")

print(f"Total sequences: {hash_table.num_sequences}")
print(f"# of unique sequences: {hash_table.num_unique_sequences}")
print(f"# of singletons: {hash_table.num_singletons}")

plt.hist(list(hash_table.counts.values()), bins=20)
plt.title("SARS-CoV-2 Genome")
plt.xlabel("Counts")
plt.ylabel("Frequency")
plt.show()
