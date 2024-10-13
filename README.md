# DNA Hash

A Python library for counting short DNA sequences for use in Bioinformatics. DNA Hash stores k-mer sequence counts by their up2bit encoding - a two-way hash that works with variable-length sequences. DNA Hash uses considerably less memory than a lookup table that stores sequences in plaintext. In addition, DNA Hash's novel autoscaling Bloom filter eliminates the need to explicitly store counts for sequences that have only been seen once.

- **Ultra-low** memory footprint
- **Embarrassingly** parallelizable
- **Open-source** and free to use commercially

> **Note:** The maximum sequence length is platform dependent. On a 64-bit machine, the max length is 31. On a 32-bit machine, the max length is 15.

> **Note:** Due to the probabilistic nature of the Bloom filter, DNA Hash may over count sequences at a bounded rate.

## References
- [1] https://github.com/JohnLonginotto/ACGTrie/blob/master/docs/UP2BIT.md.
- [2] P. Melsted et al. (2011). Efficient counting of k-mers in DNA sequences using a bloom filter.
- [3] S. Deorowicz et al. (2015). KMC 2: fast and resource-frugal k-mer counting.
- [4] A. DalPino. (2021). OkBloomer, a novel autoscaling Bloom Filter [[link](https://github.com/andrewdalpino/PyBloomer)].