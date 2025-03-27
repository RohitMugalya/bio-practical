from collections import Counter

def find_frequent_kmer(dna_string, k):

    k_mers = [dna_string[i:i+k] for i in range(len(dna_string) - k + 1)]

    k_mer_count = Counter(k_mers)

    most_frequent = k_mer_count.most_common(1)

    return most_frequent

dna = "AGCTAGCTAGCTAGGCTA"
k = 3
most_frequent_kmer = find_frequent_kmer(dna, k)
print("Most frequent k-mer:", most_frequent_kmer)
