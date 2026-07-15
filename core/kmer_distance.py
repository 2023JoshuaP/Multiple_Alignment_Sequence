"""
Cálculo de las distancias basado en k-mers.
Corresponde a la fase inicial de MUSCLE.
"""

import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from collections import Counter

def _get_kmers(seq: str, k: int = 3) -> Counter:
    kmers = Counter()
    for i in range(len(seq) - k + 1):
        kmer = seq[i:i+k]
        kmers[kmer] += 1
    return kmers

def kmer_distance(seq1: str, seq2: str, k: int = 3) -> float:
    if len(seq1) < k and len(seq2) < k:
        return 1.0

    kmers1 = _get_kmers(seq1, k)
    kmers2 = _get_kmers(seq2, k)

    shared_kmers = 0
    for kmer, count_1 in kmers1.items():
        if kmer in kmers2:
            shared_kmers += min(count_1, kmers2[kmer])

    total_1 = sum(kmers1.values())
    total_2 = sum(kmers2.values())

    max_possible_shared = min(total_1, total_2)

    if max_possible_shared == 0:
        return 1.0

    fraction_shared = shared_kmers / max_possible_shared

    return 1.0 - fraction_shared

def build_kmer_matrix(sequences: list[str], k: int = 3) -> list[list[float]]:
    n = len(sequences)
    matrix = [[0.0] * n for _ in range(n)]

    for i in range(n):
        for j in range(i + 1, n):
            d = kmer_distance(sequences[i], sequences[j], k)
            matrix[i][j] = d
            matrix[j][i] = d

    return matrix

if __name__ == "__main__":
    from core.distance_matrix import print_distance_matrix

    seqs = [
        "MKTAYIAKQRQISFV",
        "MKTAYIAKQRQISFVKS",
        "MKTAIAKQRQISFV",
        "GGGGCCCCTTTTAAAA",
    ]
    ids = ["Seq1", "Seq2", "Seq3", "Seq4"]

    print("=== Distancias K-mer (k=3) ===")
    matrix = build_kmer_matrix(seqs)
    print_distance_matrix(matrix, ids)
