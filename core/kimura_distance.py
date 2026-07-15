"""
Distancia de Kimura para las secuencias alineadas.
Utilizado en MUSCLE para refinar el Árbol con distancias
evolutivas más precisas.
"""

import sys
import os
import math
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

def pairwise_identity(aligned_seq1: str, aligned_seq2: str) -> float:
    if len(aligned_seq1) != len(aligned_seq2):
        raise ValueError("The sequences must be aligned")

    matches = 0
    valid_length = 0

    for a, b in zip(aligned_seq1, aligned_seq2):
        if a == '-' and b == '-':
            continue
        valid_length += 1
        if a == b:
            matches += 1

    if valid_length == 0:
        return 0.0

    return matches / valid_length

def kimura_distance(aligned_seq1: str, aligned_seq2: str) -> float :
    identity = pairwise_identity(aligned_seq1, aligned_seq2)
    p = 1.0 - identity

    arg = 1.0 - p - 0.2 * (p ** 2)
    if arg <= 0.0:
        return 5.0

    return -math.log(arg)

def build_kimura_matrix(aligned_sequences: list[str]) -> list[list[float]]:
    n = len(aligned_sequences)
    matrix = [[0.0] * n for _ in range(n)]

    for i in range(n):
        for j in range(i + 1, n):
            d = kimura_distance(aligned_sequences[i], aligned_sequences[j])
            matrix[i][j] = d
            matrix[j][i] = d

    return matrix

if __name__ == "__main__":
    from core.distance_matrix import print_distance_matrix

    seqs = [
        "MKTAYIAKQRQISFV--",
        "MKTAYIAKQRQISFVKS",
        "MKTA-IAKQRQISFV--",
        "-GGGGCCCCTTTTAAAA",
    ]
    ids = ["Seq1", "Seq2", "Seq3", "Seq4"]

    print("=== Distancias de Kimura (Stage 2) ===")
    matrix = build_kimura_matrix(seqs)
    print_distance_matrix(matrix, ids)
