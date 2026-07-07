"""
Algoritmo Progresivo ClustalW
"""

import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from core.distance_matrix import build_distance_matrix, print_distance_matrix
from core.neighbor_joining import build_nj_tree, print_nj_tree
from core.profile_alignment import progressive_alignment, print_alignment

class ClustalW:
    def __init__(self, match: int = 1, mismatch: int = -1, gap: int = -2):
        self.match = match
        self.mismatch = mismatch
        self.gap = gap

    def align(self, sequences: list[str], ids: list[str] = None) -> tuple:
        n = len(sequences)
        if ids is None:
            ids = [f"seq{i+1}" for i in range(n)]

        print("Step 1: Calculate distance matrix")
        dist_matrix = build_distance_matrix(sequences, self.match, self.mismatch, self.gap)

        print("Step 2: Building Guide Tree with Neighbor Joining")
        guide_tree = build_nj_tree(dist_matrix, ids)

        print("Step 3: Progressive Alignment with ClustalW")
        msa_profile, seq_order = progressive_alignment(sequences, guide_tree, self.match, self.mismatch, self.gap)

        aligned_ids = [ids[i] for i in seq_order]
        return msa_profile, aligned_ids, guide_tree, dist_matrix

if __name__ == "__main__":
    seqs = [
        "MKTAYIAKQRQISFV",
        "MKTAYIAKQRQISFVKS",
        "MKTAIAKQRQISFV",
        "GGGGCCCCTTTTAAAA",
    ]
    ids = ["Seq1", "Seq2", "Seq3", "Seq4"]

    print("=== Ejecución de ClustalW ===\n")
    clustal = ClustalW()
    msa, aligned_ids, tree, matrix = clustal.align(seqs, ids)

    print("\n--- Resultados ---")
    print("\n1. Matriz de Distancias:")
    print_distance_matrix(matrix, ids)

    print("\n2. Árbol Guía (Neighbor-Joining):")
    print_nj_tree(tree)

    print("\n3. Alineamiento Múltiple Final:")
    print_alignment(msa, aligned_ids)
