"""
Algoritmo Iterativo MUSCLE
"""

import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from core import kimura_distance
from core.kmer_distance import build_kmer_matrix
from core.upgma import UPGMANode, build_upgma_tree, get_leaf_order
from core.profile_alignment import progressive_alignment, align_profiles, print_alignment
from core.kimura_distance import build_kimura_matrix
from core.sp_score import sp_score

class MUSCLE:
    def __init__(self, match: int = 1, mismatch: int = -1, gap: int = -2, max_iters: int = 2):
        self.match = match
        self.mismatch = mismatch
        self.gap = gap
        self.max_iters = max_iters

    def align(self, sequences: list[str], ids: list[str] = None) -> tuple:
        n = len(sequences)
        if ids is None:
            ids = [f"seq{i+1}" for i in range(n)]

        print("STEP 1: Draft Progressive Alignment")
        kmer_matrix = build_kmer_matrix(sequences)
        tree1 = build_upgma_tree(kmer_matrix, ids)
        msa1, order1 = progressive_alignment(sequences, tree1, self.match, self.mismatch, self.gap)

        msa1_original_order = [""] * n
        for i, seq_idx in enumerate(order1):
            msa1_original_order[seq_idx] = msa1[i]

        print("STEP 2: Improved Progressive Alignment")
        kimura_matrix = build_kimura_matrix(msa1_original_order)
        tree2 = build_upgma_tree(kimura_matrix, ids)
        msa2, order2 = progressive_alignment(sequences, tree2, self.match, self.mismatch, self.gap)

        print("STEP 3: Refinement")
        current_msa = msa2
        current_order = order2
        current_score = sp_score(current_msa, self.match, self.mismatch, self.gap)

        internal_nodes = []
        def _collect_internals(node):
            if not node.is_leaf():
                internal_nodes.append(node)
                _collect_internals(node.left)
                _collect_internals(node.right)

        _collect_internals(tree2)

        for iteration in range(self.max_iters):
            improved = False
            print(f"Iteration {iteration+1}")

            for node in internal_nodes:
                if node == tree2:
                    continue

                profile1_seq_indices = set(get_leaf_order(node))
                all_indices = set(range(n))
                profile2_seq_indices = all_indices - profile1_seq_indices

                p1_raw = []
                p2_raw = []

                for i, seq_idx in enumerate(current_order):
                    if seq_idx in profile1_seq_indices:
                        p1_raw.append(current_msa[i])
                    else:
                        p2_raw.append(current_msa[i])

                def _remove_empty_cols(profile):
                    if not profile:
                        return []
                    keep_cols = []
                    length = len(profile[0])
                    
                    for col in range(length):
                        if any(seq[col] != '-' for seq in profile):
                            keep_cols.append(col)
                    return ["".join(seq[c] for c in keep_cols) for seq in profile]

                p1_clean = _remove_empty_cols(p1_raw)
                p2_clean = _remove_empty_cols(p2_raw)

                new_msa_raw = align_profiles(p1_clean, p2_clean, self.match, self.mismatch, self.gap)
                new_score = sp_score(new_msa_raw, self.match, self.mismatch, self.gap)

                if new_score > current_score:
                    print(f"- Better solution found by splitting node {node.label} - Score: {current_score:.1f} -> {new_score:.1f}")
                    current_score = new_score
                    current_msa = new_msa_raw
                    
                    new_order = []
                    for seq_idx in current_order:
                        if seq_idx in profile1_seq_indices:
                            new_order.append(seq_idx)
                    for seq_idx in current_order:
                        if seq_idx in profile2_seq_indices:
                            new_order.append(seq_idx)

                    current_order = new_order
                    improved = True

            if not improved:
                print("  Convergence reached (no improvements)")
                break

        aligned_ids = [ids[i] for i in current_order]
        return current_msa, aligned_ids, current_score

if __name__ == "__main__":
    seqs = [
        "MKTAYIAKQRQISFV",
        "MKTAYIAKQRQISFVKS",
        "MKTAIAKQRQISFV",
        "GGGGCCCCTTTTAAAA",
    ]
    ids = ["Seq1", "Seq2", "Seq3", "Seq4"]

    muscle = MUSCLE(max_iters=3)
    final_msa, aligned_ids, score = muscle.align(seqs, ids)

    print("\n=== Resultado Final de MUSCLE ===")
    print(f"SP-Score Final: {score:.1f}")
    print_alignment(final_msa, aligned_ids)