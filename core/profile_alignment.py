"""
Modificación de Needle-Wunsch para la alineación de Perfiles
en vez de secuencias individuales, y función de alineamiento
para recorrer el árbol guía para el MSA final.
"""

import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from core.upgma import UPGMANode

def sequence_to_profile(seq: str) -> list[str]:
    return [seq]

def _get_column(profile: list[str], index: int) -> list[str]:
    return [seq[index] for seq in profile]

def _column_score(col1: list[str], col2: list[str], match: int, mismatch: int, gap: int) -> float:
    total = 0
    count = 0
    for a in col1:
        for b in col2:
            if a == '-' and b == '-':
                score = 0
            elif a == '-' or b == '-':
                score = gap
            elif a == b:
                score = match
            else:
                score = mismatch

            total += score
            count += 1

    return total / count if count > 0 else 0.0

def align_profiles(profile1: list[str], profile2: list[str], match: int = 1, mismatch: int = -1, gap: int = -2) -> list[str]:
    n = len(profile1[0])
    m = len(profile2[0])
    
    score_matrix = [[0.0] * (m + 1) for _ in range(n + 1)]
    for i in range(n + 1):
        score_matrix[i][0] = i * gap
    for j in range(m + 1):
        score_matrix[0][j] = j * gap

    cols1 = [_get_column(profile1, i) for i in range(n)]
    cols2 = [_get_column(profile2, j) for j in range(m)]

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            col_score = _column_score(cols1[i - 1], cols2[j - 1], match, mismatch, gap)
            score_diagonal = score_matrix[i - 1][j - 1] + col_score
            score_up = score_matrix[i - 1][j] + gap
            score_left = score_matrix[i][j - 1] + gap
            score_matrix[i][j] = max(score_diagonal, score_up, score_left)

    aligned_cols_1 = []
    aligned_cols_2 = []
    i, j = n, m

    gap_col_for_profile1 = ['-'] * len(profile1)
    gap_col_for_profile2 = ['-'] * len(profile2)

    while i > 0 or j > 0:
        if i > 0 and j > 0:
            col_score = _column_score(cols1[i - 1], cols2[j - 1], match, mismatch, gap)
            score_diagonal = score_matrix[i - 1][j - 1] + col_score
            score_up = score_matrix[i - 1][j] + gap
            score_left = score_matrix[i][j - 1] + gap

            if score_matrix[i][j] == score_diagonal:
                aligned_cols_1.append(cols1[i - 1])
                aligned_cols_2.append(cols2[j - 1])
                i -= 1
                j -= 1
            elif score_matrix[i][j] == score_up:
                aligned_cols_1.append(cols1[i - 1])
                aligned_cols_2.append(gap_col_for_profile2)
                i -= 1
            else:
                aligned_cols_1.append(gap_col_for_profile1)
                aligned_cols_2.append(cols2[j - 1])
                j -= 1
        elif i > 0:
            aligned_cols_1.append(cols1[i - 1])
            aligned_cols_2.append(gap_col_for_profile2)
            i -= 1
        else:
            aligned_cols_1.append(gap_col_for_profile1)
            aligned_cols_2.append(cols2[j - 1])
            j -= 1

    aligned_cols_1.reverse()
    aligned_cols_2.reverse()

    n_seqs_1 = len(profile1)
    n_seqs_2 = len(profile2)
    new_profile = []

    for seq_idx in range(n_seqs_1):
        new_seq = "".join(col[seq_idx] for col in aligned_cols_1)
        new_profile.append(new_seq)

    for seq_idx in range(n_seqs_2):
        new_seq = "".join(col[seq_idx] for col in aligned_cols_2)
        new_profile.append(new_seq)

    return new_profile

def progressive_alignment(sequences: list[str], tree: UPGMANode, match: int = 1, mismatch: int = -1, gap: int = -2) -> tuple:
    def _align_node(node: UPGMANode) -> tuple:
        if node.is_leaf():
            return sequence_to_profile(sequences[node.seq_index]), [node.seq_index]

        left_profile, left_order = _align_node(node.left)
        right_profile, right_order = _align_node(node.right)
        merged_profile = align_profiles(left_profile, right_profile, match, mismatch, gap)
        merged_order = left_order + right_order
        return merged_profile, merged_order

    return _align_node(tree)

def print_alignment(profile: list[str], labels: list[str] = None) -> None:
    if labels is None:
        labels = [f"seq{i+1}" for i in range(len(profile))]
    width = max(len(lbl) for lbl in labels) + 2
    for lbl, seq in zip(labels, profile):
        print(f"{lbl:<{width}}{seq}")

if __name__ == "__main__":
    from core.distance_matrix import build_distance_matrix, print_distance_matrix
    from core.upgma import build_upgma_tree, print_tree

    seqs = [
        "MKTAYIAKQRQISFV",
        "MKTAYIAKQRQISFVKS",
        "MKTAIAKQRQISFV",
        "GGGGCCCCTTTTAAAA",
    ]
    ids = ["seq1", "seq2", "seq3", "seq4"]

    dist_matrix = build_distance_matrix(seqs)
    print("Distance Matrix:")
    print_distance_matrix(dist_matrix, ids)

    tree = build_upgma_tree(dist_matrix, ids)
    print("\nGuide Tree:")
    print_tree(tree)

    msa, seq_order = progressive_alignment(seqs, tree)
    labels_in_order = [ids[i] for i in seq_order]
    print("\nMSA final:")
    print_alignment(msa, labels_in_order)
