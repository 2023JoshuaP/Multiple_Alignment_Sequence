"""
Algoritmo Needle-Wunsch para el alineamiento global,
sera pieza base para los algoritmos implementados
"""

from dataclasses import dataclass

@dataclass
class AlignmentResult:
    aligned_seq1: str
    aligned_seq2: str
    score: int

def needle_wunsch(seq1: str, seq2: str, match: int, mismatch: int, gap:int) -> AlignmentResult:
    n = len(seq1)
    m = len(seq2)
    score_matrix = [[0] * (m + 1) for _ in range(n + 1)]

    for i in range(n + 1):
        score_matrix[i][0] = i * gap
    for j in range(m + 1):
        score_matrix[0][j] = j * gap

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            score_diagonal = score_matrix[i - 1][j - 1] + (match if seq1[i - 1] == seq2[j - 1] else mismatch)
            score_up = score_matrix[i - 1][j] + gap
            score_left = score_matrix[i][j - 1] + gap
            score_matrix[i][j] = max(score_diagonal, score_up, score_left)

    # Backtracking process
    aligned_seq1_chars = []
    aligned_seq2_chars = []
    i, j = n, m

    while i > 0 or j > 0:
        if i > 0 and j > 0:
            score_diagonal = score_matrix[i - 1][j - 1] + (match if seq1[i - 1] == seq2[j - 1] else mismatch)
            score_up = score_matrix[i - 1][j] + gap
            score_left = score_matrix[i][j - 1] + gap

            if score_matrix[i][j] == score_diagonal:
                aligned_seq1_chars.append(seq1[i - 1])
                aligned_seq2_chars.append(seq2[j - 1])
                i -= 1
                j -= 1
            elif score_matrix[i][j] == score_up:
                aligned_seq1_chars.append(seq1[i - 1])
                aligned_seq2_chars.append('-')
                i -= 1
            else:
                aligned_seq1_chars.append('-')
                aligned_seq2_chars.append(seq2[j - 1])
                j -= 1
        elif i > 0:
            aligned_seq1_chars.append(seq1[i - 1])
            aligned_seq2_chars.append('-')
            i -= 1
        else:
            aligned_seq1_chars.append('-')
            aligned_seq2_chars.append(seq2[j - 1])
            j -= 1

    aligned_seq1 = "".join(reversed(aligned_seq1_chars))
    aligned_seq2 = "".join(reversed(aligned_seq2_chars))

    return AlignmentResult(aligned_seq1, aligned_seq2, score_matrix[n][m])
