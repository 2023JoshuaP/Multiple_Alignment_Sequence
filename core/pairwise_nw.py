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

def needle_wunsch(seq1: str, seq2: str, match: int = 1, mismatch: int = -1, gap: int = -2, band_width: int = 400) -> AlignmentResult:
    n = len(seq1)
    m = len(seq2)
    
    if n < 500 and m < 500:
        band_width = max(n, m) + 1

    prev_score = [float('-inf')] * (m + 1)
    curr_score = [float('-inf')] * (m + 1)
    
    # Backtrace de n x m usando bytearray (1 byte por celda, memory-safe)
    # 0: None, 1: Diagonal, 2: Up, 3: Left
    backtrace = [bytearray(m + 1) for _ in range(n + 1)]

    for j in range(m + 1):
        if j > band_width: break
        prev_score[j] = j * gap
        backtrace[0][j] = 3 # Left

    for i in range(1, n + 1):
        for j in range(m + 1): curr_score[j] = float('-inf')
        
        if i <= band_width:
            curr_score[0] = i * gap
            backtrace[i][0] = 2 # Up
            
        center_j = int(i * (m / n)) if n > 0 else i
        start_j = max(1, center_j - band_width)
        end_j = min(m, center_j + band_width)

        for j in range(start_j, end_j + 1):
            s_diag = prev_score[j - 1] + (match if seq1[i - 1] == seq2[j - 1] else mismatch)
            s_up = prev_score[j] + gap
            s_left = curr_score[j - 1] + gap
            
            best_s = max(s_diag, s_up, s_left)
            curr_score[j] = best_s
            
            if best_s == s_diag:
                backtrace[i][j] = 1
            elif best_s == s_up:
                backtrace[i][j] = 2
            else:
                backtrace[i][j] = 3
                
        prev_score, curr_score = curr_score, prev_score

    final_score = prev_score[m]
    
    aligned_seq1_chars = []
    aligned_seq2_chars = []
    i, j = n, m

    while i > 0 or j > 0:
        direction = backtrace[i][j]
        if direction == 1 and i > 0 and j > 0:
            aligned_seq1_chars.append(seq1[i - 1])
            aligned_seq2_chars.append(seq2[j - 1])
            i -= 1
            j -= 1
        elif direction == 2 and i > 0:
            aligned_seq1_chars.append(seq1[i - 1])
            aligned_seq2_chars.append('-')
            i -= 1
        elif direction == 3 and j > 0:
            aligned_seq1_chars.append('-')
            aligned_seq2_chars.append(seq2[j - 1])
            j -= 1
        else:
            if i > 0 and j > 0:
                aligned_seq1_chars.append(seq1[i - 1])
                aligned_seq2_chars.append(seq2[j - 1])
                i -= 1
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

    return AlignmentResult(aligned_seq1, aligned_seq2, int(final_score))
