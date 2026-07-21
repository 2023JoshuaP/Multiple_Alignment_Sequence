"""
Matriz de distancias entre todas las secuencias
a partir del alineamiento Needle-Wunsch
"""

import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from core.pairwise_nw import needle_wunsch, AlignmentResult

def percent_identify(aligned_seq1: str, aligned_seq2: str) -> float:
    if len(aligned_seq1) != len(aligned_seq2):
        raise ValueError("The alignment sequences must have the same length.")

    length = len(aligned_seq1)
    if length == 0:
        return 0.0

    matches = sum(1 for a, b in zip(aligned_seq1, aligned_seq2) if a == b and a != '-')

    return matches / length

def pairwise_distance(seq1: str, seq2: str, match: int = 1, mismatch: int = -1, gap: int = -2) -> float:
    result: AlignmentResult = needle_wunsch(seq1, seq2, match, mismatch, gap)
    identity = percent_identify(result.aligned_seq1, result.aligned_seq2)
    
    return 1.0 - identity

def _pairwise_distance_wrapper(args):
    i, j, seq1, seq2, match, mismatch, gap = args
    return i, j, pairwise_distance(seq1, seq2, match, mismatch, gap)

def build_distance_matrix(sequences: list[str], match: int = 1, mismatch: int = -1, gap: int = -2) -> list[list[float]]:
    n = len(sequences)
    matrix = [[0.0] * n for _ in range(n)]
    
    tasks = []
    for i in range(n):
        for j in range(i + 1, n):
            tasks.append((i, j, sequences[i], sequences[j], match, mismatch, gap))
            
    # Para pocos cruces, no vale la pena el overhead de multiprocesamiento
    if len(tasks) < 5:
        for args in tasks:
            i, j, d = _pairwise_distance_wrapper(args)
            matrix[i][j] = d
            matrix[j][i] = d
    else:
        import concurrent.futures
        import multiprocessing
        cores = max(1, multiprocessing.cpu_count() - 1)
        print(f"      -> Computando matriz de distancias usando {cores} núcleos en paralelo...")
        with concurrent.futures.ProcessPoolExecutor(max_workers=cores) as executor:
            for i, j, d in executor.map(_pairwise_distance_wrapper, tasks):
                matrix[i][j] = d
                matrix[j][i] = d

    return matrix

def print_distance_matrix(matrix: list[list[float]], labels: list[str] = None) -> None:
    n = len(matrix)
    if labels is None:
        labels = [str(i) for i in range(n)]
    
    width = max(len(lbl) for lbl in labels) + 2

    header = " " * width + "".join(f"{lbl:>8}" for lbl in labels)
    print(header)
    for i in range(n):
        row = f"{labels[i]:<{width}}" + "".join(f"{matrix[i][j]:8.3f}" for j in range(n))
        print(row)

if __name__ == "__main__":
    seqs = [
        "MKTAYIAKQRQISFV",
        "MKTAYIAKQRQISFVKS",
        "MKTAIAKQRQISFV",
        "GGGGCCCCTTTTAAAA",
    ]
    ids = ["seq1", "seq2", "seq3", "seq4"]

    dist_matrix = build_distance_matrix(seqs)
    print_distance_matrix(dist_matrix, ids)
