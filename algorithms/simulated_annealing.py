"""
Algoritmo Simulated Annealing
"""

import sys
import os
import math
import random
import time

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from core.sp_score import sp_score
from core.sa_moves import get_random_neighbor
from core.profile_alignment import print_alignment
from algorithms.clustalw import ClustalW

class SimulatedAnnealing:
    def __init__(self, initial_temp: float = 100.0, cooling_rate: float = 0.95, max_iters: int = 1000, match: int = 1, mismatch: int = -1, gap: int = -2):
        self.initial_temp = initial_temp
        self.cooling_rate = cooling_rate
        self.max_iters = max_iters
        self.match = match
        self.mismatch = mismatch
        self.gap = gap

    def align(self, sequences: list[str], ids: list[str] = None) -> tuple:
        n = len(sequences)
        if ids is None:
            ids = [f"seq{i+1}" for i in range(n)]

        print("STEP 1: Generating initial state (Draft ClustalW)")
        draf_algo = ClustalW(self.match, self.mismatch, self.gap)
        current_msa, aligned_ids, tree, _ = draf_algo.align(sequences, ids)

        current_score = sp_score(current_msa, self.match, self.mismatch, self.gap)
        print(f"\nInitial Score: {current_score:.1f}")

        T = self.initial_temp
        best_msa = current_msa
        best_score = current_score

        print(f"STEP 2: Starting Simulated Annealing (T_0 = {T}, alpha = {self.cooling_rate})")
        for iteration in range(self.max_iters):
            new_msa = get_random_neighbor(current_msa)
            new_score = sp_score(new_msa, self.match, self.mismatch, self.gap)
            delta_score = new_score - current_score

            accept = False
            if delta_score >= 0:
                accept = True
            else:
                probability = math.exp(delta_score / T) if T > 0 else 0
                if random.random() < probability:
                    accept = True

            if accept:
                current_msa = new_msa
                current_score = new_score

                if current_score > best_score:
                    best_score = current_score
                    best_msa = current_msa

            T = T * self.cooling_rate

            if (iteration + 1) % 100 == 0:
                print(f"Iteration {iteration+1:4d}/{self.max_iters} | Temp: {T:.4f} | Current Score: {current_score:.1f} | Best Score: {best_score:.1f}")
            if T < 1e-8:
                print(f"Frozen temperature in iteration {iteration+1}. End of Annealing")
                break

        print("STEP 3: Terminating Algorithm")
        return best_msa, aligned_ids, best_score, tree

if __name__ == "__main__":
    seqs = [
        "MKTAYIAKQRQISFV",
        "MKTAYIAKQRQISFVKS",
        "MKTAIAKQRQISFV",
        "GGGGCCCCTTTTAAAA",
    ]
    ids = ["Seq1", "Seq2", "Seq3", "Seq4"]

    sa = SimulatedAnnealing(initial_temp=50.0, cooling_rate=0.9, max_iters=500)
    final_msa, aligned_ids, score, _ = sa.align(seqs, ids)
    
    print("\n=== Resultado Final de SA ===")
    print(f"SP-Score Final: {score:.1f}")
    print_alignment(final_msa, aligned_ids)
