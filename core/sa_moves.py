"""
Operadores de mutación para Simulated Annealing,
son funciones para alterar estocásticamente un
alineamiento múltiple
"""

import random
import copy

def get_random_neighbor(alignment: list[str]) -> list[str]:
    n_sequences = len(alignment)
    if n_sequences == 0:
        return alignment

    length = len(alignment[0])
    r = random.random()
    new_alignment = list(alignment)

    """
    Se decide que mutacion va a utilizar
    - 80% probabilidad: Intercambia gap por caracter (shift local)
    - 10% probabilidad: Insertar columna de gaps globalmente
    - 10% probabilidad: Eliminar columna de gaps vacia (solo si existe)
    """
    if r < 0.8 and length > 1:
        print("Mutation 1: GAP SWAP (shift local)")
        seq_idx = random.randint(0, n_sequences - 1)
        seq = list(new_alignment[seq_idx])
        gap_indices = [i for i, char in enumerate(seq) if char == '-']

        if gap_indices:
            g_idx = random.choice(gap_indices)
            directions = []
            if g_idx > 0 and seq[g_idx - 1] != '-':
                directions.append(-1)
            if g_idx < length - 1 and seq[g_idx + 1] != '-':
                directions.append(1)
            if directions:
                d = random.choice(directions)
                seq[g_idx], seq[g_idx + d] = seq[g_idx + d], seq[g_idx]
                new_alignment[seq_idx] = "".join(seq)

    elif r < 0.9:
        print("\nMutation 2: INSERT COLUMN GAPS")
        col_idx = random.randint(0, length)
        for i in range(n_sequences):
            seq = new_alignment[i]
            new_alignment[i] = seq[:col_idx] + "-" + seq[col_idx:]

    else:
        print("\nMutation 3: DELETE COLUMN GAPS")
        empty_cols = []
        for col in range(length):
            if all(seq[col] == '-' for seq in new_alignment):
                empty_cols.append(col)

        if empty_cols:
            col_remove = random.choice(empty_cols)
            for i in range(n_sequences):
                seq = new_alignment[i]
                new_alignment[i] = seq[:col_remove] + seq[col_remove + 1:]

    return new_alignment
