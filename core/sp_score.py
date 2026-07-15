"""
Calcula el Sum-of-Pairs score. Funcion objetivo para
decidir si un refinamiento mejoró el alineamiento o no
"""

def sp_score(alignment: list[str], match: int = 1, mismatch: int = -1, gap: int = -2) -> float:
    if not alignment:
        return 0.0

    n_seqs = len(alignment)
    if n_seqs < 2:
        return 0.0

    length = len(alignment[0])
    total_score = 0.0

    for col in range(length):
        column_chars = [seq[col] for seq in alignment]

        for i in range(n_seqs):
            for j in range(i + 1, n_seqs):
                a = column_chars[i]
                b = column_chars[j]

                if a == '-' and b == '-':
                    continue
                elif a == '-' or b == '-':
                    total_score += gap
                elif a == b:
                    total_score += match
                else:
                    total_score += mismatch

    return total_score

if __name__ == "__main__":
    # Prueba del SP-Score
    msa_good = [
        "ACTG",
        "ACTG",
        "ACTG"
    ]
    # Tres secuencias idénticas de 4 letras.
    # Pares por columna: 3 (1-2, 1-3, 2-3). Total: 3 * 4 * match(1) = 12

    msa_bad = [
        "ACTG",
        "AC-G",
        "A-TG"
    ]

    print("=== Sum-of-Pairs (SP) Score ===")
    print(f"Puntaje MSA Bueno: {sp_score(msa_good)}")
    print(f"Puntaje MSA Malo:  {sp_score(msa_bad)}")
