"""
Parser de archivos FASTA para la lectura
de las secuencias correspondientes
"""

import os

class FastaRecord:
    __slots__ = ("id", "description", "sequence")

    def __init__(self, header: str, sequence: str):
        header = header.strip()
        self.description = header
        self.id = header.split()[0] if header else ""
        self.sequence = sequence

    def __len__(self):
        return len(self.sequence)

    def __repr__(self):
        preview = self.sequence[:30] + ("..." if len(self.sequence) > 30 else "")
        return f"FastaRecord(id={self.id!r}, len={len(self.sequence)}, seq={preview!r})"

def parse_fasta(filepath: str) -> list[FastaRecord]:
    if not os.path.isfile(filepath):
        raise FileNotFoundError(f"No found the FASTA file: {filepath}")

    records: list[FastaRecord] = []
    current_header = None
    current_seq_chunks: list[str] = []

    with open(filepath, "r", encoding="utf-8") as f:
        for raw_line in f:
            line = raw_line.strip()
            
            if not line:
                continue
            if line.startswith(";"):
                continue

            if line.startswith(">"):
                if current_header is not None:
                    records.append(FastaRecord(current_header, "".join(current_seq_chunks)))
                current_header = line[1:]
                current_seq_chunks = []
            else:
                if current_header is None:
                    raise ValueError("This file don't FASTA.")

                cleaned = "".join(line.split()).upper()
                current_seq_chunks.append(cleaned)
        
        if current_header is not None:
            records.append(FastaRecord(current_header, "".join(current_seq_chunks)))

    if not records:
        raise ValueError("No sequences in FASTA")

    return records

def write_fasta(records: list[FastaRecord], filepath: str, line_width: int = 60) -> None:
    with open(filepath, "w", encoding="utf-8") as f:
        for record in records:
            f.write(f">{record.description}\n")
            seq = record.sequence
            for i in range(0, len(seq), line_width):
                f.write(seq[i:i + line_width] + "\n")

def get_sequence(records: list[FastaRecord]) -> list[str]:
    return [r.sequence for r in records]

def get_ids(records: list[FastaRecord]) -> list[str]:
    return [r.id for r in records]

if __name__ == "__main__":
    import sys

    if len(sys.argv) != 2:
        print("Use: python <file> <fasta file>")
        sys.exit(1)

    recs = parse_fasta(sys.argv[1])
    print(f"Read {len(recs)} sequences:\n")
    for r in recs:
        print(r)
