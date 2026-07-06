"""
Algoritmo UPGMA para la construcción del Árbol Guía
usado en los algoritmos progresivo de MSA
"""

import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from dataclasses import dataclass, field

@dataclass
class UPGMANode:
    left: "UPGMANode" = None
    right: "UPGMANode" = None
    height: float = 0.0
    seq_index: int = None
    size: int = 1
    label: str = ""

    def is_leaf(self) -> bool:
        return self.left is None and self.right is None

    def __repr__(self):
        if self.is_leaf():
            return f"Leaf({self.label})"
        return f"Node({self.label}, h={self.height:.3f})"

def _find_closest_pair(active_distances: dict) -> tuple:
    best_pair = None
    best_dist = float("inf")
    for (i, j), d in active_distances.items():
        if d < best_dist:
            best_dist = d
            best_pair = (i, j)
    return best_pair[0], best_pair[1], best_dist

def build_upgma_tree(distance_matrix: list[list[float]], labels: list[str] = None) -> UPGMANode:
    n = len(distance_matrix)
    if n < 2:
        raise ValueError("At least two sequences are required to construct a UPGMA tree.")
    for row in distance_matrix:
        if len(row) != n:
            raise ValueError("The distance_matrix must have cuadratic n x n")

    if labels is None:
        labels = [str(i) for i in range(n)]

    clusters: dict = {
        i: UPGMANode(seq_index=i, height=0.0, size=1, label=labels[i])
        for i in range(n)
    }

    active_distances: dict = {}
    for i in range(n):
        for j in range(i + 1, n):
            active_distances[(i, j)] = distance_matrix[i][j]
    
    next_id = n

    while len(clusters) > 1:
        i, j, d_ij = _find_closest_pair(active_distances)
        node_i = clusters[i]
        node_j = clusters[j]

        new_height = d_ij / 2.0
        new_size = node_i.size + node_j.size
        new_node = UPGMANode(left=node_i, right=node_j, height=new_height, size=new_size, label=f"({node_i.label}, {node_j.label})",)

        remaining_ids = [cid for cid in clusters.keys() if cid not in (i, j)]
        for k in remaining_ids:
            d_ik = _get_distance(active_distances, i, k)
            d_jk = _get_distance(active_distances, j, k)
            new_d = (node_i.size * d_ik + node_j.size * d_jk) / new_size
            _set_distance(active_distances, next_id, k, new_d)

        active_distances = {pair: dist for pair, dist in active_distances.items() if i not in pair and j not in pair}

        del clusters[i]
        del clusters[j]
        clusters[next_id] = new_node

        next_id += 1

    root_id = next(iter(clusters))
    return clusters[root_id]

def _get_distance(active_distances: dict, a: int, b: int) -> float:
    key = (a, b) if a < b else (b, a)
    return active_distances[key]

def _set_distance(active_distances: dict, a: int, b: int, value: float) -> None:
    key = (a, b) if a < b else (b, a)
    active_distances[key] = value

def print_tree(node: UPGMANode, indent: str = "", is_last: bool = True) -> None:
    connector = "└── " if is_last else "├── "
    if node.is_leaf():
        print(f"{indent}{connector}{node.label} (leaf)")
    else:
        print(f"{indent}{connector}[height={node.height:.3f}]")
        new_indent = indent + ("   " if is_last else "|  ")
        print_tree(node.left, new_indent, is_last=False)
        print_tree(node.right, new_indent, is_last=True)

def get_leaf_order(node: UPGMANode) -> list:
    if node.is_leaf():
        return [node.seq_index]
    return get_leaf_order(node.left) + get_leaf_order(node.right)

if __name__ == "__main__":
    from core.distance_matrix import build_distance_matrix, print_distance_matrix

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

    print("\nGuide Tree (UPGMA)")
    tree = build_upgma_tree(dist_matrix, ids)
    print_tree(tree)

    print("\nOrder leaf (post-order):", get_leaf_order(tree))
