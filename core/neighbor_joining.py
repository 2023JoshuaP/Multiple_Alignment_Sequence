"""
Neighbor-Joining para la construcción del
Árbol Guía en algoritmos progrsivos (ClustalW)
"""

import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from dataclasses import dataclass, field

@dataclass
class NJNode:
    left: "NJNode" = None
    right: "NJNode" = None
    height: float = 0.0
    seq_index: int = None
    size: int = 1
    label: str = ""

    def is_leaf(self) -> bool:
        return self.left is None and self.right is None

    def __repr__(self):
        if self.is_leaf():
            return f"Leaf({self.label})"
        return f"Node({self.label}, branch_length={self.height:.3f})"

def _get_distance(active_distances: dict, a: int, b: int) -> float:
    key = (a, b) if a < b else (b, a)
    return active_distances[key]

def _set_distance(active_distances: dict, a: int, b: int, value: float) -> None:
    key = (a, b) if a < b else (b, a)
    active_distances[key] = value

def build_nj_tree(distance_matrix: list[list[float]], labels: list[str] = None) -> NJNode:
    n = len(distance_matrix)
    if n < 2:
        raise ValueError("At least two sequences are required to construct an NJ Tree")
    
    for row in distance_matrix:
        if len(row) != n:
            raise ValueError("The distance_matrix mut be quadratic")

    if labels is None:
        labels = [str[i] for i in range(n)]

    clusters: dict = {
        i: NJNode(seq_index=i, height=0.0, size=1, label=labels[i])
        for i in range(n)
    }

    active_distances: dict = {}
    for i in range(n):
        for j in range(i + 1, n):
            active_distances[(i, j)] = distance_matrix[i][j]

    next_id = n

    while len(clusters) > 2:
        r = len(clusters)
        active_ids = list(clusters.keys())

        u = {}
        for i in active_ids:
            sum_d = 0.0
            for j in active_ids:
                if i != j:
                    sum_d += _get_distance(active_distances, i, j)
            u[i] = sum_d / (r - 2)

        best_pair = None
        min_q = float('inf')

        for i_idx in range(len(active_ids)):
            for j_idx in range(i_idx + 1, len(active_ids)):
                i = active_ids[i_idx]
                j = active_ids[j_idx]
                d_ij = _get_distance(active_distances, i, j)
                q_ij = d_ij - u[i] - u[j]

                if q_ij < min_q:
                    min_q = q_ij
                    best_pair = (i, j, d_ij)

        i, j, d_ij = best_pair
        node_i = clusters[i]
        node_j = clusters[j]

        branch_length_i = 0.5 * d_ij + 0.5 * (u[i] - u[j])
        branch_length_j = d_ij - branch_length_i

        node_i.height = branch_length_i
        node_j.height = branch_length_j

        new_size = node_i.size + node_j.size
        new_node = NJNode(left=node_i, right=node_j, size=new_size, label=f"({node_i.label}, {node_j.label}")

        remaining_ids = [cid for cid in active_ids if cid not in (i, j)]
        for k in remaining_ids:
            d_ik = _get_distance(active_distances, i, k)
            d_jk = _get_distance(active_distances, j, k)
            new_d = 0.5 * (d_ik + d_jk - d_ij)
            _set_distance(active_distances, next_id, k, new_d)

        active_distances = {pair: dist for pair, dist in active_distances.items() if i not in pair and j not in pair}

        del clusters[i]
        del clusters[j]
        clusters[next_id] = new_node
        next_id += 1

    if len(clusters) == 2:
        i, j = list(clusters.keys())
        d_ij = _get_distance(active_distances, i, j)

        node_i = clusters[i]
        node_j = clusters[j]

        branch_length_i = d_ij / 2.0
        branch_length_j = d_ij / 2.0

        node_i.height = branch_length_i
        node_j.height = branch_length_j

        root_node = NJNode(left=node_i, right=node_j, size=node_i.size + node_j.size, label=f"{node_i.label}, {node_j.label}")
        return root_node
    else:
        return list(clusters.values())[0]

def print_nj_tree(node: NJNode, indent: str = "", is_last: bool = True) -> None:
    connector = "└── " if is_last else "├── "
    if node.is_leaf():
        print(f"{indent}{connector}{node.label} (branch={node.height:.3f})")
    else:
        print(f"{indent}{connector}[branch={node.height:.3f}]")
        new_indent = indent + ("   " if is_last else "|  ")
        print_nj_tree(node.left, new_indent, is_last=False)
        print_nj_tree(node.right, new_indent, is_last=True)

def get_leaf_order(node: NJNode) -> list:
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

    print("\nGuide Tree (Neighbor-Joining)")
    tree = build_nj_tree(dist_matrix, ids)
    print_nj_tree(tree)

    print("\nOrder leaf (post-order):", get_leaf_order(tree))
