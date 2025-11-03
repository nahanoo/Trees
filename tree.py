import pandas as pd
import os.path as path
import numpy as np 
import warnings
from scipy.linalg import expm

# Module level constant with one-time excutation (memory efficient)
# All ascii encoding
_DNA_MAP = np.full(256, -1, dtype=np.int8)
for i, ch in enumerate(b"ACGTacgt"):
    _DNA_MAP[ch] = i % 4
_ID4 = np.eye(4, dtype=np.float32)  # float32 is fine for CLVs

def dna_onehot_4xN(seq: str) -> np.ndarray:
    if not seq:
        return np.zeros((4, 0), dtype=np.float32)
    b = np.frombuffer(seq.encode('ascii', 'strict'), dtype=np.uint8)
    idx = _DNA_MAP[b]              # -1 for non-ACGT
    # Map unknowns to zero rows (all zeros)
    ok = idx >= 0
    onehot = np.zeros((b.size, 4), dtype=np.float32)
    if ok.any():
        onehot[ok] = _ID4[idx[ok]]
    return onehot.T  # (4, N)

class Node:
    def __init__(self, name):
        """Initialize a node with a name."""
        self.name = name
        # Only for leaf nodes
        self.sequence = None
        self.branch_length = None
        self.parent = None
        # List of child nodes
        self.children = []
        self.seq_matrix = None
        self.log_scales = None 

    def get_probability_matrix(self, Q: np.ndarray) -> np.ndarray:
        """
        Returns (4, N) conditional likelihood-like matrix at this node.
        For leaves: one-hot. For internal nodes: product over children of (P(t) @ child_matrix).
        Caches result in self.seq_matrix.
        """
        # Already computed?
        if self.seq_matrix is not None:
            return self.seq_matrix

        # Leaf case: build one-hot from sequence
        if not self.children:
            if self.sequence is None:
                raise ValueError(f"Leaf node {self.name} has no sequence.")
            M = dna_onehot_4xN(self.sequence)

            # scaling to avoid too-small number
            s = M.sum(axis=0, keepdims=True)          # (1, N)
            s[s == 0] = 1.0
            M /= s
            self.seq_matrix = M
            self.log_scales = np.log(s, dtype=np.float32).squeeze()  # (N,)

            return self.seq_matrix

        # Internal node: ensure all children are computed, then combine
        if Q.shape != (4, 4):
            raise ValueError(f"Q must be (4,4), got {Q.shape}.")

        evom_list = []
        for child in self.children:

            child_mat = child.get_probability_matrix(Q)  # recursion ensures child seq_matrix

            # Sanity to (4,N)
            if child_mat.shape[0] != 4 and child_mat.shape[1] == 4:
                warnings.warn(f"{child.name} has (N,4); transposing to (4,N).", RuntimeWarning)
                child_mat = child_mat.T
            if child_mat.shape[0] != 4:
                raise ValueError(f"{child.name} matrix must be (4,N) or (N,4).")

            # Evolve child up the edge
            t = float(child.branch_length)
            P = expm(Q * t).astype(np.float32, copy=False)  # (4,4)
            evom_list.append(P @ child_mat)                 # (4,N)

        # Check all N agree
        N = evom_list[0].shape[1]
        if not all(m.shape == (4, N) for m in evom_list):
            shapes = [m.shape for m in evom_list]
            raise ValueError(f"Children shapes mismatch: {shapes}")

        # product and scaling
        M = np.prod(np.stack(evom_list, axis=0), axis=0)  # (k, N)
        s = M.sum(axis=0, keepdims=True)  # (1, N)
        s[s == 0] = 1.0
        M /= s
        self.seq_matrix = M

        # Elementwise product across children (independent)
        self.log_scales = np.log(s, dtype=np.float32).squeeze()  # (N,)

        # Add previous ones
        for child in self.children:
            if child.log_scales is not None:
                self.log_scales += child.log_scales

        return self.seq_matrix


class Tree:
    def __init__(self, tree_f, seq_f, branch_length_f):
        """Read tree structure, sequences, and branch lengths from files."""
        self.tree = pd.read_csv(tree_f, names=["parent", "child"])
        self.branch_length = pd.read_csv(
            branch_length_f,
            header=None,
        ).transpose()
        self.branch_length.columns = ["branch_length"]
        self.sequences = pd.read_csv(seq_f, names=["name", "sequence"], sep=" ")
        self.populate_tree()

    def populate_tree(self):
        # Always treat names as strings
        df = self.tree.astype(str).copy()

        names = pd.unique(df[["parent", "child"]].values.ravel())
        # Create all nodes
        self.nodes = {name: Node(name) for name in names}

        # Sequences as strings too
        self.sequences = self.sequences.astype({"name": str, "sequence": str})

        # Add sequences to leaf nodes
        seq_map = dict(zip(self.sequences["name"], self.sequences["sequence"]))
        for name, seq in seq_map.items():
            if name in self.nodes:
                self.nodes[name].sequence = seq

        # Attach branch lengths and add children to parents
        for i, row in df.reset_index(drop=True).iterrows():
            p = self.nodes[row["parent"]]  # strings now
            c = self.nodes[row["child"]]
            c.parent = p
            c.branch_length = float(self.branch_length.iloc[i]["branch_length"])
            p.children.append(c)

    def tips(self):
        """Get all tip (leaf) nodes."""
        return [n for n in self.nodes.values() if not n.children]

    def find_root(self):
        """Find the root node of the tree."""
        parents = set(self.tree["parent"].astype(str))
        children = set(self.tree["child"].astype(str))
        roots = list(parents - children)
        return self.nodes[roots[0]]

    def tree_loglikelihood(self, Q):
        root = self.find_root()
        Lr = root.get_probability_matrix(Q=Q)

        pi_vec = np.full(4, 0.25, dtype=np.float64)
        site_like = pi_vec @ Lr
        log_scales = root.log_scales if root.log_scales is not None else 0.0
        root_ll = float(np.sum(np.log(site_like) + log_scales))
        return root_ll

if __name__ == "__main__":
    def parse_test_data():
        tree_f = path.join("dataset", "ENSG00000013016_EHD3_NT.table.dat")
        seq_f = path.join("dataset", "ENSG00000013016_EHD3_NT.msa.dat")
        branch_length_f = path.join("dataset", "ENSG00000013016_EHD3_NT.branchlength.dat")
        return Tree(tree_f, seq_f, branch_length_f)

    T = parse_test_data()
