import pandas as pd
import os.path as path


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


def parse_test_data():
    tree_f = path.join("dataset", "ENSG00000013016_EHD3_NT.table.dat")
    seq_f = path.join("dataset", "ENSG00000013016_EHD3_NT.msa.dat")
    branch_length_f = path.join("dataset", "ENSG00000013016_EHD3_NT.branchlength.dat")
    return Tree(tree_f, seq_f, branch_length_f)


T = parse_test_data()
