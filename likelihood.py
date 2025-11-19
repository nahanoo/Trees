import numpy as np
import scipy.linalg
import math
from tree import parse_test_data


class LikelihoodCalculator:
    def __init__(self, u=0.3):
        """
        Initialize the likelihood calculator with substitution rate.
        """
        self.u = u
        # JC69 Q matrix for 4 nucleotides (A, C, G, T)
        self.Q = np.array(
            [[-3 * u, u, u, u], [u, -3 * u, u, u], [u, u, -3 * u, u], [u, u, u, -3 * u]]
        )
        # Base frequencies
        self.pi = np.array([0.25, 0.25, 0.25, 0.25])
        # Nucleotide to index mapping
        self.nuc_to_idx = {"A": 0, "C": 1, "G": 2, "T": 3}

    def get_transition_prob(self, branch_length):
        """
        Compute transition probability matrix P(t) = exp(Q*t).
        """
        return scipy.linalg.expm(self.Q * branch_length)

    def sequence_to_vector(self, nucleotide):
        """
        Convert a nucleotide to a probability vector.
        """
        vec = np.zeros(4)
        if nucleotide in self.nuc_to_idx:
            vec[self.nuc_to_idx[nucleotide]] = 1.0
        return vec

    def compute_conditional_likelihood(self, node, site_idx):
        """
        Recursively compute conditional likelihood at a node for a given site.
        """
        # Base case: leaf node
        if not node.children:  # if node has no children it's a leaf
            if node.sequence and site_idx < len(node.sequence):
                return self.sequence_to_vector(node.sequence[site_idx])
            else:
                # Missing data -
                return np.ones(4)

        # Recursive case: internal node
        # Start with all ones
        cond_like = np.ones(
            4
        )  # before seeing any children all ancestral states are equally possible

        for child in node.children:
            # Get conditional likelihood from child
            child_like = self.compute_conditional_likelihood(child, site_idx)

            # Get transition probability matrix for this branch
            P = self.get_transition_prob(child.branch_length)

            # Compute contribution from this child
            # For each ancestral state i (A, C, G or T), what's the probability?
            child_contrib = np.matmul(
                P, child_like
            )  # P1 × [1,0,0,0] #multiply matrix P by vector child_like, vec1

            # Multiply contributions from all children
            cond_like *= child_contrib  # vec1 × vec2 (element-wise)

        return cond_like

    def compute_site_likelihood(self, tree, site_idx):
        """
        Compute likelihood for a single site in the alignment.
        """
        root = tree.find_root()

        # Get conditional likelihood at root
        root_like = self.compute_conditional_likelihood(root, site_idx)

        # Multiply by base frequencies and sum
        site_like = np.matmul(root_like, self.pi)

        return math.log(site_like) if site_like > 0 else -float("inf")

    def compute_tree_likelihood(self, tree):
        """
        Compute total log-likelihood for the entire tree and alignment.
        """
        # Get alignment length from first tip
        tips = tree.tips()
        if not tips or not tips[0].sequence:
            return 0.0

        alignment_length = len(tips[0].sequence)

        # Sum log-likelihoods over all sites
        total_log_like = 0.0
        for site_idx in range(alignment_length):
            site_log_like = self.compute_site_likelihood(tree, site_idx)
            total_log_like += site_log_like

        return total_log_like


calc = LikelihoodCalculator(u=0.3)

# Compute likelihood for the tree
T = parse_test_data()
total_log_likelihood = calc.compute_tree_likelihood(T)
print(f"Total log-likelihood: {total_log_likelihood}")

# You can also compute likelihood for a single site
site_log_like = calc.compute_site_likelihood(T, site_idx=0)
