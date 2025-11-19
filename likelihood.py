import numpy as np
import pandas as pd
import math
from tree import parse_test_data
import time
from functools import lru_cache


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

    @lru_cache(maxsize=None)
    def get_transition_prob(self, branch_length: float) -> np.ndarray:
        """
        JC69 closed-form transition matrix for a branch of length t.

        P_ii = 1/4 + 3/4 * e
        P_ij = 1/4 - 1/4 * e  (i != j)
        where e = exp(-4 * u * t)
        """
        u = self.u
        t = branch_length

        # scalar
        e = math.exp(-4.0 * u * t)

        p_same = 0.25 + 0.75 * e  # diagonal entries
        p_diff = 0.25 - 0.25 * e  # off-diagonal entries

        # build 4x4 matrix: start with all p_diff, then fix diagonal
        P = np.full((4, 4), p_diff, dtype=float)
        np.fill_diagonal(P, p_same)

        return P

    def sequence_to_vector(self, nucleotide):
        """
        Convert a nucleotide to a probability vector.
        """
        vec = np.zeros(4)
        if nucleotide in self.nuc_to_idx:
            vec[self.nuc_to_idx[nucleotide]] = 1.0
        return vec

    def compute_conditional_likelihood(self, node, site_idx):
        # Leaf
        if not node.children:
            if node.sequence and site_idx < len(node.sequence):
                return self.sequence_to_vector(node.sequence[site_idx])
            else:
                return np.ones(4)

        cond_like = None

        for child in node.children:
            child_like = self.compute_conditional_likelihood(child, site_idx)
            P = self.get_transition_prob(child.branch_length)
            child_contrib = P @ child_like

            if cond_like is None:
                # first child: just take its contribution
                cond_like = child_contrib
            else:
                # subsequent children: multiply in
                cond_like *= child_contrib  # in-place

        return cond_like

    def compute_site_likelihood(self, tree, site_idx, root):
        """
        Compute likelihood for a single site in the alignment.
        """
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
            site_log_like = self.compute_site_likelihood(
                tree, site_idx, tree.find_root()
            )
            total_log_like += site_log_like

        return total_log_like


df = pd.DataFrame(columns=["run_id", "run_time", "type"])
for i in range(100):
    start = time.perf_counter()
    calc = LikelihoodCalculator(u=0.3)
    # Compute likelihood for the tree
    T = parse_test_data()
    total_log_likelihood = calc.compute_tree_likelihood(T)
    end = time.perf_counter()
    df.loc[i] = [i + 1, end - start, "optimized"]
df.to_csv("optimized_likelihood_times.csv", index=False)
