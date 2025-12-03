import os
import sys
import unittest
import numpy as np

# Add src to path
sys.path.append(os.path.join(os.path.dirname(__file__), "..", "src"))

from features import calculate_gc_content, calculate_rscu_dict
from hmm import HMM


class TestFeatures(unittest.TestCase):
    def test_gc_content_fraction(self):
        self.assertEqual(calculate_gc_content("GCGC"), 1.0)
        self.assertEqual(calculate_gc_content("ATAT"), 0.0)
        self.assertEqual(calculate_gc_content("ATGC"), 0.5)

    def test_rscu_dict_basic(self):
        seq = "ATG"  # only Met
        rscu = calculate_rscu_dict(seq)
        self.assertAlmostEqual(rscu["ATG"], 1.0)

        seq = "TTTTTT"  # two Phe codons, all TTT
        rscu = calculate_rscu_dict(seq)
        self.assertAlmostEqual(rscu["TTT"], 2.0)
        self.assertAlmostEqual(rscu["TTC"], 0.0)


class TestHMM(unittest.TestCase):
    def test_initialization(self):
        hmm = HMM(n_states=3)
        self.assertEqual(hmm.n_states, 3)
        self.assertTrue(np.allclose(np.sum(hmm.A, axis=1), 1.0))
        self.assertTrue(np.allclose(np.sum(hmm.pi), 1.0))

    def test_viterbi_shape(self):
        hmm = HMM(n_states=3)
        hmm.n_features = 2
        hmm.means = np.zeros((3, 2))
        hmm.covars = np.stack([np.eye(2) * 0.1 for _ in range(3)])

        # simple observations without context channel
        observations = [(np.zeros(2), None) for _ in range(5)]
        path = hmm.viterbi(observations)

        self.assertEqual(len(path), 5)
        self.assertTrue(np.all(path >= 0))
        self.assertTrue(np.all(path < 3))


if __name__ == "__main__":
    unittest.main()
