import os
import sys

# Add src to the path for local imports
sys.path.append(os.path.join(os.path.dirname(__file__), "..", "src"))

from features import FeaturePipeline
from hmm import HMM


def main():
    # two toy sequences
    seqs = ["ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG", "ATGAAATTTGGGCCCAAATGCTAG"]
    headers = [">g1", ">g2"]

    fp = FeaturePipeline(n_pca_components=2)
    cont_pca, ctx = fp.fit_transform(seqs, headers)
    obs = [(cont_pca[i], ctx[i]) for i in range(len(seqs))]
    h = HMM(n_states=3, n_features=cont_pca.shape[1], n_context_features=ctx.shape[1])

    # no supervised init -> baum_welch should initialize internally
    ll = h.baum_welch_train(obs, max_iters=10, verbose=True)
    print("LL:", ll)
    path = h.viterbi(obs)
    print("Viterbi:", path)


if __name__ == "__main__":
    main()
