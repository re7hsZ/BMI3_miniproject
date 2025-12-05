# Software Documentation: Foreign Gene Detection with HMM

## 1. Major Task
The goal of this software is to detect horizontally transferred (foreign) genes in the *Legionella pneumophila* genome. It implements a Hidden Markov Model (HMM) that utilizes compositional features of DNA sequences—specifically GC content, Relative Synonymous Codon Usage (RSCU), and Amino Acid composition—to classify genes into three states: Host, Ameliorated (intermediate), and Foreign.

## 2. Input/Output
**Input**:
- A FASTA file containing the coding sequences (CDS) of genes to be analyzed.
- (Optional) FASTA files of known Host and Foreign genes for training/initialization.

**Output**:
- A Tab-Separated Values (TSV) file containing:
    - `GeneID`: The identifier of the gene.
    - `State`: The predicted state (0=Host, 1=Ameliorated, 2=Foreign).
    - `Prob_Host`, `Prob_Ameliorated`, `Prob_Foreign`: Posterior probabilities for each state.

## 3. Key Algorithms
The core algorithm is a 3-state Hidden Markov Model (HMM).
- **States**:
    - **Host (0)**: Represents native genes with typical codon usage of the host.
    - **Ameliorated (1)**: Represents foreign genes that have started to adapt to the host's codon usage.
    - **Foreign (2)**: Represents recently acquired genes with distinct codon usage.
- **Emissions**: Modeled as Multivariate Gaussian distributions over the feature space (GC, GC3, RSCU, AA).
- **Transitions**: A transition matrix that favors self-transitions to model the biological reality that HGT events often transfer contiguous blocks of genes ("islands").
- **Decoding**: The Viterbi algorithm is used to determine the most likely sequence of states for the entire genome.
- **Training**: The model parameters (means, covariances) are initialized using supervised learning on reference datasets (Host vs Foreign) and can be further refined.

## 4. Benchmark
The tool includes a benchmarking script (`src/benchmark.py`) that evaluates performance using simulated or real datasets.
- **Metrics**: Confusion Matrix, ROC AUC, Precision, Recall, F1-score.
- **Visualization**:
    - **ROC Curve**: Shows the trade-off between True Positive Rate and False Positive Rate.
    - **PCA Plot**: Visualizes the gene features in 2D space, colored by predicted state, to demonstrate the separation of Host and Foreign genes.

## 5. Setting up Running Environment
**Requirements**:
- Python 3.7+
- Packages: `numpy`, `pandas`, `biopython`, `scikit-learn`, `matplotlib`, `seaborn`

**Installation**:
```bash
pip install numpy pandas biopython scikit-learn matplotlib seaborn
```

## 6. Test Running Instructions (real data by distance)
- Prepare datasets:
  ```bash
  python data/scripts/host_cds_filtered.py
  python data/scripts/donor_cds_filtered.py --distance near
  python data/scripts/donor_cds_filtered.py --distance moderate
  python data/scripts/donor_cds_filtered.py --distance distant
  ```
- **Train** (example: near; swap distance to moderate/distant as needed):
  ```bash
  python src/main.py train \
      --host data/processed/host/host_core_train.fasta \
      --foreign data/processed/near/donor/foreign_train.fasta \
      --output results/real/near/hmm_model.pkl \
      --bw_iters 10 \
      --context_weight 1.2 \
      --foreign_prior_floor 0.4
  ```
- **Predict**:
  ```bash
  python src/main.py predict \
      --input data/processed/near/host_test_with_hgt.fasta \
      --model results/real/near/hmm_model.pkl \
      --output results/real/near/predictions.tsv \
      --foreign_threshold 0.2 \
      --foreign_bias 0.2
  ```
- **Benchmark (with truth)**:
  ```bash
  python src/benchmark.py \
      --predictions results/real/near/predictions.tsv \
      --truth data/processed/near/hgt_truth.tsv \
      --output results/real/near \
      --fasta data/processed/near/host_test_with_hgt.fasta
  ```

## 7. Parameters and Explanation
- `--host`: Path to FASTA file with known host genes. Used to estimate emission parameters for State 0.
- `--foreign`: Path to FASTA file with known foreign genes. Used to estimate emission parameters for State 2.
- `--input`: Path to the target genome FASTA file.
- `--model`: Path to save/load the pickled HMM object.
- `--output`: Path for the results file.

## 8. README
See `README.md` in the root directory for a quick start guide.
