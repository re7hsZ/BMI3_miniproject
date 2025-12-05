# HGT Detection with Hidden Markov Models

Three-state HMM pipeline for detecting horizontally transferred genes in bacterial genomes, using composition and context features on real datasets (near / moderate / distant donors).

## Overview
- **States:** Host, Ameliorated, Foreign
- **Emissions:** Gaussian (GC/GC3, RSCU, AA composition) + Bernoulli (GC-shift flag, mobility/effector keywords)
- **Training:** Supervised initialization with optional Baum-Welch refinement
- **Outputs:** Per-gene probabilities, foreign flag, and benchmark plots (confusion matrix, ROC, heatmap, barplot)

## Performance

| Distance | Accuracy | Foreign Recall | Foreign Precision |
|----------|----------|----------------|-------------------|
| Near     | 96%      | 95%            | 100%              |
| Moderate | 87%      | 88%            | 96%               |
| Distant  | 95%      | 95%            | 100%              |

## Installation
```bash
pip install -r requirements.txt
```

## Project Structure
```
miniproject/
├── src/                         # Core code
│   ├── features.py              # Feature extraction (GC, RSCU, AA, context)
│   ├── hmm.py                   # HMM implementation (train, Viterbi, Forward-Backward)
│   ├── simulator.py             # Synthetic genome simulation
│   ├── benchmark.py             # Evaluation metrics and plots
│   └── main.py                  # CLI for train/predict/simulate
├── data/
│   ├── raw/                     # Raw CDS files (near/moderate/distant donors)
│   ├── processed/               # Cleaned data by distance
│   │   ├── host/                # Host genome CDS (train/test split)
│   │   ├── near/donor/          # Near-distance donor CDS
│   │   ├── moderate/donor/      # Moderate-distance donor CDS
│   │   └── distant/donor/       # Distant donor CDS
│   └── scripts/                 # Data preprocessing scripts
│       ├── host_cds_filtered.py
│       ├── donor_cds_filtered.py
│       └── generate_hgt_test.py # Generate test set with complete foreign genes
├── results/real/                # Model outputs by distance
│   ├── near/
│   ├── moderate/
│   └── distant/
├── frontend/                    # Flask dashboard (optional)
├── tests/                       # Unit tests
└── docs/documentation.md        # Technical documentation
```

## Usage

### Step 1: Data Preparation (run once)

```bash
# Filter host CDS
python data/scripts/host_cds_filtered.py

# Filter donor CDS for each distance
python data/scripts/donor_cds_filtered.py --distance near
python data/scripts/donor_cds_filtered.py --distance moderate
python data/scripts/donor_cds_filtered.py --distance distant

# Generate test sets with complete foreign genes
python data/scripts/generate_hgt_test.py --distance near
python data/scripts/generate_hgt_test.py --distance moderate
python data/scripts/generate_hgt_test.py --distance distant
```

### Step 2: Train Model

```bash
python src/main.py train \
    --host data/processed/host/host_core_train.fasta \
    --foreign data/processed/distant/donor/foreign_train.fasta \
    --output results/real/distant/hmm_model.pkl
```

**Training options:**
- `--bw_iters N` - Baum-Welch iterations (default: 0, disabled to preserve supervised signal)
- `--context_weight W` - Weight for context features (default: 1.2)
- `--foreign_prior_floor F` - Minimum prior for foreign state (default: 0.4)

### Step 3: Predict

```bash
python src/main.py predict \
    --input data/processed/distant/host_test_with_hgt.fasta \
    --model results/real/distant/hmm_model.pkl \
    --output results/real/distant/predictions.tsv
```

**Prediction options:**
- `--foreign_threshold T` - Threshold for ForeignFlag (default: 0.05)
- `--foreign_bias B` - Additive bias to Prob_Foreign (default: 0.0)
- `--foreign_log_boost L` - Log-space boost to foreign (default: 0.0)

### Step 4: Benchmark

```bash
python src/benchmark.py \
    --predictions results/real/distant/predictions.tsv \
    --truth data/processed/distant/hgt_truth.tsv \
    --output results/real/distant
```

**Outputs:** `confusion_matrix.png`, `roc_curve.png`, `heatmap.png`, `barplot.png`, `confusion_matrix_report.txt`

## Quick Run (All Distances)

```bash
# Near
python src/main.py train --host data/processed/host/host_core_train.fasta --foreign data/processed/near/donor/foreign_train.fasta --output results/real/near/hmm_model.pkl
python src/main.py predict --input data/processed/near/host_test_with_hgt.fasta --model results/real/near/hmm_model.pkl --output results/real/near/predictions.tsv
python src/benchmark.py --predictions results/real/near/predictions.tsv --truth data/processed/near/hgt_truth.tsv --output results/real/near

# Moderate
python src/main.py train --host data/processed/host/host_core_train.fasta --foreign data/processed/moderate/donor/foreign_train.fasta --output results/real/moderate/hmm_model.pkl
python src/main.py predict --input data/processed/moderate/host_test_with_hgt.fasta --model results/real/moderate/hmm_model.pkl --output results/real/moderate/predictions.tsv
python src/benchmark.py --predictions results/real/moderate/predictions.tsv --truth data/processed/moderate/hgt_truth.tsv --output results/real/moderate

# Distant
python src/main.py train --host data/processed/host/host_core_train.fasta --foreign data/processed/distant/donor/foreign_train.fasta --output results/real/distant/hmm_model.pkl
python src/main.py predict --input data/processed/distant/host_test_with_hgt.fasta --model results/real/distant/hmm_model.pkl --output results/real/distant/predictions.tsv
python src/benchmark.py --predictions results/real/distant/predictions.tsv --truth data/processed/distant/hgt_truth.tsv --output results/real/distant
```

## Frontend Dashboard (Optional)

```bash
cd frontend && python app.py
```
Open `http://127.0.0.1:5000` → Select model → Upload FASTA → Run Analysis

## Output Format

`predictions.tsv`:
| Column | Description |
|--------|-------------|
| GeneID | Gene identifier |
| State | Predicted state (0=Host, 1=Ameliorated, 2=Foreign) |
| Prob_Host | Posterior probability of Host state |
| Prob_Ameliorated | Posterior probability of Ameliorated state |
| Prob_Foreign | Posterior probability of Foreign state |
| ForeignFlag | 1 if Prob_Foreign ≥ threshold, else 0 |

## Customization

- **Transition matrix:** `src/main.py` → `_set_priors_and_transitions()`
- **Feature weights:** `src/hmm.py` → `context_weight` parameter
- **Emission distributions:** `src/hmm.py` → `_log_emission_prob()`
- **Simulation parameters:** `src/simulator.py` → `PHYLO_SCENARIOS`

## Limitations

- Requires labeled host/foreign training sets for supervised initialization
- Context features depend on informative FASTA headers
- Moderate-distance detection is more challenging due to overlapping feature distributions
- Ameliorated genes (foreign genes that have adapted to host GC) remain difficult to detect

## Tests

```bash
python -m pytest tests/test_core.py -v
```
