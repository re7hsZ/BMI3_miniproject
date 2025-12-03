# HGT Detection with Hidden Markov Models

Three-state HMM pipeline for detecting horizontally transferred genes in bacterial genomes, using composition and context features with both synthetic evaluation data and real datasets.

## Overview
- States: Host, Ameliorated, Foreign
- Emissions: Gaussian (GC/GC3, RSCU, AA composition) + Bernoulli (GC-shift flag, mobility/effector keywords)
- Transition bias: Host -> Ameliorated -> Foreign with high self-transition
- Outputs: per-gene probabilities, foreign flag, and benchmark plots (confusion matrix, ROC, heatmap, barplot)

## Installation
```bash
pip install -r requirements.txt
cd miniproject
```

## Quick Start (synthetic evaluation)
```bash
python run_demo.py
```
Inputs: `data/synthetic_eval/demo/`  
Outputs: `results/synthetic_eval/demo/latest/` (model, predictions, plots)

## Project Structure
```
miniproject/
├── src/                         # Core code (latest, no suffix)
│   ├── features.py
│   ├── hmm.py
│   ├── simulator.py
│   ├── benchmark.py
│   └── main.py                  # CLI for train/predict/simulate
├── data/
│   ├── raw/                     # Real CDS inputs
│   ├── processed/               # Cleaned real data (host/donor, truth labels)
│   └── synthetic_eval/demo/     # Synthetic evaluation FASTA files
├── results/
│   ├── synthetic_eval/demo/latest/  # Synthetic model/predictions/plots
│   └── real/                    # Real run outputs
├── docs/documentation.md        # Technical overview
├── tests/                       # Unit/sanity tests
└── run_demo.py                  # Synthetic pipeline
```

## Usage

### Synthetic (manual steps)
```bash
python src/main.py train \
    --host data/synthetic_eval/demo/train_host.fasta \
    --foreign data/synthetic_eval/demo/train_foreign.fasta \
    --output results/synthetic_eval/demo/latest/model.pkl

python src/main.py predict \
    --input data/synthetic_eval/demo/test_genome.fasta \
    --model results/synthetic_eval/demo/latest/model.pkl \
    --output results/synthetic_eval/demo/latest/predictions.tsv

python src/benchmark.py \
    --predictions results/synthetic_eval/demo/latest/predictions.tsv \
    --output results/synthetic_eval/demo/latest \
    --fasta data/synthetic_eval/demo/test_genome.fasta
```

### Real data
Prep (once):
```bash
python data/scripts/host_cds_filtered.py
python data/scripts/donor_cds_filtered.py   # writes data/processed/hgt_truth.tsv
```

Run:
```bash
python src/main.py train \
    --host data/processed/host/host_core_train.fasta \
    --foreign data/processed/donor/donor_core_train.fasta \
    --output results/real/hmm_model.pkl

python src/main.py predict \
    --input data/processed/host_test_with_hgt.fasta \
    --model results/real/hmm_model.pkl \
    --output results/real/predictions.tsv \
    --foreign_threshold 0.5

python src/benchmark.py \
    --predictions results/real/predictions.tsv \
    --truth data/processed/hgt_truth.tsv \
    --output results/real \
    --fasta data/processed/host_test_with_hgt.fasta
```

Notes:
- Train only on original genomes (no HGT tags); `main.py` checks headers.
- Predictions include `ForeignFlag` (1 if `Prob_Foreign` >= threshold).
- `simulate` supports `--phylo_scenario {near,mid,far}` to choose GC targets matching biologically plausible donors.

## Customization
- Transition/context weights: `src/hmm.py`
- Feature thresholds/PCA components: `src/features.py`
- Simulation GC targets/islands: `src/simulator.py`

## Output Format
`predictions.tsv`: GeneID, State (0/1/2), Prob_Host, Prob_Ameliorated, Prob_Foreign, ForeignFlag

## Limitations
- Depends on labeled host/foreign sets for good initialization.
- Context features rely on informative headers.
- Highly ameliorated genes remain challenging.
