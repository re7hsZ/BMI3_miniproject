"""
Generate proper HGT test set with complete foreign genes (not fragments).
This replaces the chimeric approach that inserts small foreign fragments into host genes.
"""
import os
import random
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

RANDOM_SEED = 42
HOST_TRAIN_FILE = "./data/processed/host/host_core_train.fasta"
HOST_TEST_FILE = "./data/processed/host/host_core_test_empty.fasta"

# Output paths
OUTPUT_DIR_TEMPLATE = "./data/processed/{distance}"

def create_proper_hgt_test(distance="distant", foreign_ratio=0.7):
    """
    Create test set by mixing:
    - Host test genes (labeled as TrueState=0)
    - Complete foreign genes from donor (labeled as TrueState=2)
    """
    foreign_train_file = f"./data/processed/{distance}/donor/foreign_train.fasta"
    foreign_test_file = f"./data/processed/{distance}/donor/foreign_test.fasta"
    output_hgt = f"./data/processed/{distance}/host_test_with_hgt_fixed.fasta"
    output_truth = f"./data/processed/{distance}/hgt_truth_fixed.tsv"
    
    random.seed(RANDOM_SEED)
    
    # Load host test genes
    host_records = list(SeqIO.parse(HOST_TEST_FILE, "fasta"))
    print(f"Loaded {len(host_records)} host test genes")
    
    # Load foreign genes (use test split for independence)
    foreign_records = list(SeqIO.parse(foreign_test_file, "fasta"))
    print(f"Loaded {len(foreign_records)} foreign test genes")
    
    # Calculate how many foreign genes to include
    n_host = len(host_records)
    n_foreign = int(n_host * foreign_ratio / (1 - foreign_ratio))
    n_foreign = min(n_foreign, len(foreign_records))
    
    # Sample foreign genes
    foreign_sample = random.sample(foreign_records, n_foreign)
    print(f"Selected {n_foreign} foreign genes for test mix")
    
    # Create output records
    output_records = []
    truth_labels = []
    
    # Add host genes (TrueState=0)
    for rec in host_records:
        output_records.append(rec)
        truth_labels.append((rec.id, 0))
    
    # Add complete foreign genes (TrueState=2) 
    for rec in foreign_sample:
        # Rename to include foreign marker
        new_rec = SeqRecord(
            rec.seq,
            id=f"{rec.id}_hgt_foreign",
            description=f"{rec.description} [hgt=yes] [origin=foreign]"
        )
        output_records.append(new_rec)
        truth_labels.append((new_rec.id, 2))
    
    # Shuffle to mix host and foreign
    combined = list(zip(output_records, truth_labels))
    random.shuffle(combined)
    output_records, truth_labels = zip(*combined)
    
    # Write output FASTA
    SeqIO.write(output_records, output_hgt, "fasta")
    print(f"Wrote {len(output_records)} records to {output_hgt}")
    
    # Write truth labels
    with open(output_truth, "w") as f:
        f.write("GeneID\tTrueState\n")
        for gene_id, state in truth_labels:
            f.write(f"{gene_id}\t{state}\n")
    print(f"Wrote truth labels to {output_truth}")
    
    # Print statistics
    n_host_final = sum(1 for _, s in truth_labels if s == 0)
    n_foreign_final = sum(1 for _, s in truth_labels if s == 2)
    print(f"\nTest set composition:")
    print(f"  Host genes: {n_host_final} ({100*n_host_final/len(truth_labels):.1f}%)")
    print(f"  Foreign genes: {n_foreign_final} ({100*n_foreign_final/len(truth_labels):.1f}%)")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--distance", default="distant", choices=["near", "moderate", "distant"])
    parser.add_argument("--foreign_ratio", type=float, default=0.7)
    args = parser.parse_args()
    create_proper_hgt_test(args.distance, args.foreign_ratio)
