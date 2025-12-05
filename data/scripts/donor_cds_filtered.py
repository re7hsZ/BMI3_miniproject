import argparse
import csv
import os
import random
from Bio import SeqIO

# Map distance label to raw donor CDS file
DONOR_RAW_MAP = {
    "near": "./data/raw/near_cds_from_genomic.fna",
    "moderate": "./data/raw/moderate_cds_from_genomic.fna",
    "distant": "./data/raw/distant_cds_from_genomic.fna",
}

HOST_TEST_FILE = "./data/processed/host/host_core_test_empty.fasta"
OUTPUT_DIR_TEMPLATE = "./data/processed/{distance}"
OUTPUT_FILTERED_TEMPLATE = OUTPUT_DIR_TEMPLATE + "/donor/foreign_filtered_cds.fasta"
OUTPUT_TRAIN_TEMPLATE = OUTPUT_DIR_TEMPLATE + "/donor/foreign_train.fasta"
OUTPUT_TEST_TEMPLATE = OUTPUT_DIR_TEMPLATE + "/donor/foreign_test.fasta"
OUTPUT_IDS_TEMPLATE = OUTPUT_DIR_TEMPLATE + "/donor/foreign_filtered_ids.txt"
OUTPUT_HOST_HGT_TEMPLATE = OUTPUT_DIR_TEMPLATE + "/host_test_with_hgt.fasta"
OUTPUT_METADATA_TEMPLATE = OUTPUT_DIR_TEMPLATE + "/hgt_insert_metadata.csv"
OUTPUT_TRUTH_TEMPLATE = OUTPUT_DIR_TEMPLATE + "/hgt_truth.tsv"

MIN_CDS_LENGTH = 150
N_EVAL_THRESHOLD = 0.05
VALID_NUCLEOTIDES = {"A", "T", "G", "C", "a", "t", "g", "c"}

TRAIN_TEST_RATIO = 3  # 3:1 split
RANDOM_SEED = 42

INSERT_LENGTH_MIN = 150  # lengthened inserts to amplify foreign signal
INSERT_LENGTH_MAX = 600
FRAME_SHIFT_RATIO = 0.15
INSERT_PROB = 0.8       # insert HGT into 80% of host test CDS
AVOID_REGION = 20       # allow more placement while avoiding edges

mobile_keywords = [
    "transposase", "insertion sequence", "IS element", "insertion element",
    "mobile element", "integrase", "recombinase", "resolvase",
    "excisionase", "invertase", "tyrosine recombinase", "serine recombinase",
    "phage", "prophage", "phage integrase", "phage tail", "tail protein",
    "capsid", "portal protein", "terminase", "lysin", "holin", "lysis",
    "head protein", "baseplate",
    "plasmid", "conjugation", "conjugal", "tra gene", "trb gene",
    "mobA", "mobB", "mobC", "relaxase", "type IV secretion", "T4SS", "Pilus",
    "genomic island", "island", "integrative", "ICE", "IME", "cargo gene",
    "att site", "attachment site", "acquired", "horizontal transfer",
    "HGT", "mobile region", "MGI", "CRISPR", "Cas"
]


def ensure_dirs(paths):
    for p in paths:
        os.makedirs(os.path.dirname(p), exist_ok=True)


def split_train_test(filtered_records, output_train, output_test):
    total_count = len(filtered_records)
    print(f"Total number of core CDS: {total_count}")

    random.seed(RANDOM_SEED)
    random.shuffle(filtered_records)

    train_count = int(total_count * TRAIN_TEST_RATIO / (TRAIN_TEST_RATIO + 1))
    test_count = total_count - train_count

    train_records = filtered_records[:train_count]
    test_records = filtered_records[train_count:]

    SeqIO.write(train_records, output_train, "fasta")
    SeqIO.write(test_records, output_test, "fasta")

    print(f"Train(3/4): {len(train_records)} CDS -> {output_train}")
    print(f"Test(1/4): {len(test_records)} CDS -> {output_test}")


def generate_insert_fragments(foreign_test_file):
    print(f"[DEBUG] generate_insert_fragments: loading {foreign_test_file}")
    foreign_records = list(SeqIO.parse(foreign_test_file, "fasta"))
    insert_fragments = []
    random.seed(RANDOM_SEED)

    for foreign_rec in foreign_records:
        foreign_seq = foreign_rec.seq
        foreign_len = len(foreign_seq)
        if foreign_len < INSERT_LENGTH_MIN:
            continue

        if random.random() < FRAME_SHIFT_RATIO:
            max_possible_len = min(INSERT_LENGTH_MAX, foreign_len)
            # avoid infinite loop when bounds are equal and divisible by 3
            attempts = 0
            insert_len = random.randint(INSERT_LENGTH_MIN, max_possible_len)
            while insert_len % 3 == 0 and attempts < 20:
                insert_len = random.randint(INSERT_LENGTH_MIN, max_possible_len)
                attempts += 1
            if insert_len % 3 == 0:
                # force a non-multiple-of-3 length
                insert_len = max(INSERT_LENGTH_MIN, min(max_possible_len, insert_len - 1))
        else:
            max_possible_len = min(INSERT_LENGTH_MAX, foreign_len)
            max_3x_len = (max_possible_len // 3) * 3
            min_3x_len = (INSERT_LENGTH_MIN // 3) * 3
            if min_3x_len < INSERT_LENGTH_MIN:
                min_3x_len += 3
            insert_len = random.randint(min_3x_len // 3, max_3x_len // 3) * 3

        max_start = foreign_len - insert_len
        start_pos = random.randint(0, max_start)
        insert_seq = foreign_seq[start_pos:start_pos + insert_len]

        insert_fragments.append({
            "seq": insert_seq,
            "source_id": foreign_rec.id,
            "insert_len": insert_len,
            "is_frame_shift": insert_len % 3 != 0,
            "start_in_foreign": start_pos,
            "foreign_original_len": foreign_len
        })

    print(f"Generated {len(insert_fragments)} insert fragments (15% frame-shift)")
    print("[DEBUG] generate_insert_fragments: done")
    return insert_fragments


def insert_hgt_into_host(host_test_file, insert_fragments):
    print(f"[DEBUG] insert_hgt_into_host: loading {host_test_file} with {len(insert_fragments)} fragments")
    host_records = list(SeqIO.parse(host_test_file, "fasta"))
    inserted_host_records = []
    metadata_list = []
    random.seed(RANDOM_SEED)

    for host_rec in host_records:
        host_id = host_rec.id
        host_seq = host_rec.seq
        host_len = len(host_seq)
        has_hgt = "no"
        reason = "randomly_not_inserted"
        foreign_source_id = ""
        insert_length = ""
        is_frame_shift = ""
        insert_position = ""
        host_final_length = host_len

        if random.random() < INSERT_PROB:
            fragment = random.choice(insert_fragments)
            insert_len = fragment["insert_len"]

            required_min_len = AVOID_REGION * 2 + insert_len
            if host_len < required_min_len:
                reason = "host_cds_too_short_for_insert"
            else:
                min_insert_pos = AVOID_REGION
                max_insert_pos = host_len - AVOID_REGION - insert_len
                if max_insert_pos < min_insert_pos:
                    reason = "insert_length_too_long_for_host"
                else:
                    insert_pos = random.randint(min_insert_pos, max_insert_pos)
                    new_host_seq = host_seq[:insert_pos] + fragment["seq"] + host_seq[insert_pos:]

                    has_hgt = "yes"
                    reason = ""
                    foreign_source_id = fragment["source_id"]
                    insert_length = insert_len
                    is_frame_shift = fragment["is_frame_shift"]
                    insert_position = insert_pos
                    host_final_length = len(new_host_seq)

                    new_host_rec = host_rec[:]
                    new_host_rec.id = f"{host_id}_hgt_inserted"
                    new_host_rec.description = (
                        f"{host_rec.description} "
                        f"[hgt=yes] [foreign_source={foreign_source_id}] "
                        f"[insert_length={insert_length}] [is_frame_shift={is_frame_shift}] "
                        f"[insert_position={insert_position}]"
                    ).strip()
                    new_host_rec.seq = new_host_seq
                    host_rec = new_host_rec

        inserted_host_records.append(host_rec)
        metadata_list.append({
            "host_id": host_id,
            "has_hgt": has_hgt,
            "reason": reason,
            "foreign_source_id": foreign_source_id,
            "insert_length": insert_length,
            "is_frame_shift": is_frame_shift,
            "insert_position": insert_position,
            "host_original_length": host_len,
            "host_final_length": host_final_length,
            "foreign_original_len": fragment["foreign_original_len"] if has_hgt == "yes" else ""
        })

    inserted_count = sum(1 for meta in metadata_list if meta["has_hgt"] == "yes")
    total_host = len(host_records)
    frame_shift_count = sum(1 for meta in metadata_list if meta["is_frame_shift"] is True)
    print(f"Inserted HGT into {inserted_count}/{total_host} host CDS")
    print(f"Frame-shift insertions: {frame_shift_count}/{inserted_count} (15% target)")
    print("[DEBUG] insert_hgt_into_host: done")
    return inserted_host_records, metadata_list


def write_metadata(metadata_list, output_file):
    headers = [
        "host_id", "has_hgt", "reason", "foreign_source_id",
        "insert_length", "is_frame_shift", "insert_position",
        "host_original_length", "host_final_length", "foreign_original_len"
    ]

    with open(output_file, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=headers)
        writer.writeheader()
        writer.writerows(metadata_list)
    print(f"Metadata saved to: {output_file}")


def write_truth_tsv(metadata_list, output_file):
    with open(output_file, "w", encoding="utf-8", newline="") as f:
        f.write("GeneID\tTrueState\n")
        for meta in metadata_list:
            gene_id = meta["host_id"] + ("_hgt_inserted" if meta["has_hgt"] == "yes" else "")
            true_state = 2 if meta["has_hgt"] == "yes" else 0
            f.write(f"{gene_id}\t{true_state}\n")
    print(f"Truth labels saved to: {output_file}")


def parse_args():
    parser = argparse.ArgumentParser(description="Filter donor CDS and generate HGT-inserted host test set.")
    parser.add_argument("--distance", choices=["near", "moderate", "distant"], default="near",
                        help="Donor evolutionary distance (selects raw CDS file and output subfolder).")
    return parser.parse_args()


def main():
    args = parse_args()
    distance = args.distance
    cds_file = DONOR_RAW_MAP.get(distance)
    if not cds_file or not os.path.exists(cds_file):
        raise FileNotFoundError(f"CDS file for distance '{distance}' not found: {cds_file}")

    output_dir = OUTPUT_DIR_TEMPLATE.format(distance=distance)
    output_filtered = OUTPUT_FILTERED_TEMPLATE.format(distance=distance)
    output_train = OUTPUT_TRAIN_TEMPLATE.format(distance=distance)
    output_test = OUTPUT_TEST_TEMPLATE.format(distance=distance)
    output_ids = OUTPUT_IDS_TEMPLATE.format(distance=distance)
    output_host_hgt = OUTPUT_HOST_HGT_TEMPLATE.format(distance=distance)
    output_metadata = OUTPUT_METADATA_TEMPLATE.format(distance=distance)
    output_truth = OUTPUT_TRUTH_TEMPLATE.format(distance=distance)

    ensure_dirs([
        output_filtered, output_train, output_test, output_ids,
        output_host_hgt, output_metadata, output_truth
    ])

    try:
        all_records = list(SeqIO.parse(cds_file, "fasta"))
        print(f"{len(all_records)} total CDS loaded from {cds_file}")

        filtered_by_length = [r for r in all_records if len(r.seq) >= MIN_CDS_LENGTH]
        print(f"{len(filtered_by_length)} CDS remain after length filtering (≥{MIN_CDS_LENGTH} bp)")

        filtered_by_seq_validity = []
        for record in filtered_by_length:
            seq = record.seq
            length = len(seq)
            if any(ch not in VALID_NUCLEOTIDES for ch in seq):
                continue
            n_ratio = seq.upper().count("N") / length if length > 0 else 0.0
            if n_ratio > N_EVAL_THRESHOLD:
                continue
            filtered_by_seq_validity.append(record)
        print(f"{len(filtered_by_seq_validity)} CDS remain after sequence validity filtering (N ratio ≤{N_EVAL_THRESHOLD*100}%)")

        filtered_by_mobile = []
        mobile_count = 0
        for record in filtered_by_seq_validity:
            description = record.description.lower()
            if any(keyword.lower() in description for keyword in mobile_keywords):
                mobile_count += 1
            else:
                filtered_by_mobile.append(record)
        print(f"{len(filtered_by_mobile)} CDS remain after removing {mobile_count} mobile element-related CDS")

        final_filtered_records = filtered_by_mobile
        final_filtered_ids = [record.id for record in final_filtered_records]

        SeqIO.write(final_filtered_records, output_filtered, "fasta")
        with open(output_ids, "w", encoding="utf-8") as f:
            f.write("\n".join(final_filtered_ids))

        split_train_test(final_filtered_records, output_train, output_test)

        insert_fragments = generate_insert_fragments(output_test)
        if not insert_fragments:
            raise ValueError("No valid insert fragments generated (check foreign test set length)")

        inserted_host_records, metadata_list = insert_hgt_into_host(HOST_TEST_FILE, insert_fragments)

        SeqIO.write(inserted_host_records, output_host_hgt, "fasta")
        print(f"Host test set with HGT saved to: {output_host_hgt}")

        write_metadata(metadata_list, output_metadata)
        write_truth_tsv(metadata_list, output_truth)

    except Exception as e:
        print(f"Error: {str(e)}")
        raise


if __name__ == "__main__":
    main()
