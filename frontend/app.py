import io
import os
import sys
import json
from pathlib import Path
import numpy as np
from flask import Flask, render_template, jsonify, request
from Bio import SeqIO
import pickle

# make src importable
ROOT = Path(__file__).resolve().parents[1]
sys.path.append(str(ROOT / "src"))
from features import extract_features  # noqa: E402

app = Flask(__name__)

# Configuration
app.config['ENV'] = 'production'  # or 'development'
app.config['DEBUG'] = True
app.config['MAX_CONTENT_LENGTH'] = 50 * 1024 * 1024  # 50 MB max upload size

MODEL_PATHS = {
    "near": ROOT / "results" / "real" / "near" / "hmm_model.pkl",
    "moderate": ROOT / "results" / "real" / "moderate" / "hmm_model.pkl",
    "distant": ROOT / "results" / "real" / "distant" / "hmm_model.pkl",
}


def load_model(distance: str):
    model_path = MODEL_PATHS.get(distance)
    if not model_path or not model_path.exists():
        raise FileNotFoundError(f"Model for distance '{distance}' not found at {model_path}")
    with open(model_path, "rb") as f:
        saved = pickle.load(f)
    hmm = saved["model"]
    scaler = saved["scaler"]
    genomic_gc = saved.get("genomic_gc", 0.5)
    return hmm, scaler, genomic_gc


# Cache for PCA reference data
_pca_cache = {}


def load_pca_reference(distance: str, scaler, genomic_gc: float):
    """Load host and donor training data for PCA visualization."""
    if distance in _pca_cache:
        return _pca_cache[distance]
    
    host_path = ROOT / "data" / "processed" / "host" / "host_core_train.fasta"
    donor_path = ROOT / "data" / "processed" / distance / "donor" / "foreign_train.fasta"
    
    host_data = []
    donor_data = []
    
    # Sample max 200 points from each to avoid slow rendering
    max_points = 200
    
    try:
        if host_path.exists():
            host_records = list(SeqIO.parse(str(host_path), "fasta"))[:max_points]
            for rec in host_records:
                feat = extract_features(str(rec.seq), rec.description, genomic_gc)
                host_data.append(feat[0])
        
        if donor_path.exists():
            donor_records = list(SeqIO.parse(str(donor_path), "fasta"))[:max_points]
            for rec in donor_records:
                feat = extract_features(str(rec.seq), rec.description, genomic_gc)
                donor_data.append(feat[0])
        
        if host_data:
            host_scaled = scaler.transform(np.array(host_data))[:, :2]
            host_points = [{"x": float(p[0]), "y": float(p[1])} for p in host_scaled]
        else:
            host_points = []
        
        if donor_data:
            donor_scaled = scaler.transform(np.array(donor_data))[:, :2]
            donor_points = [{"x": float(p[0]), "y": float(p[1])} for p in donor_scaled]
        else:
            donor_points = []
        
        _pca_cache[distance] = (host_points, donor_points)
        return host_points, donor_points
    
    except Exception as e:
        print(f"Warning: Could not load PCA reference data: {e}")
        return [], []


def parse_fasta_text(text: str):
    handle = io.StringIO(text)
    records = list(SeqIO.parse(handle, "fasta"))
    if not records:
        raise ValueError("No valid FASTA records found.")
    return records


def run_prediction(fasta_text: str, distance: str, foreign_threshold: float = 0.5):
    hmm, scaler, genomic_gc = load_model(distance)
    records = parse_fasta_text(fasta_text)

    headers = [r.description for r in records]
    seqs = [str(r.seq) for r in records]

    feats = [extract_features(s, h, genomic_gc) for s, h in zip(seqs, headers)]
    comp = np.array([f[0] for f in feats])
    ctx = np.array([f[1] for f in feats])
    comp_scaled = scaler.transform(comp)

    observations = [(comp_scaled[i], ctx[i]) for i in range(len(comp_scaled))]
    path = hmm.viterbi(observations)
    post = hmm.posterior(observations)
    foreign_probs = post[:, 2]
    foreign_flags = foreign_probs >= foreign_threshold

    # Count predictions
    total_genes = len(records)
    foreign_count = int(np.sum(foreign_flags))
    host_count = total_genes - foreign_count

    # basic XY using first two dimensions (or pad zeros)
    if comp_scaled.shape[1] < 2:
        xys = np.hstack([comp_scaled, np.zeros((comp_scaled.shape[0], 2 - comp_scaled.shape[1]))])
    else:
        xys = comp_scaled[:, :2]

    # Load host and donor reference points for PCA visualization
    host_points, donor_points = load_pca_reference(distance, scaler, genomic_gc)

    # Separate input points by predicted class for visualization
    input_host_points = [{"x": float(xys[i, 0]), "y": float(xys[i, 1])} 
                         for i in range(len(records)) if not foreign_flags[i]]
    input_foreign_points = [{"x": float(xys[i, 0]), "y": float(xys[i, 1])} 
                            for i in range(len(records)) if foreign_flags[i]]

    pca_data = {
        "host": host_points,
        "donor": donor_points,
        "inputHost": input_host_points,
        "inputForeign": input_foreign_points,
    }

    bar_data = [{"gene": records[i].id, "prob": float(foreign_probs[i])} for i in range(len(records))]
    heatmap_data = foreign_probs.tolist()

    result = {
        "totalGenes": total_genes,
        "foreignCount": foreign_count,
        "hostCount": host_count,
        "foreignPercent": round(100 * foreign_count / total_genes, 1) if total_genes > 0 else 0,
        "hostPercent": round(100 * host_count / total_genes, 1) if total_genes > 0 else 0,
        "pcaData": pca_data,
        "heatmapData": heatmap_data,
        "barData": bar_data,
        "details": [
            {
                "GeneID": records[i].id,
                "State": int(path[i]),
                "Prob_Host": float(post[i, 0]),
                "Prob_Ameliorated": float(post[i, 1]),
                "Prob_Foreign": float(post[i, 2]),
                "ForeignFlag": int(foreign_flags[i]),
            }
            for i in range(len(records))
        ],
    }
    return result


@app.route('/')
def index():
    """Render the main page."""
    return render_template('index.html')


@app.route('/api/v1/predict', methods=['POST'])
def predict_api():
    try:
        distance = request.form.get("distance", "near")
        foreign_threshold = float(request.form.get("foreign_threshold", 0.5))

        fasta_text = ""
        if "file" in request.files and request.files["file"].filename:
            fasta_text = request.files["file"].read().decode("utf-8", errors="ignore")
        elif request.form.get("fasta"):
            fasta_text = request.form.get("fasta")
        else:
            return jsonify({"error": "No FASTA provided"}), 400

        result = run_prediction(fasta_text, distance, foreign_threshold=foreign_threshold)
        return jsonify(result)
    except Exception as e:
        return jsonify({"error": str(e)}), 400


if __name__ == '__main__':
    print("Starting HGT-Detector Server...")
    print("Open http://127.0.0.1:5000 in your browser.")
    app.run(port=5000)
