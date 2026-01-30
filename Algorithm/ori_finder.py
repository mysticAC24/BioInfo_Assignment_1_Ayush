# ori_finder.py
from Bio import SeqIO
import math

# ---------- Utilities ----------

def is_low_complexity(motif):
    """
    Filters out low-complexity k-mers (e.g. AAAA, ATAT).
    """
    return len(set(motif)) <= 2

def background_model(seq):
    """
    Computes background nucleotide frequencies.
    """
    length = len(seq)
    return {
        "A": seq.count("A") / length,
        "C": seq.count("C") / length,
        "G": seq.count("G") / length,
        "T": seq.count("T") / length
    }

# ---------- GC Skew ----------

def compute_gc_skew(seq):
    """
    Computes cumulative GC skew.
    """
    skew = []
    value = 0
    for base in seq:
        if base == "G":
            value += 1
        elif base == "C":
            value -= 1
        skew.append(value)
    return skew

def find_ori_by_gc_skew(seq):
    """
    Approximates ORI position as minimum GC skew.
    """
    skew = compute_gc_skew(seq)
    return skew.index(min(skew))

# ---------- PWM Construction ----------

def build_pwm(instances, k):
    """
    Builds a position weight matrix with pseudocounts.
    """
    pwm = []

    for i in range(k):
        column = {"A": 1, "C": 1, "G": 1, "T": 1}  # pseudocounts

        for motif in instances:
            column[motif[i]] += 1

        total = sum(column.values())
        for base in column:
            column[base] /= total

        pwm.append(column)

    return pwm

def score_with_pwm(seq, pwm, bg):
    """
    Scores a sequence using a PWM against background frequencies.
    """
    score = 0
    k = len(pwm)

    for i in range(len(seq) - k + 1):
        window = seq[i:i+k]
        llr = 0

        for j, base in enumerate(window):
            llr += math.log(pwm[j][base] / bg[base])

        score += llr

    return score

# ---------- MEME-style k-mer Search ----------

def meme_like_search(seq, k):
    """
    MEME-inspired motif discovery for a fixed k.
    """
    bg = background_model(seq)
    best_score = float("-inf")
    best_instances = None

    for i in range(len(seq) - k + 1):
        seed = seq[i:i+k]

        if is_low_complexity(seed):
            continue

        instances = [seed]

        for j in range(len(seq) - k + 1):
            window = seq[j:j+k]
            if not is_low_complexity(window):
                instances.append(window)

        pwm = build_pwm(instances, k)
        score = score_with_pwm(seq, pwm, bg)

        if score > best_score:
            best_score = score
            best_instances = instances

    return best_instances, best_score

def meme_k_selection(seq, k_min=6, k_max=15):
    """
    Selects the best k based on PWM log-likelihood score.
    """
    best_k = None
    best_score = float("-inf")
    best_motifs = None

    for k in range(k_min, k_max + 1):
        motifs, score = meme_like_search(seq, k)

        if motifs and score > best_score:
            best_k = k
            best_score = score
            best_motifs = motifs

    return best_k, best_score, best_motifs

# ---------- Final ORI Finder ----------

def find_ori_meme_style(fasta_file, region_size=500):
    """
    Final ORI finder using GC skew + MEME-style motif discovery.
    """
    record = SeqIO.read(fasta_file, "fasta")
    seq = str(record.seq).upper()

    # ---- GC skew localization ----
    ori_center = find_ori_by_gc_skew(seq)

    start = max(0, ori_center - region_size // 2)
    end = min(len(seq), ori_center + region_size // 2)
    region = seq[start:end]

    # ---- Motif discovery ----
    k, score, motifs = meme_k_selection(region)

    return {
        "ori_position": ori_center,
        "ori_region": (start, end),
        "ori_sequence": region,      # REQUIRED by plasmid_builder
        "best_k": k
    }
