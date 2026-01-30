from Bio import SeqIO
import os

from ori_finder import find_ori_meme_style
from .restriction_sites import load_restriction_sites

# ---------- CONSTANTS ----------

MARKER_DIR = "../Antibiotic_Marker/"
RESTRICTION_SITES = load_restriction_sites()

# ---------- DESIGN FILE PARSER ----------

def parse_design_file(design_file):
    mcs_sites = []
    genes = []

    with open(design_file) as f:
        for line in f:
            if not line.strip():
                continue

            part, value = line.strip().split(",")
            part = part.lower()
            value = value.strip()

            if "site" in part:
                mcs_sites.append(value)
            else:
                genes.append(value)

    return mcs_sites, genes

# ---------- MARKER LOADER ----------

def load_marker_sequence(marker_name):
    path = os.path.join(MARKER_DIR, f"{marker_name}.fa")

    if not os.path.exists(path):
        print(f"[WARNING] Marker not found: {marker_name}")
        return ""

    record = SeqIO.read(path, "fasta")
    return str(record.seq).upper()

# ---------- RESTRICTION SITE CLEANUP ----------

def remove_restriction_sites(seq):
    for site in RESTRICTION_SITES.values():
        seq = seq.replace(site, "")
    return seq

# ---------- PLASMID ASSEMBLY ----------

def assemble_plasmid(input_fasta, design_file):

    # ---- ORI DISCOVERY ----
    ori_info = find_ori_meme_style(input_fasta)

    print("\n[ORI DISCOVERY]")
    print(f"ORI Position : {ori_info['ori_position']}")
    print(f"ORI Region   : {ori_info['ori_region']}")
    print(f"Best k       : {ori_info['best_k']}")

    # ---- EXTRACT ORI SEQUENCE ----
    start, end = ori_info["ori_region"]
    record = SeqIO.read(input_fasta, "fasta")
    plasmid_seq = str(record.seq).upper()[start:end]

    # ---- DESIGN FILE ----
    mcs_sites, genes = parse_design_file(design_file)

    # ---- ADD GENES ----
    for gene in genes:
        plasmid_seq += load_marker_sequence(gene)

    # ---- REMOVE INTERNAL RESTRICTION SITES ----
    plasmid_seq = remove_restriction_sites(plasmid_seq)

    # ---- ADD MCS (LAST) ----
    for enzyme in mcs_sites:
        if enzyme not in RESTRICTION_SITES:
            print(f"[WARNING] Unknown restriction enzyme: {enzyme}")
            continue
        plasmid_seq += RESTRICTION_SITES[enzyme]

    return plasmid_seq
