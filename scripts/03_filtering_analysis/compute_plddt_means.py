import argparse
import numpy as np
import pandas as pd
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("--mode", choices=["binder", "nonbinder"], required=True)
parser.add_argument("--base", type=str, required=True, help="Path to outputs folder e.g. /path/to/outputs")
parser.add_argument("--output", type=str, required=False, help="Directory to save output CSV (defaults to base/mode/)")
args = parser.parse_args()

base = Path(args.base) / args.mode
output_dir = Path(args.output) if args.output else base
output_csv = output_dir / f"plddt_means_{args.mode}.csv"

all_rows = []

for chunk in sorted(base.glob("chunk_*")):
    print(f"Processing {chunk.name}...")

    anchor_df = pd.read_csv(chunk / "Multiple_Anchors_input.tsv", sep="\t")
    alphafold_df = pd.read_csv(chunk / "alphafold" / "alphafold_input_file.tsv", sep="\t")

    # Build lookups
    anchor_lookup = {}
    for _, row in anchor_df.iterrows():
        anchors = [int(a) - 1 for a in str(row["anchors"]).split(";")]
        anchor_lookup[row["id"]] = anchors

    length_lookup = {}
    for _, row in alphafold_df.iterrows():
        mhc_seq, pep_seq = row["target_chainseq"].split("/")
        id_name = row["targetid"].split("/")[0]
        length_lookup[id_name] = (len(mhc_seq), len(pep_seq))

    # Loop over ID folders
    for id_dir in sorted((chunk / "alphafold").iterdir()):
        if not id_dir.is_dir():
            continue

        id_name = id_dir.name

        if id_name not in length_lookup:
            print(f"  Skipping {id_name} — not in alphafold input tsv")
            continue
        if id_name not in anchor_lookup:
            print(f"  Skipping {id_name} — not in Multiple_Anchors_input.tsv")
            continue

        mhc_len, pep_len = length_lookup[id_name]
        anchors = anchor_lookup[id_name]

        plddt_file = id_dir / f"{id_name}_model_1_model_2_ptm_plddt.npy"
        if not plddt_file.exists():
            print(f"  Skipping {id_name} — plddt file not found")
            continue

        plddt = np.load(plddt_file)

        pep_plddt = plddt[mhc_len : mhc_len + pep_len]
        pep_mean = pep_plddt.mean()

        anchor_indices = [mhc_len + a for a in anchors]
        anchor_mean = plddt[anchor_indices].mean()

        parts = id_name.split("_")
        allele = parts[0]
        peptide = parts[1]
        struct_idx = int(parts[2])

        all_rows.append({
            "id": id_name,
            "allele": allele,
            "peptide": peptide,
            "struct_idx": struct_idx,
            "chunk": chunk.name,
            "pep_mean_plddt": round(pep_mean, 3),
            "anchor_mean_plddt": round(anchor_mean, 3),
        })

df_out = pd.DataFrame(all_rows)
df_out.to_csv(output_csv, index=False)
print(f"\nDone! {len(df_out)} structures written to {output_csv}")