import pandas as pd

def clean_key(allele_key: str) -> str:
    if allele_key is None:
        return "None"
    mapping = str.maketrans({"*": "", ":": "", " ": "", "/": "_", "-": ""})
    return allele_key.translate(mapping).upper()

mhc_db = pd.read_csv("/projects/scc/MPG/MGMN/scc_mgmn_soeding/dir.project/hasmig/data/raw/mhc1_encodings.csv")

# clean the key column to match your parquet format
mhc_db["key_clean"] = mhc_db["key"].apply(clean_key)


df = pd.read_csv("/user/hasmig.aintablian01/u26864/.project/dir.project/hasmig/data/pmgen_outputs/binder/plddt_means_binder_best.csv")

print(f"Total entries: {len(df)}")

df = df[df["pep_mean_plddt"] < 80]
# df = df[df["anchor_mean_plddt"] < 80]
df = df[df["peptide"].str.len() > 8]

print(f"Filtered entries: {len(df)}")

# merge
merged = df.merge(mhc_db[["key_clean", "mhc_sequence"]],
                  left_on="allele", right_on="key_clean", how="left")

print(f"Matched: {merged['mhc_sequence'].notna().sum()} / {len(merged)}")

tsv = pd.DataFrame({
    "peptide":  merged["peptide"],
    "mhc_seq":  merged["mhc_sequence"],
    "mhc_type": 1,
    "anchors":  "",
    "id":       merged["allele"] + "_" + merged["peptide"],
})
tsv.to_csv("pmgen_input_binder.tsv", sep="\t", index=False)

