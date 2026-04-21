import pandas as pd

def clean_key(allele_key: str) -> str:
    if allele_key is None:
        return "None"
    mapping = str.maketrans({"*": "", ":": "", " ": "", "/": "_", "-": ""})
    return allele_key.translate(mapping).upper()

df = pd.read_parquet("/user/hasmig.aintablian01/u26864/.project/dir.project/hasmig/binder_sampled.parquet")
mhc_db = pd.read_csv("/projects/scc/MPG/MGMN/scc_mgmn_soeding/dir.project/hasmig/data/raw/mhc1_encodings.csv")

# clean the key column to match your parquet format
mhc_db["key_clean"] = mhc_db["key"].apply(clean_key)

# check they now match
print(df["allele"].head(5))
print(mhc_db["key_clean"].head(5))

# merge
merged = df.merge(mhc_db[["key_clean", "mhc_sequence"]],
                  left_on="allele", right_on="key_clean", how="left")

print(f"Matched: {merged['mhc_sequence'].notna().sum()} / {len(merged)}")

tsv = pd.DataFrame({
    "peptide":  merged["long_mer"],
    "mhc_seq":  merged["mhc_sequence"],
    "mhc_type": 1,
    "anchors":  "",
    "id":       merged["allele"] + "_" + merged["long_mer"],
})

tsv.to_csv("pmgen_input_non_binder.tsv", sep="\t", index=False)