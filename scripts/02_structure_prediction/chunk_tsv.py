import pandas as pd
from pathlib import Path

CHUNK_SIZE = 4000

files = {
    "/projects/scc/MPG/MGMN/scc_mgmn_soeding/dir.project/hasmig/data/pmgen_input/pmgen_input_binder.tsv": "/projects/scc/MPG/MGMN/scc_mgmn_soeding/dir.project/hasmig/data/pmgen_input_chunks/binder_chunks",
    "/projects/scc/MPG/MGMN/scc_mgmn_soeding/dir.project/hasmig/data/pmgen_input/pmgen_input_non_binder.tsv": "/projects/scc/MPG/MGMN/scc_mgmn_soeding/dir.project/hasmig/data/pmgen_input_chunks/non_binder_chunks",
}

for input_file, output_dir in files.items():
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(input_file, sep="\t")
    
    total_chunks = (len(df) + CHUNK_SIZE - 1) // CHUNK_SIZE  # ceiling division
    print(f"{input_file}: {len(df)} rows → {total_chunks} chunks")
    
    for i, start in enumerate(range(0, len(df), CHUNK_SIZE)):
        chunk = df.iloc[start:start + CHUNK_SIZE]
        out_path = f"{output_dir}/chunk_{i:03d}.tsv"
        chunk.to_csv(out_path, sep="\t", index=False)
        print(f"  Saved {out_path} ({len(chunk)} rows)")