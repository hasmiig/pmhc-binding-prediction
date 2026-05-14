"""
dataset.py
==========
PyTorch Geometric dataset for PMTopo GNN baseline.

Reads split parquets from filter_map_train_prep.py.
For each row, parses the PDB file, builds a residue-level graph
with SE(3)-invariant features, and returns a PyG Data object.

Graph construction:
  - Nodes  : all residues in the pMHC complex (MHC + peptide)
  - Edges  : Cα-Cα distance cutoff (default 8Å)
  - Node features : amino acid one-hot (20-dim)
                    + is_peptide flag (1-dim)
                    + is_anchor flag (1-dim)
                    = 22-dim total
  - Edge features : Cα-Cα distance encoded as RBF (16-dim)
                    + backbone angle cosine (1-dim)
                    = 17-dim total

Labels (both from parquet, already computed in QC step):
  - pep_mean_plddt   : mean pLDDT over peptide residues
  - anchor_mean_plddt: mean pLDDT over anchor residues

Usage:
    from dataset import PMHCDataset
    from torch_geometric.loader import DataLoader

    train_ds = PMHCDataset("path/to/fold_1/train.parquet")
    loader   = DataLoader(train_ds, batch_size=32, shuffle=True)

    for batch in loader:
        # batch.x          : node features  [N_nodes, 22]
        # batch.edge_index : edges          [2, N_edges]
        # batch.edge_attr  : edge features  [N_edges, 17]
        # batch.y          : labels         [N_graphs, 2]
        #                    col 0 = pep_mean_plddt
        #                    col 1 = anchor_mean_plddt
        # batch.batch      : graph id per node [N_nodes]
        pass
"""

import logging
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from torch_geometric.data import Data, Dataset

log = logging.getLogger(__name__)

# ─────────────────────────────────────────────────────────────────
# CONSTANTS
# ─────────────────────────────────────────────────────────────────

# Standard amino acid ordering — same as prepare_pmhc_data.py
THREE_TO_ONE = {
    'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
    'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
    'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
    'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y',
}

AA_ORDER = list("ACDEFGHIKLMNPQRSTVWY")  # 20 standard AAs
AA_TO_IDX = {aa: i for i, aa in enumerate(AA_ORDER)}

# Default graph construction parameters
DEFAULT_CUTOFF_ANGSTROM = 8.0   # Cα-Cα distance cutoff for edges
DEFAULT_RBF_BINS        = 16    # number of RBF centres for distance encoding
DEFAULT_RBF_MIN         = 2.0   # min distance (Å) for RBF
DEFAULT_RBF_MAX         = 20.0  # max distance (Å) for RBF


# ══════════════════════════════════════════════════════════════════
# PDB PARSING  (adapted from prepare_pmhc_data.py)
# ══════════════════════════════════════════════════════════════════

def find_pdb_file(pdb_dir: str | Path) -> Path:
    """
    Find the .pdb file inside a PMGen output directory.
    PMGen names it: {dir_name}_model_1_model_2_ptm.pdb
    """
    pdb_dir = Path(pdb_dir)
    if not pdb_dir.exists():
        raise FileNotFoundError(f"PDB directory not found: {pdb_dir}")

    # Try the expected PMGen naming convention first
    expected = pdb_dir / f"{pdb_dir.name}_model_1_model_2_ptm.pdb"
    if expected.exists():
        return expected

    # Fall back to any .pdb file in the directory
    pdbs = list(pdb_dir.glob("*.pdb"))
    if not pdbs:
        raise FileNotFoundError(f"No .pdb file found in: {pdb_dir}")
    return pdbs[0]


def parse_pdb(pdb_path: Path) -> dict | None:
    """
    Parse a PDB file and return per-residue data.

    Returns a dict with keys:
        residues : list of dicts, each with:
            chain_id  : str
            res_num   : str
            aa_one    : str  (single-letter AA)
            ca_coords : np.ndarray shape (3,)  — Cα xyz
            n_coords  : np.ndarray shape (3,)  — N  xyz
            c_coords  : np.ndarray shape (3,)  — C  xyz

    Returns None if no valid residues found.
    """
    residue_order  = {}  # (chain, res_num, resname) -> insertion order
    atom_data      = {}  # (chain, res_num, resname) -> {atom_name: [x,y,z]}

    with open(pdb_path, "r") as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue

            atom_name = line[12:16].strip()
            resname   = line[17:20].strip()
            chain_id  = line[21]
            res_num   = line[22:27].strip()

            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
            except ValueError:
                continue

            if resname not in THREE_TO_ONE:
                continue
            if atom_name not in ("N", "CA", "C", "O"):
                continue

            key = (chain_id, res_num, resname)
            if key not in residue_order:
                residue_order[key] = len(residue_order)
            atom_data.setdefault(key, {})[atom_name] = np.array([x, y, z])

    if not atom_data:
        return None

    # Build ordered list of residues that have all backbone atoms
    residues = []
    for key in sorted(residue_order, key=lambda k: residue_order[k]):
        chain_id, res_num, resname = key
        atoms = atom_data[key]
        if all(a in atoms for a in ("N", "CA", "C")):
            residues.append({
                "chain_id" : chain_id,
                "res_num"  : res_num,
                "aa_one"   : THREE_TO_ONE[resname],
                "ca_coords": atoms["CA"],
                "n_coords" : atoms["N"],
                "c_coords" : atoms["C"],
            })

    return {"residues": residues} if residues else None


# ══════════════════════════════════════════════════════════════════
# FEATURE ENGINEERING (SE(3)-invariant)
# ══════════════════════════════════════════════════════════════════

def aa_one_hot(aa: str) -> np.ndarray:
    """20-dim one-hot for amino acid. Unknown AAs get all-zeros."""
    vec = np.zeros(20, dtype=np.float32)
    idx = AA_TO_IDX.get(aa, None)
    if idx is not None:
        vec[idx] = 1.0
    return vec


def rbf_encode(distances: np.ndarray,
               n_bins: int = DEFAULT_RBF_BINS,
               d_min: float = DEFAULT_RBF_MIN,
               d_max: float = DEFAULT_RBF_MAX) -> np.ndarray:
    """
    Encode scalar distances as radial basis functions.

    Places n_bins Gaussian centres evenly between d_min and d_max.
    Output shape: (len(distances), n_bins).

    This is SE(3)-invariant because distance is a scalar that
    doesn't change under rotation or translation.
    """
    centres = np.linspace(d_min, d_max, n_bins)          # (n_bins,)
    sigma   = (d_max - d_min) / n_bins                   # width of each Gaussian
    dists   = distances[:, None]                          # (N, 1)
    rbf     = np.exp(-((dists - centres) ** 2) / (2 * sigma ** 2))
    return rbf.astype(np.float32)                         # (N, n_bins)


def backbone_angle_cosine(r_i: dict, r_j: dict) -> float:
    """
    Compute the cosine of the angle between the backbone vectors
    of two residues. Uses N->CA->C as the backbone direction.

    This is SE(3)-invariant: it's a dot product of unit vectors,
    which is the same regardless of global orientation.
    """
    def backbone_vec(r):
        v = r["c_coords"] - r["n_coords"]
        norm = np.linalg.norm(v)
        return v / norm if norm > 1e-6 else v

    vi = backbone_vec(r_i)
    vj = backbone_vec(r_j)
    return float(np.clip(np.dot(vi, vj), -1.0, 1.0))


# ══════════════════════════════════════════════════════════════════
# GRAPH BUILDER
# ══════════════════════════════════════════════════════════════════

def build_graph(residues: list[dict],
                peptide_seq: str,
                mhc_seq: str,
                anchor_positions: list[int],
                cutoff: float = DEFAULT_CUTOFF_ANGSTROM,
                rbf_bins: int = DEFAULT_RBF_BINS) -> Data | None:
    """
    Build a PyTorch Geometric Data object from a list of residues.

    Node features (22-dim per node):
        [0:20]  amino acid one-hot
        [20]    is_peptide  (1 if peptide residue, 0 if MHC)
        [21]    is_anchor   (1 if anchor position, 0 otherwise)

    Edge features (17-dim per edge):
        [0:16]  RBF-encoded Cα-Cα distance
        [16]    backbone angle cosine

    anchor_positions: 1-based positions within the peptide
                      (as stored in your dataset).
    """
    n = len(residues)
    if n == 0:
        return None

    # ── Identify peptide residues ──────────────────────────────
    # In PMGen PDBs, MHC chain comes first, peptide chain last.
    # We use the known lengths to assign flags.
    mhc_len = len(mhc_seq)
    pep_len = len(peptide_seq)

    # Convert anchor positions (1-based peptide indices) to
    # global residue indices (0-based in the full residue list)
    anchor_global = set()
    for pos in anchor_positions:
        global_idx = mhc_len + (pos - 1)  # pos is 1-based
        if 0 <= global_idx < n:
            anchor_global.add(global_idx)

    # ── Node features ──────────────────────────────────────────
    node_feats = []
    for i, res in enumerate(residues):
        oh        = aa_one_hot(res["aa_one"])             # (20,)
        is_pep    = float(i >= mhc_len)                   # 1 scalar
        is_anchor = float(i in anchor_global)             # 1 scalar
        feat      = np.concatenate([oh, [is_pep, is_anchor]])  # (22,)
        node_feats.append(feat)

    x = torch.tensor(np.stack(node_feats), dtype=torch.float)  # (N, 22)

    # ── Cα coordinates for distance computation ────────────────
    ca_coords = np.stack([r["ca_coords"] for r in residues])    # (N, 3)

    # ── Edge construction: Cα-Cα distance cutoff ──────────────
    src_list, dst_list = [], []
    dist_list          = []
    angle_list         = []

    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            d = np.linalg.norm(ca_coords[i] - ca_coords[j])
            if d <= cutoff:
                src_list.append(i)
                dst_list.append(j)
                dist_list.append(d)
                angle_list.append(backbone_angle_cosine(residues[i], residues[j]))

    if not src_list:
        log.warning("No edges found — check cutoff or structure")
        return None

    # ── Edge features ──────────────────────────────────────────
    dists       = np.array(dist_list, dtype=np.float32)          # (E,)
    rbf         = rbf_encode(dists, n_bins=rbf_bins)             # (E, 16)
    angles      = np.array(angle_list, dtype=np.float32)[:, None]  # (E, 1)
    edge_attr   = np.concatenate([rbf, angles], axis=1)          # (E, 17)

    edge_index  = torch.tensor([src_list, dst_list], dtype=torch.long)   # (2, E)
    edge_attr_t = torch.tensor(edge_attr, dtype=torch.float)             # (E, 17)

    return Data(x=x, edge_index=edge_index, edge_attr=edge_attr_t)


# ══════════════════════════════════════════════════════════════════
# DATASET CLASS
# ══════════════════════════════════════════════════════════════════

class PMHCDataset(Dataset):
    """
    PyTorch Geometric Dataset for pMHC structures.

    Each item is one pMHC complex represented as a graph, with
    two regression targets:
        y[:, 0] = pep_mean_plddt
        y[:, 1] = anchor_mean_plddt

    Args:
        parquet_path : path to a split parquet from filter_map_train_prep.py
        cutoff       : Cα-Cα distance cutoff in Å (default 8.0)
        rbf_bins     : number of RBF bins for distance encoding (default 16)
        verbose      : if True, log every parse error

    The parquet is expected to have columns:
        pdb_path, long_mer, mhc_sequence,
        pep_mean_plddt, anchor_mean_plddt,
        (optionally) anchors — semicolon-separated 1-based positions
    """

    def __init__(self,
                 parquet_path: str | Path,
                 cutoff: float = DEFAULT_CUTOFF_ANGSTROM,
                 rbf_bins: int = DEFAULT_RBF_BINS,
                 verbose: bool = False):
        super().__init__()
        self.parquet_path = Path(parquet_path)
        self.cutoff       = cutoff
        self.rbf_bins     = rbf_bins
        self.verbose      = verbose

        log.info(f"Loading parquet: {self.parquet_path}")
        self._df = pd.read_parquet(self.parquet_path)
        log.info(f"  {len(self._df):,} rows loaded")

        # Pre-filter rows where pdb_path is missing
        before = len(self._df)
        self._df = self._df[self._df["pdb_path"].notna()].reset_index(drop=True)
        after = len(self._df)
        if before != after:
            log.warning(f"  Dropped {before - after} rows with missing pdb_path")

    def len(self) -> int:
        return len(self._df)

    def get(self, idx: int) -> Data | None:
        """
        Load one pMHC structure and return a PyG Data object.
        Returns None if parsing fails (DataLoader will skip None items).
        """
        row = self._df.iloc[idx]

        # ── Locate and parse PDB ───────────────────────────────
        try:
            pdb_file = find_pdb_file(row["pdb_path"])
            parsed   = parse_pdb(pdb_file)
        except (FileNotFoundError, Exception) as e:
            if self.verbose:
                log.warning(f"  [{idx}] PDB error: {e}")
            return None

        if parsed is None:
            if self.verbose:
                log.warning(f"  [{idx}] PDB parsed but no valid residues")
            return None

        residues = parsed["residues"]

        # ── Anchor positions ───────────────────────────────────
        # 'anchors' column is semicolon-separated 1-based positions
        # e.g. "2;9"  ->  [2, 9]
        raw_anchors = row.get("anchors", "")
        if pd.isna(raw_anchors) or raw_anchors == "":
            anchor_positions = []
        else:
            try:
                anchor_positions = [int(p) for p in str(raw_anchors).split(";")]
            except ValueError:
                anchor_positions = []

        # ── Build graph ────────────────────────────────────────
        data = build_graph(
            residues         = residues,
            peptide_seq      = str(row["long_mer"]),
            mhc_seq          = str(row["mhc_sequence"]),
            anchor_positions = anchor_positions,
            cutoff           = self.cutoff,
            rbf_bins         = self.rbf_bins,
        )

        if data is None:
            return None

        # ── Labels ─────────────────────────────────────────────
        pep_plddt    = float(row["pep_mean_plddt"])
        anchor_plddt = float(row["anchor_mean_plddt"])
        data.y = torch.tensor([[pep_plddt, anchor_plddt]], dtype=torch.float)

        # ── Metadata (useful for evaluation) ──────────────────
        data.allele  = str(row.get("allele", ""))
        data.peptide = str(row["long_mer"])
        data.pdb_path = str(row["pdb_path"])

        return data


# ══════════════════════════════════════════════════════════════════
# QUICK SANITY CHECK
# ══════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    import sys
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s | %(levelname)s | %(message)s")

    if len(sys.argv) < 2:
        print("Usage: python dataset.py path/to/split.parquet")
        sys.exit(1)

    ds = PMHCDataset(sys.argv[1], verbose=True)
    print(f"\nDataset size: {len(ds)}")

    # Try loading first 5 items
    ok, fail = 0, 0
    for i in range(min(5, len(ds))):
        d = ds.get(i)
        if d is not None:
            print(f"  [{i}] nodes={d.x.shape[0]}, "
                  f"edges={d.edge_index.shape[1]}, "
                  f"x={d.x.shape}, "
                  f"edge_attr={d.edge_attr.shape}, "
                  f"y={d.y}")
            ok += 1
        else:
            print(f"  [{i}] FAILED to parse")
            fail += 1

    print(f"\nResult: {ok} OK, {fail} failed out of {min(5, len(ds))} tried")