"""
lr_model.py
===========
Logistic Regression baseline with attention pooling for pMHC binding prediction.

Architecture per feature (pLDDT and PAE):
    1. Initial embedding: [value, mhc_flag, pep_flag] -> (L, 3)
    2. Linear projection: (L, 3) -> (L, d)
    3. Sinusoidal positional encoding added: (L, d)
    4. RoPE applied to pairs: (L, d)
    5. Attention pooling: (L, d) -> (d,)

Both pLDDT and PAE pooled vectors are concatenated -> (2d,)
Then passed through logistic regression -> probability.

Ablation flags:
    --layer_norm          : add LayerNorm after attention pooling
    --regularization      : none (default), l1, or l2
    --reg_lambda          : regularization strength (default: 1e-4)

Folder naming convention:
    thr_{threshold}_baseline
    thr_{threshold}_layer_norm
    thr_{threshold}_l2_{reg_lambda}
    thr_{threshold}_layer_norm_l2_{reg_lambda}
    etc.

Usage:
    python lr_model.py \
        --splits_root /path/to/combined_dendrogram_splits/ \
        --output_dir  /path/to/output/ \
        --thresholds 50 60 70 80 \
        --folds 5 \
        --layer_norm \
        --regularization l2 \
        --reg_lambda 1e-4
"""

import argparse
import math
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
from torch.utils.data import Dataset, DataLoader
from sklearn.metrics import (roc_auc_score, average_precision_score,
                             f1_score, precision_score, recall_score, roc_curve)
from pathlib import Path

# ── constants ──────────────────────────────────────────────────────────────────
D       = 8    # embedding dimension
MAX_LEN = 200  # max pMHC complex length (verified from data: max=199)
BATCH   = 64
EPOCHS  = 50
LR      = 1e-3


# ── positional encodings ───────────────────────────────────────────────────────

def sinusoidal_pe(seq_len: int, d: int) -> torch.Tensor:
    pe = torch.zeros(seq_len, d)
    position = torch.arange(0, seq_len, dtype=torch.float).unsqueeze(1)
    div_term = torch.exp(
        torch.arange(0, d, 2, dtype=torch.float) * (-math.log(10000.0) / d)
    )
    pe[:, 0::2] = torch.sin(position * div_term)
    pe[:, 1::2] = torch.cos(position * div_term)
    return pe


def apply_rope(x: torch.Tensor) -> torch.Tensor:
    batch, seq_len, d = x.shape
    position = torch.arange(0, seq_len, dtype=torch.float, device=x.device)
    div_term = torch.exp(
        torch.arange(0, d, 2, dtype=torch.float, device=x.device) * (-math.log(10000.0) / d)
    )
    theta     = position.unsqueeze(1) * div_term.unsqueeze(0)
    cos_theta = torch.cos(theta)
    sin_theta = torch.sin(theta)

    x_even = x[:, :, 0::2]
    x_odd  = x[:, :, 1::2]

    x_out = torch.zeros_like(x)
    x_out[:, :, 0::2] = x_even * cos_theta - x_odd * sin_theta
    x_out[:, :, 1::2] = x_odd  * cos_theta + x_even * sin_theta
    return x_out


# ── model ──────────────────────────────────────────────────────────────────────

class FeatureEncoder(nn.Module):
    def __init__(self, d: int = D, max_len: int = MAX_LEN, layer_norm: bool = False):
        super().__init__()
        self.d          = d
        self.max_len    = max_len
        self.proj       = nn.Linear(3, d)
        self.attn_w     = nn.Linear(d, 1, bias=False)
        self.layer_norm = nn.LayerNorm(d) if layer_norm else None
        pe = sinusoidal_pe(max_len, d)
        self.register_buffer("pe", pe)

    def forward(self, x: torch.Tensor, mask: torch.Tensor) -> torch.Tensor:
        seq_len = x.shape[1]
        h = self.proj(x)
        h = h + self.pe[:seq_len].unsqueeze(0)
        h = apply_rope(h)
        scores  = self.attn_w(h).squeeze(-1)
        scores  = scores.masked_fill(mask == 0, float("-inf"))
        weights = torch.softmax(scores, dim=-1)
        pooled  = (weights.unsqueeze(-1) * h).sum(dim=1)
        if self.layer_norm is not None:
            pooled = self.layer_norm(pooled)
        return pooled


class PMHCLogisticRegression(nn.Module):
    def __init__(self, d: int = D, max_len: int = MAX_LEN, layer_norm: bool = False):
        super().__init__()
        self.plddt_encoder = FeatureEncoder(d=d, max_len=max_len, layer_norm=layer_norm)
        self.pae_encoder   = FeatureEncoder(d=d, max_len=max_len, layer_norm=layer_norm)
        self.classifier    = nn.Linear(2 * d, 1)

    def forward(self, plddt_x, pae_x, mask):
        h_plddt = self.plddt_encoder(plddt_x, mask)
        h_pae   = self.pae_encoder(pae_x, mask)
        h       = torch.cat([h_plddt, h_pae], dim=-1)
        return torch.sigmoid(self.classifier(h)).squeeze(-1)


# ── dataset ────────────────────────────────────────────────────────────────────

class PMHCDataset(Dataset):
    def __init__(self, parquet_path: str, max_len: int = MAX_LEN):
        self.max_len = max_len
        df = pd.read_parquet(parquet_path)

        df["pep_len"] = df["long_mer"].str.len()
        df["mhc_len"] = df["mhc_sequence"].str.len()

        before = len(df)
        df = df.dropna(subset=["plddt_vector", "pae_vector"]).reset_index(drop=True)
        after = len(df)
        if before != after:
            print(f"  Dropped {before - after:,} rows with missing vectors")
        print(f"  Loaded {len(df):,} rows from {parquet_path}")

        self.plddt_vectors = df["plddt_vector"].tolist()
        self.pae_vectors   = df["pae_vector"].tolist()
        self.mhc_lens      = df["mhc_len"].tolist()
        self.pep_lens      = df["pep_len"].tolist()
        self.labels        = df["label"].tolist()

    def __len__(self):
        return len(self.labels)

    def __getitem__(self, idx):
        plddt   = np.array(self.plddt_vectors[idx], dtype=np.float32)
        pae     = np.array(self.pae_vectors[idx],   dtype=np.float32)
        mhc_len = int(self.mhc_lens[idx])
        pep_len = int(self.pep_lens[idx])
        seq_len = mhc_len + pep_len

        plddt = plddt[:seq_len]
        pae   = pae[:seq_len]

        mhc_flag = np.array([1.0] * mhc_len + [0.0] * pep_len, dtype=np.float32)
        pep_flag = np.array([0.0] * mhc_len + [1.0] * pep_len, dtype=np.float32)

        plddt_x = np.stack([plddt, mhc_flag, pep_flag], axis=-1)
        pae_x   = np.stack([pae,   mhc_flag, pep_flag], axis=-1)

        plddt_x[:, 0] = plddt_x[:, 0] / 100.0
        pae_x[:, 0]   = -pae_x[:, 0] / 31.75

        pad_len = self.max_len - seq_len
        plddt_x = np.pad(plddt_x, ((0, pad_len), (0, 0)), mode="constant")
        pae_x   = np.pad(pae_x,   ((0, pad_len), (0, 0)), mode="constant")
        mask    = np.array([1.0] * seq_len + [0.0] * pad_len, dtype=np.float32)

        return (
            torch.tensor(plddt_x, dtype=torch.float32),
            torch.tensor(pae_x,   dtype=torch.float32),
            torch.tensor(mask,    dtype=torch.float32),
            torch.tensor(float(self.labels[idx]), dtype=torch.float32),
        )


# ── evaluation ─────────────────────────────────────────────────────────────────

def evaluate(model, loader, device):
    model.eval()
    all_preds  = []
    all_labels = []
    total_loss = 0.0
    criterion  = nn.BCELoss()

    with torch.no_grad():
        for plddt_x, pae_x, mask, labels in loader:
            plddt_x = plddt_x.to(device)
            pae_x   = pae_x.to(device)
            mask    = mask.to(device)
            labels  = labels.to(device)
            preds   = model(plddt_x, pae_x, mask)
            loss    = criterion(preds, labels)
            total_loss += loss.item()
            all_preds.extend(preds.cpu().numpy())
            all_labels.extend(labels.cpu().numpy())

    all_preds    = np.array(all_preds)
    all_labels   = np.array(all_labels)
    preds_binary = (all_preds >= 0.5).astype(int)
    has_both     = len(np.unique(all_labels)) > 1

    auc  = float(roc_auc_score(all_labels, all_preds))           if has_both else float("nan")
    ap   = float(average_precision_score(all_labels, all_preds)) if has_both else float("nan")
    f1   = float(f1_score(all_labels,   preds_binary, zero_division=0))
    prec = float(precision_score(all_labels, preds_binary, zero_division=0))
    rec  = float(recall_score(all_labels,  preds_binary, zero_division=0))

    if has_both:
        fpr_arr, tpr_arr, _ = roc_curve(all_labels, all_preds)
        mask_fpr = fpr_arr <= 0.10
        tpr_at_fpr10 = float(tpr_arr[mask_fpr][-1]) if mask_fpr.any() else 0.0
    else:
        tpr_at_fpr10 = float("nan")

    return {
        "loss":         total_loss / len(loader),
        "auc":          auc,
        "ap":           ap,
        "f1":           f1,
        "precision":    prec,
        "recall":       rec,
        "tpr_at_fpr10": tpr_at_fpr10,
    }


# ── helpers ────────────────────────────────────────────────────────────────────

def make_run_tag(thr: str, args) -> str:
    """Build folder name tag from threshold and ablation config."""
    tag = f"thr_{thr}"
    if args.layer_norm:
        tag += "_layer_norm"
    if args.regularization != "none":
        tag += f"_{args.regularization}_{args.reg_lambda}"
    if not args.layer_norm and args.regularization == "none":
        tag += "_baseline"
    return tag


def make_model(args, device):
    return PMHCLogisticRegression(
        d=D, max_len=MAX_LEN, layer_norm=args.layer_norm
    ).to(device)


def make_optimizer(model, args):
    if args.regularization == "l2":
        return torch.optim.Adam(model.parameters(), lr=LR, weight_decay=args.reg_lambda)
    return torch.optim.Adam(model.parameters(), lr=LR)


def compute_loss(preds, labels, model, criterion, args):
    loss = criterion(preds, labels)
    if args.regularization == "l1":
        l1 = sum(p.abs().sum() for p in model.parameters())
        loss = loss + args.reg_lambda * l1
    return loss


def metrics_row(thr, fold, eval_set, metrics, args):
    return {
        "threshold":      thr,
        "fold":           fold,
        "eval_set":       eval_set,
        "layer_norm":     args.layer_norm,
        "regularization": args.regularization,
        "reg_lambda":     args.reg_lambda if args.regularization != "none" else None,
        "auc":            metrics["auc"],
        "ap":             metrics["ap"],
        "f1":             metrics["f1"],
        "precision":      metrics["precision"],
        "recall":         metrics["recall"],
        "tpr_at_fpr10":   metrics["tpr_at_fpr10"],
    }


# ── main ───────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--splits_root",    type=str, required=True)
    parser.add_argument("--output_dir",     type=str, required=True)
    parser.add_argument("--thresholds",     nargs="+", default=["50", "60", "70", "80"])
    parser.add_argument("--folds",          type=int, default=5)
    parser.add_argument("--test_only",      action="store_true",
                        help="Skip CV folds, only retrain final model and evaluate on test")
    parser.add_argument("--layer_norm",     action="store_true",
                        help="Add LayerNorm after attention pooling")
    parser.add_argument("--regularization", choices=["none", "l1", "l2"], default="none")
    parser.add_argument("--reg_lambda",     type=float, default=1e-4)
    args = parser.parse_args()

    splits_root = Path(args.splits_root)
    output_dir  = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"Using device: {device}")
    print(f"Config: layer_norm={args.layer_norm}, "
          f"regularization={args.regularization}, "
          f"reg_lambda={args.reg_lambda}")

    n_params = sum(p.numel() for p in PMHCLogisticRegression(
        d=D, max_len=MAX_LEN, layer_norm=args.layer_norm).parameters())
    print(f"Model parameters: {n_params:,}")

    criterion   = nn.BCELoss()
    all_metrics = []

    for thr in args.thresholds:
        thr_root  = splits_root / thr / "splits"
        test_path = thr_root / "test" / "test.parquet"
        run_tag   = make_run_tag(thr, args)
        run_dir   = output_dir / run_tag
        run_dir.mkdir(parents=True, exist_ok=True)

        # ── cross-validation folds ────────────────────────────────────────────
        if not args.test_only:
            for fold in range(args.folds):
                fold_dir   = thr_root / f"fold_{fold}"
                train_path = fold_dir / "train.parquet"
                val_path   = fold_dir / "val.parquet"

                if not train_path.exists() or not val_path.exists():
                    print(f"WARNING: missing files for thr={thr} fold={fold}, skipping")
                    continue

                print(f"\n{'='*60}")
                print(f"Run={run_tag}  Fold={fold}")
                print(f"{'='*60}")

                train_ds     = PMHCDataset(str(train_path))
                val_ds       = PMHCDataset(str(val_path))
                train_loader = DataLoader(train_ds, batch_size=BATCH, shuffle=True,  num_workers=4)
                val_loader   = DataLoader(val_ds,   batch_size=BATCH, shuffle=False, num_workers=4)

                model     = make_model(args, device)
                optimizer = make_optimizer(model, args)

                fold_output_dir = run_dir / f"fold_{fold}"
                fold_output_dir.mkdir(parents=True, exist_ok=True)

                best_val_auc     = 0.0
                best_val_metrics = None

                for epoch in range(1, EPOCHS + 1):
                    model.train()
                    for plddt_x, pae_x, mask, labels in train_loader:
                        plddt_x = plddt_x.to(device)
                        pae_x   = pae_x.to(device)
                        mask    = mask.to(device)
                        labels  = labels.to(device)
                        optimizer.zero_grad()
                        preds = model(plddt_x, pae_x, mask)
                        loss  = compute_loss(preds, labels, model, criterion, args)
                        loss.backward()
                        optimizer.step()

                    val_metrics = evaluate(model, val_loader, device)
                    print(f"  Epoch {epoch:3d}/{EPOCHS} | "
                          f"val loss={val_metrics['loss']:.4f} "
                          f"auc={val_metrics['auc']:.4f} "
                          f"ap={val_metrics['ap']:.4f} "
                          f"f1={val_metrics['f1']:.4f}")

                    if val_metrics["auc"] > best_val_auc:
                        best_val_auc     = val_metrics["auc"]
                        best_val_metrics = val_metrics.copy()
                        torch.save(model.state_dict(), fold_output_dir / "best_model.pt")
                        print(f"    -> Saved best model (val AUC={best_val_auc:.4f})")

                all_metrics.append(metrics_row(thr, fold, "val", best_val_metrics, args))

        # ── final model: retrain on all folds combined ────────────────────────
        if test_path.exists():
            print(f"\n{'='*60}")
            print(f"Run={run_tag} -- Retraining final model on all folds combined")
            print(f"{'='*60}")

            all_parquet_paths = []
            for fold in range(args.folds):
                fold_dir = thr_root / f"fold_{fold}"
                for part in ("train", "val"):
                    p = fold_dir / f"{part}.parquet"
                    if p.exists():
                        all_parquet_paths.append(str(p))

            combined_df = pd.concat(
                [pd.read_parquet(p) for p in all_parquet_paths],
                ignore_index=True
            ).drop_duplicates(subset=["allele", "long_mer", "struct_idx"])
            print(f"  Combined dataset: {len(combined_df):,} rows")

            combined_path = run_dir / "combined_train.parquet"
            combined_df.to_parquet(combined_path, index=False)

            combined_ds     = PMHCDataset(str(combined_path))
            combined_loader = DataLoader(combined_ds, batch_size=BATCH,
                                         shuffle=True, num_workers=4)

            final_model     = make_model(args, device)
            final_optimizer = make_optimizer(final_model, args)
            final_output_dir = run_dir / "final_model"
            final_output_dir.mkdir(parents=True, exist_ok=True)

            for epoch in range(1, EPOCHS + 1):
                final_model.train()
                for plddt_x, pae_x, mask, labels in combined_loader:
                    plddt_x = plddt_x.to(device)
                    pae_x   = pae_x.to(device)
                    mask    = mask.to(device)
                    labels  = labels.to(device)
                    final_optimizer.zero_grad()
                    preds = final_model(plddt_x, pae_x, mask)
                    loss  = compute_loss(preds, labels, final_model, criterion, args)
                    loss.backward()
                    final_optimizer.step()
                if epoch % 10 == 0:
                    print(f"  Final model epoch {epoch}/{EPOCHS}")

            torch.save(final_model.state_dict(), final_output_dir / "final_model.pt")

            test_ds     = PMHCDataset(str(test_path))
            test_loader = DataLoader(test_ds, batch_size=BATCH, shuffle=False, num_workers=4)
            test_metrics = evaluate(final_model, test_loader, device)

            print(f"  Test | auc={test_metrics['auc']:.4f} "
                  f"ap={test_metrics['ap']:.4f} "
                  f"f1={test_metrics['f1']:.4f} "
                  f"precision={test_metrics['precision']:.4f} "
                  f"recall={test_metrics['recall']:.4f} "
                  f"tpr@fpr10={test_metrics['tpr_at_fpr10']:.4f}")

            all_metrics.append(metrics_row(thr, "test", "test", test_metrics, args))

            combined_path.unlink()

    # ── save metrics ──────────────────────────────────────────────────────────
    metrics_df   = pd.DataFrame(all_metrics)
    metrics_path = output_dir / "all_metrics.csv"

    if metrics_path.exists():
        existing   = pd.read_csv(metrics_path)
        metrics_df = pd.concat([existing, metrics_df], ignore_index=True)
        metrics_df = metrics_df.drop_duplicates(
            subset=["threshold", "fold", "eval_set", "layer_norm", "regularization"]
        )

    metrics_df.to_csv(metrics_path, index=False)
    print(f"\nMetrics saved to {metrics_path}  ({len(metrics_df):,} rows total)")

    val_df = metrics_df[metrics_df["eval_set"] == "val"]
    if not val_df.empty:
        print("\n=== Summary (median val AUC per threshold / config) ===")
        summary = val_df.groupby(
            ["threshold", "layer_norm", "regularization"]
        )["auc"].median().round(4)
        print(summary.to_string())