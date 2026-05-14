"""
train.py
========
Training script for  GNN baseline.

Trains the model to predict pep_mean_plddt and
anchor_mean_plddt from pMHC structure graphs.

Supports both HLA and anchor split modes from
filter_map_train_prep.py.

Usage (single fold):
    python train.py \\
        --splits_dir /cbscratch/.../splits/hla \\
        --split_mode hla \\
        --fold 1 \\
        --output_dir /cbscratch/.../results/baseline_gnn

Usage (SLURM array, all 5 folds):
    #SBATCH --array=1-5
    python train.py \\
        --splits_dir /cbscratch/.../splits/hla \\
        --split_mode hla \\
        --fold $SLURM_ARRAY_TASK_ID \\
        --output_dir /cbscratch/.../results/baseline_gnn

Output (per fold):
    output_dir/
        fold_1/
            best_model.pt       ← best checkpoint by val loss
            metrics.json        ← train/val loss per epoch + test metrics
            config.json         ← all hyperparameters used
"""

import json
import logging
import argparse
from pathlib import Path

import torch
import torch.nn as nn
from torch_geometric.loader import DataLoader
from scipy.stats import spearmanr

from dataset import PMHCDataset
from model import PMTopoGNN

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s",
)
log = logging.getLogger(__name__)


# ══════════════════════════════════════════════════════════════════
# UTILITIES
# ══════════════════════════════════════════════════════════════════

def collate_skip_none(batch):
    """
    Custom collate function for DataLoader.
    Skips None items returned by dataset.get() when PDB parsing fails.
    """
    batch = [b for b in batch if b is not None]
    if not batch:
        return None
    from torch_geometric.data import Batch
    return Batch.from_data_list(batch)


def run_epoch(model, loader, optimizer, device, train=True):
    """
    Run one epoch of training or evaluation.

    Args:
        model     : GNN
        loader    : DataLoader
        optimizer : torch optimizer (None if train=False)
        device    : torch device
        train     : if True, update weights; if False, eval only

    Returns:
        mean_loss       : float, mean MSE loss over all batches
        all_preds       : list of predicted values (for Spearman)
        all_labels      : list of true values (for Spearman)
    """
    model.train(train)
    criterion = nn.MSELoss()

    total_loss = 0.0
    n_batches  = 0
    all_preds  = []
    all_labels = []

    for batch in loader:
        if batch is None:
            continue

        batch = batch.to(device)

        # Extract peptide mask from node features
        # x[:, 20] is the is_peptide flag set in dataset.py
        is_peptide_mask = batch.x[:, 20].bool()

        # Forward pass
        pred = model(
            x               = batch.x,
            edge_index      = batch.edge_index,
            edge_attr       = batch.edge_attr,
            batch           = batch.batch,
            is_peptide_mask = is_peptide_mask,
        )   # [batch_size, 2]

        # Labels: [batch_size, 2]
        # y was stored as [1, 2] per graph; DataLoader stacks them
        labels = batch.y.view(-1, 2)

        loss = criterion(pred, labels)

        if train:
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()

        total_loss += loss.item()
        n_batches  += 1

        # Collect for Spearman (detach from graph, move to cpu)
        all_preds.append(pred.detach().cpu())
        all_labels.append(labels.detach().cpu())

    mean_loss  = total_loss / max(n_batches, 1)
    all_preds  = torch.cat(all_preds,  dim=0)   # [N, 2]
    all_labels = torch.cat(all_labels, dim=0)   # [N, 2]

    return mean_loss, all_preds, all_labels


def compute_spearman(preds, labels):
    """
    Compute Spearman correlation for pep_plddt and anchor_plddt separately.

    Args:
        preds  : torch.Tensor [N, 2]
        labels : torch.Tensor [N, 2]

    Returns:
        dict with 'pep_spearman' and 'anchor_spearman'
    """
    pep_r, _    = spearmanr(preds[:, 0].numpy(), labels[:, 0].numpy())
    anchor_r, _ = spearmanr(preds[:, 1].numpy(), labels[:, 1].numpy())
    return {
        "pep_spearman"   : float(pep_r),
        "anchor_spearman": float(anchor_r),
    }


# ══════════════════════════════════════════════════════════════════
# MAIN TRAINING FUNCTION
# ══════════════════════════════════════════════════════════════════

def train_fold(args):
    """Train and evaluate for one fold."""

    # ── Paths ─────────────────────────────────────────────────
    splits_dir = Path(args.splits_dir)
    fold_out   = Path(args.output_dir) / f"fold_{args.fold}"
    fold_out.mkdir(parents=True, exist_ok=True)

    if args.split_mode == "hla":
        train_path = splits_dir / f"fold_{args.fold}" / "train.parquet"
        val_path   = splits_dir / f"fold_{args.fold}" / "val.parquet"
        test_path  = splits_dir / "test.parquet"
    elif args.split_mode == "anchor":
        train_path = splits_dir / f"fold_{args.fold}" / "train.parquet"
        val_path   = splits_dir / f"fold_{args.fold}" / "val.parquet"
        test_path  = splits_dir / f"fold_{args.fold}" / "test.parquet"
    else:
        raise ValueError(f"Unknown split_mode: {args.split_mode}")

    log.info("=" * 60)
    log.info(f"Fold {args.fold} | split_mode={args.split_mode}")
    log.info(f"  train : {train_path}")
    log.info(f"  val   : {val_path}")
    log.info(f"  test  : {test_path}")
    log.info(f"  output: {fold_out}")
    log.info("=" * 60)

    # ── Save config ────────────────────────────────────────────
    config = vars(args)
    with open(fold_out / "config.json", "w") as f:
        json.dump(config, f, indent=2)

    # ── Device ────────────────────────────────────────────────
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    log.info(f"Using device: {device}")

    # ── Datasets ──────────────────────────────────────────────
    log.info("Loading datasets...")
    train_ds = PMHCDataset(train_path)
    val_ds   = PMHCDataset(val_path)
    test_ds  = PMHCDataset(test_path)

    log.info(f"  train: {len(train_ds):,} structures")
    log.info(f"  val  : {len(val_ds):,} structures")
    log.info(f"  test : {len(test_ds):,} structures")

    # ── DataLoaders ───────────────────────────────────────────
    train_loader = DataLoader(
        train_ds,
        batch_size  = args.batch_size,
        shuffle     = True,
        num_workers = args.num_workers,
        collate_fn  = collate_skip_none,
    )
    val_loader = DataLoader(
        val_ds,
        batch_size  = args.batch_size,
        shuffle     = False,
        num_workers = args.num_workers,
        collate_fn  = collate_skip_none,
    )
    test_loader = DataLoader(
        test_ds,
        batch_size  = args.batch_size,
        shuffle     = False,
        num_workers = args.num_workers,
        collate_fn  = collate_skip_none,
    )

    # ── Model ─────────────────────────────────────────────────
    model = PMTopoGNN(
        node_in_dim = 22,
        edge_in_dim = 17,
        hidden_dim  = args.hidden_dim,
        num_layers  = args.num_layers,
        heads       = args.heads,
        dropout     = args.dropout,
    ).to(device)

    n_params = sum(p.numel() for p in model.parameters())
    log.info(f"Model parameters: {n_params:,}")

    # ── Optimiser + scheduler ─────────────────────────────────
    optimizer = torch.optim.Adam(model.parameters(), lr=args.lr)

    # Reduce learning rate by half if val loss plateaus for 5 epochs
    scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
        optimizer, mode="min", factor=0.5, patience=5, verbose=True
    )

    # ── Training loop ─────────────────────────────────────────
    best_val_loss  = float("inf")
    patience_count = 0
    history        = []

    for epoch in range(1, args.epochs + 1):

        # Train
        train_loss, _, _ = run_epoch(
            model, train_loader, optimizer, device, train=True
        )

        # Validate
        with torch.no_grad():
            val_loss, val_preds, val_labels = run_epoch(
                model, val_loader, None, device, train=False
            )

        val_spearman = compute_spearman(val_preds, val_labels)
        scheduler.step(val_loss)

        log.info(
            f"Epoch {epoch:3d} | "
            f"train_loss={train_loss:.4f} | "
            f"val_loss={val_loss:.4f} | "
            f"pep_r={val_spearman['pep_spearman']:.3f} | "
            f"anchor_r={val_spearman['anchor_spearman']:.3f}"
        )

        # Record history
        history.append({
            "epoch"          : epoch,
            "train_loss"     : train_loss,
            "val_loss"       : val_loss,
            "pep_spearman"   : val_spearman["pep_spearman"],
            "anchor_spearman": val_spearman["anchor_spearman"],
        })

        # Save best model
        if val_loss < best_val_loss:
            best_val_loss = val_loss
            patience_count = 0
            torch.save(model.state_dict(), fold_out / "best_model.pt")
            log.info(f"  → new best val_loss={best_val_loss:.4f}, model saved")
        else:
            patience_count += 1
            if patience_count >= args.patience:
                log.info(
                    f"Early stopping: val_loss has not improved "
                    f"for {args.patience} epochs"
                )
                break

    # ── Test evaluation ───────────────────────────────────────
    log.info("Loading best model for test evaluation...")
    model.load_state_dict(torch.load(fold_out / "best_model.pt",
                                     map_location=device,
                                     weights_only=True))

    with torch.no_grad():
        test_loss, test_preds, test_labels = run_epoch(
            model, test_loader, None, device, train=False
        )

    test_spearman = compute_spearman(test_preds, test_labels)

    log.info("=" * 60)
    log.info("TEST RESULTS")
    log.info(f"  test_loss       : {test_loss:.4f}")
    log.info(f"  pep_spearman    : {test_spearman['pep_spearman']:.4f}")
    log.info(f"  anchor_spearman : {test_spearman['anchor_spearman']:.4f}")
    log.info("=" * 60)

    # ── Save metrics ──────────────────────────────────────────
    metrics = {
        "fold"           : args.fold,
        "best_val_loss"  : best_val_loss,
        "test_loss"      : test_loss,
        "test_spearman"  : test_spearman,
        "history"        : history,
    }
    with open(fold_out / "metrics.json", "w") as f:
        json.dump(metrics, f, indent=2)

    log.info(f"Results saved to {fold_out}")
    return metrics


# ══════════════════════════════════════════════════════════════════
# ARGUMENT PARSING
# ══════════════════════════════════════════════════════════════════

def parse_args():
    parser = argparse.ArgumentParser(
        description="Train PMTopo GNN baseline for pLDDT prediction",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # ── Data ──────────────────────────────────────────────────
    parser.add_argument(
        "--splits_dir", required=True,
        help="Path to splits directory (splits/hla or splits/anchor)"
    )
    parser.add_argument(
        "--split_mode", required=True, choices=["hla", "anchor"],
        help="hla: global test set + k-fold CV | anchor: fold-level test"
    )
    parser.add_argument(
        "--fold", type=int, required=True,
        help="Which fold to train (1-5)"
    )
    parser.add_argument(
        "--output_dir", required=True,
        help="Directory to save checkpoints and metrics"
    )

    # ── Model ─────────────────────────────────────────────────
    parser.add_argument("--hidden_dim", type=int,   default=64)
    parser.add_argument("--num_layers", type=int,   default=3)
    parser.add_argument("--heads",      type=int,   default=4)
    parser.add_argument("--dropout",    type=float, default=0.1)

    # ── Training ──────────────────────────────────────────────
    parser.add_argument("--lr",          type=float, default=1e-3)
    parser.add_argument("--epochs",      type=int,   default=100)
    parser.add_argument("--batch_size",  type=int,   default=32)
    parser.add_argument("--patience",    type=int,   default=10,
                        help="Early stopping patience (epochs)")
    parser.add_argument("--num_workers", type=int,   default=2,
                        help="DataLoader worker processes")

    return parser.parse_args()


# ══════════════════════════════════════════════════════════════════
# ENTRY POINT
# ══════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    args = parse_args()
    train_fold(args)