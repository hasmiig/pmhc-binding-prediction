"""
model.py
========
GNN baseline model for PMTopo pLDDT prediction.

Architecture:
    - 3 × GATConv layers (graph attention)
    - each layer updates node representations using
      attention-weighted messages from neighbors
    - after message passing, pool over peptide nodes only
    - linear head outputs two predictions:
        [pep_mean_plddt, anchor_mean_plddt]

Why GATConv:
    Not all neighbors contribute equally to a residue's pLDDT.
    Anchor residues contacting MHC pockets matter more than
    distant residues. Attention lets the model learn these
    weights automatically.

Input (from dataset.py):
    x          : [N_nodes, 22]   node features
    edge_index : [2, N_edges]    graph connectivity
    edge_attr  : [N_edges, 17]   edge features
    batch      : [N_nodes]       which graph each node belongs to
                                 (used by DataLoader for batching)

Output:
    pred : [batch_size, 2]   predicted pep and anchor mean pLDDT

Usage:
    from model import PMTopoGNN

    model = PMTopoGNN(
        node_in_dim  = 22,   # must match dataset.py
        edge_in_dim  = 17,   # must match dataset.py
        hidden_dim   = 64,
        num_layers   = 3,
        heads        = 4,
        dropout      = 0.1,
    )

    pred = model(batch.x, batch.edge_index, batch.edge_attr,
                 batch.batch, batch.is_peptide_node)
    # pred shape: [batch_size, 2]
"""

import torch
import torch.nn as nn
from torch_geometric.nn import GATConv, global_mean_pool


class PMTopoGNN(nn.Module):
    """
    Graph Attention Network for pMHC pLDDT prediction.

    Args:
        node_in_dim  : input node feature dimension (22 from dataset.py)
        edge_in_dim  : input edge feature dimension (17 from dataset.py)
        hidden_dim   : hidden dimension for GNN layers
        num_layers   : number of GATConv layers (default 3)
        heads        : number of attention heads per GAT layer (default 4)
        dropout      : dropout probability (default 0.1)
    """

    def __init__(self,
                 node_in_dim: int = 22,
                 edge_in_dim: int = 17,
                 hidden_dim: int  = 64,
                 num_layers: int  = 3,
                 heads: int       = 4,
                 dropout: float   = 0.1):
        super().__init__()

        self.num_layers = num_layers
        self.dropout    = nn.Dropout(dropout)

        # ── Input projection ──────────────────────────────────
        # Projects raw 22-dim node features to hidden_dim.
        # This lets the GNN layers all work in the same space.
        self.node_proj = nn.Linear(node_in_dim, hidden_dim)

        # ── Edge feature projection ────────────────────────────
        # GATConv with edge_dim projects edge features internally.
        # We project edges to hidden_dim for consistency.
        self.edge_proj = nn.Linear(edge_in_dim, hidden_dim)

        # ── GATConv layers ─────────────────────────────────────
        # Each layer: node receives attention-weighted messages
        # from neighbors, incorporating edge features.
        #
        # concat=True means the `heads` attention heads are
        # concatenated → output dim = hidden_dim * heads.
        # We then project back to hidden_dim after each layer
        # so dimensions stay consistent across layers.
        self.gat_layers  = nn.ModuleList()
        self.proj_layers = nn.ModuleList()
        self.norm_layers = nn.ModuleList()

        for i in range(num_layers):
            in_dim = hidden_dim  # input is always hidden_dim (after proj)

            self.gat_layers.append(
                GATConv(
                    in_channels  = in_dim,
                    out_channels = hidden_dim,
                    heads        = heads,
                    edge_dim     = hidden_dim,   # edge features fed into attention
                    concat       = True,         # concatenate heads
                    dropout      = dropout,
                )
            )
            # Project concatenated heads back to hidden_dim
            self.proj_layers.append(
                nn.Linear(hidden_dim * heads, hidden_dim)
            )
            # Layer norm stabilises training
            self.norm_layers.append(
                nn.LayerNorm(hidden_dim)
            )

        # ── Output head ────────────────────────────────────────
        # After pooling over peptide nodes, we have a hidden_dim
        # vector per structure. Two linear layers predict:
        #   output[0] = pep_mean_plddt
        #   output[1] = anchor_mean_plddt
        self.output_head = nn.Sequential(
            nn.Linear(hidden_dim, hidden_dim // 2),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim // 2, 2),        # 2 targets
        )

    def forward(self,
                x: torch.Tensor,
                edge_index: torch.Tensor,
                edge_attr: torch.Tensor,
                batch: torch.Tensor,
                is_peptide_mask: torch.Tensor) -> torch.Tensor:
        """
        Forward pass.

        Args:
            x               : [N_nodes, 22]   node features
            edge_index      : [2, N_edges]     graph connectivity
            edge_attr       : [N_edges, 17]    edge features
            batch           : [N_nodes]        graph id per node
                                               (from PyG DataLoader)
            is_peptide_mask : [N_nodes]        bool, True for peptide nodes

        Returns:
            pred : [batch_size, 2]  predicted pep and anchor pLDDT
        """

        # ── Project inputs to hidden_dim ───────────────────────
        x         = self.node_proj(x)           # [N, hidden_dim]
        edge_attr = self.edge_proj(edge_attr)   # [E, hidden_dim]

        # ── Message passing ────────────────────────────────────
        for gat, proj, norm in zip(
                self.gat_layers, self.proj_layers, self.norm_layers):

            x_new = gat(x, edge_index, edge_attr=edge_attr)
            # x_new shape: [N, hidden_dim * heads]

            x_new = proj(x_new)        # [N, hidden_dim]
            x_new = norm(x_new)        # layer normalisation
            x_new = torch.relu(x_new)  # non-linearity
            x_new = self.dropout(x_new)

            # Residual connection: add input to output
            # This helps gradients flow and stabilises training
            x = x + x_new             # [N, hidden_dim]

        # ── Pool over peptide nodes only ───────────────────────
        # We only care about the peptide's representation,
        # not the whole complex. So we mask out MHC nodes
        # before pooling.
        #
        # is_peptide_mask is True for peptide residues.
        # We use global_mean_pool restricted to those nodes.

        pep_x     = x[is_peptide_mask]              # [N_pep, hidden_dim]
        pep_batch = batch[is_peptide_mask]           # [N_pep]

        # global_mean_pool averages all peptide nodes per graph
        pooled = global_mean_pool(pep_x, pep_batch) # [batch_size, hidden_dim]

        # ── Predict pLDDT ──────────────────────────────────────
        pred = self.output_head(pooled)              # [batch_size, 2]

        return pred


# ══════════════════════════════════════════════════════════════════
# QUICK SANITY CHECK
# ══════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    import torch

    print("Running model sanity check...")

    # Simulate a batch of 2 structures
    # Structure 1: 189 nodes (180 MHC + 9 peptide)
    # Structure 2: 192 nodes (183 MHC + 9 peptide)
    N1, N2   = 189, 192
    N        = N1 + N2       # total nodes in batch
    E        = 2400          # approximate edges

    # Fake node features
    x = torch.randn(N, 22)

    # Fake edge index (random connectivity for testing)
    edge_index = torch.randint(0, N, (2, E))

    # Fake edge features
    edge_attr = torch.randn(E, 17)

    # Batch vector: first N1 nodes belong to graph 0,
    # next N2 nodes belong to graph 1
    batch = torch.cat([
        torch.zeros(N1, dtype=torch.long),
        torch.ones(N2, dtype=torch.long),
    ])

    # Peptide mask: last 9 nodes of each structure are peptide
    is_peptide = torch.zeros(N, dtype=torch.bool)
    is_peptide[N1-9 : N1]     = True   # peptide of structure 1
    is_peptide[N1+N2-9 : ]    = True   # peptide of structure 2

    # Build model
    model = PMTopoGNN(
        node_in_dim = 22,
        edge_in_dim = 17,
        hidden_dim  = 64,
        num_layers  = 3,
        heads       = 4,
        dropout     = 0.1,
    )

    print(f"Model parameters: "
          f"{sum(p.numel() for p in model.parameters()):,}")

    # Forward pass
    pred = model(x, edge_index, edge_attr, batch, is_peptide)

    print(f"Input:  {N} nodes, {E} edges, batch_size=2")
    print(f"Output: {pred.shape}  (should be [2, 2])")
    print(f"Predictions:\n  structure 0: pep={pred[0,0]:.2f}, anchor={pred[0,1]:.2f}")
    print(f"               structure 1: pep={pred[1,0]:.2f}, anchor={pred[1,1]:.2f}")
    print("Sanity check passed.")