"""
gnn.py - Graph Neural Network training (GAT + GraphSAGE)
=========================================================
- Detects PyTorch + torch-geometric availability gracefully.
- When unavailable, NumpyFallback runs HistGradientBoostingClassifier
  with the same train/test interface.
- All per-fold scaling is fit on the training fold only.
- validate_no_leakage is called at the start of every fold.
- Per-fold alpha (GNN/bypass gate) values are persisted to
  results/gnn_alpha_per_fold.csv.
- All temporary array modifications use _with_modified_matrices for safe
  restoration via try/finally.

Public:
    run_gnn_pipeline(arch, dataset, matrices, valid_species, valid_genera,
                     cfg, seed) -> dict
"""
from __future__ import annotations

import logging
from copy import deepcopy
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.preprocessing import StandardScaler

from build_labels import DataLeakageError, validate_no_leakage
from evaluation import metrics as compute_metrics
from model import scale_train_test, pair_feature_matrix, get_device

log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Optional torch import
# ---------------------------------------------------------------------------
try:
    import torch
    import torch.nn as nn
    import torch.nn.functional as F
    from torch_geometric.nn import GATConv, SAGEConv
    HAS_TORCH = True
except ImportError:
    torch = None  # type: ignore
    HAS_TORCH = False


# ---------------------------------------------------------------------------
# Modified-matrices context helper
# ---------------------------------------------------------------------------
def _with_modified_matrices(matrices: dict, modifications: dict, fn):
    """Temporarily overwrite matrices[k] in place, run fn, then restore."""
    saved: dict = {}
    for k in modifications:
        v = matrices.get(k)
        saved[k] = v.copy() if hasattr(v, "copy") else v
    try:
        for k, v in modifications.items():
            target = matrices[k]
            if hasattr(target, "shape") and hasattr(v, "shape"):
                target[:] = v
            else:
                matrices[k] = v
        return fn()
    finally:
        for k, v in saved.items():
            target = matrices.get(k)
            if hasattr(target, "shape") and hasattr(v, "shape"):
                target[:] = v
            else:
                matrices[k] = v


# ---------------------------------------------------------------------------
# Numpy fallback model
# ---------------------------------------------------------------------------
class NumpyFallback:
    """Classical-ML fallback when PyTorch is unavailable."""

    def __init__(self, cfg: dict, seed: int):
        self.clf = HistGradientBoostingClassifier(
            max_iter=300, max_depth=6, learning_rate=0.06,
            random_state=seed)
        self.scaler = StandardScaler()
        self._fitted = False

    def fit(self, X_train: np.ndarray, y_train: np.ndarray):
        X = self.scaler.fit_transform(X_train).astype(np.float32)
        X = np.where(np.isfinite(X), X, 0.0)
        self.clf.fit(X, y_train)
        self._fitted = True

    def predict_proba(self, X_test: np.ndarray) -> np.ndarray:
        assert self._fitted, "Must call fit before predict_proba"
        X = self.scaler.transform(X_test).astype(np.float32)
        X = np.where(np.isfinite(X), X, 0.0)
        p = self.clf.predict_proba(X)[:, 1]
        return np.clip(np.where(np.isfinite(p), p, 0.5), 0.0, 1.0).astype(np.float32)


# ---------------------------------------------------------------------------
# Torch model definitions
# ---------------------------------------------------------------------------
if HAS_TORCH:
    class _GNNModel(nn.Module):
        def __init__(self, arch: str, in_phage: int, in_host: int,
                     hidden: int, out: int, dropout: float,
                     n_layers: int, gat_heads: int, edge_dim: int):
            super().__init__()
            self.arch = arch
            self.in_phage = in_phage
            self.in_host = in_host
            # project both node types into a common hidden dim
            self.p_proj = nn.Linear(in_phage, hidden)
            self.h_proj = nn.Linear(in_host, hidden)
            self.edge_mlp = (nn.Sequential(nn.Linear(edge_dim, hidden), nn.ReLU(),
                                            nn.Linear(hidden, hidden))
                             if edge_dim > 0 else None)

            self.convs = nn.ModuleList()
            for _ in range(n_layers):
                if arch == "GAT":
                    self.convs.append(GATConv(hidden, hidden // max(1, gat_heads),
                                              heads=gat_heads, dropout=dropout,
                                              add_self_loops=True))
                else:
                    self.convs.append(SAGEConv(hidden, hidden))
            self.dropout = nn.Dropout(dropout)
            self.head = nn.Sequential(
                nn.Linear(hidden * 2 + (hidden if edge_dim > 0 else 0), hidden),
                nn.ReLU(),
                nn.Dropout(dropout),
                nn.Linear(hidden, 1),
            )
            self.bypass = nn.Sequential(
                nn.Linear(in_phage + in_host + edge_dim, hidden),
                nn.ReLU(),
                nn.Dropout(dropout),
                nn.Linear(hidden, 1),
            )
            self.alpha = nn.Parameter(torch.tensor(0.0))  # sigmoid -> 0.5

        def encode(self, x_p, x_h, edge_index):
            n_p = x_p.size(0)
            x = torch.cat([self.p_proj(x_p), self.h_proj(x_h)], dim=0)
            for conv in self.convs:
                x = conv(x, edge_index)
                x = F.elu(x)
                x = self.dropout(x)
            return x, n_p

        def forward(self, x_p, x_h, edge_index, src, dst, edge_feats_dst):
            x_all, n_p = self.encode(x_p, x_h, edge_index)
            z_p = x_all[src]
            z_h = x_all[dst + n_p]
            if self.edge_mlp is not None and edge_feats_dst is not None:
                e = self.edge_mlp(edge_feats_dst)
                logits_gnn = self.head(torch.cat([z_p, z_h, e], dim=1)).squeeze(-1)
            else:
                logits_gnn = self.head(torch.cat([z_p, z_h], dim=1)).squeeze(-1)
            # bypass branch operates on raw concat features
            bp = torch.cat([x_p[src], x_h[dst], edge_feats_dst], dim=1) \
                if edge_feats_dst is not None else torch.cat([x_p[src], x_h[dst]], dim=1)
            logits_bp = self.bypass(bp).squeeze(-1)
            a = torch.sigmoid(self.alpha)
            return a * logits_gnn + (1 - a) * logits_bp


# ---------------------------------------------------------------------------
# GAT graph builder (directed edges, no manual self-loop dup)
# ---------------------------------------------------------------------------
def _build_train_edges(dataset: pd.DataFrame, train_mask: np.ndarray,
                       edge_feats_scaled: np.ndarray, n_phages: int):
    """Return (edge_index [2, E], edge_attr [E, F]) for positive training edges only."""
    assert HAS_TORCH
    pm = (dataset["label"].values == 1) & train_mask
    if not pm.any():
        return (torch.zeros((2, 0), dtype=torch.long),
                torch.zeros((0, edge_feats_scaled.shape[1]), dtype=torch.float32))
    src = torch.tensor(dataset.loc[pm, "phage_idx"].astype(int).values, dtype=torch.long)
    dst = torch.tensor(dataset.loc[pm, "host_idx"].astype(int).values + n_phages, dtype=torch.long)
    ef = torch.tensor(edge_feats_scaled[pm], dtype=torch.float32)
    ei = torch.stack([src, dst], dim=0)
    return ei, ef


# ---------------------------------------------------------------------------
# Per-fold training
# ---------------------------------------------------------------------------
def _run_fold_torch(arch: str, train_mask: np.ndarray, test_mask: np.ndarray,
                    dataset: pd.DataFrame, matrices: dict, cfg: dict,
                    seed: int, val_frac: float = 0.15) -> tuple[np.ndarray, np.ndarray, dict]:
    """Train one fold using torch.  Returns (probs_on_test, y_test, fold_meta)."""
    assert HAS_TORCH
    device = get_device() or torch.device("cpu")
    torch.manual_seed(seed)

    train_df = dataset[train_mask].reset_index(drop=True)
    test_df = dataset[test_mask].reset_index(drop=True)
    try:
        validate_no_leakage(train_df, test_df, context=f"{arch} fold")
    except DataLeakageError as e:
        log.error(f"  Leakage in fold ({arch}): {e}")
        return np.array([]), np.array([]), {"skipped": True, "reason": "leakage"}

    # Per-fold scaler on edge features (train only)
    E_raw = matrices["EDGE_FEATS_RAW"]
    e_train = E_raw[train_mask]
    e_test = E_raw[test_mask]
    scaler = StandardScaler()
    if e_train.shape[1]:
        scaler.fit(e_train)
        e_train_s = scaler.transform(e_train).astype(np.float32)
        e_full_s = scaler.transform(E_raw).astype(np.float32)
        e_test_s = scaler.transform(e_test).astype(np.float32)
    else:
        e_train_s = e_train.astype(np.float32)
        e_full_s = E_raw.astype(np.float32)
        e_test_s = e_test.astype(np.float32)
    e_full_s = np.where(np.isfinite(e_full_s), e_full_s, 0.0).astype(np.float32)
    e_test_s = np.where(np.isfinite(e_test_s), e_test_s, 0.0).astype(np.float32)

    # Node feature scaling (fit on phages/hosts that appear in training only)
    P = matrices["PHAGE_BASE"]
    H = matrices["HOST_BASE"]
    n_phages = matrices["n_phages"]
    n_hosts = matrices["n_hosts"]

    train_p = train_df["phage_idx"].dropna().astype(int).unique()
    train_h = train_df["host_idx"].dropna().astype(int).unique()
    P_scaler = StandardScaler().fit(P[train_p]) if len(train_p) > 1 else StandardScaler().fit(P)
    H_scaler = StandardScaler().fit(H[train_h]) if len(train_h) > 1 else StandardScaler().fit(H)
    P_s = P_scaler.transform(P).astype(np.float32)
    H_s = H_scaler.transform(H).astype(np.float32)
    P_s = np.where(np.isfinite(P_s), P_s, 0.0)
    H_s = np.where(np.isfinite(H_s), H_s, 0.0)

    edge_index, edge_attr = _build_train_edges(dataset, train_mask, e_full_s, n_phages)

    # Build inputs
    x_p = torch.tensor(P_s, dtype=torch.float32, device=device)
    x_h = torch.tensor(H_s, dtype=torch.float32, device=device)
    edge_index = edge_index.to(device)

    hidden = int(cfg["models"]["hidden_dim"])
    out_dim = int(cfg["models"]["out_dim"])
    dropout = float(cfg["models"]["dropout"])
    lr = float(cfg["models"]["learning_rate"])
    wd = float(cfg["models"]["weight_decay"])
    epochs = int(cfg["models"]["epochs"])
    patience = int(cfg["models"]["patience"])
    gat_heads = int(cfg["models"]["gat_heads"])
    n_layers = int(cfg["models"]["gat_layers"])

    model = _GNNModel(arch=arch,
                      in_phage=P.shape[1], in_host=H.shape[1],
                      hidden=hidden, out=out_dim, dropout=dropout,
                      n_layers=n_layers, gat_heads=gat_heads,
                      edge_dim=e_full_s.shape[1]).to(device)
    opt = torch.optim.Adam(model.parameters(), lr=lr, weight_decay=wd)

    # Pair tensors
    src_all = torch.tensor(dataset["phage_idx"].astype(int).values, dtype=torch.long, device=device)
    dst_all = torch.tensor(dataset["host_idx"].astype(int).values, dtype=torch.long, device=device)
    y_all = torch.tensor(dataset["label"].astype(float).values, dtype=torch.float32, device=device)
    e_all = torch.tensor(e_full_s, dtype=torch.float32, device=device)

    train_idx = np.where(train_mask)[0]
    test_idx = np.where(test_mask)[0]
    rng = np.random.default_rng(seed)
    rng.shuffle(train_idx)
    n_val = max(1, int(len(train_idx) * val_frac))
    val_idx = train_idx[:n_val]
    train_idx_inner = train_idx[n_val:]
    train_idx_inner_t = torch.tensor(train_idx_inner, dtype=torch.long, device=device)
    val_idx_t = torch.tensor(val_idx, dtype=torch.long, device=device)
    test_idx_t = torch.tensor(test_idx, dtype=torch.long, device=device)

    bce = nn.BCEWithLogitsLoss()
    best_val_auc = -np.inf
    best_state = None
    bad = 0
    epochs_trained = 0

    for ep in range(epochs):
        model.train()
        opt.zero_grad()
        logits = model(x_p, x_h, edge_index,
                       src_all[train_idx_inner_t],
                       dst_all[train_idx_inner_t],
                       e_all[train_idx_inner_t])
        loss = bce(logits, y_all[train_idx_inner_t])
        loss.backward()
        opt.step()
        epochs_trained = ep + 1

        if (ep + 1) % 2 == 0 or ep == 0:
            model.eval()
            with torch.no_grad():
                v_logits = model(x_p, x_h, edge_index,
                                 src_all[val_idx_t], dst_all[val_idx_t],
                                 e_all[val_idx_t])
                v_probs = torch.sigmoid(v_logits).cpu().numpy()
                y_val = y_all[val_idx_t].cpu().numpy()
            try:
                from sklearn.metrics import roc_auc_score
                v_auc = float(roc_auc_score(y_val, v_probs))
            except Exception:
                v_auc = 0.5
            if v_auc > best_val_auc + 1e-4:
                best_val_auc = v_auc
                best_state = deepcopy(model.state_dict())
                bad = 0
            else:
                bad += 1
                if bad >= patience:
                    break

    if best_state is not None:
        model.load_state_dict(best_state)
    model.eval()
    with torch.no_grad():
        t_logits = model(x_p, x_h, edge_index,
                         src_all[test_idx_t], dst_all[test_idx_t],
                         e_all[test_idx_t])
        probs = torch.sigmoid(t_logits).cpu().numpy().astype(np.float32)
    y_true = y_all[test_idx_t].cpu().numpy().astype(np.float32)
    alpha_val = float(torch.sigmoid(model.alpha).item())
    if alpha_val > 0.6:
        cat = "GNN-dominant"
    elif alpha_val < 0.4:
        cat = "bypass-dominant"
    else:
        cat = "balanced"
    log.info(f"    alpha={alpha_val:.4f} ({cat})")
    fold_meta = {
        "alpha": alpha_val,
        "best_val_auc": float(best_val_auc) if np.isfinite(best_val_auc) else 0.5,
        "epochs_trained": int(epochs_trained),
        "skipped": False,
    }
    return probs, y_true, fold_meta


def _run_fold_fallback(train_mask: np.ndarray, test_mask: np.ndarray,
                       dataset: pd.DataFrame, matrices: dict, cfg: dict,
                       seed: int) -> tuple[np.ndarray, np.ndarray, dict]:
    """Classical fallback fold using NumpyFallback."""
    train_df = dataset[train_mask].reset_index(drop=True)
    test_df = dataset[test_mask].reset_index(drop=True)
    try:
        validate_no_leakage(train_df, test_df, context="GNN-fallback fold")
    except DataLeakageError as e:
        log.error(f"  Leakage in fallback fold: {e}")
        return np.array([]), np.array([]), {"skipped": True, "reason": "leakage"}

    X = pair_feature_matrix(dataset, matrices, use_sanitized_names=False)
    y = dataset["label"].astype(float).values
    fb = NumpyFallback(cfg, seed)
    fb.fit(X[train_mask], y[train_mask])
    probs = fb.predict_proba(X[test_mask])
    return probs.astype(np.float32), y[test_mask].astype(np.float32), {
        "alpha": float("nan"), "best_val_auc": float("nan"),
        "epochs_trained": 0, "skipped": False, "fallback": True,
    }


def _run_fold(arch: str, train_mask: np.ndarray, test_mask: np.ndarray,
              dataset: pd.DataFrame, matrices: dict, cfg: dict, seed: int):
    if HAS_TORCH:
        try:
            return _run_fold_torch(arch, train_mask, test_mask, dataset,
                                   matrices, cfg, seed)
        except Exception as e:
            log.exception(f"  Torch fold failed ({arch}): {e}; using NumpyFallback")
    return _run_fold_fallback(train_mask, test_mask, dataset, matrices, cfg, seed)


# ---------------------------------------------------------------------------
# Pipeline runner
# ---------------------------------------------------------------------------
def _loso_iter(dataset: pd.DataFrame, valid_species: list[str]):
    for sp in valid_species:
        test_mask = (dataset["host"] == sp).values
        train_mask = ~test_mask
        yield sp, train_mask, test_mask


def _logo_iter(dataset: pd.DataFrame, valid_genera: list[str]):
    for gen in valid_genera:
        test_mask = (dataset["host_genus"] == gen).values
        train_mask = ~test_mask
        yield gen, train_mask, test_mask


def run_gnn_pipeline(arch: str, dataset: pd.DataFrame, matrices: dict,
                     valid_species: list[str], valid_genera: list[str],
                     cfg: dict, seed: int) -> dict:
    results_dir = Path(cfg["paths"]["results_dir"])
    results_dir.mkdir(parents=True, exist_ok=True)

    loso_rows = []
    pooled_probs = []
    pooled_y = []
    pred_records = []
    alpha_records = []
    loso_path = results_dir / f"gnn_{arch.lower()}_loso.csv"

    for i, (sp, tr, te) in enumerate(_loso_iter(dataset, valid_species)):
        probs, y_true, meta = _run_fold(arch, tr, te, dataset, matrices, cfg, seed)
        if meta.get("skipped") or len(probs) == 0:
            continue
        m = compute_metrics(y_true, probs)
        loso_rows.append({"species": sp, "n_test": int(te.sum()),
                          **{k: float(v) for k, v in m.items()
                             if isinstance(v, (int, float))},
                          "alpha": meta.get("alpha", float("nan")),
                          "best_val_auc": meta.get("best_val_auc", float("nan")),
                          "epochs_trained": meta.get("epochs_trained", 0)})
        pooled_probs.append(probs)
        pooled_y.append(y_true)
        alpha_records.append({"arch": arch, "species": sp,
                              "alpha": meta.get("alpha", float("nan"))})
        # incremental persistence
        pd.DataFrame(loso_rows).to_csv(loso_path, index=False)
        # capture per-pair predictions
        idx = np.where(te)[0]
        for j, ix in enumerate(idx):
            pred_records.append({"species": sp,
                                 "phage": dataset.iloc[ix]["phage"],
                                 "host": dataset.iloc[ix]["host"],
                                 "y_true": float(y_true[j]),
                                 "y_prob": float(probs[j])})

    loso_df = pd.DataFrame(loso_rows)
    loso_mean = float(loso_df["roc_auc"].mean()) if len(loso_df) else 0.0
    loso_std = float(loso_df["roc_auc"].std()) if len(loso_df) > 1 else 0.0
    if pooled_probs:
        pooled_probs_v = np.concatenate(pooled_probs)
        pooled_y_v = np.concatenate(pooled_y)
        try:
            from sklearn.metrics import roc_auc_score
            loso_pooled = float(roc_auc_score(pooled_y_v, pooled_probs_v))
        except Exception:
            loso_pooled = 0.5
    else:
        pooled_probs_v = np.array([])
        pooled_y_v = np.array([])
        loso_pooled = 0.5

    pd.DataFrame(alpha_records).to_csv(
        results_dir / f"gnn_{arch.lower()}_alpha_per_fold.csv", index=False)

    # LOGO loop (use logo_validation module for consistency)
    try:
        from logo_validation import run_logo
        logo_df = run_logo(arch, dataset, matrices, valid_genera, cfg, seed,
                          fold_fn=lambda tr, te: _run_fold(arch, tr, te, dataset, matrices, cfg, seed))
    except Exception as e:
        log.warning(f"  LOGO failed for {arch}: {e}")
        logo_df = pd.DataFrame()
    logo_mean = float(logo_df["roc_auc"].mean()) if "roc_auc" in logo_df else 0.0

    # Unseen-strain MC
    try:
        from strain_sim import run_unseen_strain
        mc_df = run_unseen_strain(arch, dataset, matrices, cfg, seed,
                                 fold_fn=lambda tr, te: _run_fold(arch, tr, te, dataset, matrices, cfg, seed))
    except Exception as e:
        log.warning(f"  Strain sim failed for {arch}: {e}")
        mc_df = pd.DataFrame()
    mc_auc = float(mc_df["auc"].mean()) if "auc" in mc_df else 0.0

    # Cocktail optimization
    try:
        from cocktail_optimizer import run_cocktail
        ctail_df = run_cocktail(arch, dataset, pd.DataFrame(pred_records),
                                cfg, seed)
        means_c = ctail_df.get("mean_coverage", pd.Series([0.0])).mean()
        pct75_c = float((ctail_df.get("coverage_geq_75", pd.Series([0])).mean()) if len(ctail_df) else 0.0)
    except Exception as e:
        log.warning(f"  Cocktail failed for {arch}: {e}")
        ctail_df = pd.DataFrame()
        means_c = 0.0
        pct75_c = 0.0

    pred_df = pd.DataFrame(pred_records)
    pred_df.to_csv(results_dir / f"gnn_{arch.lower()}_predictions.csv", index=False)

    return {
        "loso_df": loso_df,
        "loso_mean": loso_mean,
        "loso_std": loso_std,
        "loso_pooled": loso_pooled,
        "logo_df": logo_df,
        "logo_mean": logo_mean,
        "mc_df": mc_df,
        "mc_auc": mc_auc,
        "ctail_df": ctail_df,
        "means_c": float(means_c),
        "pct75_c": float(pct75_c),
        "pred_df": pred_df,
        "all_proba": pooled_probs_v,
        "y_true": pooled_y_v,
        "alpha_values": [r["alpha"] for r in alpha_records],
        "mean_auc": loso_mean,
        "per_species_auc": loso_df.set_index("species")["roc_auc"] if "roc_auc" in loso_df else pd.Series(),
    }
