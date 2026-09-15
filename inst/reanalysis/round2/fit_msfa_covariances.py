#!/usr/bin/env python
"""Fit the shared plus zone-specific MSFA covariance model from sufficient statistics.

The Gaussian negative log-likelihood uses each zone's covariance and exact row count,
so every spatially adjusted cell contributes without constructing per-cell factor lists.
"""
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
import torch

torch.set_num_threads(4)
torch.set_default_dtype(torch.float64)


def positive_factors(matrix, rank):
    vals, vecs = np.linalg.eigh((matrix + matrix.T) / 2)
    order = np.argsort(vals)[::-1][:rank]
    vals = np.maximum(vals[order], 0)
    return vecs[:, order] * np.sqrt(vals)[None, :]


def partial_corr(sigma, ridge=1e-6):
    theta = np.linalg.inv(sigma + np.eye(sigma.shape[0]) * ridge)
    den = np.sqrt(np.outer(np.diag(theta), np.diag(theta)))
    pcor = -theta / den
    np.fill_diagonal(pcor, 1.0)
    return np.clip((pcor + pcor.T) / 2, -1, 1)


def fit_one(zones, rank, steps=350, seed=2026):
    torch.manual_seed(seed)
    covs = np.stack([np.asarray(z["cov"], dtype=float) for z in zones])
    ns = np.array([int(z["n"]) for z in zones], dtype=float)
    g = covs.shape[1]
    if covs.shape != (len(zones), g, g):
        raise ValueError("Each covariance must be square and share cell-type order")
    if np.any(ns < g + 1):
        raise ValueError("Each zone needs more rows than variables")
    if not np.all(np.isfinite(covs)):
        raise ValueError("Non-finite covariance")
    scales = np.sqrt(np.maximum(np.mean(np.diagonal(covs, axis1=1, axis2=2), axis=0), 1e-8))
    corrs = covs / scales[None, :, None] / scales[None, None, :]
    pooled = np.average(corrs, axis=0, weights=ns)
    phi0 = positive_factors(pooled * 0.65, rank)
    lam0 = np.stack([positive_factors((c - phi0 @ phi0.T) * 0.55, rank) for c in corrs])
    psi0 = np.maximum(np.diagonal(corrs, axis1=1, axis2=2) -
                      np.sum(phi0**2, axis=1)[None, :] -
                      np.sum(lam0**2, axis=2), 0.05)
    phi = torch.nn.Parameter(torch.tensor(np.ascontiguousarray(phi0)))
    lam = torch.nn.Parameter(torch.tensor(np.ascontiguousarray(lam0)))
    logpsi = torch.nn.Parameter(torch.tensor(np.ascontiguousarray(np.log(psi0))))
    target = torch.tensor(corrs)
    weights = torch.tensor(ns / ns.sum())
    opt = torch.optim.LBFGS([phi, lam, logpsi], lr=0.6, max_iter=steps,
                            tolerance_grad=1e-8, tolerance_change=1e-10,
                            line_search_fn="strong_wolfe")
    state = {"calls": 0}

    def closure():
        opt.zero_grad()
        psi = torch.nn.functional.softplus(logpsi) + 1e-5
        shared = phi @ phi.T
        sigma = shared[None, :, :] + lam @ lam.transpose(1, 2) + torch.diag_embed(psi)
        chol = torch.linalg.cholesky(sigma)
        logdet = 2 * torch.log(torch.diagonal(chol, dim1=1, dim2=2)).sum(1)
        solved = torch.cholesky_solve(target, chol)
        trace = torch.diagonal(solved, dim1=1, dim2=2).sum(1)
        loss = ((logdet + trace) * weights).sum() / 2
        loss = loss + 0.004 * phi.square().mean() + 0.012 * lam.square().mean()
        loss.backward()
        state["calls"] += 1
        return loss

    final = float(opt.step(closure).detach())
    with torch.no_grad():
        ph = phi.numpy().copy()
        la = lam.numpy().copy()
        ps = (torch.nn.functional.softplus(logpsi) + 1e-5).numpy().copy()
    shared = (ph @ ph.T) * scales[:, None] * scales[None, :]
    full = np.stack([(ph @ ph.T + la[q] @ la[q].T + np.diag(ps[q])) *
                     scales[:, None] * scales[None, :] for q in range(len(zones))])
    return {"shared": shared, "full": full, "phi": ph * scales[:, None],
            "lambda": la * scales[None, :, None], "psi": ps * scales[None, :]**2,
            "loss": final, "calls": state["calls"]}


def run(manifest_path, output_path, rank, steps, seed):
    manifest = json.loads(Path(manifest_path).read_text(encoding="utf-8"))
    labels = manifest["cell_types"]
    zones = []
    for entry in manifest["zones"]:
        cov = pd.read_csv(entry["cov"], index_col=0).loc[labels, labels].to_numpy()
        zones.append({"name": entry["name"], "n": entry["n"], "cov": cov})
    diag = np.stack([np.diag(z["cov"]) for z in zones])
    active = np.flatnonzero(np.max(diag, axis=0) > 1e-10)
    if len(active) < 2:
        raise ValueError("At least two nonconstant cell types are needed")
    modeled = [{"name": z["name"], "n": z["n"],
                "cov": z["cov"][np.ix_(active, active)]} for z in zones]
    result = fit_one(modeled, min(rank, len(active) - 1), steps, seed)
    output = Path(output_path)
    output.mkdir(parents=True, exist_ok=True)

    def embed_cov(small):
        full = np.zeros((len(labels), len(labels)))
        full[np.ix_(active, active)] = small
        return full

    def embed_pcor(small):
        full = np.eye(len(labels))
        full[np.ix_(active, active)] = partial_corr(small)
        return full

    for q, zone in enumerate(zones):
        tag = zone["name"].replace(" ", "_")
        cov = result["full"][q]
        pd.DataFrame(embed_cov(cov), index=labels, columns=labels).to_csv(output / f"cov_{tag}.csv")
        pd.DataFrame(embed_pcor(cov), index=labels, columns=labels).to_csv(output / f"pcor_{tag}.csv")
    pd.DataFrame(embed_cov(result["shared"]), index=labels, columns=labels).to_csv(output / "cov_Shared.csv")
    pd.DataFrame(embed_pcor(result["shared"]), index=labels, columns=labels).to_csv(output / "pcor_Shared.csv")
    (output / "factor_fit.json").write_text(json.dumps({
        "method": "MSFA covariance maximum likelihood from sufficient statistics",
        "rank_shared": min(rank, len(active) - 1), "rank_zone": min(rank, len(active) - 1),
        "inactive_cell_types": [labels[i] for i in range(len(labels)) if i not in active],
        "n": {z["name"]: z["n"] for z in zones},
        "loss": result["loss"], "optimizer_calls": result["calls"], "seed": seed,
    }, indent=2), encoding="utf-8")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--rank", type=int, required=True)
    parser.add_argument("--steps", type=int, default=350)
    parser.add_argument("--seed", type=int, default=2026)
    args = parser.parse_args()
    run(args.manifest, args.output, args.rank, args.steps, args.seed)
