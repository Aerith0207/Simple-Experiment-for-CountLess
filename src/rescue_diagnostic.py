"""
Count-Less CL-CU vs CL-MU: 100-packet flow cross-layer rescue analysis.

Run from the project root (with SplitMix64):
    python src/rescue_diagnostic.py

This is a software-mechanism diagnostic, NOT an official paper reproduction.
"""

import csv
from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from experiment import CountLess, make_hashes

# ---------------- User-controlled parameters ----------------
MEMORY_BYTES = 8192
FLOW_COUNTS = [5000, 10000, 20000]
SEEDS = [42, 43, 44, 45, 46]
TARGET_SIZE = 100
THRESHOLD = 0.10   # The paper's strict mouse flow survival condition
PROJECT_ROOT = Path(__file__).resolve().parents[1]
OUTPUT = PROJECT_ROOT / "results" / "04_rescue_analysis" / (
    "rescue_" + datetime.now().strftime("%Y%m%d_%H%M%S")
)
OUTPUT.mkdir(parents=True, exist_ok=True)


def make_trace(n_flows, seed):
    """Use the SAME size mix and packet shuffle as the previous stress test."""
    rng = np.random.default_rng(seed)
    n1 = int(n_flows * 0.80)
    n10 = int(n_flows * 0.15)
    n100 = int(n_flows * 0.04)
    n1000 = n_flows - n1 - n10 - n100
    truth = np.array(
        [1] * n1 + [10] * n10 + [100] * n100 + [1000] * n1000,
        dtype=np.int64,
    )
    rng.shuffle(truth)
    packets = np.repeat(np.arange(n_flows, dtype=np.int32), truth)
    rng.shuffle(packets)
    return packets, truth


def check_hash_distribution():
    """Catch accidental use of the old linear hash on consecutive flow IDs."""
    sketch = CountLess(MEMORY_BYTES, make_hashes(42 + 2026), "CL-MU")
    occupied = len({sketch.positions(i)[0] for i in range(5000)})
    print(f"Hash check (5000 IDs, L1 occupied positions): {occupied}")
    if not (2500 <= occupied <= 3600):
        raise RuntimeError(
            "Unexpected hash distribution. Verify experiment.py uses the "
            "previously tested SplitMix64 hash, NOT the old linear hash."
        )


def layer_diagnostics(sketch, flow_ids, truth):
    """For a group of flows, compare each valid layer against the same threshold."""
    layer_values = np.array([
        [sketch.reg[j][sketch.positions(int(fid))[j]] for j in range(3)]
        for fid in flow_ids
    ], dtype=np.int64)
    limits = np.asarray(sketch.limits, dtype=np.int64)
    valid = layer_values < limits[None, :]

    sizes = truth[flow_ids]
    layer_error = (layer_values - sizes[:, None]) / sizes[:, None]
    layer_good = valid & (np.abs(layer_error) < THRESHOLD)

    final_values = np.asarray([sketch.query(int(fid)) for fid in flow_ids], dtype=np.int64)
    if not np.all(final_values >= sizes):
        raise AssertionError("Underestimation occurred; inspect overflow/query behavior")

    # query() is the minimum over non-saturated layers
    valid_min = np.min(np.where(valid, layer_values, np.iinfo(np.int64).max), axis=1)
    if not np.array_equal(valid_min, final_values):
        raise AssertionError("Layerwise query reconstruction disagrees with CountLess.query()")

    final_good = (final_values - sizes) / sizes < THRESHOLD
    assert np.array_equal(final_good, layer_good.any(axis=1))

    rescued = ~layer_good[:, 0] & final_good
    l1_failed = ~layer_good[:, 0]
    return {
        "values": layer_values,
        "valid": valid,
        "layer_good": layer_good,
        "final_values": final_values,
        "final_good": final_good,
        "rescue": rescued,
        "l1_failed": l1_failed,
        "metrics": {
            "L1_FSR_pct": 100 * layer_good[:, 0].mean(),
            "L2_FSR_pct": 100 * layer_good[:, 1].mean(),
            "L3_FSR_pct": 100 * layer_good[:, 2].mean(),
            "Final_FSR_pct": 100 * final_good.mean(),
            "Rescue_pct": 100 * rescued.mean(),  # of ALL 100-packet flows
            "Rescue_given_L1_fail_pct": 100 * rescued.sum() / max(1, l1_failed.sum()),
            "All_layers_fail_pct": 100 * (~layer_good.any(axis=1)).mean(),
            "ARE_100": np.mean((final_values - sizes) / sizes),
        },
    }


def save_csv(path, rows):
    with path.open("w", newline="", encoding="utf-8-sig") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def run():
    check_hash_distribution()
    rows = []
    cases = []
    for n_flows in FLOW_COUNTS:
        for seed in SEEDS:
            packets, truth = make_trace(n_flows, seed)
            target_ids = np.flatnonzero(truth == TARGET_SIZE)
            hashes = make_hashes(seed + 2026)
            per_method = {}
            for method in ("CL-CU", "CL-MU"):
                sketch = CountLess(MEMORY_BYTES, hashes, method)
                assert sketch.used_bits <= MEMORY_BYTES * 8
                for flow in packets:
                    sketch.update(int(flow))
                diag = layer_diagnostics(sketch, target_ids, truth)
                per_method[method] = diag
                rows.append({
                    "n_flows": n_flows,
                    "n_packets": len(packets),
                    "seed": seed,
                    "method": method,
                    "n_100_flows": len(target_ids),
                    **diag["metrics"],
                })

            cu = per_method["CL-CU"]
            mu = per_method["CL-MU"]
            cu_good_mu_bad = cu["final_good"] & ~mu["final_good"]
            mu_good_cu_bad = mu["final_good"] & ~cu["final_good"]
            print(
                f"Flows={n_flows:5d}, seed={seed}, "
                f"100-pkt FSR CU={100 * cu['final_good'].mean():5.1f}%, "
                f"MU={100 * mu['final_good'].mean():5.1f}%, "
                f"CU-only good={cu_good_mu_bad.sum():3d}, "
                f"MU-only good={mu_good_cu_bad.sum():3d}"
            )
            if n_flows == max(FLOW_COUNTS):
                # Save examples where CU succeeds but MU fails, for diagnosis.
                bad_indices = np.flatnonzero(cu_good_mu_bad)[:25]
                for k in bad_indices:
                    cases.append({
                        "seed": seed,
                        "flow_id": int(target_ids[k]),
                        "true_size": TARGET_SIZE,
                        "CU_est": int(cu["final_values"][k]),
                        "MU_est": int(mu["final_values"][k]),
                        **{f"CU_L{j+1}": int(cu["values"][k, j]) for j in range(3)},
                        **{f"MU_L{j+1}": int(mu["values"][k, j]) for j in range(3)},
                        "CU_rescued": bool(cu["rescue"][k]),
                        "MU_rescued": bool(mu["rescue"][k]),
                    })

    save_csv(OUTPUT / "rescue_details.csv", rows)
    if cases:
        save_csv(OUTPUT / "cu_good_mu_bad_examples.csv", cases)

    metrics = ["L1_FSR_pct", "L2_FSR_pct", "L3_FSR_pct", "Final_FSR_pct",
               "Rescue_pct", "Rescue_given_L1_fail_pct", "All_layers_fail_pct", "ARE_100"]
    summary = []
    for n in FLOW_COUNTS:
        for method in ("CL-CU", "CL-MU"):
            group = [r for r in rows if r["n_flows"] == n and r["method"] == method]
            row = {"n_flows": n, "method": method, "n_seeds": len(group)}
            for metric in metrics:
                a = np.asarray([r[metric] for r in group], dtype=float)
                row[metric + "_mean"] = float(a.mean())
                row[metric + "_std"] = float(a.std(ddof=1)) if len(a) > 1 else 0.0
            summary.append(row)
    save_csv(OUTPUT / "rescue_summary.csv", summary)

    for metric, label in [("Final_FSR_pct", "100-packet Final FSR (%)"),
                          ("Rescue_pct", "100-packet Cross-layer Rescue (% of group)")]:
        fig, ax = plt.subplots(figsize=(8, 4.8))
        for method in ("CL-CU", "CL-MU"):
            grp = [r for r in summary if r["method"] == method]
            x = [r["n_flows"] for r in grp]
            y = [r[metric + "_mean"] for r in grp]
            sd = [r[metric + "_std"] for r in grp]
            ax.errorbar(x, y, yerr=sd, marker="o", capsize=3, label=method)
        ax.set_xlabel("Number of distinct flows")
        ax.set_ylabel(label)
        ax.set_ylim(bottom=0, top=105 if metric == "Final_FSR_pct" else None)
        ax.grid(alpha=0.25)
        ax.legend()
        fig.tight_layout()
        fig.savefig(OUTPUT / f"{metric}_comparison.png", dpi=170)
        plt.close(fig)

    print("\n=== Five-seed mean at highest load ===")
    for r in summary:
        if r["n_flows"] == max(FLOW_COUNTS):
            print(
                f"{r['method']}: L1={r['L1_FSR_pct_mean']:.2f}%, "
                f"Rescue={r['Rescue_pct_mean']:.2f}%, "
                f"FinalFSR={r['Final_FSR_pct_mean']:.2f}%, "
                f"ARE100={r['ARE_100_mean']:.4f}"
            )
    print("\nSaved output to:", OUTPUT)


if __name__ == "__main__":
    run()
