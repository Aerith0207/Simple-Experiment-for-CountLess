
import csv
from datetime import datetime
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from experiment import CountMin, CountLess, make_hashes


# =====================================
# 1. 压力测试参数
# =====================================

MEMORY_BYTES = 8192

FLOW_COUNTS = [1000, 2000, 5000, 10000, 20000]
SEEDS = [42, 43, 44, 45, 46]

METHODS = ["CM", "CM-CU", "CL-CM", "CL-CU", "CL-MU"]

# 每次运行创建独立结果文件夹，避免覆盖原结果
PROJECT_ROOT = Path(__file__).resolve().parents[1]
OUT = PROJECT_ROOT / "results" / "03_stress_test" / (
    "stress_" + datetime.now().strftime("%Y%m%d_%H%M%S")
)
OUT.mkdir(parents=True, exist_ok=True)


# =====================================
# 2. 生成固定流大小组成的数据
# =====================================

def generate_trace(n_flows, seed):
    rng = np.random.default_rng(seed)

    n1 = int(n_flows * 0.80)
    n10 = int(n_flows * 0.15)
    n100 = int(n_flows * 0.04)
    n1000 = n_flows - n1 - n10 - n100

    # truth[i] = 第 i 条 Flow 的真实包数
    truth = np.array(
        [1] * n1 +
        [10] * n10 +
        [100] * n100 +
        [1000] * n1000,
        dtype=np.int64
    )

    # 随机分配 Flow Size，避免大流集中在小 ID
    rng.shuffle(truth)

    # 按真实包数构造逐包到达序列
    packets = np.repeat(
        np.arange(n_flows, dtype=np.int32),
        truth
    )

    # 随机打乱包到达顺序
    rng.shuffle(packets)

    return packets, truth


# =====================================
# 3. 评价指标
# =====================================

def evaluate(sketch, truth):
    n = len(truth)
    estimates = np.array(
        [sketch.query(i) for i in range(n)], dtype=np.int64
    )
    if np.any(estimates < truth):
        raise AssertionError("发现低估，请检查饱和与查询逻辑")

    errors = np.abs(estimates - truth)
    relative = errors / truth
    mouse = truth <= 255
    elephant = truth > 255

    result = {
        "AAE": float(errors.mean()),
        "ARE": float(relative.mean()),
        "Mouse_FSR": float(100 * np.mean(relative[mouse] < 0.10)),
        "Elephant_FSR": float(100 * np.mean(relative[elephant] < 0.01)),
    }

    for size in [1, 10, 100, 1000]:
        selected = truth == size
        if not np.any(selected):
            result[f"ARE_size_{size}"] = float("nan")
            result[f"FSR_size_{size}"] = float("nan")
            continue
        threshold = 0.10 if size <= 255 else 0.01
        result[f"ARE_size_{size}"] = float(relative[selected].mean())
        result[f"FSR_size_{size}"] = float(
            100 * np.mean(relative[selected] < threshold)
        )

    return result


# =====================================
# 4. 底层哈希碰撞比例
# =====================================

def layer1_collision_rate(sketch, n_flows):
    positions = [
        sketch.positions(i)[0]
        for i in range(n_flows)
    ]

    counts = np.bincount(
        positions,
        minlength=sketch.widths[0]
    )

    # 统计映射到非独占位置的 Flow 比例
    colliding_flows = counts[counts > 1].sum()

    return 100 * colliding_flows / n_flows


# =====================================
# 5. 实验主体
# =====================================

def run():
    rows = []

    for n_flows in FLOW_COUNTS:
        print(f"\n======= {n_flows} Flows =======")

        for seed in SEEDS:
            packets, truth = generate_trace(n_flows, seed)
            hashes = make_hashes(seed + 2026)

            sketches = {
                "CM": CountMin(MEMORY_BYTES, hashes),
                "CM-CU": CountMin(
                    MEMORY_BYTES, hashes, conservative=True
                ),
                "CL-CM": CountLess(
                    MEMORY_BYTES, hashes, "CL-CM"
                ),
                "CL-CU": CountLess(
                    MEMORY_BYTES, hashes, "CL-CU"
                ),
                "CL-MU": CountLess(
                    MEMORY_BYTES, hashes, "CL-MU"
                ),
            }

            collision_rate = layer1_collision_rate(
                sketches["CL-MU"], n_flows
            )

            for method, sketch in sketches.items():
                assert sketch.used_bits <= MEMORY_BYTES * 8

                for flow in packets:
                    sketch.update(int(flow))

                metrics = evaluate(sketch, truth)

                rows.append({
                    "n_flows": n_flows,
                    "n_packets": len(packets),
                    "seed": seed,
                    "method": method,
                    **metrics,
                    "L1_collision_pct": (
                        collision_rate
                        if method.startswith("CL-")
                        else float("nan")
                    ),
                    "L1_saturation_pct": (
                        sketch.saturation_l1()
                        if isinstance(sketch, CountLess)
                        else float("nan")
                    ),
                })

            print(
                f"Seed={seed}  "
                f"Packets={len(packets)}  "
                f"L1 collision={collision_rate:.2f}%"
            )

    # 保存逐次实验结果
    with (OUT / "stress_results.csv").open(
        "w", newline="", encoding="utf-8-sig"
    ) as f:
        writer = csv.DictWriter(
            f, fieldnames=list(rows[0].keys())
        )
        writer.writeheader()
        writer.writerows(rows)

    return rows


# =====================================
# 6. 汇总五个随机种子的均值与标准差
# =====================================

def summarize(rows):
    summary = []

    for n in FLOW_COUNTS:
        for method in METHODS:
            group = [
                r for r in rows
                if r["n_flows"] == n
                and r["method"] == method
            ]

            row = {
                "n_flows": n,
                "method": method,
            }

            metrics = [
                "ARE", "AAE", "Mouse_FSR", "Elephant_FSR",
                "L1_collision_pct", "L1_saturation_pct",
                "ARE_size_1", "ARE_size_10", "ARE_size_100", "ARE_size_1000",
                "FSR_size_1", "FSR_size_10", "FSR_size_100", "FSR_size_1000",
            ]
            for metric in metrics:
                values = np.array(
                    [r[metric] for r in group],
                    dtype=float
                )

                # CM 的 Layer 1 指标为 NaN
                if np.all(np.isnan(values)):
                    row[metric + "_mean"] = float("nan")
                    row[metric + "_std"] = float("nan")
                else:
                    row[metric + "_mean"] = float(
                        np.nanmean(values)
                    )
                    row[metric + "_std"] = float(
                        np.nanstd(values, ddof=1)
                    )

            summary.append(row)

    with (OUT / "stress_summary.csv").open(
        "w", newline="", encoding="utf-8-sig"
    ) as f:
        writer = csv.DictWriter(
            f, fieldnames=list(summary[0].keys())
        )
        writer.writeheader()
        writer.writerows(summary)

    return summary


# =====================================
# 7. 生成趋势图
# =====================================

def plot_summary(summary):
    for metric in ["ARE", "Mouse_FSR", "Elephant_FSR"]:
        plt.figure(figsize=(8, 5))

        for method in METHODS:
            selected = [
                r for r in summary
                if r["method"] == method
            ]

            x = [r["n_flows"] for r in selected]
            y = [r[metric + "_mean"] for r in selected]

            plt.plot(x, y, marker="o", label=method)

        plt.xlabel("Number of distinct flows")
        plt.ylabel(metric)
        plt.title(f"Count-Less Stress Test: {metric}")

        if metric == "ARE":
            plt.yscale("symlog", linthresh=0.01)
        else:
            plt.ylim(0, 105)

        plt.grid(alpha=0.25)
        plt.legend()
        plt.tight_layout()

        plt.savefig(
            OUT / f"stress_{metric}.png",
            dpi=180
        )
        plt.close()


def plot_mouse_groups(summary):
    selected = [r for r in summary if r["method"] == "CL-MU"]
    plt.figure(figsize=(9, 5))
    for size in [1, 10, 100]:
        x = [r["n_flows"] for r in selected]
        y = [r[f"FSR_size_{size}_mean"] for r in selected]
        plt.plot(x, y, marker="o", label=f"{size}-packet flows")
    plt.xlabel("Number of distinct flows")
    plt.ylabel("FSR (%)")
    plt.title("CL-MU: Mouse FSR by Flow Size")
    plt.ylim(0, 105)
    plt.grid(alpha=0.25)
    plt.legend()
    plt.tight_layout()
    plt.savefig(OUT / "mouse_group_fsr.png", dpi=180)
    plt.close()


# =====================================
# 8. 主程序
# =====================================

if __name__ == "__main__":
    rows = run()
    summary = summarize(rows)
    plot_summary(summary)
    plot_mouse_groups(summary)

    print("\n实验完成，结果目录：", OUT)

    print("\nCL-MU Summary:")
    for row in summary:
        if row["method"] == "CL-MU":
            print(
                f"Flows={row['n_flows']:5d}, "
                f"ARE={row['ARE_mean']:.4f}, "
                f"MouseFSR={row['Mouse_FSR_mean']:.2f}%, "
                f"ElephantFSR={row['Elephant_FSR_mean']:.2f}%"
            )
