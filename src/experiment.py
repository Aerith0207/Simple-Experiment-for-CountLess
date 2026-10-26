
import csv
import time
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from collections import Counter



# ======================================
# 1. 实验参数
# ======================================

N_PACKETS = 100_000       # 数据包总数
N_FLOWS = 5_000           # 不同 Flow 的数量

MEMORY_BYTES = 65536
ZIPFS = [2.0]
SEEDS = [42, 43, 44, 45, 46]

PRIME = 2_147_483_647
PROJECT_ROOT = Path(__file__).resolve().parents[1]
OUTPUT_DIR = PROJECT_ROOT / "results" / "01_core_experiment" / "generated"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


# ======================================
# 2. 可复现的哈希函数
# ======================================

MASK64 = (1 << 64) - 1


def make_hashes(seed):
    rng = np.random.default_rng(seed)

    # 三层使用独立随机 salt
    return [
        int(x)
        for x in rng.integers(
            0, 2**64,
            size=3,
            dtype=np.uint64
        )
    ]


def hash_index(flow, params, width):
    # SplitMix64：将连续 Flow ID 混合成分散的哈希值
    x = (int(flow) ^ int(params)) & MASK64

    x = (x + 0x9E3779B97F4A7C15) & MASK64
    x = ((x ^ (x >> 30)) * 0xBF58476D1CE4E5B9) & MASK64
    x = ((x ^ (x >> 27)) * 0x94D049BB133111EB) & MASK64
    x = x ^ (x >> 31)

    return x % width

def test_hash_distribution():
    print("\n=== Hash Distribution Check ===")

    for seed in SEEDS:
        hashes = make_hashes(seed + 2026)

        sketch = CountLess(
            MEMORY_BYTES,
            hashes,
            "CL-MU"
        )

        positions = [
            sketch.positions(flow)[0]
            for flow in range(N_FLOWS)
        ]

        counts = Counter(positions)

        unique_positions = len(counts)

        isolated = sum(
            counts[pos] == 1
            for pos in positions
        )

        print(
            f"Seed={seed}: "
            f"Occupied={unique_positions}, "
            f"Isolated={isolated / N_FLOWS:.2%}"
        )


# ======================================
# 3. 标准 Count-Min 和 CM-CU
# ======================================

class CountMin:
    def __init__(self, memory_bytes, hashes, conservative=False):
        self.hashes = hashes
        self.conservative = conservative
        self.width = memory_bytes * 8 // (3 * 32)

        assert self.width > 0
        self.reg = [[0] * self.width for _ in range(3)]
        self.used_bits = 3 * 32 * self.width

    def positions(self, flow):
        return [
            hash_index(flow, self.hashes[j], self.width)
            for j in range(3)
        ]

    def update(self, flow):
        pos = self.positions(flow)
        values = [
            self.reg[j][pos[j]]
            for j in range(3)
        ]

        if self.conservative:
            minimum = min(values)
            for j in range(3):
                if values[j] == minimum:
                    self.reg[j][pos[j]] += 1 #最小值加1
        else:
            for j in range(3):
                self.reg[j][pos[j]] += 1

    def query(self, flow):
        pos = self.positions(flow)
        return min(
            self.reg[j][pos[j]]
            for j in range(3)
        )


# ======================================
# 4. Count-Less：CL-CM / CL-CU / CL-MU
# ======================================

class CountLess:
    def __init__(self, memory_bytes, hashes, mode="CL-MU"):
        self.hashes = hashes
        self.mode = mode

        # d=3, r=4
        # 内存：16w*8 + 4w*16 + w*32 = 224w bit
        w = memory_bytes * 8 // 224

        assert w > 0
        self.widths = [16 * w, 4 * w, w]
        self.limits = [255, 65535, 2**32 - 1]

        self.reg = [
            [0] * self.widths[j]
            for j in range(3)
        ]
        self.used_bits = 224 * w

    def positions(self, flow):
        return [
            hash_index(flow, self.hashes[j], self.widths[j])
            for j in range(3)
        ]

    def update(self, flow):
        pos = self.positions(flow)

        # CL-CM：对全部未饱和层更新
        if self.mode == "CL-CM":
            for j in range(3):
                if self.reg[j][pos[j]] < self.limits[j]:
                    self.reg[j][pos[j]] += 1

        # CL-CU：预先寻找全局最小值
        elif self.mode == "CL-CU":
            valid = [
                (j, self.reg[j][pos[j]])
                for j in range(3)
                if self.reg[j][pos[j]] < self.limits[j]
            ]

            if valid:
                minimum = min(value for _, value in valid)
                for j, value in valid:
                    if value == minimum:
                        self.reg[j][pos[j]] += 1

        # CL-MU：按层顺序更新临时最小值
        elif self.mode == "CL-MU":
            valmin = float("inf")

            for j in range(3):
                value = self.reg[j][pos[j]]

                # 饱和层不可继续更新
                if value >= self.limits[j]:
                    continue

                if value < valmin:
                    value += 1
                    self.reg[j][pos[j]] = value
                    valmin = value

        else:
            raise ValueError("Unknown update mode")

    def query(self, flow):
        pos = self.positions(flow)

        # 查询时排除已经饱和的层
        valid = [
            self.reg[j][pos[j]]
            for j in range(3)
            if self.reg[j][pos[j]] < self.limits[j]
        ]

        return min(valid) if valid else self.limits[-1]

    def saturation_l1(self):
        values = self.reg[0]
        return 100 * sum(v == 255 for v in values) / len(values)


# ======================================
# 5. 正确性测试
# ======================================

def unit_tests():
    hashes = make_hashes(123)
    flow = 777

    def check_mu(before, expected):
        sketch = CountLess(8192, hashes, "CL-MU")
        pos = sketch.positions(flow)

        for j in range(3):
            sketch.reg[j][pos[j]] = before[j]

        sketch.update(flow)

        actual = tuple(
            sketch.reg[j][pos[j]]
            for j in range(3)
        )

        assert actual == expected, (actual, expected)

    check_mu((15, 18, 12), (16, 18, 13))
    check_mu((105, 100, 80), (106, 101, 81))
    check_mu((255, 7, 9), (255, 8, 9))

    cu = CountLess(8192, hashes, "CL-CU")
    pos = cu.positions(flow)

    for j, value in enumerate((15, 18, 12)):
        cu.reg[j][pos[j]] = value

    cu.update(flow)

    assert tuple(
        cu.reg[j][pos[j]] for j in range(3)
    ) == (15, 18, 13)

    # 单一 Flow 且未饱和，应精确计数
    mu = CountLess(8192, hashes, "CL-MU")
    for _ in range(20):
        mu.update(flow)

    assert mu.query(flow) == 20

    print("PASS: 所有基础单元测试通过")


# ======================================
# 6. 生成 Zipf 偏斜流量
# ======================================

def make_trace(alpha, seed):
    rng = np.random.default_rng(seed)

    ranks = np.arange(1, N_FLOWS + 1, dtype=float)
    probs = ranks ** (-alpha)
    probs /= probs.sum()

    # 每条 Flow 至少出现一次，用于固定 Cardinality
    truth = np.ones(N_FLOWS, dtype=np.int64)

    truth += rng.multinomial(
        N_PACKETS - N_FLOWS, probs
    )

    packets = np.repeat(
        np.arange(N_FLOWS, dtype=np.int32),
        truth
    )

    rng.shuffle(packets)

    return packets, truth


# ======================================
# 7. 误差指标
# ======================================

def evaluate(sketch, truth):
    estimates = np.array(
        [sketch.query(i) for i in range(N_FLOWS)],
        dtype=np.int64
    )
    assert np.all(estimates >= truth), \
    "ERROR: 存在 Flow 被低估，请检查算法实现"

    absolute_error = np.abs(estimates - truth)
    relative_error = absolute_error / truth

    mouse = truth <= 255
    elephant = truth > 255

    return {
        "AAE": float(np.mean(absolute_error)),
        "ARE": float(np.mean(relative_error)),
        "Mouse_FSR_pct": float(
            100 * np.mean(relative_error[mouse] < 0.10)
        ) if mouse.any() else float("nan"),
        "Elephant_FSR_pct": float(
            100 * np.mean(relative_error[elephant] < 0.01)
        ) if elephant.any() else float("nan"),
        "Elephant_flows": int(elephant.sum()),
        "Elephant_packet_share": float(
            truth[elephant].sum() / truth.sum()
        ),
    }


# ======================================
# 8. 运行实验
# ======================================

def run_experiments():
    results = []

    for alpha in ZIPFS:
        for seed in SEEDS:
            packets, truth = make_trace(alpha, seed)
            mouse = truth <= 255
            elephant = truth > 255

            print(f"Total Packets: {truth.sum()}")
            print(f"Total Flows: {len(truth)}")
            print(f"Single-packet Flows: {(truth == 1).sum()}")
            print(f"Mouse Flows: {mouse.sum()}")
            print(f"Elephant Flows: {elephant.sum()}")
            print(f"Max Flow Size: {truth.max()}")
            hashes = make_hashes(seed + 2026)

            print(f"\nZipf={alpha}, Seed={seed}")

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

            for name, sketch in sketches.items():
                assert sketch.used_bits <= MEMORY_BYTES * 8

                start = time.perf_counter()

                for flow in packets:
                    sketch.update(int(flow))

                elapsed = time.perf_counter() - start
                metrics = evaluate(sketch, truth)

                saturation = (
                    sketch.saturation_l1()
                    if isinstance(sketch, CountLess)
                    else float("nan")
                )

                row = {
                    "zipf": alpha,
                    "seed": seed,
                    "method": name,
                    "memory_bits": sketch.used_bits,
                    **metrics,
                    "L1_saturated_pct": saturation,
                    "throughput_Mpps": N_PACKETS / elapsed / 1e6,
                }

                results.append(row)

                print(
                    f"{name:6s} "
                    f"AAE={row['AAE']:.3f} "
                    f"ARE={row['ARE']:.3f} "
                    f"MouseFSR={row['Mouse_FSR_pct']:.2f}% "
                    f"ElephantFSR={row['Elephant_FSR_pct']:.2f}%"
                )

    output_csv = OUTPUT_DIR / (
    f"results_{MEMORY_BYTES // 1024}KB.csv"
    )

    with output_csv.open(
        "w", newline="", encoding="utf-8-sig"
    ) as file:
        writer = csv.DictWriter(
            file, fieldnames=list(results[0].keys())
        )
        writer.writeheader()
        writer.writerows(results)

    return results


# ======================================
# 9. 画图
# ======================================

def plot_results(results):
    methods = ["CM", "CM-CU", "CL-CM", "CL-CU", "CL-MU"]

    for metric in ["ARE", "AAE"]:
        plt.figure(figsize=(8, 5))

        for method in methods:
            ys = []

            for alpha in ZIPFS:
                values = [
                    r[metric]
                    for r in results
                    if r["method"] == method
                    and r["zipf"] == alpha
                ]
                ys.append(np.mean(values))

            plt.plot(
                ZIPFS, ys, marker="o", label=method
            )

        plt.xlabel("Zipf skewness")
        plt.ylabel(metric + " (lower is better)")
        plt.title("Frequency Estimation: " + metric)
        plt.yscale("symlog", linthresh=1)
        plt.xticks(ZIPFS)
        plt.legend()
        plt.grid(alpha=0.25)
        plt.tight_layout()

        filepath = OUTPUT_DIR / f"{metric.lower()}_vs_zipf.png"
        plt.savefig(filepath, dpi=180)
        plt.close()

    print(f"\n实验结果已保存到 {OUTPUT_DIR}")


# ======================================
# 10. 主程序
# ======================================

if __name__ == "__main__":
    unit_tests()
    test_hash_distribution()
    experiment_results = run_experiments()
    plot_results(experiment_results)

