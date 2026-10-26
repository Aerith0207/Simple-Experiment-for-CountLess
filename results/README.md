# 结果索引

## 01_core_experiment

- `per_seed/results_8KB_splitmix64.csv`：8 KB、Zipf 1.2/2.0/3.0、5 种子的正式 SplitMix64 结果。
- `per_seed/results_16KB.csv`、`results_32KB.csv`、`results_64KB.csv`：Zipf=2.0 的内存扩展实验。
- `figures/latest_64KB_*_vs_zipf.png`：最后一次 64 KB、Zipf=2.0 单点运行生成的图，不建议作为 Zipf 趋势图使用。
- `generated/`：今后运行 `src/experiment.py` 的输出位置。

## 02_memory_sensitivity

- `memory_summary.csv`：8/16/32/64 KB 的五种子均值与样本标准差。
- `memory_AAE.png`、`memory_ARE.png`、`memory_Mouse_FSR.png`：内存敏感性汇总图。

这一目录是组会中“基础验证：不同内存下的 ARE”的首选数据源。

## 03_stress_test

- `20261025_172617_initial/`：初版压力测试，包含总体 ARE、Mouse FSR、Elephant FSR。
- `20261025_181944_flow_size_breakdown/`：完整版，额外包含 1/10/100/1000-packet Flow 的分组 ARE 和 FSR。

如果只选一组进行分析，优先使用 `20261025_181944_flow_size_breakdown/`。

## 04_rescue_analysis

- `20261025_232226/rescue_summary.csv`：CL-CU/CL-MU 的 L1/L2/L3 FSR、Final FSR、Rescue Rate 和 ARE_100。
- `rescue_details.csv`：逐种子数据。
- `cu_good_mu_bad_examples.csv`：CU 达标而 MU 失败的具体 Flow 案例。
- `Final_FSR_pct_comparison.png`、`Rescue_pct_comparison.png`：可用于汇报的汇总图。

## 99_legacy

- `results_old_linear_hash.csv`：使用旧线性哈希的历史调试结果。仅用于说明“哈希方案导致种子敏感”，不应与正式 SplitMix64 结果混用。

