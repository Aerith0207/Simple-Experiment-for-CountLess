# Count-Less 软件复现实验

本项目在 Python 环境中复现 Count-Less 的核心计数与更新机制，重点研究 Pyramid 结构、Minimum Update（MU）与 Conservative Update（CU）的精度差异。这是软件级机制验证，不是 Count-Less 论文全部实验的严格复现。

## 项目导航

| 路径 | 内容 |
| --- | --- |
| `src/experiment.py` | CM、CM-CU、CL-CM、CL-CU、CL-MU 的核心实现；Zipf 与内存实验 |
| `src/stress_test.py` | 固定 8 KB 下的流数增长压力测试和 Flow Size 分组 |
| `src/rescue_diagnostic.py` | 100-packet Flow 的逐层 FSR 和 Cross-layer Rescue 诊断 |
| `docs/实验设计与结果.md` | 可检索的实验设计、数据与研究边界 |
| `docs/实验总结.docx` | 正式 Word 版实验报告 |
| `docs/README.md` | 文档目录与阅读顺序 |
| `results/README.md` | 结果目录、正式数据与历史调试数据的区分 |

## 实验主线

1. 验证五种算法的基本实现，并检查哈希随机性。
2. 在 Zipf 合成流量下对比基础精度。
3. 改变内存预算，观察 CL-MU 与 CL-CU 的差距。
4. 固定 8 KB 增加整体测量负载，寻找精度退化场景。
5. 按 1、10、100、1000-packet Flow 分组，定位总体 ARE 与 FSR 的来源。
6. 针对 100-packet Flow 分析 CU/MU 差异和跨层挽救机制。
7. 下一步进行 Update Noise Decomposition，量化 MU 相对 CU 的额外更新。

## 快速运行

在项目根目录中执行：

```powershell
.\.venv\Scripts\python.exe src\experiment.py
.\.venv\Scripts\python.exe src\stress_test.py
.\.venv\Scripts\python.exe src\rescue_diagnostic.py
```

如需重建环境：

```powershell
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
```

新生成结果会自动进入对应的 `results` 子目录，不依赖当前命令行所在路径。

## 当前最重要的结论

- 在当前合成流量下，ARE 稳定表现为 `CL-CU < CL-MU < CL-CM < CM-CU < CM`。
- CL-MU 随内存增加明显改善；8–64 KB 下 ARE 从 0.7521 降至 0.0152。
- 20,000 Flows 时，单包流约贡献 CL-MU 总体 ARE 的 98.2%，总体指标不能代替 Flow Size 分组。
- 100-packet Flow 上，20,000 Flows 时 CL-CU 与 CL-MU 的 Final FSR 分别为 91.70% 和 64.63%。
- CL-MU 的跨层挽救率更高主要是因为第一层失败更多，不能据此认定其挽救机制更优。

## 解读边界

- 数据为自生成合成流量，尚未使用真实 CAIDA 轨迹。
- 尚未实现真正的 FCM 和 Elastic 基线。
- CL-CU 是软件理想参考，不代表可直接满足单遍流水线约束。
- 压力测试中 Flow 数与总包数同时增加，因此属于整体负载增长，不能单独归因于 Cardinality。

