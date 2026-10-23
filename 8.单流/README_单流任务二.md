# 单流任务二：400 MHz CPU/UE 分侧仿真

本目录是原仓库的独立副本，新增的单流任务二入口位于根目录：

```matlab
run_single_stream_task2                    % smoke：6 个固定样本，检查全流程并生成全部图
run_single_stream_task2('Mode','full')     % 400 个冻结 UE 位置的完整仿真
```

也可指定新的结果目录，避免覆盖已有结果：

```matlab
run_single_stream_task2('Mode','full','OutputDir','results/single_stream_full');
```

仿真固定为**单流**，吞吐量按 **400 MHz** 换算。CPU 与 UE 的结果写入不同的 `CPU/`、`UE/` 子目录；CPU 侧的部分反馈方案仅使用 R=3 个参考 AP 的方向或位置消息，因此信道估计/消息维度不随协作 AP 数 B 增长。UE 侧以本地观测 CSI 计算共相位，并给出“全反馈-完美 CSI（上界）”以展示完美信道信息的性能增益。

吞吐量与反馈开销图的基线包含：

- 全反馈-Type-I（UE-CSI）：每 AP 选择一条过采样 DFT 波束，并进行 QPSK AP 间相位对齐。
- 全反馈-Type-II（UE-CSI）：每 AP 选择两条 DFT 波束，反馈 3 bit 幅度和 2 bit 内部相位后进行 QPSK AP 间相位对齐。
- Docomo 近场 LoS（UE-CSI，仅 UE 图）：本地 4 倍过采样 DFT、QPSK `beta`、AP 逻辑栅格 DFT 和三轮交替更新；它采用完整 UE CSI 选码，故不与 CPU 几何反馈链路直接比较。
- 部分反馈-DFT-index（R=3）
- 部分反馈-HAQSE（R=3）
- 部分反馈-位置反馈（R=3）
- 全反馈-完美 CSI（SVD 性能上界）

Docomo 的原始 `simulate_proposed` 对 4 个 AP 全枚举 `4^(B-1)` 个 QPSK `beta` 组合。为适配本任务的 5–10 个协作 AP，本副本保留其“`beta` / AP-DFT / 本地波束”交替优化结构，但将该步骤改为固定首 AP 的 QPSK 坐标更新。每次更新的反馈比特数由本地波束索引、`beta` 相位和 AP-DFT 索引共同统计。

统一绘图函数为 `plot_single_stream_task2.m`。所有图均无标题、使用中文图例和较大字体，覆盖：定位误差与 SNR/天线数/AP 数；吞吐量与 SNR/天线数/莱斯因子；以及比特数与 AP 数。定位误差图仅显示能够输出位置估计的三种部分反馈方案；Type-I、Type-II 与 Docomo 均不虚构定位误差曲线。

输出中的 `single_stream_trials.csv` 是逐样本数据，`single_stream_summary.csv` 为分组均值，`experiment_contract.json` 记录仿真口径。
