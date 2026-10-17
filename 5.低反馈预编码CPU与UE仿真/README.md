# 低反馈预编码正式仿真：UPA 与 ULA 的 CPU/UE 链路

本目录收录完整正式实验：32×32 UPA 和 64 阵元 ULA 两种 AP 阵列；每种阵列均包含 CPU 端几何重构共相位、UE 端测量 CSI 计算并反馈共相位两条链路。四组结果均为 400 个共享 UE 位置、每组 140 种设置、每组 56,000 条逐样本记录，含 LOS 与 LOS+NLOS 两类信道。

实验总结及必要公式见[正式仿真总结报告](正式仿真总结报告.md)。UPA 的链路实现约定见[实现细节](实现细节.md)，ULA 的独立结果说明见[ULA64 仿真结果说明](ULA64_仿真结果说明.md)。

## 四组正式结果

|阵列|共相位端|正式结果目录|核查|
|---|---|---|---|
|UPA 32×32|CPU|`results/full400_UPA_CPU_20261017/`|目录内 `delivery_audit.json`|
|UPA 32×32|UE|`results/full400_UPA_UE_20261017/`|目录内 `delivery_audit.json`|
|ULA 64|CPU|`results/full400_ULA_CPU_20261017/`|`results/ULA_full_audit.json`|
|ULA 64|UE|`results/full400_ULA_UE_20261017/`|`results/ULA_full_audit.json`|

每组目录含逐样本 CSV、汇总 CSV、反馈预算、运行元数据、日志，以及吞吐量、波束成形增益和反馈开销的 PNG/MATLAB FIG 图。UPA 和 ULA 的 CPU/UE 对照表分别为 `results/UPA_CPU_UE_B8_OBS30_comparison.csv` 和 `results/ULA_CPU_UE_B8_OBS30_comparison.csv`。

## 复现入口

使用 MATLAB R2024a，在本目录运行对应入口。正式结果已随分支提供；如需重跑，请另设新的 `OutputDir`，避免覆盖已保存的正式数据。

```matlab
verify_upa_cpu;
verify_upa_ue;
verify_ula;

run_upa_cpu('Mode','full');
run_upa_ue('Mode','full');
run_ula_cpu('Mode','full');
run_ula_ue('Mode','full');

audit_ula_results;
```

四组实验共享协议 v2.0 的输入文件与随机样本。CPU 接口只接收公共部署配置和已量化方向消息；UE 链路仅额外反馈各协作 AP 相对参考 AP 的 4 bit 共相位，不把完整 CSI 发给 CPU。
