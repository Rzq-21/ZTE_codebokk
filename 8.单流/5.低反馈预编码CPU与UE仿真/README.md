# 低反馈预编码正式仿真：UPA 与 ULA 的 CPU/UE 链路

本目录保留协议 v2.0 的七方法正式实验：32×32 UPA 和 64 阵元 ULA 两种 AP 阵列，每种均包含 CPU 几何共相位和 UE 测量共相位两条链路。原四组结果各为 400 位置、140 配置、56,000 条记录，覆盖 LOS 与 LOS+NLOS。

[正式仿真总结报告](正式仿真总结报告.md)已更新为包含 POS R3/R5 的九方法主实验，每组180配置、72000条展示记录。新增方法的实现、测试和九方法完整结果统一见 [6.位置反馈方案](../6.位置反馈方案/README.md)；本页下面的入口与表格仍对应原七方法版本。

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
