# 位置反馈扩展结果

本轮CPU与UE统一采用相位敏感度坐标bit分配，正式400个共享位置，原七方法加POS R3/R5。

- `位置反馈仿真结果报告.md`：中文结果、设置、诊断和实验图。
- `selected_B8_OBS30.csv`：8 AP、观测30 dB的全部方法及相对Type-I的反馈节省和吞吐量差异。
- `equal_bit_comparison.csv`：相同实际消息bit下的位置反馈与HAQSE方向反馈比较。
- `one_percent_thresholds_all.csv`：相对未量化位置参考的1%平均吞吐量损失门槛，仅在已扫描预算内判断。
- `weighted_localization_summary.csv`：观测能量加权定位的独立消融，不改变主实验定位器。
- `runtime_benchmark.csv`：本机热缓存、单计算线程、31次重复的分项计算耗时。
- `complete_audit.json`：全量记录、功率、信息边界、CPU/UE波束及共享基线审计。
- `source_and_input_SHA256.json`：实现、公共输入和原有结果文件的来源校验。

`ULA/`、`UPA/`各含配置、分块、位置主实验和补充实验CSV。其`CPU/`、`UE/`子目录提供合并原七方法后的9方法完整结果，每组180配置、72000条展示记录，以及PNG和MATLAB FIG。原数据按面板复用，共享基线不能视为独立样本。

`preflight/`保存回归核查，`smoke_ULA/`、`smoke_UPA/`为三个原共享样本的小规模接口核查，`pilot_ULA/`、`pilot_UPA/`为各24个独立新位置的预实验。它们不计入正式400点统计。

ULA主坐标24 bit按12/12分配，UPA主坐标36 bit按13/14/9分配；UE链路另加4(B−1) bit相位。所有主图吞吐量采用原固定rho及20 MHz带宽，未扣除反馈资源；主BF图sqrt(mean(G))与预算功率BF图mean(G)分别标注。

完整运行结束并通过审计后，才生成`complete_audit.json`及最终中文报告。复现入口见上两级目录的`位置反馈扩展说明.md`。
