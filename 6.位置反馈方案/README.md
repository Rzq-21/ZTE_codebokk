# 位置反馈方案

本目录集中整理位置反馈 POS R3/R5 的实现、测试、原始结果和复现材料。CPU、UE 两条共相位链路统一按相对 AP1 路径差敏感度分配坐标 bit；主实验 ULA 为24 bit，UPA 为36 bit。原正式报告中的四组主图和结果表已并入本方法，见[九方法主实验总结](../5.低反馈预编码CPU与UE仿真/正式仿真总结报告.md)。

## 从哪里开始

- [算法与运行说明](位置反馈扩展说明.md)：消息格式、定位、分配与 CPU/UE 处理流程。
- [完整位置反馈实验报告](results/position_phase_v1_20261018/位置反馈仿真结果报告.md)：主实验、预算扫描、未量化参考、加权定位变体、诊断和耗时。
- [结果文件索引](results/position_phase_v1_20261018/README.md)：逐样本数据、汇总、PNG、MATLAB FIG 及审计。
- [完整结果审计](results/position_phase_v1_20261018/complete_audit.json)：每点400样本、CPU/UE位置与波束一致、功率和原七方法回归检查。

## 实现和测试文件

|位置|用途|
|---|---|
|`+posfb/config.m`、`allocate.m`|公共坐标域、固定设计网格、相位敏感度与整数 bit 分配|
|`+posfb/encode.m`、`decode.m`|UE坐标编码与CPU纯整数消息解码；解码边界不接受真实位置或观测CSI|
|`+posfb/prediction.m`|由解码位置生成全部协作AP的最终DFT波束|
|`+posfb/cpu_phase.m`、`ue_phase.m`、`phases.m`|CPU几何模型相位、UE观测相位及配对仿真|
|`+posfb/weighted_position.m`|独立的观测能量加权定位变体，不替换主实验无权定位|
|`+posfb/result_row.m`、`columns.m`、`summarize.m`|新增结果字段、逐样本记录和统计|
|`verify_position_extension.m`|旧存档回归、坐标消息边界、量化误差界、功率与两分支一致性|
|`run_position_extension.m`|smoke、独立pilot、完整400位置实验，以及按20样本分块保存|
|`finalize_position_extension.m`|复用旧七方法数据，合并POS R3/R5并输出九方法主图和预算图|
|`benchmark_position_extension.m`|预热后31次分项本机耗时测量|
|`analyze_position_extension.py`|全量审计、实际等bit对照、门槛统计、汇总与报告生成|

`+upa/`、`+ula/`、`lib/`和`shared/`是从原第5目录复制的冻结依赖和输入，用于保证本目录可以单独运行。`results/full400_{ULA,UPA}_{CPU,UE}_20261017/results_trials.csv`为原七方法逐样本对照快照，仅供回归和合并；来源仍是第5目录的正式数据。

## 已完成的实验

主实验覆盖 ULA/UPA、CPU/UE、LOS/LOS+NLOS、R3/R5；观测 SNR 为20–40 dB，协作AP数为5–10。两种阵列均使用原400个位置、散射系数和噪声种子。

每个阵列/相位组包含原七方法加POS R3/R5，共180配置、72000条展示记录；四组合计新增位置主实验64000条。主图包括吞吐量随观测SNR、吞吐量随AP数、BF gain随观测SNR和反馈开销随AP数，共16张。补充实验156800条，另有8张预算图；24张图均附可编辑FIG。

独立pilot每阵列24个新位置；预算扫描与方向精度扫描固定B=8、观测30 dB。加权定位、未量化位置参考和计算耗时均独立标识，不混入固定主预算的九方法主图。主BF图为`sqrt(mean(G))`，预算功率图为`mean(G)`。

## 独立复现

将整个本目录复制到新位置，使用 MATLAB R2024a 进入该目录。所有 MATLAB 计算依赖、公共输入和旧对照CSV均在本目录内；不需要把第5目录加入搜索路径。为重跑选择新的输出目录。

```matlab
verify_position_extension('results/new_run/preflight');
run_position_extension('ULA','Mode','smoke','OutputDir','results/new_run/smoke_ULA');
run_position_extension('UPA','Mode','smoke','OutputDir','results/new_run/smoke_UPA');
run_position_extension('ULA','Mode','pilot','OutputDir','results/new_run/pilot_ULA');
run_position_extension('UPA','Mode','pilot','OutputDir','results/new_run/pilot_UPA');
run_position_extension('ULA','Mode','full','OutputDir','results/new_run/ULA');
run_position_extension('UPA','Mode','full','OutputDir','results/new_run/UPA');
finalize_position_extension('ULA','results/new_run/ULA');
finalize_position_extension('UPA','results/new_run/UPA');
benchmark_position_extension('results/new_run');
```

随后在本目录运行，Python需安装numpy、pandas：

```text
python analyze_position_extension.py results/new_run
```

已保存的完整结果位于`results/position_phase_v1_20261018/`。smoke结果用于与后续完整结果对齐，因此新输出目录中应先生成对应smoke数据，再运行全量分析。

## 来源与口径

目录整理不重新抽样、不改变信道和预编码实现。`source_and_input_SHA256_原始执行.json`保留初次运行校验；`source_and_input_SHA256.json`记录当前独立目录的源码与输入校验；`目录整理核查.json`记录冻结依赖、对照数据和原始结果的一致性检查。

本目录通过`.gitattributes`保留文件原始字节，避免跨平台检出时的换行转换影响SHA256核查。

坐标、方向与相位均按上行净消息字段计费；协议头、CRC、训练和CPU至AP下发不在计费范围内。吞吐量采用原固定rho及20 MHz带宽，未扣除反馈资源。Type-I/II在CPU面板中仍使用UE测量相位；混合信道是原K=10 dB的LOS主导六径场景。
