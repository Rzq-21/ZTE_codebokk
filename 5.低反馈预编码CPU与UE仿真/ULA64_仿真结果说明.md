# ULA 64天线：CPU端与UE端低反馈仿真结果

完成日期：2026-10-17。结果遵循共享仿真协议 v2.0，与UPA结果使用同一批400个UE位置、散射系数和逐AP噪声种子。

## 实验设置

- AP：64阵元ULA，半波长间距，7 GHz；DFT码本512个方向（8倍过采样），每个AP索引9 bit。
- UE：4阵元ULA。UE高度固定为0 m，位置在50 m×50 m场地内的20–30 m方形区域随机取400点。
- 部署：冻结的`ULA64_G50_H6_ORIGINAL`，每边两台核心AP，最小水平中心间距6 m；活动AP顺序为`[1,2,3,6,8,4,5,7,9,10]`。
- 方向反馈：3参考AP或5参考AP。HAQSE每AP反馈14 bit圆周量化空间频率；DFT每AP反馈一个9 bit码本索引。其余AP方向由CPU按位置估计结果推断。
- 共相位：CPU链路用定位结果构造几何LOS信道并连续计算相位，不反馈相位；UE链路用带噪实测CSI计算相对共相位，并对其余`B−1`个AP各反馈4 bit。
- 信道：逐阵元精确球面波LOS；LOS+NLOS为LOS与6条共享单跳散射路径的Rician混合，K=10 dB。
- 性能：等AP功率、总发射功率固定。观测SNR为20–40 dB；吞吐量按固定10 dB数据工作点、20 MHz带宽换算。每点取400样本均值，不加置信区间。
- 主实验方法：Ideal SVD、Type-II、Type-I、HAQSE/DFT各自R=3和R=5；协作AP数B=5–10。

## 主要数值：B=8、观测SNR=30 dB

吞吐量单位为Mbps；BF增益为`sqrt(mean(G))`。Type-I/II是CPU与UE两份结果共享的基线，二者逐样本完全一致。

|信道|方法|CPU吞吐量|UE吞吐量|CPU/UE反馈bit|CPU/UE位置RMSE|CPU/UE全AP对相位RMSE|
|---|---|---:|---:|---:|---:|---:|
|LOS|Type-I|100.64|100.64|100 / 100|—|8.61° / 8.61°|
|LOS|HAQSE R=3|100.34|100.65|42 / 70|2.25 mm / 2.25 mm|19.71° / 9.15°|
|LOS|DFT R=3|89.42|100.62|27 / 55|4.03 cm / 4.03 cm|90.09° / 9.10°|
|LOS+NLOS|Type-I|97.88|97.88|100 / 100|—|8.57° / 8.57°|
|LOS+NLOS|HAQSE R=3|94.98|98.09|42 / 70|8.37 mm / 8.37 mm|63.23° / 8.89°|
|LOS+NLOS|DFT R=3|86.90|98.07|27 / 55|4.05 cm / 4.05 cm|90.56° / 8.91°|

在这组设置下，UE端DFT R=3以55 bit总反馈（比Type-I的100 bit少45 bit）达到与Type-I几乎相同的吞吐量。HAQSE R=3 UE链路为70 bit，也接近Type-I性能。两种方法中，R=3方向码本及相位反馈的总开销均明显低于Type-I。

CPU端和UE端的低反馈方案使用相同方向消息，所以位置RMSE和实际方向波束完全一致；性能差异主要来自相位来源。LOS下DFT R=3的位置RMSE约4 cm，CPU据此重构的相对相位RMSE约90°，UE利用实测CSI后约9°。因此CPU DFT吞吐量明显下降是当前几何相位恢复链路的实际敏感性表现，不能把它误读为UE端DFT方向码本本身效果差。HAQSE方向更精细，CPU LOS结果更接近UE端；但在LOS+NLOS下，CPU的几何LOS模型没有还原散射造成的残余相位，故其相位RMSE和吞吐量仍显著落后于UE链路。

## 输出文件

- CPU正式结果：`results/full400_ULA_CPU_20261017/`
- UE正式结果：`results/full400_ULA_UE_20261017/`
- 两端在B=8、SNR=30 dB的对照表：`results/ULA_CPU_UE_B8_OBS30_comparison.csv`
- 完整性审计：`results/ULA_full_audit.json`。检查结果为两端各56,000条记录、140种配置；共享基线接收功率最大差为0；CPU/UE低反馈方向码、波束索引和位置差异均为0；两端各生成4幅正式图。
- CPU图：吞吐量-SNR、吞吐量-AP数、BF增益-SNR、反馈开销-AP数。
- UE图：同样四幅；文件名前缀为`E1_UE`、`E2_UE`、`E3_UE`。

## 复现入口

在MATLAB R2024a中切换到项目根目录后运行：

```matlab
verify_ula;                 % 271项单样本、几何和反馈口径核验
run_ula_cpu('Mode','full'); % 400位置CPU相位链路
run_ula_ue('Mode','full');  % 同一400位置UE相位链路
audit_ula_results;          % 全量结果一致性与B=8对照表
```

正式结果已保存在指定目录；重新运行时请另选`OutputDir`，以免覆盖已完成数据。
