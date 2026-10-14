# 频率扰动与位置估计的关系

本文件夹是《2.对比方案》的物理模型对照组，不修改原有代码。

## 运行

在 MATLAB 中运行：

```matlab
cd('E:\github\codebook_ZTE_task2\3.频率扰动与位置估计的关系（物理模型对照组）');
frequency_perturbation_position_error
```

脚本会生成：

- `frequency_position_error_results.mat`
- `frequency_error_vs_position_error_scatter.png`
- `frequency_error_vs_position_error_trend.png`
- `per_ap_dft_frequency_error.png`

## 对照逻辑

1. 完全沿用《2.对比方案》的 AP 位置、UE 位置、AP 朝向、16×16 UPA、用户 4 天线、7 GHz、4 倍 DFT 过采样和参考 AP 1:3。
2. 不添加信道噪声，也不重新估计信道。
3. 从真实几何关系计算参考 AP 的真实二维空间频率 `f=[f_y,f_z]`。
4. 用无噪声 LOS 信道和原对比代码的 DFT 相关性搜索得到实际 DFT 索引；同时计算最近 DFT 网格点作为诊断基线。
5. 图 1 在真实频率上加入指定 RMS 的随机扰动；扰动后的频率仍然先经过 DFT-index 网格量化，再按照 DFT-index 的方向反演和最小二乘射线交汇求位置。
6. 图 2 不加入人为扰动，而是固定同一个真实空间频率，分别使用 1、2、4、8、16 倍 DFT 过采样进行最近网格量化，比较不同量化分辨率下的频率 RMS 误差和位置估计偏差；同时显示连续真实频率完美估计的零误差基准。
7. 图 1 的横轴使用 DFT-index 量化后的有效空间频率 RMS 误差，图 2 的横轴使用真实频率量化 RMS 误差，纵轴均为用户位置欧氏误差，单位分别为 cycles/element 和 m。原始人为加入的扰动误差也保存在 `results.input_frequency_error_rms` 中。

图 2 中每个方向的量化比特数为 `log2(16*O)`：O=1、2、4、8、16 分别对应 4、5、6、7、8 bit；当前 O=4 的点与无噪声 DFT-index 基线相对应。

脚本保留了原 `solve_position_from_dft_index.m` 的索引拆分方式：线性 DFT 索引的内层索引解释为局部 y 频率，外层索引解释为局部 z 频率。原物理模型中两个索引跨度都为 64，因此该约定不会引入额外参数变化，但可以直接暴露原 DFT-index 索引映射本身造成的位置误差。

如果无噪声 DFT 网格基线已经产生较大的位置误差，说明 DFT-index 曲线较差主要来自频率网格量化和位置误差传播；如果基线位置误差很小，则应继续排查原对比代码中的预编码、相位配置或信道重构环节。
