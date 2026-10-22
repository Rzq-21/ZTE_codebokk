# 全方法双流双极化最小验证

固定B=8、观测SNR=30dB。`confirm`使用原冻结样本库前32个位置，`smoke`使用399/400号位置，二者不重叠。不同方法和流数共享信道与标准噪声。两种阵列、两种信道，每个位置16种设计和2种流数，确认4096条、预检查256条。没有进行新一轮参数扫描。

## 阅读报告

- [Markdown扩展讲解](双流双极化全方法最小验证与扩展说明.md)
- [Word扩展讲解](双流双极化全方法最小验证与扩展说明.docx)

两份报告由同一组内容和CSV生成，讲解信道维度、空间波束、2×Ns系数矩阵、各方法消息、Type-II双波束复合、功率与反馈、结果及适用边界。Type-I/II沿用仓库简化定义，不声称完整标准码本符合性。Type-II不是独立4维波束子空间扩展。

## 数据与核查

|文件|内容|
|---|---|
|`confirm/trials.csv`|4096条逐样本结果，包含实际整数空间和系数消息|
|`confirm/summary.csv`|每配置32位置的均值、最弱接收模式和反馈bit|
|`smoke/trials.csv`|256条预检查记录，附与既有POS存档的回归|
|`paired_trials.csv`|同位置同方法的单流与双流配对速率差|
|`paired_summary.csv`|配对均值、正增量个数、最差增量和实际bit|
|`legacy_reduction_audit.json`|Type-I/II在单模式单极化退化条件下的8项一致性检查|
|`complete_audit.json`|独立整数解码、功率、接收秩、容量关系及POS回归|
|`source_and_input_SHA256.json`|本次计算、分析和报告源码与冻结输入的字节校验|
|`document_QA.json`|最终Word结构、页数与逐页视觉检查记录|

SVD等功率参考和注水容量放松每AP功率约束；容量行不对应指定发射矩阵，其逐流SINR和接收模式字段为空。所有可实施方法每AP功率1/B，总功率1。UE消息不重复计算旧共相位，流数预先约定，不含RI或协议头开销。反馈未从吞吐量中扣除。

## 复现

MATLAB R2024a进入第7目录，保留同仓库第6目录冻结依赖及既有`+dmimo`模块：

```matlab
verify_minimal_all_methods;
run_minimal_all_methods('smoke');
run_minimal_all_methods('confirm');
```

然后用含numpy/pandas/python-docx/Pillow的Python环境运行：

```text
python analyze_minimal_all_methods.py
python build_minimal_reports.py
```

当前runner对固定最小输出重新计算，不改变既有400位置数据。smoke的POS回归需要`results/dual_v1_20261022/{ULA,UPA}/trials.csv`。Word已在Windows Word 16中导出并逐页检查；页面图和中间PDF仅用于内部排版检查，不作为发布附件。
