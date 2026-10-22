"""Independent CSV/message audit, paired comparisons and Chinese report. No plots."""
from pathlib import Path
import hashlib
import json
import math
import os
import sys
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
OUT = ROOT / 'results' / 'dual_v1_20261022'
BASE = ROOT.parent / '6.位置反馈方案'
KEYS = ['array', 'stage', 'scenario', 'channel', 'B', 'observation_SNR_dB', 'R', 'streams', 'method', 'qA', 'qP']


def freeze_sources():
    files = list(ROOT.glob('*.m')) + list(ROOT.glob('*.py')) + list(ROOT.glob('+dmimo/*.m'))
    for pattern in ['+posfb/*.m', '+upa/*.m', '+ula/*.m', 'lib/*.m', 'shared/*']:
        files.extend(BASE.glob(pattern))
    entries = [dict(path=os.path.relpath(p, ROOT).replace('\\', '/'), bytes=p.stat().st_size,
                    sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in sorted(files) if p.is_file()]
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / 'source_and_input_SHA256.json').write_text(json.dumps(dict(version='DUAL_V1_20261022', files=entries), ensure_ascii=False, indent=2), encoding='utf-8')
    print('Frozen source/input hashes:', len(entries))


def mdtable(frame, columns, labels=None):
    labels = labels or columns
    rows = ['|' + '|'.join(labels) + '|', '|' + '|'.join(['---'] * len(columns)) + '|']
    for values in frame[columns].itertuples(index=False, name=None):
        out = []
        for v in values:
            if isinstance(v, (float, np.floating)):
                out.append('—' if not np.isfinite(v) else f'{v:.3f}')
            else:
                out.append(str(v))
        rows.append('|' + '|'.join(out) + '|')
    return '\n'.join(rows)


def paired(frame, fixed, left_filter, right_filter, merge_keys, label):
    left = frame.loc[left_filter].copy()
    right = frame.loc[right_filter].copy()
    t = left.merge(right, on=merge_keys, suffixes=('_left', '_right'), validate='one_to_one')
    assert len(t) == len(left), (label, len(t), len(left))
    records = []
    for key, g in t.groupby(fixed, dropna=False):
        if not isinstance(key, tuple):
            key = (key,)
        d = g.MMSE_Mbps_left.to_numpy() - g.MMSE_Mbps_right.to_numpy()
        se = float(np.std(d, ddof=1) / np.sqrt(len(d)))
        r = dict(zip(fixed, key))
        r.update(comparison=label, samples=len(g), left_mean_Mbps=g.MMSE_Mbps_left.mean(),
                 right_mean_Mbps=g.MMSE_Mbps_right.mean(), delta_Mbps=d.mean(),
                 delta_CI95_low=d.mean() - 1.96 * se, delta_CI95_high=d.mean() + 1.96 * se,
                 percent_of_right=100 * d.mean() / g.MMSE_Mbps_right.mean(), positive_fraction=float(np.mean(d > 0)))
        records.append(r)
    return pd.DataFrame(records)


def audit_one(array):
    t = pd.read_csv(OUT / array / 'trials.csv')
    assert len(t) == 148800 and not t.duplicated(KEYS + ['sample_id']).any()
    sizes = t.groupby(KEYS, dropna=False).size()
    assert sizes.eq(400).all() and set(t.sample_id) == set(range(1, 401))
    real = ~t.method.str.startswith('CAP_')
    assert np.max(np.abs(t.loc[real, 'precoder_norm_squared'] - 1)) < 1e-10
    assert t.loc[real, 'per_AP_power_error'].max() < 1e-10
    assert (t.MMSE_SE <= t.joint_SE + 1e-8).all()
    assert (t.joint_SE <= t.rank_capacity + 1e-7).all()
    assert np.max(np.abs(t.loc[t.streams.eq(1), 'MMSE_SE'] - t.loc[t.streams.eq(1), 'joint_SE'])) < 1e-8
    finite_ue = real & t.qA.gt(0)
    expected_coeff = 2 * t.B * t.streams * (t.qA + t.qP) - t.streams * t.qP
    assert t.loc[finite_ue, 'coefficient_bits'].eq(expected_coeff[finite_ue]).all()
    assert t.loc[finite_ue, 'feedback_bits'].eq(t.loc[finite_ue, 'direction_bits'] + expected_coeff[finite_ue]).all()
    # Replay every transmitted coefficient message independently of MATLAB.
    max_norm_error = 0.0
    count = 0
    for row in t.loc[finite_ue].itertuples(index=False):
        a = np.asarray(json.loads(row.coefficient_amplitude_indices))
        if row.streams == 1:
            a = a.reshape(-1, 1)
        p = np.asarray(json.loads(row.coefficient_phase_indices)).ravel()
        assert a.shape == (2 * row.B, row.streams) and p.size == 2 * row.B * row.streams - row.streams
        assert np.all(a == np.round(a)) and np.all((a >= 0) & (a < 2 ** row.qA))
        assert np.all(p == np.round(p)) and np.all((p >= 0) & (p < 2 ** row.qP))
        amp = a / (2 ** row.qA - 1)
        assert np.all(np.linalg.norm(amp.reshape(row.B, 2, row.streams), axis=(1, 2)) > 0)
        phase = np.zeros((2 * row.B, row.streams))
        # MATLAB logical-index order is column-major.
        mask = np.ones(phase.shape, dtype=bool)
        mask[np.arange(row.streams), np.arange(row.streams)] = False
        flat = phase.ravel(order='F')
        flat[mask.ravel(order='F')] = p
        phase = flat.reshape(phase.shape, order='F')
        C = amp * np.exp(2j * np.pi * phase / 2 ** row.qP)
        C = C.reshape(row.B, 2, row.streams)
        C /= np.linalg.norm(C, axis=(1, 2))[:, None, None] * np.sqrt(row.B)
        max_norm_error = max(max_norm_error, abs(np.sum(np.abs(C) ** 2) - 1))
        count += 1
    # Coordinate decoding, allocations and CPU/UE identical position and beams.
    cfg = json.loads((OUT / array / 'experiment_contract.json').read_text(encoding='utf-8'))['position']
    bits = np.asarray(cfg['main_allocation']).ravel()
    lo = np.asarray(cfg['lo']).ravel()
    hi = np.asarray(cfg['hi']).ravel()
    coords = t[t.method.eq('POS_CPU')]
    for row in coords.itertuples(index=False):
        indices = np.asarray(json.loads(row.coordinate_indices)).ravel()
        pos = np.asarray(json.loads(row.decoded_position)).ravel()
        assert np.all((indices >= 0) & (indices < 2 ** bits))
        assert np.max(np.abs(pos[:len(bits)] - (lo + (indices + .5) * (hi - lo) / 2 ** bits))) < 1e-10
    p = t[t.method.isin(['POS_CPU', 'POS_UE_Q']) & (t.stage != 'COEFFICIENT')]
    group = ['stage', 'scenario', 'channel', 'sample_id', 'B', 'observation_SNR_dB', 'R', 'streams']
    assert p.groupby(group).decoded_position.nunique().eq(1).all()
    assert p.groupby(group).beam_indices.nunique().eq(1).all()
    audit = dict(array=array, records=len(t), case_count=len(sizes), samples_per_case=400,
                 coefficient_messages_replayed=count, coefficient_replay_max_norm_error=max_norm_error,
                 coordinate_messages_replayed=len(coords), position_and_beam_equality=True,
                 overflow_records=int(t.quantizer_overflow.sum()),
                 maximum_joint_capacity_excess=float((t.joint_SE - t.rank_capacity)[real].max()))
    return t, audit


def main():
    tables, audits = [], []
    for array in ['ULA', 'UPA']:
        t, a = audit_one(array)
        tables.append(t)
        audits.append(a)
        print('AUDITED', array, a['records'], flush=True)
    t = pd.concat(tables, ignore_index=True)
    main_t = t[t.stage.eq('MAIN')]
    merge = ['array', 'scenario', 'channel', 'B', 'observation_SNR_dB', 'R', 'method', 'qA', 'qP', 'sample_id']
    fixed = ['array', 'channel', 'B', 'observation_SNR_dB', 'R', 'method', 'qA', 'qP']
    gains = paired(main_t, fixed, main_t.streams.eq(2), main_t.streams.eq(1), merge, 'two_stream_minus_one')
    gains.to_csv(OUT / 'dual_stream_gain_paired.csv', index=False)
    # Repeat full direction comparator for each R, maintaining one-to-one sample pairing.
    compare = []
    for R in [3, 5]:
        p = main_t[(main_t.R.eq(R) & main_t.method.eq('POS_UE_Q')) | main_t.method.eq('FULLDIR_UE_Q')].copy()
        p['comparison_R'] = R
        merge2 = ['array', 'scenario', 'channel', 'B', 'observation_SNR_dB', 'streams', 'qA', 'qP', 'sample_id', 'comparison_R']
        compare.append(paired(p, ['array', 'channel', 'B', 'observation_SNR_dB', 'streams', 'comparison_R'],
                              p.method.eq('POS_UE_Q'), p.method.eq('FULLDIR_UE_Q'), merge2, 'position_minus_full_direction'))
    comparison = pd.concat(compare, ignore_index=True)
    comparison.to_csv(OUT / 'UE_vs_full_direction_paired.csv', index=False)
    mismatch = t[(t.stage.eq('MISMATCH') | t.stage.eq('MAIN')) & t.B.eq(8) & t.observation_SNR_dB.eq(30)]
    effects = []
    for scenario in ['tilt20', 'phase15', 'combined']:
        p = mismatch[mismatch.scenario.isin(['nominal', scenario])]
        merge3 = ['array', 'channel', 'R', 'streams', 'method', 'qA', 'qP', 'sample_id']
        p = p[p.method.ne('POS_CPU_CALIBRATED')]
        result = paired(p, ['array', 'channel', 'R', 'streams', 'method', 'qA', 'qP'],
                        p.scenario.eq(scenario), p.scenario.eq('nominal'), merge3, scenario + '_minus_nominal')
        effects.append(result)
    pd.concat(effects).to_csv(OUT / 'mismatch_paired.csv', index=False)
    # Independent pilot coverage and recorded runtime assertions.
    for array in ['ULA', 'UPA']:
        a = json.loads((OUT / f'pilot_verified_{array}' / 'audit.json').read_text())
        assert a['samples'] == 24 and a['records'] == 4896 and a['status'] == 'complete'
    manifest = json.loads((OUT / 'source_and_input_SHA256.json').read_text(encoding='utf-8'))
    for entry in manifest['files']:
        file = ROOT / entry['path']
        assert hashlib.sha256(file.read_bytes()).hexdigest() == entry['sha256'], file
    complete = dict(status='passed', version='DUAL_V1_20261022', main_and_supplement_records=len(t),
                    legacy_rank_records=3200, independent_pilot_records=9792, arrays=audits,
                    frozen_source_and_input_files_verified=len(manifest['files']),
                    paired_CI='normal approximation mean +/- 1.96*std(paired delta)/sqrt(400); conditional on frozen topology and polar model',
                    net_feedback_resource_subtraction=False)
    (OUT / 'complete_audit.json').write_text(json.dumps(complete, ensure_ascii=False, indent=2), encoding='utf-8')

    selected = main_t[main_t.B.eq(8) & main_t.observation_SNR_dB.eq(30)]
    s = selected.groupby(['array', 'channel', 'R', 'streams', 'method'], as_index=False).agg(
        MMSE_Mbps=('MMSE_Mbps', 'mean'), joint_Mbps=('joint_Mbps', 'mean'), feedback_bits=('feedback_bits', 'mean'))
    s.to_csv(OUT / 'selected_B8_OBS30.csv', index=False)
    rank = pd.read_csv(OUT / 'legacy_rank_trials.csv')
    rank['equal_R2_gain_Mbps'] = 20 * (rank.beam_R2_equal_SE - rank.beam_R1_SE)
    r = rank.groupby(['array', 'channel', 'R'], as_index=False).agg(
        sigma2_ratio_median=('full_sigma2_ratio', 'median'), energy_retention=('beam_energy_retention', 'mean'),
        ideal_two_stream_gain_Mbps=('equal_R2_gain_Mbps', 'mean'))
    r.to_csv(OUT / 'legacy_rank_summary.csv', index=False)
    gs = gains[gains.B.eq(8) & gains.observation_SNR_dB.eq(30) & gains.R.eq(3) & gains.method.isin(['POS_CPU', 'POS_UE_Q'])]
    cs = comparison[comparison.B.eq(8) & comparison.observation_SNR_dB.eq(30) & comparison.streams.eq(2) & comparison.comparison_R.eq(3)]
    mis_s = mismatch[(mismatch.R.eq(3)) & (mismatch.streams.eq(2)) & mismatch.method.isin(['POS_CPU', 'POS_UE_Q', 'POS_CPU_CALIBRATED'])]
    mis_s = mis_s.groupby(['array', 'channel', 'scenario', 'method'], as_index=False).agg(MMSE_Mbps=('MMSE_Mbps', 'mean'))
    diag = pd.read_csv(OUT / 'cpu_oracle_diagnostics.csv')
    assert len(diag) == 9600
    assert diag.groupby(['array', 'channel', 'streams', 'variant']).size().eq(400).all()
    daudit = json.loads((OUT / 'cpu_diagnostic_audit.json').read_text())
    assert daudit['status'] == 'passed' and daudit['maximum_main_regression_Mbps'] < 1e-8
    dsummary = diag.groupby(['array', 'channel', 'streams', 'variant'], as_index=False).agg(MMSE_Mbps=('MMSE_Mbps', 'mean'))
    dsummary.to_csv(OUT / 'cpu_oracle_summary.csv', index=False)
    complete['cpu_diagnostic_records'] = len(diag)
    complete['cpu_main_regression_maximum_Mbps'] = daudit['maximum_main_regression_Mbps']
    (OUT / 'complete_audit.json').write_text(json.dumps(complete, ensure_ascii=False, indent=2), encoding='utf-8')
    text = ['# 双极化与双流位置反馈仿真报告', '', '完成日期：2026-10-22。四步实验均已执行；完整模型与消息协议见[目录说明](../../README.md)。', '',
            '## 1. 先检查原单极化信道的多流潜力', '',
            '原4接收空间位置、8协作AP的信道并非整体秩1。下表使用原400位置、原单极化真信道和存档POS解码位置；第二模态比为σ₂²/σ₁²。位置波束子空间的等功率理想双流参考在全部样本上高于单流。该参考允许在波束子空间内自由分配AP功率且不计反馈，不能直接视作可实施双流收益。', '',
            mdtable(r, list(r.columns)), '',
            '## 2. 双极化单流控制与3. UE双流实现', '',
            '每个原空间阵元改为两个±45°极化端口，UE为4位置×2极化。以下单/双流采用相同硬件、总功率1和每AP功率1/B。位置仍生成所有空间波束；多流系数由2×Ns矩阵承担，不能把旧单个AP共相位简单复制到两个流。', '',
            'B=8、观测30dB、R3定点表（全方向对照和容量上界为R0）：', '',
            mdtable(s[s.R.isin([0, 3])], list(s.columns)), '',
            '双流减单流的逐样本配对结果，单位Mbps；95%区间为配对差值的正态近似均值区间：', '',
            mdtable(gs, ['array', 'channel', 'method', 'delta_Mbps', 'delta_CI95_low', 'delta_CI95_high', 'positive_fraction']), '',
            '双流POS UE R3与全方向UE采用相同系数量化的配对差值：', '',
            mdtable(cs, ['array', 'channel', 'left_mean_Mbps', 'right_mean_Mbps', 'delta_Mbps', 'delta_CI95_low', 'delta_CI95_high', 'percent_of_right']), '',
            '系数幅度3bit、相位4bit。B=8双流系数216bit；ULA的POS总反馈240bit、全方向288bit，节省16.67%；UPA为252/328bit，节省23.17%。CPU只反馈24/36bit位置，但其系数依赖公开极化与LOS模型。以上均未扣除导频和反馈资源占用，不是净吞吐量。连续系数没有有限反馈bit；全方向对照为本项目构造，不等同完整NR Type-I。', '',
            '## 4. CPU已知极化与失配试验', '',
            '固定B=8、观测30dB，分别将UE绕y轴倾斜20°、加入逐物理AP第二极化端口RMS15°相位偏差、叠加二者。UE通过有噪观测重新设计系数，CPU默认仍使用标称极化。CALIBRATED对照额外知道真实倾角/端口校准值，但共享同一解码位置与波束；未计该额外信息成本，作为诊断。', '',
            mdtable(mis_s, list(mis_s.columns)), '',
            '相位失配是发射端口的列酉变换，本身不改变真实信道容量；若CPU性能下降而UE保持稳定，原因是系数设计没有匹配实际端口相位。倾角失配也可能改变真实信道的可用空间模式，速率不保证随倾角单调下降，不能把其全部变化归因于算法误差。逐样本配对数据见`mismatch_paired.csv`。', '',
            '### CPU标称极化下残差的追加诊断', '',
            'ULA在标称LOS下仍出现CPU与UE差距，追加R3定点实验：保持主实验估计位置生成的波束不变，只替换系数设计所用的信息。CPU_estimated_p重算主实验；ORACLE_true_p_same_beams用真实位置但仍是平面LOS模型；ORACLE_true_effective_H_same_beams用真实有效信道。后两项均没有可实施有限反馈预算，也不是一般最优上界。', '',
            mdtable(dsummary, list(dsummary.columns)), '',
            '## 曲线与量化开销', '']
    for array in ['ULA', 'UPA']:
        text += [f'### {array}', '', f'![观测SNR]({array}/rate_SNR.png)', '', f'![AP数]({array}/rate_B.png)', '',
                 f'![极化失配]({array}/polarization_mismatch.png)', '', f'![反馈开销]({array}/feedback_rate.png)', '']
    text += ['## 审计、范围与后续判断', '',
             '两阵列各148800条正式/补充记录，共297600条；另有原信道秩诊断3200条、独立24位置pilot共9792条及CPU残差诊断9600条。正式实验每配置400位置。54项预检查及全量独立CSV审计通过，整数系数消息逐条解码复核，CPU/UE位置与波束一致，无坐标域溢出。总功率、每AP功率、MMSE≤联合解码≤对应秩容量上界均通过核查。', '',
             '该结果支持在受控双极化模型下扩展位置反馈到双流；位置承担公共空间波束信息，UE系数承担极化混合、跨AP相干关系与流间干扰控制。CPU方案是否足够，应结合上表在目标极化校准条件下的实际差距判断。继续降低UE双流开销，需要结构化的极化/子空间系数码本；不能沿用旧单流52/64bit来声称同样的压缩率。', '',
             '统计结论条件于冻结AP部署、6个散射体及一次冻结极化耦合/端口误差实例。这里没有覆盖完整38.901、实测天线方向图、互耦、宽带频率选择性、多用户调度、移动跟踪、独立部署重复和极化误差强度扫描，不能据此给出一般系统性能保证。容量上界放松了每AP功率约束；可实施系数为SVD初始化后逐AP归一的启发式方案。', '',
             '文件：`legacy_rank_trials.csv`、`ULA/UPA/trials.csv`为逐样本数据；`dual_stream_gain_paired.csv`、`UE_vs_full_direction_paired.csv`、`mismatch_paired.csv`为配对比较；`complete_audit.json`与`source_and_input_SHA256.json`为审计和复现校验。每张PNG附同名MATLAB FIG。']
    report = '\n'.join(text) + '\n'
    (OUT / '双极化多流仿真报告.md').write_text(report, encoding='utf-8')
    # Append a compact result section to the original main report, preserving earlier experiments.
    original = ROOT.parent / '5.低反馈预编码CPU与UE仿真' / '正式仿真总结报告.md'
    marker = '\n## 10. 双极化与双流扩展（2026-10-22）'
    old = original.read_text(encoding='utf-8').split(marker)[0].rstrip()
    compact = s[(s.R.eq(3)) & s.method.isin(['POS_CPU', 'POS_UE_Q'])].copy()
    append = marker + '\n\n原单极化信道秩诊断、双极化单流控制、UE双流实现及CPU极化失配测试已完成。两个阵列沿用原400位置、原AP空间阵元；每个空间位置增加两个极化端口，UE为8端口。下表比较同一双极化硬件、同一总功率1、每AP功率1/B下的单/双流，B=8、观测30dB、R3。不能与前文单极化绝对速率直接相减，解释为多流收益。\n\n'
    append += mdtable(compact, list(compact.columns))
    append += '\n\n双流UE反馈系数为216bit，再加ULA 24bit或UPA 36bit坐标，总量240/252bit；同系数量化的全方向对照为288/328bit，节省16.67%/23.17%。这与前文单流52/64bit属于不同消息格式。CPU位置消息仍为24/36bit，其性能依赖极化、姿态和端口相位模型。\n\n完整四阶段结果、失配/校准诊断、SNR/AP曲线、量化扫描及全量审计见[双极化多流报告](../7.双极化多流扩展/results/dual_v1_20261022/双极化多流仿真报告.md)，代码与复现见[第7目录](../7.双极化多流扩展/README.md)。新增信道为受控Jones极化扩展，不声称完整3GPP或实测验证；主速率采用线性MMSE，另存联合解码速率。\n'
    original.write_text(old + '\n' + append, encoding='utf-8')
    print('COMPLETE: 297600 trial records audited; report and original report updated.', flush=True)


if __name__ == '__main__':
    if '--freeze' in sys.argv:
        freeze_sources()
    else:
        main()
