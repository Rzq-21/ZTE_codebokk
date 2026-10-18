"""Audit the completed MATLAB extension and write a source-grounded Chinese report."""
from pathlib import Path
import argparse
import hashlib
import json
import platform
import numpy as np
import pandas as pd


def markdown(frame, formats=None):
    formats = formats or {}
    def cell(name, value):
        if pd.isna(value):
            return '—'
        if name in formats:
            return format(value, formats[name])
        return str(value).replace('|', '/')
    lines = ['|' + '|'.join(map(str, frame.columns)) + '|', '|' + '|'.join(['---'] * len(frame.columns)) + '|']
    lines += ['|' + '|'.join(cell(k, v) for k, v in row.items()) + '|' for row in frame.to_dict('records')]
    return '\n'.join(lines)


def sha256(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def analyze(base):
    root = Path(__file__).resolve().parent
    base = base.resolve()
    selected = []
    thresholds = []
    audit = {'arrays': {}, 'complete': False}
    weighted = []
    equal_bits = []
    for array in ['ULA', 'UPA']:
        directory = base / array
        new = pd.read_csv(directory / 'position_main_trials.csv')
        supp = pd.read_csv(directory / 'supplement_trials.csv')
        assert len(new) == 32000
        smoke = pd.read_csv(base / f'smoke_{array}' / 'position_main_trials.csv').set_index(['case_id', 'sample_id']).sort_index()
        aligned_smoke = new.set_index(['case_id', 'sample_id']).loc[smoke.index]
        smoke_difference = float(np.max(np.abs(smoke.received_power - aligned_smoke.received_power)))
        assert smoke_difference < 1e-12
        expected_supp = 76800 if array == 'ULA' else 80000
        assert len(supp) == expected_supp
        for frame in [new, supp]:
            assert not frame.duplicated(['case_id', 'sample_id']).any()
            assert np.isfinite(frame.received_power).all()
            assert np.isfinite(frame.throughput_Mbps).all()
            assert frame.power_gain.between(0, 1 + 1e-8).all()
            assert np.max(np.abs(frame.precoder_norm_squared - 1)) < 1e-10
            assert (frame.groupby('case_id').size() == 400).all()
            keys = ['experiment_id', 'method_id', 'B', 'observation_SNR_dB', 'channel_model', 'coordinate_budget', 'q_frequency', 'sample_id']
            # NaN values are stable strings in this join; continuous references are explicit.
            left = frame[frame.panel_branch.eq('CPU')].copy()
            right = frame[frame.panel_branch.eq('UE')].copy()
            def join_key(data):
                return data[keys].fillna('NA').astype(str).agg('|'.join, axis=1)
            left.index = join_key(left)
            right.index = join_key(right)
            assert left.index.is_unique and right.index.is_unique
            right = right.loc[left.index]
            assert (left.actual_beam_indices == right.actual_beam_indices).all()
            assert np.array_equal(left[['estimated_x_m', 'estimated_y_m', 'estimated_z_m']].to_numpy(), right[['estimated_x_m', 'estimated_y_m', 'estimated_z_m']].to_numpy())
            finite = np.isfinite(left.feedback_bits.to_numpy())
            assert np.allclose(right.feedback_bits.to_numpy()[finite] - left.feedback_bits.to_numpy()[finite], 4 * (left.B.to_numpy()[finite] - 1))
        summaries = {}
        main_by_branch = {}
        for branch in ['CPU', 'UE']:
            t = pd.read_csv(directory / branch / 'results_trials.csv')
            s = pd.read_csv(directory / branch / 'results_summary.csv')
            assert len(t) == 72000 and len(s) == 180
            assert (s.samples == 400).all() and s.method_id.nunique() == 9
            old = pd.read_csv(root / 'results' / f'full400_{array}_{branch}_20261017' / 'results_trials.csv')
            archived = t[~t.method_id.str.startswith('POS')].set_index(['case_id', 'sample_id']).sort_index()
            old = old.set_index(['case_id', 'sample_id']).sort_index()
            for col in ['received_power', 'power_gain', 'throughput_Mbps', 'feedback_bits']:
                assert np.allclose(archived[col], old[col], rtol=1e-12, atol=1e-12, equal_nan=True)
            main_by_branch[branch] = t
            summaries[branch] = s
            pick = s[(s.B == 8) & (s.observation_SNR_dB == 30)].copy()
            pick['array_group'] = array
            pick['branch'] = branch
            losses = pd.read_csv(directory / branch / 'paired_vs_TypeI.csv')
            pick = pick.merge(losses[['case_id', 'feedback_saving_pct', 'throughput_loss_pct', 'paired_delta_CI95_low_Mbps', 'paired_delta_CI95_high_Mbps']], on='case_id', how='left', validate='one_to_one')
            selected.append(pick)
        shared_keys = ['sample_id', 'channel_model', 'B', 'observation_SNR_dB', 'method_id']
        c = main_by_branch['CPU'].query("method_id in ['SVD','TypeI','TypeII']").set_index(shared_keys).sort_index()
        u = main_by_branch['UE'].query("method_id in ['SVD','TypeI','TypeII']").set_index(shared_keys).sort_index()
        shared_diff = float(np.max(np.abs(c.received_power - u.received_power)))
        assert shared_diff == 0
        th = pd.read_csv(directory / 'one_percent_thresholds.csv')
        th.insert(0, 'array', array)
        thresholds.append(th)
        ss = pd.read_csv(directory / 'supplement_summary.csv')
        position_sweep = ss[ss.experiment_id.eq('BUDGET') & ss.method_id.isin(['POS_R3', 'POS_R5'])].copy()
        frequency_sweep = ss[ss.experiment_id.eq('FREQUENCY')].copy()
        matched = position_sweep.merge(frequency_sweep, on=['channel_model', 'panel_branch', 'R', 'feedback_bits'], suffixes=('_position', '_frequency'), validate='one_to_one')
        eq = matched[['channel_model', 'panel_branch', 'R', 'feedback_bits', 'coordinate_budget_position', 'q_frequency_frequency', 'mean_throughput_Mbps_position', 'mean_throughput_Mbps_frequency']].copy()
        eq['array'] = array
        eq['position_throughput_improvement_pct'] = 100 * (eq.mean_throughput_Mbps_position / eq.mean_throughput_Mbps_frequency - 1)
        equal_bits.append(eq)
        w = ss[ss.experiment_id.eq('WEIGHTED')].copy()
        w.insert(0, 'array_group', array)
        weighted.append(w)
        audit['arrays'][array] = {
            'main_new_records': len(new), 'supplement_records': len(supp),
            'merged_records_per_branch': 72000, 'cases_per_branch': 180,
            'CPU_UE_coordinate_and_beam_equality': True,
            'archived_seven_methods_preserved': True,
            'shared_baseline_power_difference': shared_diff,
            'smoke_to_final_power_difference': smoke_difference,
            'main_overflows': int(new.quantizer_overflow.sum()),
            'supplement_overflows': int(supp.quantizer_overflow.sum()),
            'maximum_power_gain': float(max(new.power_gain.max(), supp.power_gain.max())),
        }
    selected = pd.concat(selected, ignore_index=True)
    thresholds = pd.concat(thresholds, ignore_index=True)
    weighted = pd.concat(weighted, ignore_index=True)
    equal_bits = pd.concat(equal_bits, ignore_index=True)
    selected.to_csv(base / 'selected_B8_OBS30.csv', index=False, encoding='utf-8-sig')
    thresholds.to_csv(base / 'one_percent_thresholds_all.csv', index=False, encoding='utf-8-sig')
    weighted.to_csv(base / 'weighted_localization_summary.csv', index=False, encoding='utf-8-sig')
    equal_bits.to_csv(base / 'equal_bit_comparison.csv', index=False, encoding='utf-8-sig')
    # Freeze the provenance of inputs, old performance data and final simulator source.
    paths = list((root / '+posfb').glob('*.m')) + [root / x for x in ['run_position_extension.m', 'verify_position_extension.m', 'finalize_position_extension.m', 'benchmark_position_extension.m', 'analyze_position_extension.py']]
    paths += [x for x in (root / 'shared').iterdir() if x.is_file()]
    paths += list((root / '+upa').glob('*.m')) + list((root / '+ula').glob('*.m')) + list((root / 'lib').glob('*.m'))
    for array in ['ULA','UPA']:
        for branch in ['CPU','UE']:
            paths.append(root / 'results' / f'full400_{array}_{branch}_20261017' / 'results_trials.csv')
    manifest = {str(path.relative_to(root)): sha256(path) for path in sorted(set(paths))}
    (base / 'source_and_input_SHA256.json').write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding='utf-8')
    audit['complete'] = True
    audit['total_new_main_records'] = 64000
    audit['total_supplement_records'] = 156800
    audit['total_merged_display_records'] = 288000
    audit['plots_PNG'] = len(list(base.glob('ULA/**/*.png'))) + len(list(base.glob('UPA/**/*.png')))
    audit['plots_FIG'] = len(list(base.glob('ULA/**/*.fig'))) + len(list(base.glob('UPA/**/*.fig')))
    assert audit['plots_PNG'] == 24 and audit['plots_FIG'] == 24
    (base / 'complete_audit.json').write_text(json.dumps(audit, ensure_ascii=False, indent=2), encoding='utf-8')

    overview = selected[selected.method_id.isin(['TypeI','HAQSE_R3','DFT_R3','POS_R3'])][['array_group','channel_model','branch','method_id','feedback_bits','mean_throughput_Mbps','throughput_loss_pct','feedback_saving_pct']].copy()
    overview.columns = ['阵列','信道','相位端','方法','反馈bit','吞吐量Mbps','相对Type-I损失%','相对Type-I节省bit%']
    pos = selected[selected.method_id.str.startswith('POS')][['array_group','channel_model','branch','method_id','feedback_bits','mean_throughput_Mbps','bf_gain_report','position_RMSE_m','all_pair_phase_RMSE_deg','throughput_loss_pct','feedback_saving_pct']].copy()
    pos.columns = ['阵列','信道','相位端','方法','反馈bit','吞吐量Mbps','BF幅度口径','位置RMSE_m','全AP对相位RMSE_deg','相对Type-I损失%','反馈节省%']
    threshold_view = thresholds.copy()
    threshold_view.columns = ['阵列','信道','相位端','R','最低已测坐标bit','总反馈bit','相对连续位置损失%','连续位置吞吐量Mbps','扫描最小bit','扫描最大bit']
    figures = []
    for array in ['ULA','UPA']:
        for branch in ['CPU','UE']:
            folder = base / array / branch
            figures += [f'### {array} {branch}', f'![吞吐量随观测SNR变化]({(folder / "E1_throughput_SNR.png").relative_to(base).as_posix()})', f'![吞吐量随AP数变化]({(folder / "E2_throughput_AP.png").relative_to(base).as_posix()})', f'![BF gain随观测SNR变化]({(folder / "E1_BF_gain_SNR.png").relative_to(base).as_posix()})', f'![反馈开销随AP数变化]({(folder / "E3_feedback_AP.png").relative_to(base).as_posix()})']
        figures += [f'### {array} 预算扫描', f'![CPU吞吐量随反馈bit变化]({(base / array / "E4_CPU_throughput_bits.png").relative_to(base).as_posix()})', f'![UE吞吐量随反馈bit变化]({(base / array / "E4_UE_throughput_bits.png").relative_to(base).as_posix()})']
    timing = pd.read_csv(base / 'runtime_benchmark.csv') if (base / 'runtime_benchmark.csv').exists() else pd.DataFrame()
    conclusions = []
    for ar in ['ULA', 'UPA']:
        parts = []
        for ch in ['LOS', 'LOS_NLOS_K10']:
            pp = selected[(selected.array_group == ar) & (selected.channel_model == ch) & (selected.method_id == 'POS_R3')].set_index('branch')
            label = '纯LOS' if ch == 'LOS' else 'LOS+NLOS'
            parts.append(f'{label}下，CPU/UE吞吐量为{pp.loc["CPU", "mean_throughput_Mbps"]:.3f}/{pp.loc["UE", "mean_throughput_Mbps"]:.3f} Mbps，相对Type-I损失为{pp.loc["CPU", "throughput_loss_pct"]:.3f}%/{pp.loc["UE", "throughput_loss_pct"]:.3f}%')
        conclusions.append(f'{ar}的R3位置反馈：' + '；'.join(parts) + '。')
    text = '\n\n'.join([
        '# 位置反馈加入CPU与UE仿真的结果',
        '本轮按用户确认，在CPU和UE两条链路统一采用相对AP1路径差敏感度分配坐标bit。位置编码在UE端使用连续HAQSE方向及原无权射线定位；CPU只接收量化坐标及公共配置。原正式部署、400位置、观测噪声和七条方法均保留。',
        '## 主要结果',
        '\n\n'.join(conclusions),
        markdown(overview, {'吞吐量Mbps':'.3f','相对Type-I损失%':'.3f','相对Type-I节省bit%':'.2f','反馈bit':'.0f'}),
        '正损失表示低于Type-I，负值表示略高于项目Type-I实现，不表示超过真实信道SVD上界。不同阵列部署、天线数及定位维度不同，吞吐量只在各自场景内比较。',
        '位置反馈主要把多AP方向压缩为低维几何消息。UE共相位链已用观测信道补偿AP间相位，故它的主要收益体现为降低消息开销。CPU混合信道分支仍用LOS模型，不恢复散射相位；其与UE的性能差异同时包含定位误差及模型残差，不能仅归因于坐标量化。R3/R5均按实测结果保留，不要求参考AP增加后性能单调改善。',
        '## 编码与信息流',
        'ULA固定24 bit，x/y各12 bit，步长均约2.930 mm。UPA固定36 bit，x/y/z分别13/14/9 bit，步长约1.465/0.732/13.672 mm。传输域为x/y 19–31 m，UPA z为−1–6 m。分配权重由B=8公共几何与7×7或7×7×5设计网格生成，随后在B=5–10、R3/R5、CPU/UE间固定。未用正式400点优化编码。',
        'CPU共相位只由解码位置的LOS几何模型生成。UE分支用同一坐标恢复同一组实际波束，再用全部活动AP的带噪信道计算4 bit相对相位；B=8额外28 bit。R5增加UE参考观测和定位计算，固定坐标格式下不增加坐标字段。所有相位都与实际量化后波束匹配。',
        '## 位置与相位诊断',
        markdown(pos, {'反馈bit':'.0f','吞吐量Mbps':'.3f','BF幅度口径':'.4f','位置RMSE_m':'.5f','全AP对相位RMSE_deg':'.2f','相对Type-I损失%':'.3f','反馈节省%':'.2f'}),
        '## 扫描范围内的1%编码损失门槛',
        markdown(threshold_view, {'最低已测坐标bit':'.0f','总反馈bit':'.0f','相对连续位置损失%':'.3f','连续位置吞吐量Mbps':'.3f','扫描最小bit':'.0f','扫描最大bit':'.0f'}),
        '门槛参照同一连续位置、同一全位置波束生成规则及同一相位分支的未量化参考。它是已扫描预算中的最低点，不是理论最小bit，也不是相对Type-I损失或每个用户的保证。HAQSE方向量化每维4–16 bit，保留ULA圆周网格、UPA端点网格。等bit比较仅限实际预算重合点。加权定位作为独立消融，主实验没有更换定位器。',
        '本轮加权消融以各参考AP带噪观测块的Frobenius范数平方为权重，再拟合空间频率残差，标识为energy_weighted_frequency；它与原探索报告使用去噪能量的实现不同，属于独立变体，不声称逐项复现原报告的加权结果。',
        '全部实际等bit匹配点见equal_bit_comparison.csv，包含相同参考数、信道和共相位分支的吞吐量。ULA24 bit与R3方向每维8 bit可直接比较；UPA36 bit与R3方向每维6 bit可直接比较。该比较保留各自实际消息解码和波束规则，体现完整信息格式差异，不单独归因为定位器优劣。',
        '## 计算耗时',
        markdown(timing, {'median_ms':'.4f','p10_ms':'.4f','p90_ms':'.4f'}) if not timing.empty else '详见runtime_benchmark.csv。',
        '耗时为本机MATLAB R2024a单计算线程、固定输入预热后31次测量的中位数及10/90分位。各操作分项测量，缓存为热状态；不等同于终端硬件时延、反馈空口时延或实时实现性能。',
        '## 统计与核查',
        '每个阵列/相位组180配置、72000条展示记录，其中新增位置方法16000条。四组合计新增主实验64000条；预算、方向精度和加权定位补充156800条。独立预实验每阵列24个新位置，不并入正式400点统计。CPU/UE共享基线按面板引用，重复展示记录不作为独立试验。',
        '预检对原七方法进行同输入回归。完整核查确认旧记录保留、每点400样本、唯一键无重复、CPU/UE坐标和实际波束一致、功率归一化及G≤1。详见preflight/preflight_audit.json和complete_audit.json。主图只显示均值；paired_vs_TypeI.csv提供配对吞吐差的常态近似95%区间及尾部分位，不在主图叠加阴影。',
        '频谱效率逐样本使用log2(1+rho*实际接收功率)，rho固定为原10/Sref，20 MHz换算名义吞吐量。主图BF为sqrt(mean(G))，预算附图为mean(G)，轴标签分别标注。反馈按上行消息字段计费，未计协议头、CRC、训练、CPU至AP下发，也未从名义吞吐量中扣除反馈资源。',
        f'正式主实验坐标溢出：ULA {audit["arrays"]["ULA"]["main_overflows"]}条、UPA {audit["arrays"]["UPA"]["main_overflows"]}条（CPU/UE分别计记录）。坐标越界按预定饱和规则处理，未删除不利样本。混合信道仍为K=10 dB的LOS主导受控六径模型；未扩展到纯遮挡NLOS、移动、硬件未知相位或实际UCI调度。',
        '## 复现与文件',
        '本目录下ULA、UPA分别包含主实验、补充实验、坐标配置和分配表；CPU、UE子目录提供9方法完整CSV、汇总、配对比较、PNG和MATLAB FIG。运行入口为run_position_extension.m，合并绘图为finalize_position_extension.m，完整审计与本报告由analyze_position_extension.py生成。source_and_input_SHA256.json记录源码、公共输入及旧数据哈希。',
        '```matlab\nverify_position_extension("results/new_run/preflight");\nrun_position_extension("ULA","Mode","smoke","OutputDir","results/new_run/smoke_ULA");\nrun_position_extension("UPA","Mode","smoke","OutputDir","results/new_run/smoke_UPA");\nrun_position_extension("ULA","Mode","full","OutputDir","results/new_run/ULA");\nrun_position_extension("UPA","Mode","full","OutputDir","results/new_run/UPA");\nfinalize_position_extension("ULA","results/new_run/ULA");\nfinalize_position_extension("UPA","results/new_run/UPA");\nbenchmark_position_extension("results/new_run");\n```',
        '## 实验图',
        '\n\n'.join(figures),
    ])
    (base / '位置反馈仿真结果报告.md').write_text(text + '\n', encoding='utf-8')
    print(json.dumps(audit, ensure_ascii=False, indent=2))
    print(str(base / '位置反馈仿真结果报告.md'))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('output', type=Path)
    analyze(parser.parse_args().output)
