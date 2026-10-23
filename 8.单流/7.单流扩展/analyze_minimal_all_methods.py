"""Independent message audit and descriptive statistics for the fixed minimal run."""
from pathlib import Path
import hashlib
import json
import os
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
OUT = ROOT / 'results/minimal_all_methods_20261022'
BASE = ROOT.parent / '6.位置反馈方案'


def main():
    t = pd.read_csv(OUT / 'confirm/trials.csv')
    expected = {'HAQSE_R3', 'HAQSE_R5', 'DFT_R3', 'DFT_R5', 'POS_R3', 'POS_R5', 'TypeI', 'TypeII', 'SVD_equal', 'CAP_rank_le_N'}
    assert set(t.method) == expected and len(t) == 4096
    keys = ['array', 'channel', 'method', 'branch', 'streams']
    assert t.groupby(keys).size().eq(32).all() and not t.duplicated(keys + ['sample_id']).any()
    assert set(t.sample_id) == set(range(1, 33))
    real = t.branch.ne('ideal')
    assert np.abs(t.loc[real, 'precoder_norm_squared'] - 1).max() < 1e-10
    assert t.loc[real, 'per_AP_power_error'].max() < 1e-10
    assert (t.MMSE_Mbps <= t.joint_Mbps + 1e-7).all()
    assert (t.joint_Mbps <= t.rank_capacity_Mbps + 1e-6).all()
    assert np.abs(t.loc[t.streams.eq(1), 'MMSE_Mbps'] - t.loc[t.streams.eq(1), 'joint_Mbps']).max() < 1e-7
    messages = 0
    coordinate_messages = 0
    for r in t.loc[real].itertuples(index=False):
        msg = json.loads(r.geometry_message)
        D = 9 if r.array == 'ULA' else 14
        if r.method.startswith('POS'):
            bits = np.asarray(msg['coordinate_bits']).ravel()
            index = np.asarray(msg['integer_codes']).ravel()
            assert bits.tolist() == ([12, 12] if r.array == 'ULA' else [13, 14, 9])
            assert np.all((index >= 0) & (index < 2 ** bits))
            lo = np.array([19, 19] if r.array == 'ULA' else [19, 19, -1])
            hi = np.array([31, 31] if r.array == 'ULA' else [31, 31, 6])
            p = np.asarray(json.loads(r.decoded_position)).ravel()
            assert np.max(np.abs(p[:len(bits)] - lo - (index + .5) * (hi - lo) / 2 ** bits)) < 1e-10
            db = int(bits.sum())
            ib = 0
            coordinate_messages += 1
        elif r.method.startswith(('HAQSE', 'DFT')):
            index = np.asarray(msg['integer_codes']).ravel()
            dim = 1 if r.array == 'ULA' else 2
            assert len(msg['reference_ap_ids']) == r.R and index.size == dim * r.R
            if r.method.startswith('HAQSE'):
                assert msg['q_frequency'] == 14
                limit, db = 2 ** 14, dim * r.R * 14
            else:
                limit, db = (512 if r.array == 'ULA' else 128), D * r.R
            assert np.all(index == np.round(index)) and np.all((index >= 0) & (index < limit))
            ib = 0
        else:
            index = np.asarray(msg['codes']).ravel()
            assert msg['kind'] == r.method
            assert np.all(index == np.round(index)) and np.all((index >= 0) & (index < (512 if r.array == 'ULA' else 128)))
            db = D * 8 * (1 if r.method == 'TypeI' else 2)
            ib = 0 if r.method == 'TypeI' else 5 * 8
            if r.method == 'TypeII':
                a, p = np.asarray(msg['amplitude_indices']), np.asarray(msg['phase_indices'])
                assert a.size == 8 and p.size == 8 and np.all((a >= 0) & (a < 8)) and np.all((p >= 0) & (p < 4))
        assert db == r.direction_bits and ib == r.internal_bits
        cb = 0
        if r.branch != 'CPU':
            cm = json.loads(r.coefficient_message)
            a = np.asarray(cm['amplitude_indices']).reshape(16, r.streams)
            ph = np.asarray(cm['phase_indices']).ravel()
            assert cm['B'] == 8 and cm['streams'] == r.streams and cm['qA'] == 3 and cm['qP'] == 4
            assert np.all(a == np.round(a)) and np.all((a >= 0) & (a < 8))
            assert ph.size == 16 * r.streams - r.streams and np.all((ph >= 0) & (ph < 16))
            phase = np.zeros((16, r.streams));mask = np.ones_like(phase, dtype=bool)
            mask[np.arange(r.streams), np.arange(r.streams)] = False
            flat = phase.ravel(order='F');flat[mask.ravel(order='F')] = ph
            C = a / 7 * np.exp(2j * np.pi * flat.reshape(phase.shape, order='F') / 16)
            C = C.reshape(8, 2, r.streams)
            lengths = np.linalg.norm(C, axis=(1, 2));assert np.all(lengths > 0)
            C /= lengths[:, None, None] * np.sqrt(8)
            assert abs(np.sum(abs(C) ** 2) - 1) < 1e-12
            cb = a.size * 3 + ph.size * 4
            messages += 1
        assert r.coefficient_bits == cb and r.feedback_bits == db + ib + cb
        sinr = np.asarray(json.loads(r.stream_SINR))
        assert abs(20 * np.log2(1 + sinr).sum() - r.MMSE_Mbps) < 1e-8
    geom = t[t.branch.isin(['CPU', 'UE'])]
    paired_keys = ['array', 'channel', 'method', 'sample_id', 'streams']
    assert geom.groupby(paired_keys).geometry_message.nunique().eq(1).all()
    assert geom.groupby(paired_keys).decoded_position.nunique().eq(1).all()
    merge = ['array', 'channel', 'method', 'branch', 'sample_id']
    pairs = t[t.streams.eq(1)].merge(t[t.streams.eq(2)], on=merge, suffixes=('_1', '_2'), validate='one_to_one')
    pairs['gain_Mbps'] = pairs.MMSE_Mbps_2 - pairs.MMSE_Mbps_1
    records = []
    for key, g in pairs.groupby(['array', 'channel', 'method', 'branch']):
        records.append(dict(zip(['array', 'channel', 'method', 'branch'], key), samples=len(g),
                            single_Mbps=g.MMSE_Mbps_1.mean(), dual_Mbps=g.MMSE_Mbps_2.mean(),
                            gain_Mbps=g.gain_Mbps.mean(), gain_percent=100 * g.gain_Mbps.mean() / g.MMSE_Mbps_1.mean(),
                            positive_count=int((g.gain_Mbps > 0).sum()), minimum_gain_Mbps=g.gain_Mbps.min(),
                            single_bits=g.feedback_bits_1.mean(), dual_bits=g.feedback_bits_2.mean(),
                            minimum_received_mode_ratio=g.received_mode_ratio_2.min()))
    result = pd.DataFrame(records)
    result.to_csv(OUT / 'paired_summary.csv', index=False)
    pairs[merge + ['MMSE_Mbps_1', 'MMSE_Mbps_2', 'gain_Mbps']].to_csv(OUT / 'paired_trials.csv', index=False)
    # Check unchanged POS against the previously published 400-position experiment.
    regression = 0.0
    for array in ['ULA', 'UPA']:
        old = pd.read_csv(ROOT / 'results/dual_v1_20261022' / array / 'trials.csv')
        old = old[old.stage.eq('MAIN') & old.B.eq(8) & old.observation_SNR_dB.eq(30) & old.method.isin(['POS_CPU', 'POS_UE_Q']) & old.sample_id.le(32)].copy()
        old['method_new'] = 'POS_R' + old.R.astype(str)
        old['branch_new'] = np.where(old.method.eq('POS_CPU'), 'CPU', 'UE')
        new = t[t.array.eq(array) & t.method.str.startswith('POS')]
        z = new.merge(old, left_on=['channel', 'sample_id', 'method', 'branch', 'streams'], right_on=['channel', 'sample_id', 'method_new', 'branch_new', 'streams'], validate='one_to_one')
        assert len(z) == len(new)
        err = np.abs(z.MMSE_Mbps_x - z.MMSE_Mbps_y).max()
        regression = max(regression, float(err));assert err < 1e-8
    dual = t[real & t.streams.eq(2)]
    assert len(dual) == 1792 and dual.received_mode_ratio.min() > 1e-10
    for name in ['smoke/audit.json', 'confirm/audit.json', 'legacy_reduction_audit.json']:
        assert json.loads((OUT / name).read_text())['status'] == 'passed'
    audit = dict(status='passed', records=len(t), samples_per_case=32, case_count=128,
                 original_method_families=9, extra_capacity_reference=True,
                 realizable_dual_records=len(dual), all_received_rank_two=True,
                 minimum_received_mode_ratio=float(dual.received_mode_ratio.min()),
                 coefficient_messages_independently_replayed=messages, coordinate_messages_replayed=coordinate_messages,
                 CPU_UE_geometry_equal=True, maximum_POS_regression_Mbps=regression,
                 total_power_max_error=float(np.abs(t.loc[real, 'precoder_norm_squared'] - 1).max()),
                 per_AP_power_max_error=float(t.loc[real, 'per_AP_power_error'].max()),
                 coordinate_overflow_records=int(t.quantizer_overflow.sum()),
                 negative_gain_cases=result[(result.positive_count < 32) & result.branch.ne('ideal')].to_dict('records'),
                 interpretation='minimal engineering feasibility; fixed 32 positions and nominal polarization; not a ranking or deployment-general validation')
    (OUT / 'complete_audit.json').write_text(json.dumps(audit, ensure_ascii=False, indent=2), encoding='utf-8')
    print(json.dumps(audit, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
