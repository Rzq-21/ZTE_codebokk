%% 任务2补充仿真：AP/UE侧CJT、O=8与稳定自适应参考AP
% 七条曲线：Ideal SVD、Type-II、Type-I、HAQSE（3参考AP）、
% DFT-index（3参考AP）、HAQSE（5参考AP）、DFT-index（5参考AP）。
%
% 物理部署采用文件夹《1.仅含HAQSE》和《2.对比方案》的原始配置：
% 8个AP的坐标、Euler朝向均与原始脚本一致；每个AP为16x16 UPA，O1=O2=8。
% 两个原始文件的物理点集相同，但AP3/AP8的编号顺序不同。这里采用《2.对比方案》
% 中与alpha_b、beta_b朝向向量一一对应的编号，避免把面板朝向配给错误的坐标。
% UE固定在《1.仅含HAQSE》的p_ue=[25,25,10]^T m。AP数扫描使用原始全局编号
% 的逐级包含集合[1,2,...,B]；不同AP数复用同一观测噪声/NLoS realization。
% 参考AP采用“长期稳定自适应”选择：先在覆盖区先验位置基于无噪声LoS构造
% K=3/K=5的参考集合；AP数增加时沿用上一步集合，只有候选集合的综合评分
% 明显更高才允许替换。这样避免每次Monte-Carlo和每个B值都跳变参考集合。
% 评分结合DFT峰值集中度与几何分散性；未选AP仍不提供完整CSI，因此AP侧链路
% 仍严格是“参考AP实测、其余AP重构”。
%
% LoS+NLoS信道采用《信道》中的Rician模型，并保留相对于固定UE至AP1参考距离的
% 大尺度衰减；LoS与i.i.d.复高斯NLoS共享同一距离衰减后按K=20合成。
% HAQSE和DFT-index均以“UE到AP”的局部方向为定位中间量；定位时统一
% 转换为“AP到UE”的射线。AP/CPU侧CJT中，参考AP采用实测块H_observed，
% 未反馈AP采用几何重构块H_model，拼成处理信道后计算局部DFT束和AP间相位；
% UE侧CJT中，UE使用完整实测信道计算接收合并和量化的AP间相位，局部束仍由
% HAQSE/DFT-index各自的位置估计产生。每张性能图均显示纯LoS/LoS+NLoS下的
% AP侧、UE侧共四种场景。
% 吞吐量采用相对于固定单AP参考增益的有效信道增益，不再除以每个AP数
% 各自的最大奇异值；因此图2保留了增加AP带来的绝对阵列增益。

clear; clc; close all;
rng(20261015, 'twister');

cfg = default_config();
method_names = {'Ideal SVD','Type-II','Type-I','HAQSE 3 AP', ...
    'DFT-index 3 AP','HAQSE 5 AP','DFT-index 5 AP'};
scenarios = scenario_list();
print_stable_reference_plan(cfg);
validate_direction_convention(cfg);
haqse_validation = validate_haqse_noiseless(cfg);

% 原始物理配置使用固定UE位置；蒙特卡洛样本只改变观测噪声与NLoS realization。
user_positions = sample_user_positions(cfg, cfg.N_MC);
result = struct();
result.cfg = cfg;
result.method_names = method_names;
result.scenarios = scenarios;
result.user_positions = user_positions;
result.data_snr_list = cfg.data_snr_list;
result.observation_snr_dB = cfg.observation_snr_dB;
result.ap_list = cfg.ap_list;
result.haqse_noiseless_validation = haqse_validation;

for is = 1:numel(scenarios)
    [result.rate_vs_snr{is}, result.diagnostic_samples_snr{is}, ...
        result.diagnostic_summary_snr{is}] = sweep_snr(cfg, scenarios(is), user_positions);
    [result.rate_vs_ap{is}, result.diagnostic_samples_ap{is}, ...
        result.diagnostic_summary_ap{is}] = sweep_ap(cfg, scenarios(is), user_positions);
end
result.feedback_bits = feedback_bits(cfg, scenarios);
for is = 1:numel(scenarios)
    print_diagnostic_summary(scenarios(is), result.diagnostic_summary_snr{is});
end

this_dir = fileparts(mfilename('fullpath'));
if isempty(this_dir), this_dir = pwd; end
save(fullfile(this_dir, 'task2_supplement_fixed_origin_6_results.mat'), 'result');

plot_ap_rate_pair(cfg.data_snr_list, result.rate_vs_snr, method_names, scenarios, ...
    '数据传输SNR（dB）', '任务2补充 图1 AP/UE侧CJT：O=8与稳定自适应参考AP', ...
    fullfile(this_dir, '任务2_fixed_origin_6_图1_SNR与吞吐量_AP侧与UE侧.png'));
plot_ap_rate_pair(cfg.ap_list, result.rate_vs_ap, method_names, scenarios, ...
    '激活AP数', '任务2补充 图2 AP/UE侧CJT：O=8、稳定自适应参考AP、固定总功率', ...
    fullfile(this_dir, '任务2_fixed_origin_6_图2_AP数与吞吐量_AP侧与UE侧.png'));
plot_ap_feedback_pair(cfg, result.feedback_bits, method_names, scenarios, ...
    fullfile(this_dir, '任务2_fixed_origin_6_图3_AP数与反馈开销_AP侧与UE侧.png'));
plot_stable_reference_deployment(cfg, ...
    fullfile(this_dir, '任务2_fixed_origin_6_图4_稳定自适应参考AP物理部署.png'));

fprintf('任务2补充仿真已完成。\n');
fprintf('原始部署：8 AP、16x16 UPA、O=8；AP坐标/朝向来自文件夹1和2，固定UE=[25,25,10]^T m。\n');
fprintf('Rician信道：K=%.0f（%.2f dB），带宽=%.0f MHz，固定观测SNR=%.1f dB。\n', ...
    cfg.rice_K, 10*log10(cfg.rice_K), cfg.bandwidth_hz/1e6, cfg.observation_snr_dB);
fprintf('图1横轴为相对于单AP参考链路的数据传输SNR；HAQSE频率量化：%d bit。\n', cfg.b_freq);
fprintf('本脚本同时计算AP侧CJT（参考实测块+未反馈AP重构块）和UE侧CJT（完整实测信道相位）。\n');
fprintf('AP侧参考AP采用长期稳定自适应选择；未选AP按参考块拟合误差降权，且不参与CPU实测CSI拼接。\n');
fprintf('AP激活顺序（全局编号）：%s；固定网络总发射功率。\n', ...
    mat2str(cfg.ap_activation_order));
fprintf('蒙特卡洛次数=%d：固定UE，仅重复观测噪声与NLoS realization。\n', cfg.N_MC);
fprintf('结果保存到：\n%s\n', this_dir);

%% ========================================================================
function cfg = default_config()
    cfg.fc = 7e9;
    cfg.c = 3e8;
    cfg.lambda_c = cfg.c/cfg.fc;
    cfg.d = cfg.lambda_c/2;
    cfg.k_c = 2*pi/cfg.lambda_c;

    % 文件夹《1.仅含HAQSE》《2.对比方案》的原始物理部署固定参数。
    cfg.B_full = 8;
    cfg.Nh = 16;
    cfg.Nv = 16;
    cfg.Nu = 4;
    cfg.O1 = 8;
    cfg.O2 = 8;
    cfg.origin_ue_position = [25;25;10];
    % 使用原始全局编号的逐级包含激活集合，B个AP时激活全局1:B。
    cfg.ap_activation_order = 1:cfg.B_full;
    [cfg.pb_full, cfg.Rb_full] = original_layout(cfg);

    % 自适应参考AP选择：覆盖区先验仅取部署规划中心，并非本次真实UE位置。
    % 各候选AP本地从导频计算一个DFT峰值集中度，CPU只融合该质量分数与
    % 已知AP坐标来选择K个参考AP，不会使用未选AP的完整CSI。
    cfg.reference_selection_prior = [25;25;0];
    cfg.reference_quality_weight = 0.40;
    cfg.reference_selection_regularization = 1e-6;
    % 只有新候选集合比继承集合高出此评分裕量时，才允许替换参考AP。
    cfg.reference_replacement_margin = 0.12;
    % 未反馈重构块的权重由参考块模型失配估计；用于抑制错误块的相消影响。
    cfg.unobserved_min_weight = 0.35;
    cfg.unobserved_error_weight = 1.50;

    % 任务2补充仿真参数。
    cfg.Q_phase = 4;
    cfg.b_freq = 12;
    cfg.b_typeII_amp = 3;
    cfg.b_typeII_inner_phase = 2;
    cfg.Qopt = 3;
    cfg.L_path_est = 1;
    % 图1扫描数据传输SNR；频率/DFT索引估计所用的观测SNR固定，二者不可混用。
    cfg.data_snr_list = -10:5:30;
    cfg.observation_snr_dB = 30;
    cfg.ap_list = 4:8;
    cfg.data_snr_for_ap_sweep_dB = 10;
    cfg.bandwidth_hz = 100e6;
    cfg.N_MC = 400;
    cfg.N_noiseless_validation = 1;
    % 保留《1》中固定UE位置，不引入新的随机用户分布。
    cfg.user_box = [cfg.origin_ue_position, cfg.origin_ue_position];
    % 保留相对距离衰减；以固定UE到AP1的距离作为单位幅度参考距离。
    cfg.pathloss_reference_m = norm(cfg.pb_full(:,1)-cfg.origin_ue_position);
    % 数据SNR以参考距离处单AP的理想阵列增益Nu*Nb为基准。该常数不随
    % 激活AP数变化，因此不会消除多AP带来的绝对信道增益。
    cfg.snr_gain_reference = cfg.Nu*cfg.Nh*cfg.Nv;
    cfg.trial_seed_base = 20262000;
    cfg.noise_seed_offset = 500000;
    % 《信道》报告规定的弱NLoS Rician参数。
    cfg.rice_K = 20;
    cfg.reference_plan = build_stable_reference_plan(cfg);
end

function scenarios = scenario_list()
    scenarios(1) = struct('name', '纯LoS + AP侧拼接重构信道CJT', ...
        'channel_mode', 1, 'link_mode', 1);
    scenarios(2) = struct('name', 'LoS+NLoS + AP侧拼接重构信道CJT', ...
        'channel_mode', 2, 'link_mode', 1);
    scenarios(3) = struct('name', '纯LoS + UE侧完整实测信道CJT', ...
        'channel_mode', 1, 'link_mode', 2);
    scenarios(4) = struct('name', 'LoS+NLoS + UE侧完整实测信道CJT', ...
        'channel_mode', 2, 'link_mode', 2);
end

function [avg_rate, diagnostic_samples, diagnostic_summary] = sweep_snr(cfg, scenario, user_positions)
    % 波束增益与数据传输SNR无关：每个位置仅做一次固定观测SNR下的反馈估计，
    % 再代入不同数据SNR的容量公式，避免横轴各点因随机重抽样而抖动。
    B = cfg.B_full;
    diagnostic_samples = collect_trial_samples(cfg, B, cfg.observation_snr_dB, scenario, user_positions);
    avg_rate = nan(7, numel(cfg.data_snr_list));
    for ix = 1:numel(cfg.data_snr_list)
        avg_rate(:,ix) = average_rate_from_gain(cfg, diagnostic_samples.link_gain, ...
            cfg.data_snr_list(ix));
    end
    diagnostic_summary = summarize_diagnostics(diagnostic_samples);
end

function [avg_rate, diagnostic_samples, diagnostic_summary] = sweep_ap(cfg, scenario, user_positions)
    avg_rate = nan(7, numel(cfg.ap_list));
    diagnostic_samples = cell(1,numel(cfg.ap_list));
    diagnostic_summary = cell(1,numel(cfg.ap_list));
    for ix = 1:numel(cfg.ap_list)
        diagnostic_samples{ix} = collect_trial_samples(cfg, cfg.ap_list(ix), ...
            cfg.observation_snr_dB, scenario, user_positions);
        avg_rate(:,ix) = average_rate_from_gain(cfg, diagnostic_samples{ix}.link_gain, ...
            cfg.data_snr_for_ap_sweep_dB);
        diagnostic_summary{ix} = summarize_diagnostics(diagnostic_samples{ix});
    end
end

function samples = collect_trial_samples(cfg, B, observation_snr_dB, scenario, user_positions)
    samples = initialize_diagnostics(cfg.N_MC);
    [pb, Rb] = active_layout(cfg, B);
    for imc = 1:cfg.N_MC
        % 相同信道模式、用户编号在不同B及AP/UE链路中使用同一种子。
        % 配合逐级包含的AP顺序，可保证较小B的信道块是较大B的前缀。
        trial_seed = cfg.trial_seed_base + 100000*scenario.channel_mode + imc;
        [gain, trial_diagnostic] = one_trial(cfg, B, pb, Rb, user_positions(:,imc), ...
            observation_snr_dB, scenario, trial_seed);
        samples.link_gain(:,imc) = gain;
        samples.normalized_gain(:,imc) = trial_diagnostic.normalized_gain;
        samples.position_error_m(:,imc) = trial_diagnostic.position_error_m;
        samples.frequency_rms_cpe(:,imc) = trial_diagnostic.frequency_rms_cpe;
        samples.processing_channel_nmse(:,imc) = trial_diagnostic.processing_channel_nmse;
        samples.oracle_gain(:,imc) = trial_diagnostic.oracle_gain;
    end
end

function avg_rate = average_rate_from_gain(cfg, link_gain, data_snr_dB)
    % link_gain=||H*W||^2/(Nu*Nb)采用固定单AP参考值，绝不随激活AP数改变。
    % 因此R=BW*log2(1+rho*link_gain)保留了增加AP带来的阵列增益。
    sample_rate = cfg.bandwidth_hz/1e6 * log2(1 + ...
        10^(data_snr_dB/10)*link_gain);
    avg_rate = mean(sample_rate, 2, 'omitnan');
end

function [gain, diagnostic] = one_trial(cfg, B, pb, Rb, p_ue, observation_snr_dB, ...
        scenario, trial_seed)
    cfg_case = cfg;
    cfg_case.B = B;
    cfg_case.Nb = cfg.Nh*cfg.Nv;
    cfg_case.N = B*cfg_case.Nb;

    rng(trial_seed, 'twister');
    [H_clean, dft_basis] = build_channel(cfg_case, p_ue, pb, Rb, scenario.channel_mode);
    % 观测噪声相对于参考距离处单位阵元信道固定，不随激活AP集合重新标定。
    noise_power = 1 / 10^(observation_snr_dB/10);
    % 独立重置噪声种子，使不同B的前缀信道块获得完全相同的观测噪声。
    rng(trial_seed + cfg.noise_seed_offset, 'twister');
    H_observed = H_clean + sqrt(noise_power/2) * ...
        (randn(size(H_clean)) + 1j*randn(size(H_clean)));
    [~, S_clean, ~] = svd(H_clean, 'econ');
    sigma1_sq = S_clean(1,1)^2;
    [~,~,V_observed] = svd(H_observed, 'econ');

    gain = nan(7,1);
    gain(1) = sigma1_sq/cfg.snr_gain_reference;
    diagnostic = empty_trial_diagnostic();

    W_typeII = typeII_precoder(V_observed, dft_basis, cfg_case);
    gain(2) = effective_link_gain(H_clean, W_typeII, cfg);
    W_typeI = typeI_precoder(V_observed, dft_basis, cfg_case);
    gain(3) = effective_link_gain(H_clean, W_typeI, cfg);

    for ii = 1:2
        K = 3 + 2*(ii-1);
        % 图片规定K=B是合法端点：此时H_complete中的所有AP均为实测块，
        % 不再含未反馈AP的几何重构块。
        if B < K, continue; end
        % 同一B、同一K在所有Monte-Carlo realization及LoS/NLoS场景中采用
        % 相同的长期参考集合，避免随机选择切换造成图2的人为跳变。
        refs = stable_reference_aps(cfg_case, B, K);
        u_true = true_local_spatial_frequency(p_ue, refs, pb, Rb, cfg_case);

        % Oracle仅作诊断：真实位置代入同一条AP/UE处理链路。K=3和K=5
        % 分别记录，以反映参考实测块占比变化。
        [W_oracle, ~, ~] = limited_feedback_precoder(H_observed, p_ue, refs, pb, Rb, ...
            cfg_case, dft_basis, scenario.link_mode, 'haqse');
        diagnostic.oracle_gain(ii) = normalized_gain(H_clean, W_oracle, sigma1_sq);

        [p_h, f_h] = estimate_position_haqse(H_observed, refs, pb, Rb, cfg_case);
        [W_h, ~, H_processing_h] = limited_feedback_precoder(H_observed, p_h, refs, pb, Rb, ...
            cfg_case, dft_basis, scenario.link_mode, 'haqse');
        gain(2 + 2*ii) = effective_link_gain(H_clean, W_h, cfg_case);
        diagnostic = record_limited_feedback_diagnostic(diagnostic, 2*ii-1, p_h, p_ue, ...
            -f_h, u_true, H_processing_h, H_clean);

        [p_d, u_d] = estimate_position_dft(H_observed, refs, pb, Rb, cfg_case, dft_basis);
        [W_d, ~, H_processing_d] = limited_feedback_precoder(H_observed, p_d, refs, pb, Rb, ...
            cfg_case, dft_basis, scenario.link_mode, 'dft');
        gain(3 + 2*ii) = effective_link_gain(H_clean, W_d, cfg_case);
        diagnostic = record_limited_feedback_diagnostic(diagnostic, 2*ii, p_d, p_ue, ...
            u_d, u_true, H_processing_d, H_clean);
    end
    % 单独保存相对于同一B理想SVD的预编码效率，仅用于诊断，不再代入吞吐量。
    diagnostic.normalized_gain = gain/max(gain(1),realmin);
end

function [W, H_model, H_processing] = limited_feedback_precoder(H_observed, p_est, refs, ...
        pb, Rb, cfg, dft_basis, link_mode, beam_mode)
    H_model = build_channel(cfg, p_est, pb, Rb, 1);
    if link_mode == 1
        % AP/CPU侧：参考AP保留实测信道，未反馈AP使用几何重构信道。
        % 几何模型只负责未反馈AP；参考AP实测块由AP本地测量/前传提供。
        % 先用参考实测块与模型块的失配估计未反馈块可靠度；失配大时降低
        % 未反馈块在SVD与最终预编码中的权重，避免错误相位造成相消。
        ap_weight = model_block_reliability(H_observed, H_model, refs, cfg);
        H_processing = H_model;
        for b = 1:cfg.B
            if ~ismember(b, refs)
                idx = (b-1)*cfg.Nb+(1:cfg.Nb);
                H_processing(:,idx) = ap_weight(b)*H_model(:,idx);
            end
        end
        for ii = 1:numel(refs)
            b = refs(ii);
            idx = (b-1)*cfg.Nb+(1:cfg.Nb);
            H_processing(:,idx) = H_observed(:,idx);
        end
        [~,~,V_processing] = svd(H_processing, 'econ');
        W = model_dft_precoder(V_processing, dft_basis, cfg, ap_weight);
    else
        % UE侧仅从完整实测信道计算接收合并向量和AP间CJT相位；局部束
        % 仍由HAQSE/DFT各自的定位结果产生，避免多种方案退化为同一Type-I。
        H_processing = H_model;
        [U_observed,~,~] = svd(H_observed, 'econ');
        if strcmp(beam_mode, 'haqse')
            % HAQSE保留连续空间频率优势，不再二次量化到DFT网格。
            local = continuous_local_beams(H_model, cfg);
        elseif strcmp(beam_mode, 'dft')
            [~,~,V_model] = svd(H_model, 'econ');
            local = local_dft_beams(V_model, dft_basis, cfg);
        else
            error('未知的UE侧局部波束模式：%s。', beam_mode);
        end
        W = apply_observed_quantized_phase(local, H_observed, U_observed(:,1), cfg);
    end
end

function ap_weight = model_block_reliability(H_observed, H_model, refs, cfg)
    % 参考AP同时拥有实测块和同一位置估计产生的模型块，可据此估计几何重构
    % 的可信度。该可信度不读取未反馈AP的实测CSI。
    reference_nmse = zeros(1,numel(refs));
    for ii = 1:numel(refs)
        b = refs(ii);
        idx = (b-1)*cfg.Nb+(1:cfg.Nb);
        reference_nmse(ii) = norm(H_model(:,idx)-H_observed(:,idx),'fro')^2 / ...
            max(norm(H_observed(:,idx),'fro')^2, realmin);
    end
    model_error = median(reference_nmse);
    reconstructed_weight = 1/(1 + cfg.unobserved_error_weight*sqrt(model_error));
    reconstructed_weight = max(cfg.unobserved_min_weight, min(1, reconstructed_weight));
    ap_weight = reconstructed_weight*ones(cfg.B,1);
    ap_weight(refs) = 1;
end

function [p_est, u_est] = estimate_position_dft(H, refs, pb, Rb, cfg, dft_basis)
    index_dft = zeros(1, numel(refs));
    for ii = 1:numel(refs)
        b = refs(ii);
        block = H(:, (b-1)*cfg.Nb+(1:cfg.Nb));
        [~,~,V] = svd(block, 'econ');
        % V(:,1)是阵列响应的共轭。用Hermitian相关得到UE→AP方向频率，
        % 与HAQSE的频率方向保持一致。
        index_dft(ii) = best_dft_index(V(:,1), dft_basis, cfg);
    end
    u_est = dft_index_to_ue_to_ap_frequency(index_dft, cfg);
    p_est = solve_ue_to_ap_rays(u_est(1,:)*cfg.lambda_c/cfg.d, ...
        u_est(2,:)*cfg.lambda_c/cfg.d, refs, pb, Rb);
end

function [p_est, f_est, f_raw] = estimate_position_haqse(H, refs, pb, Rb, cfg)
    f_est = zeros(2, numel(refs));
    f_raw = zeros(2, numel(refs));
    for ii = 1:numel(refs)
        b = refs(ii);
        block = H(:, (b-1)*cfg.Nb+(1:cfg.Nb));
        tensor = reshape(block.', [cfg.Nv, cfg.Nh, cfg.Nu]);
        f = MMV_CF_mQSE2_SumOfRatios(tensor, cfg.L_path_est, 1, 1, cfg.Qopt);
        f_raw(:,ii) = f(1,1:2).';
        % HAQSE输出定义为f=-d/lambda*[q_y,q_z]，单位是cycles/element；
        % 先记录原始输出，再在同一规范[-0.5,0.5]内按统一b_freq量化。
        f_est(:,ii) = endpoint_quantize(f_raw(:,ii), cfg.b_freq);
    end
    qy = max(-1, min(1, -f_est(1,:)*cfg.lambda_c/cfg.d));
    qz = max(-1, min(1, -f_est(2,:)*cfg.lambda_c/cfg.d));
    p_est = solve_ue_to_ap_rays(qy, qz, refs, pb, Rb);
end

function f = endpoint_quantize(f, bits)
    levels = 2^bits;
    f = -0.5 + round((f+0.5)*(levels-1))/(levels-1);
    f = max(-0.5, min(0.5, f));
end

function p = solve_index_rays(index_dft, refs, pb, Rb, cfg)
    u = dft_index_to_ue_to_ap_frequency(index_dft, cfg);
    p = solve_ue_to_ap_rays(u(1,:)*cfg.lambda_c/cfg.d, ...
        u(2,:)*cfg.lambda_c/cfg.d, refs, pb, Rb);
end

function u = dft_index_to_ue_to_ap_frequency(index_dft, cfg)
    % scores的行对应Nv*O1（阵列局部y），列对应Nh*O2（阵列局部z）。
    Ny = cfg.Nv*cfg.O1;
    Nz = cfg.Nh*cfg.O2;
    ky = mod(index_dft-1, Ny);
    kz = floor((index_dft-1)/Ny);
    ky(ky/Ny > 0.5) = ky(ky/Ny > 0.5)-Ny;
    kz(kz/Nz > 0.5) = kz(kz/Nz > 0.5)-Nz;
    % u=d/lambda*q，单位cycles/element；DFT索引给出的正负号与UE→AP一致。
    u = [ky/Ny; kz/Nz];
end

function p = solve_ue_to_ap_rays(qy, qz, refs, pb, Rb)
    % HAQSE/DFT反馈的qy、qz统一定义为UE→AP局部方向余弦。当前物理部署
    % 中所有面板均面向服务区中心、UE均在面板正前方，因此UE→AP的局部
    % 法向分量恒为负。定位射线必须取其相反方向，即AP→UE。
    A = zeros(3);
    rhs = zeros(3,1);
    for ii = 1:numel(refs)
        qx_ue_to_ap = -sqrt(max(0, 1-qy(ii)^2-qz(ii)^2));
        q_ue_to_ap = [qx_ue_to_ap;qy(ii);qz(ii)];
        g_ap_to_ue = -Rb(:,:,refs(ii))*q_ue_to_ap;
        g_ap_to_ue = g_ap_to_ue/norm(g_ap_to_ue);
        P = eye(3)-g_ap_to_ue*g_ap_to_ue.';
        A = A+P;
        rhs = rhs+P*pb(:,refs(ii));
    end
    p = A\rhs;
end

function W = typeI_precoder(V, dft_basis, cfg)
    local = local_dft_beams(V, dft_basis, cfg);
    % Type-I由UE基于完整实测信道选择单束，并反馈量化相对相位。
    W = apply_quantized_model_phase(local, V, cfg);
end

function W = model_dft_precoder(V_model, dft_basis, cfg, ap_weight)
    % AP/CPU链路：局部束和相干相位均由拼接信道产生，因而不需要UE
    % 额外反馈AP间相位字段。
    local = local_dft_beams(V_model, dft_basis, cfg);
    W = apply_model_phase(local, V_model, cfg, ap_weight);
end

function local = continuous_local_beams(H_model, cfg)
    % 每个AP模型块的主右奇异向量即连续方向束，不经过DFT网格量化。
    local = zeros(cfg.Nb, cfg.B);
    for b = 1:cfg.B
        idx = (b-1)*cfg.Nb+(1:cfg.Nb);
        [~,~,V_block] = svd(H_model(:,idx), 'econ');
        local(:,b) = V_block(:,1)/norm(V_block(:,1));
    end
end

function W = apply_observed_quantized_phase(local, H_observed, u_observed, cfg)
    % UE只反馈各AP相对于第一个AP的Q_phase比特有效信道相位。
    gains = zeros(cfg.B,1);
    for b = 1:cfg.B
        idx = (b-1)*cfg.Nb+(1:cfg.Nb);
        gains(b) = u_observed' * H_observed(:,idx) * local(:,b);
    end
    W = apply_quantized_relative_phase(local, gains, cfg);
end

function local = local_dft_beams(V, dft_basis, cfg)
    local = zeros(cfg.Nb, cfg.B);
    for b = 1:cfg.B
        idx = (b-1)*cfg.Nb+(1:cfg.Nb);
        beam_index = best_dft_index(V(idx,1), dft_basis, cfg);
        local(:,b) = dft_beam(beam_index, dft_basis);
    end
end

function W = apply_model_phase(local, V, cfg, ap_weight)
    W = zeros(cfg.N,1);
    for b = 1:cfg.B
        idx = (b-1)*cfg.Nb+(1:cfg.Nb);
        coefficient = local(:,b)'*V(idx,1);
        W(idx) = ap_weight(b)*coefficient/max(abs(coefficient),realmin)*local(:,b);
    end
    W = W/norm(W);
end

function W = apply_quantized_model_phase(local, V, cfg)
    gains = zeros(cfg.B,1);
    for b = 1:cfg.B
        idx = (b-1)*cfg.Nb+(1:cfg.Nb);
        gains(b) = V(idx,1)'*local(:,b);
    end
    W = apply_quantized_relative_phase(local, gains, cfg);
end

function W = apply_quantized_relative_phase(local, gains, cfg)
    % 第一个AP是相位参考，因整体公共相位不影响单流增益，故只反馈其余(B-1)个
    % 相对相位。该实现与feedback_bits中的(B-1)*Q_phase严格一致。
    W = zeros(cfg.N,1);
    phase_cb = exp(1j*2*pi*(0:2^cfg.Q_phase-1)/2^cfg.Q_phase);
    for b = 1:cfg.B
        idx = (b-1)*cfg.Nb+(1:cfg.Nb);
        if b == 1
            W(idx) = local(:,b);
        else
            desired = exp(1j*(angle(gains(1))-angle(gains(b))));
            [~,k] = min(abs(phase_cb-desired));
            W(idx) = phase_cb(k)*local(:,b);
        end
    end
    W = W/norm(W);
end

function W = typeII_precoder(V, dft_basis, cfg)
    % 单流Type-II：每AP选两个最强DFT束，并量化第二束的相对幅度和相位。
    local = zeros(cfg.Nb, cfg.B);
    gains = zeros(cfg.B,1);
    amp_cb = linspace(0,1,2^cfg.b_typeII_amp);
    phase_cb = exp(1j*2*pi*(0:2^cfg.b_typeII_inner_phase-1)/2^cfg.b_typeII_inner_phase);
    for b = 1:cfg.B
        idx = (b-1)*cfg.Nb+(1:cfg.Nb);
        v = V(idx,1);
        selected = top_dft_indices(v, dft_basis, cfg, 2);
        F = [dft_beam(selected(1),dft_basis), dft_beam(selected(2),dft_basis)];
        coeff = F\v;
        [~,r] = max(abs(coeff));
        other = 3-r;
        ratio = abs(coeff(other))/max(abs(coeff(r)), realmin);
        [~,ai] = min(abs(amp_cb-ratio));
        desired = exp(1j*(angle(coeff(other))-angle(coeff(r))));
        [~,phase_index] = min(abs(phase_cb-desired));
        c = zeros(2,1);
        c(r) = 1;
        c(other) = amp_cb(ai)*phase_cb(phase_index);
        local(:,b) = F*c;
        local(:,b) = local(:,b)/norm(local(:,b));
        gains(b) = v'*local(:,b);
    end
    W = apply_quantized_relative_phase(local, gains, cfg);
end

function ratio = normalized_gain(H, W, sigma1_sq)
    ratio = norm(H*W)^2/sigma1_sq;
end

function value = effective_link_gain(H, W, cfg)
    % 固定参考增益Nu*Nb不随AP数变化；该量可直接用于跨AP数吞吐量比较。
    value = norm(H*W)^2/cfg.snr_gain_reference;
end

function diagnostic = empty_trial_diagnostic()
    % 行顺序固定为：HAQSE-3、DFT-3、HAQSE-5、DFT-5。
    diagnostic.position_error_m = nan(4,1);
    diagnostic.frequency_rms_cpe = nan(4,1);
    diagnostic.processing_channel_nmse = nan(4,1);
    diagnostic.oracle_gain = nan(2,1); % K=3与K=5
    diagnostic.normalized_gain = nan(7,1);
end

function samples = initialize_diagnostics(sample_count)
    samples.method_names = {'HAQSE 3 AP','DFT-index 3 AP','HAQSE 5 AP','DFT-index 5 AP'};
    samples.link_gain = nan(7,sample_count);
    samples.normalized_gain = nan(7,sample_count);
    samples.position_error_m = nan(4,sample_count);
    samples.frequency_rms_cpe = nan(4,sample_count);
    samples.processing_channel_nmse = nan(4,sample_count);
    samples.oracle_gain = nan(2,sample_count);
end

function diagnostic = record_limited_feedback_diagnostic(diagnostic, row, p_est, p_true, ...
        u_est, u_true, H_processing, H_clean)
    diagnostic.position_error_m(row) = norm(p_est-p_true);
    % u=d/lambda*q是归一化空间频率，单位cycles/element；将HAQSE的-f与
    % DFT-index的u均转换到该同一规范后，才可直接比较RMS误差。
    diagnostic.frequency_rms_cpe(row) = sqrt(mean(abs(u_est(:)-u_true(:)).^2));
    % AP侧对应H_complete的NMSE；UE侧对应产生局部束的几何模型NMSE。
    diagnostic.processing_channel_nmse(row) = norm(H_processing-H_clean,'fro')^2 / ...
        max(norm(H_clean,'fro')^2, realmin);
end

function u = true_local_spatial_frequency(p_ue, refs, pb, Rb, cfg)
    u = zeros(2,numel(refs));
    for ii = 1:numel(refs)
        b = refs(ii);
        q_ue_to_ap = Rb(:,:,b)'*(pb(:,b)-p_ue);
        q_ue_to_ap = q_ue_to_ap/norm(q_ue_to_ap);
        u(:,ii) = cfg.d/cfg.lambda_c*q_ue_to_ap(2:3);
    end
end

function summary = summarize_diagnostics(samples)
    summary.method_names = samples.method_names;
    position_count = sum(~isnan(samples.position_error_m),2);
    gain_count = sum(~isnan(samples.link_gain),2);
    summary.mean_position_error_m = mean(samples.position_error_m,2,'omitnan');
    summary.std_position_error_m = std(samples.position_error_m,0,2,'omitnan');
    summary.ci95_position_error_m = 1.96*summary.std_position_error_m ./ max(sqrt(position_count),1);
    summary.mean_frequency_rms_cpe = mean(samples.frequency_rms_cpe,2,'omitnan');
    summary.mean_processing_channel_nmse = mean(samples.processing_channel_nmse,2,'omitnan');
    summary.mean_link_gain = mean(samples.link_gain,2,'omitnan');
    summary.mean_normalized_gain = mean(samples.normalized_gain,2,'omitnan');
    summary.std_gain = std(samples.link_gain,0,2,'omitnan');
    summary.ci95_gain = 1.96*summary.std_gain ./ max(sqrt(gain_count),1);
    summary.mean_oracle_gain = mean(samples.oracle_gain,2,'omitnan');
end

function [H, dft_basis] = build_channel(cfg, p_ue, pb, Rb, channel_mode)
    p_lcs = array_positions(cfg);
    ue_lcs = zeros(3,cfg.Nu);
    ue_lcs(1,:) = ((0:cfg.Nu-1)-(cfg.Nu-1)/2)*cfg.d;
    H = zeros(cfg.Nu,cfg.N);
    for b = 1:cfg.B
        idx = (b-1)*cfg.Nb+(1:cfg.Nb);
        H_los = geometric_los_block(p_ue, pb(:,b), Rb(:,:,b), p_lcs, ue_lcs, cfg);
        % 保留相对于pathloss_reference_m的距离衰减，不再逐AP单独归一化。
        large_scale_amp = cfg.pathloss_reference_m/norm(pb(:,b)-p_ue);
        if channel_mode == 1
            H(:,idx) = H_los;
        else
            % NLoS与LoS共享同一个大尺度距离衰减，再按Rician K因子合成。
            H_nlos = large_scale_amp * ...
                (randn(cfg.Nu,cfg.Nb)+1j*randn(cfg.Nu,cfg.Nb))/sqrt(2);
            H(:,idx) = sqrt(cfg.rice_K/(cfg.rice_K+1))*H_los + ...
                sqrt(1/(cfg.rice_K+1))*H_nlos;
        end
    end
    dft_basis = make_dft_basis(cfg);
end

function H = geometric_los_block(p_ue, p_ap, R_ap, p_lcs, ue_lcs, cfg)
    d_vec = p_ap-p_ue;
    range = norm(d_vec);
    d_hat = d_vec/range;
    ap_dir = R_ap'*d_hat;
    phase_center = cfg.pathloss_reference_m*exp(-1j*cfg.k_c*range)/range;
    ap_phase = exp(-1j*cfg.k_c*(ap_dir.'*p_lcs));
    H = zeros(cfg.Nu,cfg.Nb);
    for m = 1:cfg.Nu
        ue_phase = exp(1j*cfg.k_c*(d_hat.'*ue_lcs(:,m)));
        H(m,:) = phase_center*ue_phase*ap_phase;
    end
end

function basis = make_dft_basis(cfg)
    % 仅保存两个一维DFT字典，避免构造完整Kron码本。
    basis.h = exp(1j*2*pi*(0:cfg.Nh-1).'*(0:cfg.O2*cfg.Nh-1)/(cfg.O2*cfg.Nh))/sqrt(cfg.Nh);
    basis.v = exp(1j*2*pi*(0:cfg.Nv-1).'*(0:cfg.O1*cfg.Nv-1)/(cfg.O1*cfg.Nv))/sqrt(cfg.Nv);
end

function indices = top_dft_indices(v, basis, cfg, count)
    X = reshape(v, cfg.Nv, cfg.Nh);
    scores = basis.v' * X * conj(basis.h);
    [~,order] = sort(abs(scores(:)), 'descend');
    indices = order(1:count);
end

function index = best_dft_index(v, basis, cfg)
    X = reshape(v, cfg.Nv, cfg.Nh);
    % 全链路统一使用Hermitian相关；输出是UE→AP局部方向对应的DFT频率。
    scores = basis.v' * X * conj(basis.h);
    [~,index] = max(abs(scores(:)));
end

function beam = dft_beam(index, basis)
    Nv_grid = size(basis.v,2);
    iv = mod(index-1, Nv_grid)+1;
    ih = floor((index-1)/Nv_grid)+1;
    beam = kron(basis.h(:,ih), basis.v(:,iv));
end

function p = array_positions(cfg)
    p = zeros(3,cfg.Nb);
    k = 0;
    for nh = 1:cfg.Nh
        for nv = 1:cfg.Nv
            k = k+1;
            p(:,k) = [0;(nv-1-(cfg.Nv-1)/2)*cfg.d; ...
                (nh-1-(cfg.Nh-1)/2)*cfg.d];
        end
    end
end

function [pb, Rb] = original_layout(cfg)
    % 原始AP坐标：文件夹1/2的物理点集完全相同。
    % 此处选择文件夹2的AP3/AP8编号，以使其与原始alpha_b、beta_b逐项配对。
    pb = [0,  30, 0,  50, -1,  0,  0,  20; ...
          0, -10,50,  50, -1, -1, 51, -10; ...
          0,   0, 0,   0,  0,  0,  0,   0];
    assert(size(pb,2) == cfg.B_full, '原始AP坐标数与B_full不一致。');

    % 原始脚本给出的面板Euler朝向（单位：度）。
    alpha_b = [0;130;-45;-125;49;43;98;90]*pi/180;
    beta_b  = [0;10;-5;2;6;-4;8;-7]*pi/180;
    gamma_b = zeros(cfg.B_full,1);
    Rb = zeros(3,3,cfg.B_full);
    for b = 1:cfg.B_full
        ca = cos(alpha_b(b)); sa = sin(alpha_b(b));
        cb = cos(beta_b(b));  sb = sin(beta_b(b));
        cg = cos(gamma_b(b)); sg = sin(gamma_b(b));
        % 与原始文件中的Euler旋转矩阵写法保持一致。
        Rb(:,:,b) = [ca*cb, ca*sb*sg-sa*cg, ca*sb*cg+sa*sg; ...
                     sa*cb, sa*sb*sg+ca*cg, sa*sb*cg-ca*sg; ...
                     -sb,   cb*sg,            cb*cg];
    end
end

function [pb, Rb] = active_layout(cfg, B)
    active = active_ap_indices(cfg, B);
    pb = cfg.pb_full(:,active);
    Rb = cfg.Rb_full(:,:,active);
end

function active = active_ap_indices(cfg, B)
    assert(B >= 4 && B <= cfg.B_full, '仅支持4至8个激活AP。');
    % 严格逐级包含：A_4属于A_5属于...属于A_8。
    active = cfg.ap_activation_order(1:B);
end

function plan = build_stable_reference_plan(cfg)
    % 用固定的覆盖区先验、无噪声LoS和O=8码本预先建立长期参考AP计划。
    % 因此该计划与每次观测噪声、NLoS realization无关，保证B扫描可比。
    plan.B = cfg.ap_list;
    plan.refs3 = cell(1,numel(cfg.ap_list));
    plan.refs5 = cell(1,numel(cfg.ap_list));
    previous3 = [];
    previous5 = [];
    for ib = 1:numel(cfg.ap_list)
        B = cfg.ap_list(ib);
        cfg_case = cfg;
        cfg_case.B = B;
        cfg_case.Nb = cfg.Nh*cfg.Nv;
        cfg_case.N = B*cfg_case.Nb;
        [pb,Rb] = active_layout(cfg_case, B);
        [H_prior,dft_basis] = build_channel(cfg_case, cfg.reference_selection_prior, pb, Rb, 1);
        quality = reference_quality_scores(H_prior, B, cfg_case, dft_basis);
        plan.refs3{ib} = choose_stable_reference_subset(quality, B, 3, pb, cfg_case, previous3);
        previous3 = plan.refs3{ib};
        if B >= 5
            plan.refs5{ib} = choose_stable_reference_subset(quality, B, 5, pb, cfg_case, previous5);
            previous5 = plan.refs5{ib};
        end
    end
end

function refs = stable_reference_aps(cfg, B, K)
    ib = find(cfg.reference_plan.B == B, 1, 'first');
    assert(~isempty(ib), '稳定参考AP计划中缺少B=%d。', B);
    if K == 3
        refs = cfg.reference_plan.refs3{ib};
    elseif K == 5
        refs = cfg.reference_plan.refs5{ib};
    else
        error('仅支持K=3或K=5。');
    end
    assert(numel(refs) == K, 'B=%d、K=%d的参考AP计划无效。', B, K);
end

function refs = choose_stable_reference_subset(quality, B, K, pb, cfg, previous_refs)
    assert(K <= B && (K == 3 || K == 5), '仅支持K=3或K=5，且K不能超过激活AP数。');
    if K == B
        refs = 1:B;
        return;
    end
    [candidate_refs,candidate_score] = best_reference_subset(quality, B, K, pb, cfg);
    if isempty(previous_refs)
        refs = candidate_refs;
        return;
    end
    inherited_score = reference_subset_score(previous_refs, quality, pb, cfg);
    % 新集合必须显著优于继承集合才替换，形成跨B的选择滞回。
    if candidate_score-inherited_score >= cfg.reference_replacement_margin
        refs = candidate_refs;
    else
        refs = previous_refs;
    end
end

function [refs,best_score] = best_reference_subset(quality, B, K, pb, cfg)
    candidates = nchoosek(1:B,K);
    score = -inf(size(candidates,1),1);
    for row = 1:size(candidates,1)
        score(row) = reference_subset_score(candidates(row,:), quality, pb, cfg);
    end
    [best_score,best] = max(score);
    refs = candidates(best,:);
end

function score = reference_subset_score(subset, quality, pb, cfg)
    % det越大代表依据覆盖区先验形成的射线交会几何越稳定、越分散。
    A = zeros(3);
    for ii = 1:numel(subset)
        b = subset(ii);
        g = cfg.reference_selection_prior-pb(:,b);
        g = g/max(norm(g),realmin);
        A = A + eye(3)-g*g.';
    end
    geometry_score = log(max(det(A + cfg.reference_selection_regularization*eye(3)), realmin));
    quality_score = sum(log(max(quality(subset),realmin)));
    score = geometry_score + cfg.reference_quality_weight*quality_score;
end

function quality = reference_quality_scores(H, B, cfg, dft_basis)
    quality = zeros(B,1);
    for b = 1:B
        idx = (b-1)*cfg.Nb+(1:cfg.Nb);
        [~,S,V] = svd(H(:,idx), 'econ');
        dominant_ratio = S(1,1)^2 / max(norm(S,'fro')^2, realmin);
        peak_concentration = dft_peak_concentration(V(:,1), dft_basis, cfg);
        quality(b) = sqrt(max(dominant_ratio,realmin) * ...
            max(peak_concentration,realmin));
    end
    quality = quality/max(max(quality),realmin);
end

function score = dft_peak_concentration(v, dft_basis, cfg)
    % 只输出标量峰值质量；不输出候选AP的DFT索引或完整信道。
    X = reshape(v, cfg.Nv, cfg.Nh);
    spectrum = dft_basis.v' * X * conj(dft_basis.h);
    power = abs(spectrum(:)).^2;
    score = max(power)/max(sum(power), realmin);
end

function positions = sample_user_positions(cfg, count)
    % 原始方案采用固定UE；不同样本只服务于观测噪声/NLoS的蒙特卡洛平均。
    positions = repmat(cfg.origin_ue_position, 1, count);
end

function bits = feedback_bits(cfg, scenarios)
    index_bits = ceil(log2(cfg.Nv*cfg.O1))+ceil(log2(cfg.Nh*cfg.O2));
    bits = nan(numel(scenarios), numel(cfg.ap_list), 7);
    for is = 1:numel(scenarios)
        for i = 1:numel(cfg.ap_list)
            B = cfg.ap_list(i);
            % Type-I/II以及UE实测相位链路均反馈(B-1)个相对相位。
            % AP侧拼接重构链路由CPU从H_complete计算相位，不计UE相位字段。
            codebook_phase_bits = (B-1)*cfg.Q_phase;
            if scenarios(is).link_mode == 1
                haqse3_bits = 3*2*cfg.b_freq;
                dft3_bits = 3*index_bits;
                haqse5_bits = ternary(B>=5,5*2*cfg.b_freq,NaN);
                dft5_bits = ternary(B>=5,5*index_bits,NaN);
            else
                % UE侧仍只用K个参考AP的HAQSE频率/DFT索引生成局部束，
                % 另反馈(B-1)个实测有效信道相对相位，故各方案不再退化为Type-I。
                haqse3_bits = 3*2*cfg.b_freq+codebook_phase_bits;
                dft3_bits = 3*index_bits+codebook_phase_bits;
                haqse5_bits = ternary(B>=5,5*2*cfg.b_freq+codebook_phase_bits,NaN);
                dft5_bits = ternary(B>=5,5*index_bits+codebook_phase_bits,NaN);
            end
            bits(is,i,:) = [0, ...
                B*(2*index_bits+cfg.b_typeII_amp+cfg.b_typeII_inner_phase)+codebook_phase_bits, ...
                B*index_bits+codebook_phase_bits, ...
                haqse3_bits, dft3_bits, haqse5_bits, dft5_bits];
        end
    end
end

function validate_direction_convention(cfg)
    % 无噪声几何闭环：用精确UE→AP方向参与定位，验证方向约定正确。
    [pb,Rb] = active_layout(cfg, cfg.B_full);
    refs = 1:3;
    rng_state = rng;
    rng(20261015, 'twister');
    errors = zeros(1,cfg.N_noiseless_validation);
    for it = 1:numel(errors)
        p_ue = cfg.origin_ue_position;
        qy = zeros(1,numel(refs));
        qz = zeros(1,numel(refs));
        for ii = 1:numel(refs)
            b = refs(ii);
            q = Rb(:,:,b)'*(pb(:,b)-p_ue);
            q = q/norm(q);
            qy(ii) = q(2);
            qz(ii) = q(3);
        end
        p_est = solve_ue_to_ap_rays(qy, qz, refs, pb, Rb);
        errors(it) = norm(p_est-p_ue);
    end
    rng(rng_state);
    max_error = max(errors);
    assert(max_error < 1e-9, '定位方向闭环验证失败，最大误差=%.3e m。', max_error);
    fprintf('无噪声方向闭环验证通过：固定UE的最大定位误差=%.3e m。\n', max_error);
end

function report = validate_haqse_noiseless(cfg)
    % 端到端无噪声自检：已知几何位置→生成LoS张量→HAQSE/DFT频率→定位。
    % 该检查能够发现外部MMV_CF函数的单位、符号或维度与本脚本约定不一致，
    % 而不仅仅验证几何射线方向是否正确。
    cfg_case = cfg;
    cfg_case.B = cfg.B_full;
    cfg_case.Nb = cfg.Nh*cfg.Nv;
    cfg_case.N = cfg_case.B*cfg_case.Nb;
    [pb,Rb] = active_layout(cfg_case, cfg_case.B);
    refs = 1:3;
    rng_state = rng;
    rng(20261016,'twister');
    count = cfg.N_noiseless_validation;
    pos_error = nan(1,count);
    freq_error_vs_dft = nan(1,count);
    raw_range_ok = true;
    raw_abs_max = 0;
    for it = 1:count
        p_ue = cfg_case.origin_ue_position;
        [H,dft_basis] = build_channel(cfg_case, p_ue, pb, Rb, 1);
        [p_h,f_h,f_raw] = estimate_position_haqse(H, refs, pb, Rb, cfg_case);
        [~,u_dft] = estimate_position_dft(H, refs, pb, Rb, cfg_case, dft_basis);
        raw_abs_max = max(raw_abs_max, max(abs(f_raw(:))));
        raw_range_ok = raw_range_ok && all(abs(f_raw(:)) <= 0.5+1e-9);
        pos_error(it) = norm(p_h-p_ue);
        freq_error_vs_dft(it) = sqrt(mean(abs(-f_h(:)-u_dft(:)).^2));
    end
    rng(rng_state);
    report.sample_count = count;
    report.mean_position_error_m = mean(pos_error);
    report.max_position_error_m = max(pos_error);
    report.mean_frequency_rms_vs_dft_cpe = mean(freq_error_vs_dft);
    report.raw_frequency_range_ok = raw_range_ok;
    report.raw_frequency_abs_max_cpe = raw_abs_max;
    fprintf(['HAQSE无噪声端到端检查：%d点，平均定位误差=%.3e m，最大误差=%.3e m，', ...
        '与DFT粗峰的统一频率RMS=%.3e cycles/element，原始|f|max=%.3f。\n'], ...
        count, report.mean_position_error_m, report.max_position_error_m, ...
        report.mean_frequency_rms_vs_dft_cpe, report.raw_frequency_abs_max_cpe);
    if ~raw_range_ok
        warning(['HAQSE原始频率超出[-0.5,0.5] cycles/element。请先核对', ...
            'MMV_CF_mQSE2_SumOfRatios的单位/符号，端点量化不能作为修复。']);
    end
end

function value = ternary(condition, yes_value, no_value)
    if condition, value = yes_value; else, value = no_value; end
end

function print_stable_reference_plan(cfg)
    fprintf('长期稳定参考AP计划（基于覆盖区先验与无噪声LoS）：\n');
    for ib = 1:numel(cfg.reference_plan.B)
        B = cfg.reference_plan.B(ib);
        refs3 = cfg.reference_plan.refs3{ib};
        fprintf('  B=%d，K=3：%s', B, mat2str(refs3));
        if B >= 5
            refs5 = cfg.reference_plan.refs5{ib};
            fprintf('；K=5：%s', mat2str(refs5));
        end
        fprintf('\n');
    end
end

function print_diagnostic_summary(scenario, summary)
    fprintf('\n%s：固定观测SNR下的诊断均值\n', scenario.name);
    fprintf('  Oracle K=3（真实位置代入同一数据流）归一化增益：%.4f\n', ...
        summary.mean_oracle_gain(1));
    fprintf('  Oracle K=5（真实位置代入同一数据流）归一化增益：%.4f\n', ...
        summary.mean_oracle_gain(2));
    for ii = 1:numel(summary.method_names)
        if ~isnan(summary.mean_position_error_m(ii))
            fprintf(['  %-15s：位置误差 %.4f m，频率RMS %.4e cycles/element，', ...
                '处理信道NMSE %.4f，固定参考增益 %.4f，相对G %.4f\n'], ...
                summary.method_names{ii}, ...
                summary.mean_position_error_m(ii), summary.mean_frequency_rms_cpe(ii), ...
                summary.mean_processing_channel_nmse(ii), summary.mean_link_gain(ii+3), ...
                summary.mean_normalized_gain(ii+3));
        end
    end
end

function plot_ap_rate_pair(x, data, method_names, scenarios, xlabel_text, title_text, destination)
    % 场景顺序固定为：纯LoS/AP侧、LoS+NLoS/AP侧、纯LoS/UE侧、LoS+NLoS/UE侧。
    fig = figure('Color','w','Position',[80,80,1360,940]);
    styles = {'k-*','b-o','c-d','r-s','m-^','g-v'};
    for is = 1:numel(scenarios)
        subplot(2,2,is); hold on;
        Y = data{is};
        for m = 1:7
            if m == 7
                plot(x,Y(m,:),'-p','Color',[0.45,0.25,0.75], ...
                    'LineWidth',1.2,'MarkerSize',5);
            else
                plot(x,Y(m,:),styles{m},'LineWidth',1.2,'MarkerSize',5);
            end
        end
        grid on;
        xlabel(xlabel_text);
        ylabel('平均吞吐量（Mbit/s）');
        title(scenarios(is).name);
        if is == 1, legend(method_names,'Location','best','FontSize',8); end
    end
    sgtitle(title_text);
    save_png(fig, destination);
end

function plot_ap_feedback_pair(cfg, bits, method_names, scenarios, destination)
    fig = figure('Color','w','Position',[80,80,1360,940]);
    for is = 1:numel(scenarios)
        subplot(2,2,is); hold on;
        Y = squeeze(bits(is,:,:));
        for m = 2:7
            plot(cfg.ap_list,Y(:,m),'-o','LineWidth',1.2,'MarkerSize',5);
        end
        grid on;
        xlabel('激活AP数');
        ylabel('用户侧反馈开销（bit）');
        title(scenarios(is).name);
        if is == 1, legend(method_names(2:7),'Location','northwest','FontSize',8); end
    end
    sgtitle(['任务2补充 图3 AP/UE侧激活AP数与反馈开销', ...
        '（反馈比特不受LoS/NLoS影响；UE侧含AP间相对相位字段）']);
    save_png(fig, destination);
end

function plot_stable_reference_deployment(cfg, destination)
    % 图4展示长期稳定参考计划：其选择仅由覆盖区先验、无噪声LoS和部署几何决定，
    % 在所有Monte-Carlo及两种信道场景中保持一致，因此可直接对应图2的比较条件。
    cfg_case = cfg;
    cfg_case.B = cfg.B_full;
    cfg_case.Nb = cfg.Nh*cfg.Nv;
    cfg_case.N = cfg_case.B*cfg_case.Nb;
    [pb,Rb] = active_layout(cfg_case, cfg_case.B);
    refs3 = stable_reference_aps(cfg, cfg.B_full, 3);
    refs5 = stable_reference_aps(cfg, cfg.B_full, 5);

    fig = figure('Color','w','Position',[100,100,1180,560]);
    subplot(1,2,1);
    draw_reference_deployment(cfg, pb, Rb, refs3, 3);
    subplot(1,2,2);
    draw_reference_deployment(cfg, pb, Rb, refs5, 5);
    sgtitle('原始物理部署与长期稳定自适应参考AP（覆盖区先验、无噪声LoS）');
    save_png(fig, destination);
end

function draw_reference_deployment(cfg, pb, Rb, refs, K)
    hold on;
    % 淡蓝色方块表示参考选择的覆盖区几何先验，中心为[25,25,0]^T m。
    patch([15,35,35,15], [15,15,35,35], [0.80,0.88,0.93], ...
        'FaceAlpha',0.45, 'EdgeColor','none', 'HandleVisibility','off');
    unselected = setdiff(1:cfg.B_full, refs);
    h_unselected = plot(pb(1,unselected), pb(2,unselected), 'o', ...
        'Color',[0.12,0.42,0.67], 'MarkerFaceColor',[0.12,0.42,0.67], ...
        'MarkerSize',7, 'LineStyle','none');
    h_selected = plot(pb(1,refs), pb(2,refs), 'o', ...
        'Color',[0.78,0.20,0.20], 'MarkerFaceColor',[0.78,0.20,0.20], ...
        'MarkerSize',8, 'LineStyle','none');
    % 原始部署中AP1/AP5/AP6以及AP3/AP7位置相邻。为使编号互不遮挡，
    % 为每个固定AP编号设置独立的显示偏移；文字只显示数字，不再重复“AP”。
    label_offset = [ 1.4,  2.0;  ... % AP1: [0,0]
                     1.2, -3.0;  ... % AP2: [30,-10]
                     2.0, -1.6;  ... % AP3: [0,50]
                     1.3,  1.8;  ... % AP4: [50,50]
                    -3.8, -2.5;  ... % AP5: [-1,-1]
                     2.1, -3.2;  ... % AP6: [0,-1]
                    -3.8,  1.8;  ... % AP7: [0,51]
                    -1.2,  2.0];     % AP8: [20,-10]
    h_arrow = gobjects(1,1);
    for b = 1:cfg.B_full
        direction = Rb(1:2,1,b);
        direction = direction/max(norm(direction), realmin);
        arrow_color = [0.12,0.42,0.67];
        if ismember(b, refs), arrow_color = [0.78,0.20,0.20]; end
        h = quiver(pb(1,b), pb(2,b), 3.6*direction(1), 3.6*direction(2), 0, ...
            'Color',arrow_color, 'LineWidth',1.35, 'MaxHeadSize',0.7);
        if b == 1, h_arrow = h; end
        text(pb(1,b)+label_offset(b,1), pb(2,b)+label_offset(b,2), sprintf('%d',b), ...
            'FontSize',10, 'Color',arrow_color, 'FontWeight','bold', ...
            'BackgroundColor','white', 'Margin',1, ...
            'HorizontalAlignment','center', 'VerticalAlignment','middle');
    end
    h_ue = plot(cfg.origin_ue_position(1), cfg.origin_ue_position(2), 'kp', ...
        'MarkerFaceColor',[1,0.78,0], 'MarkerSize',11);
    plot(cfg.reference_selection_prior(1), cfg.reference_selection_prior(2), 'kx', ...
        'MarkerSize',9, 'LineWidth',1.4, 'HandleVisibility','off');
    grid on; box on; axis equal;
    xlim([-8,58]); ylim([-15,58]);
    xlabel('x（m）'); ylabel('y（m）');
    title(sprintf('K=%d：稳定参考AP={%s}', K, strtrim(sprintf('%d ',refs))));
    legend([h_unselected,h_selected,h_ue,h_arrow], ...
        {'未选AP','稳定参考AP','固定UE','面板法向'}, ...
        'Location','southoutside','Orientation','horizontal','FontSize',9);
end

function save_png(fig, destination)
    tmp = [tempname(tempdir),'.png'];
    cleanup = onCleanup(@() delete_if_present(tmp)); %#ok<NASGU>
    exportgraphics(fig,tmp,'Resolution',160,'BackgroundColor','white');
    movefile(tmp,destination,'f');
end

function delete_if_present(file)
    if exist(file,'file'), delete(file); end
end
