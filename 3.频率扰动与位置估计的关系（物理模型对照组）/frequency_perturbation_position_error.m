%% 频率扰动与位置误差关系
% DFT-index 定位方案的物理模型对照组。
%
% 本脚本保持 2.对比方案/mc_beam_gain_vs_snr.m 中的几何模型和天线配置不变。
% 不添加信道噪声，而是在两个局部空间频率上施加扰动，
% 然后将无扰动和有扰动的频率都经过相同的 DFT-index 量化、方向反演，
% 最后计算位置误差。
%
% 本脚本中的空间频率约定：
%   f(:,1) = 局部 y 方向空间频率，单位为周/阵元
%   f(:,2) = 局部 z 方向空间频率，单位为周/阵元
% 对于原始模型，f = [q_y, q_z] * d/lambda_c，其中 q 是 AP 局部坐标系下
% 从 AP 指向 UE 的方向向量。

clear; clc; close all;
rng(20, 'twister');

cfg = default_config();
[pb, Rb, p_lcs, p_gcs, R_ue, p_ue_lcs, dft_hv, H_total_clean] = ...
    build_physical_model(cfg);

% 尽管本对照组进行位置估计时只需要 AP 中心位置和朝向，
% 这里仍完整构建整个几何模型，以明确保留原始天线布局和物理模型参数。
cfg.pb = pb;
cfg.Rb = Rb;
cfg.p_lcs = p_lcs;
cfg.p_gcs = p_gcs;
cfg.R_ue = R_ue;
cfg.p_ue_lcs = p_ue_lcs;
cfg.dft_hv = dft_hv;
cfg.H_total_clean = H_total_clean;

index_known = cfg.index_known(:).';
K = numel(index_known);

% 根据原始 AP/UE 几何关系计算真实局部空间频率。
f_true = true_local_frequency(pb, Rb, cfg.p_ue, index_known, ...
    cfg.d, cfg.lambda_c);

% 记录两种无噪声 DFT-index 基线：
%   1) 将真实频率量化到最近频率网格；
%   2) 使用原始 DFT 相关性搜索得到的实际索引。
% 第二种基线用于计算最终的位置误差。
[f_dft_grid, index_dft_grid] = quantize_frequency_to_dft_index(f_true, cfg);
index_dft = select_dft_index_from_clean_channel(H_total_clean, ...
    dft_hv, cfg, index_known);
f_dft = frequency_from_dft_index(index_dft, cfg);

% 分别使用连续真实频率和 DFT-index 量化频率进行位置估计。
[p_from_true_f, ~] = solve_position_from_frequency(f_true, pb, Rb, ...
    index_known, cfg.lambda_c, cfg.d);
[p_from_dft_index, ~] = solve_position_from_dft_index_exact(index_dft, ...
    pb, Rb, index_known, cfg.Nh, cfg.Nv, cfg.O1, cfg.O2, ...
    cfg.lambda_c, cfg.d);

f_error_dft = rms_all(f_dft - f_true);
f_error_dft_grid = rms_all(f_dft_grid - f_true);
position_error_true_f = norm(p_from_true_f - cfg.p_ue);
position_error_dft = norm(p_from_dft_index - cfg.p_ue);

fprintf('===== 物理模型 DFT-index 对照组 =====\n');
fprintf('AP总数 = %d，参考AP = %s，阵列规模 = %d x %d，载波频率 = %.3f GHz\n', ...
    cfg.B, mat2str(index_known), cfg.Nh, cfg.Nv, cfg.fc/1e9);
fprintf('真实用户位置                 = [% .6f, % .6f, % .6f]^T m\n', cfg.p_ue);
fprintf('使用真实连续频率估计的位置   = [% .6f, % .6f, % .6f]^T m\n', p_from_true_f);
fprintf('真实连续频率对应的位置误差   = %.6e m\n', position_error_true_f);
fprintf('DFT网格步长 (y,z)            = [%.6e, %.6e] 周/阵元\n', ...
    1/(cfg.Nh*cfg.O2), 1/(cfg.Nv*cfg.O1));
fprintf('最近网格频率RMS误差           = %.6e 周/阵元\n', ...
    f_error_dft_grid);
fprintf('DFT-index频率RMS误差          = %.6e 周/阵元\n', f_error_dft);
fprintf('DFT-index位置估计误差         = %.6e m\n\n', position_error_dft);

fprintf('参考AP真实空间频率 [f_y, f_z]：\n');
disp(f_true);
fprintf('参考AP的DFT空间频率 [f_y, f_z]：\n');
disp(f_dft);
fprintf('参考AP的DFT索引：\n');
disp(index_dft);
fprintf('最近网格DFT索引（仅用于诊断）：\n');
disp(index_dft_grid);

% 设置受控的频率扰动等级。对于每次试验，delta_f 表示所有参考 AP、
% 两个空间频率方向合并计算得到的 RMS 误差。
cfg.delta_f_min = 1e-5;
cfg.delta_f_max = 2 * max(1/(cfg.Nh*cfg.O2), 1/(cfg.Nv*cfg.O1));
cfg.N_delta = 25;
cfg.N_MC = 200;
delta_f_list = [0; logspace(log10(cfg.delta_f_min), ...
    log10(cfg.delta_f_max), cfg.N_delta-1).'];

% frequency_error_rms 是经过 DFT-index 量化后的有效频率误差。
% input_frequency_error_rms 是经过 DFT-index 量化前人为加入的频率扰动误差。
% 位置估计始终使用量化后的 DFT-index 频率。
frequency_error_rms = zeros(numel(delta_f_list), cfg.N_MC);
input_frequency_error_rms = zeros(numel(delta_f_list), cfg.N_MC);
position_error = zeros(numel(delta_f_list), cfg.N_MC);
position_estimates = zeros(3, numel(delta_f_list), cfg.N_MC);

for id = 1:numel(delta_f_list)
    delta_f = delta_f_list(id);
    for imc = 1:cfg.N_MC
        if delta_f == 0
            f_target = f_true;
            % 零扰动基线使用实际无噪声 DFT-index 搜索结果。
            f_indexed = f_dft;
            index_current = index_dft;
        else
            % 施加 RMS 值严格为 delta_f 的各向同性扰动。
            % 整个循环中不添加信道噪声或测量噪声。
            e = randn(K, 2);
            e = e / rms_all(e);
            f_target = f_true + delta_f * e;
            % 人为加入的扰动也必须经过 DFT-index 方案，
            % 这是本物理模型对照组的核心操作。
            [f_indexed, index_current] = ...
                quantize_frequency_to_dft_index(f_target, cfg);
        end

        [p_est, ~] = solve_position_from_dft_index_exact(index_current, ...
            pb, Rb, index_known, cfg.Nh, cfg.Nv, cfg.O1, cfg.O2, ...
            cfg.lambda_c, cfg.d);
        input_frequency_error_rms(id, imc) = rms_all(f_target - f_true);
        frequency_error_rms(id, imc) = rms_all(f_indexed - f_true);
        position_error(id, imc) = norm(p_est - cfg.p_ue);
        position_estimates(:, id, imc) = p_est;
    end
end

% 计算绘图所需的统计量，并保存结果，便于后续与 DFT-index 波束增益曲线比较。
mean_position_error = mean(position_error, 2);
median_position_error = median(position_error, 2);
p95_position_error = zeros(numel(delta_f_list), 1);
for id = 1:numel(delta_f_list)
    p95_position_error(id) = percentile_linear(position_error(id, :), 95);
end

% 不同 DFT 量化分辨率下的真实频率位置偏差对照。
% 该对照不加入随机频率扰动，只对同一个真实 f_true 使用不同 DFT 网格量化。
resolution_O_list = [1, 2, 4, 8, 16];
resolution_bits_per_dim = log2(cfg.Nh * resolution_O_list);
resolution_grid_step = 1 ./ (cfg.Nh * resolution_O_list);
resolution_frequency_error_rms = zeros(size(resolution_O_list));
resolution_position_error = zeros(size(resolution_O_list));
resolution_f_quantized = cell(size(resolution_O_list));
resolution_index = cell(size(resolution_O_list));

for ir = 1:numel(resolution_O_list)
    O_resolution = resolution_O_list(ir);
    [f_quantized, index_quantized] = ...
        quantize_frequency_to_dft_index(f_true, cfg, ...
        O_resolution, O_resolution);
    p_quantized = solve_position_from_dft_index_exact(index_quantized, ...
        pb, Rb, index_known, cfg.Nh, cfg.Nv, ...
        O_resolution, O_resolution, cfg.lambda_c, cfg.d);

    resolution_f_quantized{ir} = f_quantized;
    resolution_index{ir} = index_quantized;
    resolution_frequency_error_rms(ir) = rms_all(f_quantized - f_true);
    resolution_position_error(ir) = norm(p_quantized - cfg.p_ue);
end

fprintf('\n===== 不同 DFT 量化分辨率对照 =====\n');
fprintf('真实空间频率完美估计的位置误差 = %.6e m\n', ...
    position_error_true_f);
fprintf('过采样倍数 | 每方向量化比特数 | 网格步长 | 频率RMS误差 | 位置估计偏差\n');
for ir = 1:numel(resolution_O_list)
    fprintf('O = %2d      | %8.0f bit       | %.6e | %.6e | %.6e m\n', ...
        resolution_O_list(ir), resolution_bits_per_dim(ir), ...
        resolution_grid_step(ir), resolution_frequency_error_rms(ir), ...
        resolution_position_error(ir));
end

results = struct();
results.cfg = cfg;
results.index_known = index_known;
results.f_true = f_true;
results.f_dft = f_dft;
results.index_dft = index_dft;
results.f_dft_grid = f_dft_grid;
results.index_dft_grid = index_dft_grid;
results.f_error_dft = f_error_dft;
results.f_error_dft_grid = f_error_dft_grid;
results.p_from_true_f = p_from_true_f;
results.p_from_dft_index = p_from_dft_index;
results.position_error_true_f = position_error_true_f;
results.position_error_dft = position_error_dft;
results.delta_f_list = delta_f_list;
results.frequency_error_rms = frequency_error_rms;
results.input_frequency_error_rms = input_frequency_error_rms;
results.position_error = position_error;
results.position_estimates = position_estimates;
results.mean_position_error = mean_position_error;
results.median_position_error = median_position_error;
results.p95_position_error = p95_position_error;
results.resolution_O_list = resolution_O_list;
results.resolution_bits_per_dim = resolution_bits_per_dim;
results.resolution_grid_step = resolution_grid_step;
results.resolution_frequency_error_rms = resolution_frequency_error_rms;
results.resolution_position_error = resolution_position_error;
results.resolution_f_quantized = resolution_f_quantized;
results.resolution_index = resolution_index;

this_dir = fileparts(mfilename('fullpath'));
if isempty(this_dir)
    this_dir = pwd;
end
save(fullfile(this_dir, 'frequency_position_error_results.mat'), 'results');

% 图 1：所有经过 DFT-index 处理的受控扰动试验结果。
% 黄色星形表示无噪声、无扰动的 DFT-index 基线。
fig1 = figure('Color', 'w');
nonzero = delta_f_list > 0;
h_trials_all = loglog(frequency_error_rms(nonzero, :), position_error(nonzero, :), ...
    '.', 'Color', [0.20 0.45 0.85], 'MarkerSize', 8);
h_trials = h_trials_all(1);  % 取一个蓝色散点作为图例代表。
hold on;
h_baseline = loglog(f_error_dft, position_error_dft, 'kp', ...
    'MarkerSize', 12, 'MarkerFaceColor', [1.0 0.75 0.1]);
grid on;
xlabel('Effective spatial-frequency error after DFT-index (cycles/element)');
ylabel('Position error (m)');
title('DFT-index control group with frequency perturbation');
legend([h_trials, h_baseline], {'Perturbed f -> DFT-index -> position', ...
    'No-noise DFT-index baseline'}, ...
    'Location', 'northwest');
save_png_robust(fig1, fullfile(this_dir, 'frequency_error_vs_position_error_scatter.png'));

% 图 2：不同 DFT 量化分辨率下，真实频率量化误差与位置偏差的关系。
% 该图不加入人为扰动；每个点对应一个不同的 DFT 过采样倍数。
fig2 = figure('Color', 'w');
h_quantization = plot(resolution_frequency_error_rms, ...
    resolution_position_error, 'r-o', 'LineWidth', 1.4, ...
    'MarkerSize', 6, 'MarkerFaceColor', 'w');
hold on;
h_dft_baseline = plot(f_error_dft, position_error_dft, 'kp', ...
    'MarkerSize', 12, 'MarkerFaceColor', [1.0 0.75 0.1]);
h_perfect = plot(0, 0, 'g*', 'MarkerSize', 12, ...
    'LineWidth', 1.4);
for ir = 1:numel(resolution_O_list)
    text(resolution_frequency_error_rms(ir), ...
        resolution_position_error(ir), ...
        sprintf('  O=%d, %d bit', resolution_O_list(ir), ...
        resolution_bits_per_dim(ir)), 'FontSize', 8);
end
text(0.02 * max(resolution_frequency_error_rms), ...
    0.05 * max(resolution_position_error), ...
    'Perfect continuous f: (0, 0)', 'Color', [0.0 0.45 0.0]);
grid on;
xlabel('True-frequency quantization RMS error (cycles/element)');
ylabel('Position estimation bias (m)');
title('Position bias versus DFT quantization resolution');
legend([h_quantization, h_dft_baseline, h_perfect], ...
    {'Nearest-grid quantization', ...
    'Actual no-noise DFT-index baseline (O=4)', ...
    'Perfect continuous-frequency estimate'}, ...
    'Location', 'northwest');
xlim([0, 1.15 * max(resolution_frequency_error_rms)]);
ylim([0, 1.15 * max(resolution_position_error)]);
save_png_robust(fig2, fullfile(this_dir, 'frequency_error_vs_position_error_trend.png'));

% 图 3：显示无扰动 DFT 基线下三个参考 AP 的频率误差，
% 用于判断是否有某个 AP 主导了位置误差。
fig3 = figure('Color', 'w');
bar(1:K, abs(f_dft - f_true));
grid on;
xlabel('Reference AP index in index\_known');
ylabel('|DFT frequency error| (cycles/element)');
title('Per-AP DFT-index frequency quantization error');
xticks(1:K);
xticklabels(string(index_known));
legend({'f_y error', 'f_z error'}, 'Location', 'northwest');
save_png_robust(fig3, fullfile(this_dir, 'per_ap_dft_frequency_error.png'));

fprintf('结果数据和图片已保存到：\n%s\n', this_dir);

%% 本地函数

function cfg = default_config()
    % 以下参数与对比方案文件夹中的 mc_beam_gain_vs_snr.m 保持一致。
    cfg.B = 8;
    cfg.index_known = 1:3;
    cfg.Nh = 16;
    cfg.Nv = 16;
    cfg.Nb = cfg.Nh * cfg.Nv;
    cfg.N = cfg.B * cfg.Nb;
    cfg.Nu = 4;
    cfg.Ns = 1;
    cfg.L_path = 1;
    cfg.fc = 7e9;
    cfg.c = 3e8;
    cfg.lambda_c = cfg.c / cfg.fc;
    cfg.d = cfg.lambda_c / 2;
    cfg.k_c = 2*pi / cfg.lambda_c;
    cfg.O1 = 4;
    cfg.O2 = 4;
    cfg.p_ue = [25; 25; 0];
end

function [pb, Rb, p_lcs, p_gcs, R_ue, p_ue_lcs, dft_hv, H_total_clean] = ...
        build_physical_model(cfg)
    B = cfg.B;
    Nh = cfg.Nh;
    Nv = cfg.Nv;
    Nb = cfg.Nb;
    d = cfg.d;

    pb = zeros(3, B);
    pb(:,1) = [0; 0; 0];
    pb(:,2) = [30; -10; 0];
    pb(:,3) = [0; 50; 0];
    pb(:,4) = [50; 50; 0];
    pb(:,5) = [-1; -1; 0];
    pb(:,6) = [0; -1; 0];
    pb(:,7) = [0; 51; 0];
    pb(:,8) = [20; -10; 0];

    alpha_b = [0; 130; -45; -125; 49; 43; 98; 90] / 180*pi;
    beta_b = [0; 10; -5; 2; 6; -4; 8; -7] / 180*pi;
    gamma_b = zeros(B, 1);

    Rb = zeros(3, 3, B);
    for b = 1:B
        ca = cos(alpha_b(b)); sa = sin(alpha_b(b));
        cb = cos(beta_b(b)); sb = sin(beta_b(b));
        cg = cos(gamma_b(b)); sg = sin(gamma_b(b));
        Rb(:,:,b) = [ca*cb, ca*sb*sg-sa*cg, ca*sb*cg+sa*sg; ...
            sa*cb, sa*sb*sg+ca*cg, sa*sb*cg-ca*sg; ...
            -sb, cb*sg, cb*cg];
    end

    p_lcs = zeros(3, Nb);
    p_gcs = zeros(3, Nb, B);
    for b = 1:B
        idx_local = 0;
        for nh = 1:Nh
            for nv = 1:Nv
                y = ((nv-1) - (Nv-1)/2) * d;
                z = ((nh-1) - (Nh-1)/2) * d;
                idx_local = idx_local + 1;
                p_lcs(:,idx_local) = [0; y; z];
                p_gcs(:,idx_local,b) = Rb(:,:,b) * p_lcs(:,idx_local) + pb(:,b);
            end
        end
    end

    R_ue = eye(3);
    ue_offset = ((0:cfg.Nu-1) - (cfg.Nu-1)/2) * d;
    p_ue_lcs = zeros(3, cfg.Nu);
    p_ue_lcs(1,:) = ue_offset;

    dft_h = exp(1j*2*pi*(0:Nh-1).'*(0:cfg.O2*Nh-1)/(cfg.O2*Nh)) / sqrt(Nh);
    dft_v = exp(1j*2*pi*(0:Nv-1).'*(0:cfg.O1*Nv-1)/(cfg.O1*Nv)) / sqrt(Nv);
    dft_hv = kron(dft_h, dft_v);

    % 根据对比方案中的几何模型生成精确的无噪声 LOS 信道。
    H_L = zeros(cfg.Nu, Nb, B);
    H_total_clean = zeros(cfg.Nu, cfg.N);
    for b = 1:B
        d_bu = pb(:,b) - cfg.p_ue;       % 用户 -> AP
        r_bu = norm(d_bu);
        d_hat = d_bu / r_bu;
        d_hat_bLCS = Rb(:,:,b).' * d_hat;
        d_hat_uLCS = R_ue.' * d_hat;
        phase_center = exp(-1j*cfg.k_c*r_bu) / r_bu;
        for m = 1:cfg.Nu
            phase_UE = exp(1j*cfg.k_c * ...
                (d_hat_uLCS.' * p_ue_lcs(:,m)));
            for n = 1:Nb
                phase_AP = exp(-1j*cfg.k_c * ...
                    (d_hat_bLCS.' * p_lcs(:,n)));
                H_L(m,n,b) = phase_center * phase_AP * phase_UE;
            end
        end
        idx_b = (b-1)*Nb + (1:Nb);
        H_total_clean(:,idx_b) = H_L(:,:,b);
    end
end

function f_true = true_local_frequency(pb, Rb, p_ue, index_known, d, lambda_c)
    K = numel(index_known);
    f_true = zeros(K, 2);
    for ii = 1:K
        b = index_known(ii);
        ap_to_ue = p_ue - pb(:,b);
        ap_to_ue = ap_to_ue / norm(ap_to_ue);
        q_local = Rb(:,:,b).' * ap_to_ue;
        f_true(ii,:) = [q_local(2), q_local(3)] * d / lambda_c;
    end
end

function [f_dft, index_dft] = quantize_frequency_to_dft_index( ...
        f_true, cfg, O1, O2)
    if nargin < 3 || isempty(O1)
        O1 = cfg.O1;
    end
    if nargin < 4 || isempty(O2)
        O2 = cfg.O2;
    end

    NhO2 = cfg.Nh * O2;
    NvO1 = cfg.Nv * O1;
    f_dft = zeros(size(f_true));
    index_dft = zeros(size(f_true, 1), 1);

    for ii = 1:size(f_true, 1)
        kh = round(f_true(ii,1) * NhO2);
        kv = round(f_true(ii,2) * NvO1);
        kh_mod = mod(kh, NhO2);
        kv_mod = mod(kv, NvO1);
        f_dft(ii,1) = kh / NhO2;
        f_dft(ii,2) = kv / NvO1;
        % 保留原 solve_position_from_dft_index.m 中的索引拆分方式：
        % 内层索引解释为局部 y 频率，外层索引解释为局部 z 频率。
        % 原始物理模型中两个索引跨度均为 64，因此这与对比代码使用的
        % 线性索引约定在数值上相同，同时保留 y/z 的物理含义。
        index_dft(ii) = kv_mod * (NhO2) + kh_mod + 1;
    end
end

function index_dft = select_dft_index_from_clean_channel(H_total_clean, ...
        dft_hv, cfg, index_known)
    % 使用对比代码中的方法，精确选择无噪声条件下的 DFT-index。
    index_dft = zeros(1, cfg.B);
    for ii = 1:numel(index_known)
        b = index_known(ii);
        idx_b = (b-1)*cfg.Nb + (1:cfg.Nb);
        H_b = H_total_clean(:,idx_b);
        [~,~,V_b] = svd(H_b, 'econ');
        v_b = V_b(:,1);
        [~, index_dft(b)] = max(abs(v_b.' * dft_hv) / ...
            max(norm(v_b), realmin));
    end
    index_dft = index_dft(index_known);
end

function f = frequency_from_dft_index(index_dft, cfg)
    % 对原 solve_position_from_dft_index.m 索引拆分方式进行反变换。
    index_dft = index_dft(:).';
    index_dfty = mod(index_dft - 1, cfg.Nh * cfg.O2);
    index_dftz = floor((index_dft - 1) / (cfg.Nv * cfg.O1));
    mask_y = index_dfty / (cfg.Nh * cfg.O2) > 0.5;
    mask_z = index_dftz / (cfg.Nv * cfg.O1) > 0.5;
    index_dfty(mask_y) = index_dfty(mask_y) - cfg.Nh * cfg.O2;
    index_dftz(mask_z) = index_dftz(mask_z) - cfg.Nv * cfg.O1;
    f = [index_dfty / (cfg.Nh * cfg.O2); ...
        index_dftz / (cfg.Nv * cfg.O1)].';
end

function [p_best, info] = solve_position_from_frequency(f, pb, Rb, ...
        index_known, lambda_c, d)
    % 根据 DFT-index 角度约定，将连续空间频率反演为方向。
    % 使用与原 solve_position_from_dft_index.m 相同的局部方向模型，
    % 但允许直接输入受到扰动的频率 f。
    K = numel(index_known);
    A = zeros(3,3);
    rhs = zeros(3,1);
    q_global = zeros(3, K);

    for ii = 1:K
        b = index_known(ii);
        q_y = f(ii,1) * lambda_c / d;
        q_z = f(ii,2) * lambda_c / d;
        q_y = max(-1, min(1, q_y));
        q_z = max(-1, min(1, q_z));
        q_x = sqrt(max(0, 1 - q_y^2 - q_z^2));
        q_local = [q_x; q_y; q_z];
        q_global(:,ii) = Rb(:,:,b) * q_local;
        q_global(:,ii) = q_global(:,ii) / norm(q_global(:,ii));

        P = eye(3) - q_global(:,ii) * q_global(:,ii).';
        A = A + P;
        rhs = rhs + P * pb(:,b);
    end

    p_best = A \ rhs;
    info.q_global = q_global;
    info.rcond_A = rcond(A);
end

function [p_best, info] = solve_position_from_dft_index_exact(index_dft, ...
        pb, Rb, index_known, Nh, Nv, O1, O2, lambda_c, d)
    % 这是原始 DFT-index 反演映射的本地实现，
    % 使本对照组代码自包含，同时不修改原始源文件。
    index_dft = index_dft(:).';
    index_dft_known = index_dft(index_known);
    B_known = numel(index_known);

    index_dfty = mod(index_dft_known - 1, Nh * O2);
    index_dftz = floor((index_dft_known - 1) / (Nv * O1));

    neg_y = index_dfty / (Nh * O2) > 0.5;
    neg_z = index_dftz / (Nv * O1) > 0.5;
    index_dfty(neg_y) = index_dfty(neg_y) - Nh * O2;
    index_dftz(neg_z) = index_dftz(neg_z) - Nv * O1;

    theta_est = acos(index_dftz / (Nv * O1) * (lambda_c / d));
    phi_est = asin(index_dfty / (Nh * O2) * (lambda_c / d) ./ sin(theta_est));
    theta_est = reshape(theta_est, B_known, 1);
    phi_est = reshape(phi_est, B_known, 1);

    q_local = [sin(theta_est).*cos(phi_est), ...
        sin(theta_est).*sin(phi_est), cos(theta_est)].';

    A = zeros(3,3);
    rhs = zeros(3,1);
    q_global = zeros(3, B_known);
    for ii = 1:B_known
        b = index_known(ii);
        q_global(:,ii) = Rb(:,:,b) * q_local(:,ii);
        q_global(:,ii) = q_global(:,ii) / norm(q_global(:,ii));
        P = eye(3) - q_global(:,ii) * q_global(:,ii).';
        A = A + P;
        rhs = rhs + P * pb(:,b);
    end

    p_best = A \ rhs;
    info.q_local = q_local;
    info.q_global = q_global;
    info.rcond_A = rcond(A);
end

function value = rms_all(x)
    value = sqrt(mean(abs(x(:)).^2));
end

function value = percentile_linear(x, pct)
    x = sort(x(:));
    if isempty(x)
        value = NaN;
        return;
    end
    if numel(x) == 1
        value = x(1);
        return;
    end
    pos = 1 + (numel(x)-1) * pct / 100;
    lo = floor(pos);
    hi = ceil(pos);
    if lo == hi
        value = x(lo);
    else
        value = x(lo) + (pos-lo) * (x(hi)-x(lo));
    end
end

function save_png_robust(fig, destination)
% 先保存到唯一的 ASCII 临时路径，再移动到指定输出路径。
% 这样可以避免非 ASCII 路径或图片查看器占用旧 PNG 文件时，
% MATLAB 的 print/saveas 保存失败。
destination = char(destination);
output_dir = fileparts(destination);
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

drawnow;
temp_file = [tempname(tempdir), '.png'];
cleanup_obj = onCleanup(@() delete_if_exists(temp_file)); %#ok<NASGU>

try
    exportgraphics(fig, temp_file, 'Resolution', 150, ...
        'BackgroundColor', 'white');
catch export_err
    % 兼容不支持 exportgraphics 的旧版本 MATLAB。
    try
        print(fig, temp_file, '-dpng', '-r150');
    catch print_err
        error('frequency_perturbation_position_error:SavePngFailed', ...
            ['无法保存图片。exportgraphics 报错：%s；print 报错：%s'], ...
            export_err.message, print_err.message);
    end
end

[ok, msg] = movefile(temp_file, destination, 'f');
if ~ok
    % 如果旧目标文件正被图片查看器占用，则保留新生成的结果，
    % 不让整个仿真中断。
    [output_dir, output_name, output_ext] = fileparts(destination);
    fallback = fullfile(output_dir, sprintf('%s_%s%s', output_name, ...
        datestr(now, 'yyyymmdd_HHMMSSFFF'), output_ext));
    [ok_fallback, fallback_msg] = movefile(temp_file, fallback, 'f');
    if ok_fallback
        warning('frequency_perturbation_position_error:SavePngFallback', ...
            ['目标 PNG 无法覆盖，已改存为：%s。原错误：%s'], ...
            fallback, msg);
    else
        error('frequency_perturbation_position_error:MovePngFailed', ...
            ['临时 PNG 已生成，但无法写入目标文件：%s (%s)；备用文件也失败：%s'], ...
            destination, msg, fallback_msg);
    end
end
end

function delete_if_exists(filename)
if exist(filename, 'file')
    delete(filename);
end
end
