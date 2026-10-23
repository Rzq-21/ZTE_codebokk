%% Monte-Carlo: equal-amplitude continuous phase versus 6-bit phase
clear;clc;close all;
rng(20);

% Monte-Carlo comparison: same geometry, independent channel-estimation noise per trial.
N_MC = 100;
gain_cont_mc = nan(N_MC,1);       % continuous phase, equal AP amplitudes
gain_quant_mc = nan(N_MC,1);      % 6-bit phase, equal AP amplitudes
gain_diff_mc = nan(N_MC,1);       % continuous minus 6-bit, linear gain
gain_diff_db_mc = nan(N_MC,1);    % continuous minus 6-bit, dB gain

for mc = 1:N_MC

%% 参数设置
B = 8;                          % UPA/基站数
index_known = 1:3;              % 用于定位的 UPA 索引

Nh = 16;                         % UPA 水平阵元数，对应 nh / 局部 z 方向
Nv = 16;                         % UPA 垂直阵元数，对应 nv / 局部 y 方向
Nb = Nh * Nv;                   % 每个 UPA 的阵元数
N = B * Nb;                     % 全局虚拟发射阵列阵元数

Nu = 4;                         % 用户 ULA 阵元数，同时作为 MMV 快拍数
Ns = 1;                         % 当前按单流选择一个 DFT 码字
L_path = 1;                     % 当前纯 LOS 只有一条空间路径

fc = 7e9;
c = 3e8;
lambda_c = c / fc;
d = lambda_c / 2;
k_c = 2*pi / lambda_c;

Qopt = 3;                       % 二维 HAQSE 细化迭代次数
Fs_spatial_y = 1;               % 空间采样归一化：cycles/element
Fs_spatial_z = 1;

% 二维过采样 DFT 码本。码字排列为 kron(dft_h,dft_v)，垂直索引变化更快。
O1 = 4;                         % 垂直方向过采样因子
O2 = 4;                         % 水平方向过采样因子

dft_h = exp(1j*2*pi*(0:Nh-1).'*(0:O2*Nh-1)/(O2*Nh)) / sqrt(Nh);% 水平
dft_v = exp(1j*2*pi*(0:Nv-1).'*(0:O1*Nv-1)/(O1*Nv)) / sqrt(Nv);% 垂直
dft_hv = kron(dft_h,dft_v);

color = lines(10);
%% UPA 中心位置
pb = zeros(3,B);
% AP分布位置很接近
pb(:,1) = [ 0;   0; 0];pb(:,2) = [30; -10; 0];pb(:,3) = [20; -10; 0];pb(:,4) = [50;  50; 0];
pb(:,5) = [-1;  -1; 0];pb(:,6) = [ 0;  -1; 0];pb(:,7) = [ 0;  51; 0];pb(:,8) = [ 0;  50; 0];

% % AP分布分散
% pb(:,1) = [ 52;  90; 0];pb(:,2) = [ 52; -40; 0];pb(:,3) = [ -2;  90; 0];pb(:,4) = [ 90;  52; 0];
% pb(:,5) = [-40;  -2; 0];pb(:,6) = [ -2; -40; 0];pb(:,7) = [-40;  52; 0];pb(:,8) = [ 90;  -2; 0];
%% UPA 局部坐标系到全局坐标系的旋转矩阵
% alpha_b = [0; 130; -45; -125] / 180 * pi;
% beta_b  = [0; 10; -5; 2] / 180 * pi;
alpha_b = [0; 130; -45; -125;    49;     43;    98;     90] / 180 * pi;
beta_b  = [0; 10; -5; 2;    6;    -4;    8;    -7] / 180 * pi;
gamma_b = zeros(B,1);

Rb = zeros(3,3,B);
for b = 1:B
    ca = cos(alpha_b(b)); sa = sin(alpha_b(b));
    cb = cos(beta_b(b));  sb = sin(beta_b(b));
    cg = cos(gamma_b(b)); sg = sin(gamma_b(b));

    Rb(:,:,b) = [ ca*cb, ca*sb*sg-sa*cg, ca*sb*cg+sa*sg; ...
        sa*cb, sa*sb*sg+ca*cg, sa*sb*cg-ca*sg; ...
        -sb,   cb*sg,          cb*cg];
end

%% UPA 阵元位置
p_lcs = zeros(3,Nb);
p_gcs = zeros(3,Nb,B);
p_virtual = zeros(3,N);

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

    idx_b = (b-1)*Nb + (1:Nb);
    p_virtual(:,idx_b) = p_gcs(:,:,b);
end

%% 一维阵元索引映射
n0 = zeros(1,N);
for b = 1:B
    for nh = 1:Nh
        for nv = 1:Nv
            n_local = (nh-1)*Nv + nv;
            n_global = (b-1)*Nb + n_local;
            n0(n_global) = n_global;
        end
    end
end

%% 用户位置和用户 ULA
p_ue = [25; 25; 10];

alpha_ue = 0 / 180 * pi;
beta_ue  = 0 / 180 * pi;
gamma_ue = 0 / 180 * pi;

ca = cos(alpha_ue); sa = sin(alpha_ue);
cb = cos(beta_ue);  sb = sin(beta_ue);
cg = cos(gamma_ue); sg = sin(gamma_ue);

R_ue = [ ca*cb, ca*sb*sg-sa*cg, ca*sb*cg+sa*sg; ...
    sa*cb, sa*sb*sg+ca*cg, sa*sb*cg-ca*sg; ...
    -sb,   cb*sg,          cb*cg];

ue_offset = ((0:Nu-1) - (Nu-1)/2) * d;
p_ue_lcs = zeros(3,Nu);
p_ue_lcs(1,:) = ue_offset;    % 用户 ULA 沿局部 x 轴排列
p_ue_gcs = R_ue * p_ue_lcs + p_ue;

%% 纯 LOS 信道
H_L = zeros(Nu,Nb,B);
H_LOS_DL = zeros(Nu,N);

d_bu = zeros(3,B);
r_GCS_bu = zeros(1,B);
d_tilde_bu = zeros(3,B);
d_tilde_bu_bLCS = zeros(3,B);
d_tilde_bu_uLCS = zeros(3,B);

for b = 1:B
    % 该方向是 UE -> AP
    d_bu(:,b) = pb(:,b) - p_ue;
    r_GCS_bu(b) = norm(d_bu(:,b));
    d_tilde_bu(:,b) = d_bu(:,b) / r_GCS_bu(b);

    % 全局坐标方向转换到局部坐标：使用旋转矩阵转置
    d_tilde_bu_bLCS(:,b) = Rb(:,:,b).' * d_tilde_bu(:,b);
    d_tilde_bu_uLCS(:,b) = R_ue.' * d_tilde_bu(:,b);

    alpha_bu = 1 / r_GCS_bu(b);
    phase_center = exp(-1j*k_c*r_GCS_bu(b));

    for m = 1:Nu
        phase_UE = exp(1j*k_c * (d_tilde_bu_uLCS(:,b).' * p_ue_lcs(:,m)));

        for n = 1:Nb
            phase_AP = exp(-1j*k_c * (d_tilde_bu_bLCS(:,b).' * p_lcs(:,n)));

            H_L(m,n,b) = alpha_bu * phase_center  * phase_AP * phase_UE;
        end
    end
end

for b = 1:B
    idx_b = (b-1)*Nb + (1:Nb);
    H_LOS_DL(:,idx_b) = H_L(:,:,b);
end

H_total_clean = H_LOS_DL;

% 给信道估计添加复高斯白噪声
SNR_H_dB = 30;        % 信道估计信噪比
signal_power = mean(abs(H_total_clean(:)).^2);% 平均信号功率
noise_power = signal_power / 10^(SNR_H_dB/10);% 根据目标 SNR 计算噪声功率
H_noise = sqrt(noise_power/2) * (randn(size(H_total_clean)) + 1j*randn(size(H_total_clean)));% 生成复高斯白噪声
H_total = H_total_clean + H_noise;% 带噪声的信道估计

%% 每个 UPA 的二维空间 HAQSE 估计
% nu_v 对应 p_lcs 的 y/nv 方向；nu_h 对应 p_lcs 的 z/nh 方向。
nu_v_hat = zeros(1,B);nu_h_hat = zeros(1,B);
q_y_hat = zeros(1,B);q_z_hat = zeros(1,B);q_x_abs_hat = zeros(1,B);
index_dft_feedback = nan(1,B);

% 实际反馈场景只对已反馈 AP 运行 HAQSE；其余 AP 的方向和信道均未知。
for b = index_known
    idx_b = (b-1)*Nb + (1:Nb);
    H_b = H_total(:,idx_b);          % Nu×Nb = 4×16

    S_mat = H_b.';                   % 外部采用 Nb×Nu 矩阵；Nb×Nu = 16×4
    S_tensor = reshape(S_mat,[Nv,Nh,Nu]); % 送入二维 HAQSE 时恢复 UPA 的二维结构

    estimated_freqs = MMV_CF_mQSE2_SumOfRatios(S_tensor,L_path,Fs_spatial_y,Fs_spatial_z,Qopt);
    
    nu_v_hat(b) = estimated_freqs(1,1);% 垂直方向空间频率
    nu_h_hat(b) = estimated_freqs(1,2);% 水平方向空间频率
    b_nu = 14;
    Q = 2^b_nu;
    nu_v_hat(b) = -0.5 + round((nu_v_hat(b) + 0.5) * (Q-1)) / (Q-1);
    nu_h_hat(b) = -0.5 + round((nu_h_hat(b) + 0.5) * (Q-1)) / (Q-1);


    % 空间频率 cycles/element -> 局部方向余弦
    q_y_hat(b) = -nu_v_hat(b) * lambda_c / d;% 用户到UPA方向在 UPA 局部 y 轴的方向余弦
    q_z_hat(b) = -nu_h_hat(b) * lambda_c / d;% 用户到UPA方向在 UPA 局部 z 轴的方向余弦

    % 数值噪声可能使 q_y^2+q_z^2 略大于 1，这里做保护
    q_y_hat(b) = max(-1,min(1,q_y_hat(b)));
    q_z_hat(b) = max(-1,min(1,q_z_hat(b)));
    q_x_abs_hat(b) = sqrt(max(0,1-q_y_hat(b)^2-q_z_hat(b)^2));

    if mc == 1
        fprintf(['UPA %d: nu_v=%.8f, nu_h=%.8f, ', 'q_y=%.8f, q_z=%.8f\n'], ...
                 b,nu_v_hat(b),nu_h_hat(b),q_y_hat(b),q_z_hat(b));

        fprintf([' true UE->AP local q_y=%.8f, q_z=%.8f, ', 'errors=(%.3e, %.3e)\n'], ...
                 d_tilde_bu_bLCS(2,b),d_tilde_bu_bLCS(3,b), ...
                 q_y_hat(b)-d_tilde_bu_bLCS(2,b), ...
                 q_z_hat(b)-d_tilde_bu_bLCS(3,b));
    end
end

%% 多 UPA 位置交汇 + 估计用户位置
[p_ue_est,best_signs,best_score,t_best, q_local_selected,q_global_selected] = ...
solve_position_front_back(pb,Rb,index_known, q_y_hat,q_z_hat,q_x_abs_hat);

position_error = norm(p_ue_est-p_ue);

if mc == 1
    fprintf('\n===== HAQSE 用户位置估计 =====\n');
    fprintf('\n真实用户位置:   [%.8f, %.8f, %.8f]^T\n',p_ue);
    fprintf('估计用户位置:   [%.8f, %.8f, %.8f]^T\n',p_ue_est);
    fprintf('位置误差: %.8e m\n',position_error);
end
%% 根据用户估计位置重构未反馈 UPA 的信道
index_unfeedback = setdiff(1:B,index_known);

% 只把已反馈 AP 的信道交给重构函数；未反馈 AP 的块在这里明确屏蔽。
H_feedback = zeros(size(H_total),'like',H_total);
for ii = 1:numel(index_known)
    b = index_known(ii);
    idx_b = (b-1)*Nb + (1:Nb);
    H_feedback(:,idx_b) = H_total(:,idx_b);
end

[H_unfeedback_est,H_model_est,nu_unfeedback,index_dft_unfeedback] = ...
    infer_unfeedback_los_channel(p_ue_est,index_unfeedback,index_known,H_feedback, ...
    pb,Rb,R_ue,p_ue_lcs,p_lcs,k_c);

% 将已反馈 AP 的实测信道和未反馈 AP 的重构信道合并，便于后续预编码/检测。
H_complete_est = H_feedback;
for ii = 1:numel(index_unfeedback)
    b = index_unfeedback(ii);
    idx_b = (b-1)*Nb + (1:Nb);
    H_complete_est(:,idx_b) = H_unfeedback_est(:,:,ii);
end

%% 重构信道 NMSE
% H_complete_est 是已反馈信道和重构信道拼接后的完整信道
H_error = H_complete_est - H_total_clean;
nmse_complete = sum(abs(H_error(:)).^2) / max(sum(abs(H_total_clean(:)).^2), realmin);% 线性 NMSE
nmse_complete_dB = 10 * log10(max(nmse_complete, realmin));% dB 形式

if mc == 1
    fprintf('\n===== 重构信道 NMSE =====\n');
    fprintf('完整信道 NMSE = %.6e\n', nmse_complete);
    fprintf('完整信道 NMSE = %.4f dB\n', nmse_complete_dB);
end

%% 最优预编码矩阵
[U,S,V] = svd(H_total_clean,'econ');
sigma = diag(S);

%% 根据估计信道设计预编码器
[U_est, S_est, V_est] = svd(H_complete_est, 'econ');
sigma_est = diag(S_est);

% DFT matrix selection
index_est_dftb = zeros(B,Ns);
for b = 1 : B
    idx_b = (b-1)*Nb + (1:Nb);
    for i = 1:Ns
        v_est_local = V_est(idx_b,1);
        [~,index_est_dftb(b,i)] = max(abs(v_est_local' * dft_hv) / norm(v_est_local));
    end
end
% 基于估计信道的连续 LS 系数：同时保留 AP 间幅度和相位
coeff_est_ls = zeros(B,1,'like',V_est);
w_est = zeros(Nb,1);
W_est_global = zeros(N,Ns);

for it2 = 1 : B
    idx_b = (it2-1)*Nb + (1:Nb);

    v_est_local = V_est(idx_b,1);
    f_est_local = dft_hv(:,index_est_dftb(it2));

    coeff_est_ls(it2) = (f_est_local' * v_est_local) / (f_est_local' * f_est_local);
    coeff_est_phase = coeff_est_ls(it2) /  max(abs(coeff_est_ls(it2)));
    w_est = coeff_est_phase * f_est_local;

    W_est_global(idx_b) = w_est; % 拼接到第i条流的全局波束向量中

end


W_est_global = W_est_global / norm(W_est_global) ;
gain_est = sqrt(norm(H_total_clean * W_est_global)^2  / sigma(1)^2);
 
   
gain_ideal = norm(H_total_clean*V(:,1))^2/sigma(1)^2;
%% 原来相位6bit量化
N_quan = 6;
coeff_list = exp(1j*pi * (0:(2^N_quan-1))/2^(N_quan-1));

w_est_Q = zeros(Nb,1);
W_est_Q_global = zeros(N,Ns);
index_est_Q_coeffb = zeros(1,B);

for it2 = 1 : B
    idx_b = (it2-1)*Nb + (1:Nb);

    v_est_local = V_est(idx_b,1);
    [~,index_est_Q_coeffb(it2)] = max(real( (v_est_local' * dft_hv(:,index_est_dftb(it2))) .* coeff_list)) ;
    w_est = coeff_list(index_est_Q_coeffb(it2)) * dft_hv(:,index_est_dftb(it2));
   
    W_est_Q_global(idx_b) = w_est; % 拼接到第i条流的全局波束向量中

end

W_est_Q_global = W_est_Q_global / norm(W_est_Q_global) ;

gain_est_Q = sqrt(norm(H_total_clean * W_est_Q_global)^2  / sigma(1)^2);

% Store this trial's gains. Both schemes use the same estimated channel and DFT beams.
gain_cont_mc(mc) = gain_est;
gain_quant_mc(mc) = gain_est_Q;
gain_diff_mc(mc) = gain_est - gain_est_Q;
gain_diff_db_mc(mc) = 10*log10(max(gain_est,realmin)) -10*log10(max(gain_est_Q,realmin));

end % Monte-Carlo loop





