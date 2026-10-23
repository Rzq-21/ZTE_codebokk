%% Monte-Carlo: normalized beam gain versus channel-estimation SNR
% Six curves:
%   SVD       - gain_ideal, unconstrained full-CSI right singular vector
%   Full CSI  - gain_FullCSI, full-CSI DFT/phase-quantized precoder
%   Type-I    - gain_typei, noisy-CSI DFT/phase-quantized precoder
%   Type-II   - gain_typeII, noisy-CSI two-beam Type-II precoder
%   HAQSE     - gain_est, HAQSE position based reconstruction
%   DFT-index - gain_est_raw, DFT-index position based reconstruction
%
% The single-trial implementation is kept in run_one_trial below so that
% every SNR point uses the same geometry and independent channel noise.

clear; clc; close all;
rng(20);

SNR_list = 20:5:40;
N_MC = 100;                         % increase for smoother curves

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
cfg.Qopt = 3;
cfg.Fs_spatial_y = 1;
cfg.Fs_spatial_z = 1;
cfg.O1 = 4;
cfg.O2 = 4;
cfg.N_quan = 6;
cfg.coeff_list = exp(1j*pi * (0:2^cfg.N_quan-1) / ...
    2^(cfg.N_quan-1));
cfg.TypeII_L = 2;
cfg.N_amp_typeII = 3;
cfg.N_phase_typeII = 2;
cfg.N_quan_typeII = 8;
cfg.amp_cb_typeII = linspace(0,1,2^cfg.N_amp_typeII);
cfg.phase_cb_typeII = exp(1j*2*pi*(0:2^cfg.N_phase_typeII-1) / ...
    2^cfg.N_phase_typeII);
cfg.coeff_list_typeII = exp(1j*pi*(0:2^cfg.N_quan_typeII-1) / ...
    2^(cfg.N_quan_typeII-1));
cfg.p_ue = [25;25;0];
cfg.alpha_ue = 0;
cfg.beta_ue = 0;
cfg.gamma_ue = 0;

[cfg.pb,cfg.Rb,cfg.p_lcs,cfg.p_gcs,cfg.R_ue,cfg.p_ue_lcs, ...
    cfg.dft_hv,H_total_clean] = build_geometry(cfg);

gain_ideal_mc = zeros(numel(SNR_list),N_MC);
gain_full_csi_mc = zeros(numel(SNR_list),N_MC);
gain_typei_mc = zeros(numel(SNR_list),N_MC);
gain_typeii_mc = zeros(numel(SNR_list),N_MC);
gain_haqse_mc = zeros(numel(SNR_list),N_MC);
gain_dft_index_mc = zeros(numel(SNR_list),N_MC);

for isnr = 1:numel(SNR_list)
    fprintf('SNR = %d dB (%d trials)\n',SNR_list(isnr),N_MC);
    for imc = 1:N_MC
        [gain_ideal_mc(isnr,imc),gain_full_csi_mc(isnr,imc), ...
            gain_typei_mc(isnr,imc),gain_typeii_mc(isnr,imc), ...
            gain_haqse_mc(isnr,imc),gain_dft_index_mc(isnr,imc)] = ...
            run_one_trial(SNR_list(isnr),cfg,H_total_clean);
    end
end

gain_ideal = mean(gain_ideal_mc,2,'omitnan');
gain_FullCSI = mean(gain_full_csi_mc,2,'omitnan');
gain_typei = mean(gain_typei_mc,2,'omitnan');
gain_typeII = mean(gain_typeii_mc,2,'omitnan');
gain_est = mean(gain_haqse_mc,2,'omitnan');
gain_est_raw = mean(gain_dft_index_mc,2,'omitnan');


gain_ideal = sqrt(gain_ideal) ;
gain_FullCSI = sqrt(gain_FullCSI) ;
gain_typei = sqrt(gain_typei);
gain_typeII = sqrt(gain_typeII);
gain_est = sqrt(gain_est);
gain_est_raw = sqrt(gain_est_raw) ;


figure('Color','w');
plot(SNR_list,gain_ideal,'k-*','LineWidth',1.5,'MarkerSize',7); hold on;
% plot(SNR_list,gain_FullCSI,'b-o','LineWidth',1.5,'MarkerSize',7);
plot(SNR_list,gain_typei,'c-d','LineWidth',1.5,'MarkerSize',7);
plot(SNR_list,gain_typeII,'b-o','LineWidth',1.5,'MarkerSize',7);
plot(SNR_list,gain_est,'r-s','LineWidth',1.5,'MarkerSize',7);
plot(SNR_list,gain_est_raw,'m-^','LineWidth',1.5,'MarkerSize',7);
grid on;
xlabel('SNR (dB)');
ylabel('Normalized beam gain');
% legend({'SVD','Full CSI','Type-I','Type-II','HAQSE','DFT-index position'},'Location','best');
legend({'SVD','Type-I','Type-II','HAQSE','DFT-index position'}, 'Location','best');
title(sprintf('Monte-Carlo normalized beam gain',N_MC));


%% ========================================================================
function [gain_ideal,gain_FullCSI,gain_typei,gain_typeII,gain_est,gain_est_raw] = ...
    run_one_trial(SNR_H_dB,cfg,H_total_clean)

B = cfg.B;
Nb = cfg.Nb;
N = cfg.N;
Nu = cfg.Nu;
Ns = cfg.Ns;
index_known = cfg.index_known;
index_unfeedback = setdiff(1:B,index_known);
dft_hv = cfg.dft_hv;
coeff_list = cfg.coeff_list;

signal_power = mean(abs(H_total_clean(:)).^2);
noise_power = signal_power / 10^(SNR_H_dB/10);
H_noise = sqrt(noise_power/2) * ...
    (randn(size(H_total_clean)) + 1j*randn(size(H_total_clean)));
H_total = H_total_clean + H_noise;

% ---------------------------------------------------------------
% HAQSE position estimate from the three known APs
% ---------------------------------------------------------------
nu_v_hat = zeros(1,B);
nu_h_hat = zeros(1,B);
q_y_hat = zeros(1,B);
q_z_hat = zeros(1,B);
q_x_abs_hat = zeros(1,B);

for b = index_known
    idx_b = (b-1)*Nb + (1:Nb);
    H_b = H_total(:,idx_b);
    S_tensor = reshape(H_b.',[cfg.Nv,cfg.Nh,Nu]);
    estimated_freqs = MMV_CF_mQSE2_SumOfRatios( ...
        S_tensor,cfg.L_path,cfg.Fs_spatial_y,cfg.Fs_spatial_z,cfg.Qopt);

    nu_v_hat(b) = estimated_freqs(1,1);
    nu_h_hat(b) = estimated_freqs(1,2);


    b_nu = 14;
    Q = 2^b_nu;

    nu_v_hat(b) = -0.5 + round((nu_v_hat(b) + 0.5) * (Q-1)) / (Q-1);

    nu_h_hat(b) = -0.5 + round((nu_h_hat(b) + 0.5) * (Q-1)) / (Q-1);



    q_y_hat(b) = max(-1,min(1,-nu_v_hat(b)*cfg.lambda_c/cfg.d));
    q_z_hat(b) = max(-1,min(1,-nu_h_hat(b)*cfg.lambda_c/cfg.d));
    q_x_abs_hat(b) = sqrt(max(0,1-q_y_hat(b)^2-q_z_hat(b)^2));
end

p_ue_est = solve_position_front_back(cfg.pb,cfg.Rb,index_known, ...
    q_y_hat,q_z_hat,q_x_abs_hat);

% ---------------------------------------------------------------
% DFT-index position estimate. Use all Nu receive antennas through
% the local SVD; v.' (ordinary transpose) restores the h_b' convention.
% ---------------------------------------------------------------
index_dft_feedback = nan(1,B);
for b = index_known
    idx_b = (b-1)*Nb + (1:Nb);
    H_b = H_total(:,idx_b);
    [~,~,V_b] = svd(H_b,'econ');
    v_b = V_b(:,1);
    [~,index_dft_feedback(b)] = max(abs(v_b.'*dft_hv)/max(norm(v_b),realmin));
end

p_ue_est_raw = solve_position_from_dft_index( ...
    index_dft_feedback,index_known,cfg.pb,cfg.Rb,cfg.p_gcs, ...
    cfg.Nh,cfg.Nv,cfg.O2,cfg.O1,cfg.lambda_c,cfg.d,dft_hv);

% ---------------------------------------------------------------
% Build the two reconstructed full channels.
% ---------------------------------------------------------------
H_feedback = zeros(size(H_total),'like',H_total);
for ii = 1:numel(index_known)
    b = index_known(ii);
    idx_b = (b-1)*Nb + (1:Nb);
    H_feedback(:,idx_b) = H_total(:,idx_b);
end

[H_unfeedback_est,~,~,~] = infer_unfeedback_los_channel( ...
    p_ue_est,index_unfeedback,index_known,H_feedback,cfg.pb,cfg.Rb, ...
    cfg.R_ue,cfg.p_ue_lcs,cfg.p_lcs,cfg.k_c);
H_complete_est = H_feedback;
for ii = 1:numel(index_unfeedback)
    b = index_unfeedback(ii);
    idx_b = (b-1)*Nb + (1:Nb);
    H_complete_est(:,idx_b) = H_unfeedback_est(:,:,ii);
end

[H_unfeedback_est_raw,~,~,~] = infer_unfeedback_los_channel( ...
    p_ue_est_raw,index_unfeedback,index_known,H_feedback,cfg.pb,cfg.Rb, ...
    cfg.R_ue,cfg.p_ue_lcs,cfg.p_lcs,cfg.k_c);
H_complete_est_raw = H_feedback;
for ii = 1:numel(index_unfeedback)
    b = index_unfeedback(ii);
    idx_b = (b-1)*Nb + (1:Nb);
    H_complete_est_raw(:,idx_b) = H_unfeedback_est_raw(:,:,ii);
end

% ---------------------------------------------------------------
% SVD reference and full-CSI DFT/phase-quantized precoder
% ---------------------------------------------------------------
[~,S,V] = svd(H_total_clean,'econ');
sigma = diag(S);
gain_ideal = norm(H_total_clean*V(:,1))^2 / sigma(1)^2;

W_global = build_dft_precoder_from_global_v(V,dft_hv,coeff_list,B,Nb,Ns);
gain_FullCSI = norm(H_total_clean*W_global)^2 / sigma(1)^2;

% ---------------------------------------------------------------
% Type-I: noisy-CSI SVD followed by one DFT beam and a quantized
% phase per AP.  The selected precoder is evaluated on the clean
% channel using the same normalization as the other curves.
% ---------------------------------------------------------------
[~,~,V_noisy] = svd(H_total,'econ');
W_typei_global = build_dft_precoder_from_global_v( ...
    V_noisy,dft_hv,coeff_list,B,Nb,Ns);
gain_typei = norm(H_total_clean*W_typei_global)^2 / sigma(1)^2;

% ---------------------------------------------------------------
% Type-II: for every AP, select two DFT beams from the same
% orthogonal group using noisy CSI. Quantize their relative amplitude
% with 3 bits and relative phase with 2 bits, then apply an 8-bit
% common phase across APs. Evaluate the result on the clean channel.
% ---------------------------------------------------------------
W_typeII = build_typeii_precoder_from_global_v( ...
    V_noisy,dft_hv,cfg.coeff_list_typeII,B,cfg.Nh,cfg.Nv, ...
    cfg.O1,cfg.O2,cfg.TypeII_L,cfg.amp_cb_typeII, ...
    cfg.phase_cb_typeII);
gain_typeII = norm(H_total_clean*W_typeII)^2 / sigma(1)^2;

% ---------------------------------------------------------------
% HAQSE reconstructed-channel precoder
% ---------------------------------------------------------------
[~,~,V_est] = svd(H_complete_est,'econ');
W_est_global = build_dft_precoder_from_global_v( ...
    V_est,dft_hv,coeff_list,B,Nb,Ns);
gain_est = norm(H_total_clean*W_est_global)^2 / sigma(1)^2;

% ---------------------------------------------------------------
% DFT-index-position reconstructed-channel precoder
% ---------------------------------------------------------------
[~,~,V_est_raw] = svd(H_complete_est_raw,'econ');
W_est_raw_global = build_dft_precoder_from_global_v( ...
    V_est_raw,dft_hv,coeff_list,B,Nb,Ns);
gain_est_raw = norm(H_total_clean*W_est_raw_global)^2 / sigma(1)^2;
end

%% ========================================================================
function W_typeII = build_typeii_precoder_from_global_v( ...
    V_noisy,dft_hv,coeff_list_typeII,B,Nh,Nv,O1,O2,TypeII_L, ...
    amp_cb_typeII,phase_cb_typeII)

Nb = Nh*Nv;
N = B*Nb;
W_typeII_local = zeros(Nb,B);
typeII_local_gain = zeros(B,1);

for b = 1:B
    idx_b = (b-1)*Nb + (1:Nb);
    v_noisy_local = V_noisy(idx_b,1);
    local_gain_all = v_noisy_local' * dft_hv;

    % Search all O1*O2 orthogonal groups. In each group retain the
    % strongest L beams, then choose the group with maximum energy.
    best_group_metric = -inf;
    selected_beams = zeros(1,TypeII_L);
    selected_gain = zeros(1,TypeII_L);
    for q_h = 0:O2-1
        for q_v = 0:O1-1
            group_beams = zeros(1,Nb);
            kk = 0;
            for n_h = 0:Nh-1
                for n_v = 0:Nv-1
                    kk = kk+1;
                    m_h = O2*n_h+q_h;
                    m_v = O1*n_v+q_v;
                    % dft_hv=kron(dft_h,dft_v): v-index changes faster.
                    group_beams(kk) = m_h*(O1*Nv)+m_v+1;
                end
            end

            group_gain = local_gain_all(group_beams);
            [~,group_order] = sort(abs(group_gain),'descend');
            candidate_beams = group_beams(group_order(1:TypeII_L));
            candidate_gain = local_gain_all(candidate_beams);
            group_metric = sum(abs(candidate_gain).^2);
            if group_metric > best_group_metric
                best_group_metric = group_metric;
                selected_beams = candidate_beams;
                selected_gain = candidate_gain;
            end
        end
    end

    % The strongest selected beam is the unit-amplitude, zero-phase
    % reference for the remaining Type-II beam coefficients.
    [~,ref_idx] = max(abs(selected_gain));
    ref_gain = selected_gain(ref_idx);
    relative_amp = abs(selected_gain)/max(abs(ref_gain),eps);
    coeff_typeII = zeros(TypeII_L,1);
    for il = 1:TypeII_L
        if il == ref_idx
            coeff_typeII(il) = 1;
        else
            [~,amp_index] = min(abs(amp_cb_typeII-relative_amp(il)));
            desired_relative_phase = exp(1j* ...
                (angle(ref_gain)-angle(selected_gain(il))));
            [~,phase_index] = min(abs( ...
                phase_cb_typeII-desired_relative_phase));
            coeff_typeII(il) = amp_cb_typeII(amp_index) * ...
                phase_cb_typeII(phase_index);
        end
    end

    local_precoder = dft_hv(:,selected_beams)*coeff_typeII;
    local_norm = norm(local_precoder);
    if local_norm <= realmin
        error('The generated local Type-II precoder is zero.');
    end
    local_precoder = local_precoder/local_norm;
    W_typeII_local(:,b) = local_precoder;
    typeII_local_gain(b) = v_noisy_local'*local_precoder;
end

% Quantize the common phase of every AP after forming its local beam.
W_typeII = zeros(N,1);
for b = 1:B
    [~,ap_phase_index] = max(real( ...
        typeII_local_gain(b).*coeff_list_typeII));
    idx_b = (b-1)*Nb + (1:Nb);
    W_typeII(idx_b) = coeff_list_typeII(ap_phase_index) * ...
        W_typeII_local(:,b);
end

W_norm = norm(W_typeII,'fro');
if W_norm <= realmin
    error('The generated Type-II precoder is zero.');
end
W_typeII = W_typeII/W_norm;
end

%% ========================================================================
function W_global = build_dft_precoder_from_global_v(V,dft_hv,coeff_list,B,Nb,Ns)
N = B*Nb;
W_global = zeros(N,Ns);
for b = 1:B
    idx_b = (b-1)*Nb + (1:Nb);
    for i = 1:Ns
        v_local = V(idx_b,i);
        [~,index_dft] = max(abs(v_local'*dft_hv)/max(norm(v_local),realmin));

        % 使用量化相位
        % local_gain = v_local' * dft_hv(:,index_dft);
        % [~,index_coeff] = max(real(local_gain .* coeff_list));
        % W_global(idx_b,i) = coeff_list(index_coeff) * dft_hv(:,index_dft);

        % 根据最小二乘法，使用连续相位
        f_local = dft_hv(:,index_dft);
        coeff_ls(b) = (f_local' * v_local) / (f_local' * f_local);
        coeff_phase = coeff_ls(b) /  max(abs(coeff_ls(b)));
        W_global(idx_b,i) = coeff_phase * f_local;
    end
end
W_norm = norm(W_global,'fro');
if W_norm <= realmin
    error('The generated DFT precoder is zero.');
end
W_global = W_global / W_norm;
end

%% ========================================================================
function [pb,Rb,p_lcs,p_gcs,R_ue,p_ue_lcs,dft_hv,H_total_clean] = ...
    build_geometry(cfg)

B = cfg.B;
Nh = cfg.Nh;
Nv = cfg.Nv;
Nb = cfg.Nb;
N = cfg.N;
Nu = cfg.Nu;
d = cfg.d;
k_c = cfg.k_c;
p_ue = cfg.p_ue;

pb = zeros(3,B);
pb(:,1) = [0;0;0];
pb(:,2) = [30;-10;0];
pb(:,3) = [0;50;0];
pb(:,4) = [50;50;0];
pb(:,5) = [-1;-1;0];
pb(:,6) = [0;-1;0];
pb(:,7) = [0;51;0];
pb(:,8) = [20;-10;0];

alpha_b = [0;130;-45;-125;49;43;98;90] / 180*pi;
beta_b = [0;10;-5;2;6;-4;8;-7] / 180*pi;
gamma_b = zeros(B,1);
Rb = zeros(3,3,B);
for b = 1:B
    ca = cos(alpha_b(b)); sa = sin(alpha_b(b));
    cb = cos(beta_b(b)); sb = sin(beta_b(b));
    cg = cos(gamma_b(b)); sg = sin(gamma_b(b));
    Rb(:,:,b) = [ca*cb, ca*sb*sg-sa*cg, ca*sb*cg+sa*sg; ...
        sa*cb, sa*sb*sg+ca*cg, sa*sb*cg-ca*sg; ...
        -sb, cb*sg, cb*cg];
end

p_lcs = zeros(3,Nb);
p_gcs = zeros(3,Nb,B);
for b = 1:B
    idx_local = 0;
    for nh = 1:Nh
        for nv = 1:Nv
            y = ((nv-1)-(Nv-1)/2)*d;
            z = ((nh-1)-(Nh-1)/2)*d;
            idx_local = idx_local + 1;
            p_lcs(:,idx_local) = [0;y;z];
            p_gcs(:,idx_local,b) = Rb(:,:,b)*p_lcs(:,idx_local)+pb(:,b);
        end
    end
end

R_ue = eye(3);
ue_offset = ((0:Nu-1)-(Nu-1)/2)*d;
p_ue_lcs = zeros(3,Nu);
p_ue_lcs(1,:) = ue_offset;

H_L = zeros(Nu,Nb,B);
H_total_clean = zeros(Nu,N);
for b = 1:B
    d_bu = pb(:,b)-p_ue;
    r_bu = norm(d_bu);
    d_hat = d_bu/r_bu;
    d_hat_bLCS = Rb(:,:,b).'*d_hat;
    d_hat_uLCS = R_ue.'*d_hat;
    phase_center = exp(-1j*k_c*r_bu)/r_bu;
    for m = 1:Nu
        phase_UE = exp(1j*k_c*(d_hat_uLCS.'*p_ue_lcs(:,m)));
        for n = 1:Nb
            phase_AP = exp(-1j*k_c*(d_hat_bLCS.'*p_lcs(:,n)));
            H_L(m,n,b) = phase_center*phase_AP*phase_UE;
        end
    end
    idx_b = (b-1)*Nb + (1:Nb);
    H_total_clean(:,idx_b) = H_L(:,:,b);
end

dft_h = exp(1j*2*pi*(0:Nh-1).'*(0:cfg.O2*Nh-1)/(cfg.O2*Nh))/sqrt(Nh);
dft_v = exp(1j*2*pi*(0:Nv-1).'*(0:cfg.O1*Nv-1)/(cfg.O1*Nv))/sqrt(Nv);
dft_hv = kron(dft_h,dft_v);
end
