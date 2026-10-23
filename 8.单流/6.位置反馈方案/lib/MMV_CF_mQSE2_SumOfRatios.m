function estimated_freqs = MMV_CF_mQSE2_SumOfRatios(S_tensor, L, Fs1, Fs2, Qopt)
% MMV_CF_mQSE2_SumOfRatios: 论文CF-mQSE2算法的2D-MMV扩展版
% 采用 "Sum of Ratios" (非相干累积) 策略
%
% 输入:
%   S_tensor        - 输入信号张量 (N1 x N2 x T)
%   L               - 信号分量个数
%   Fs1, Fs2        - 两个维度上的采样率
%   Qopt            - 精细估计的迭代次数 (通常为2或3)
%
% 输出:
%   estimated_freqs - 估计出的 L x 2 频率矩阵

    [N1, N2, T] = size(S_tensor);
    n1_vec = (0:N1-1)';
    n2_vec = (0:N2-1);

    % --- 参数设置 ---
    q1 = min(1/sqrt(3*N1*N2), 0.3232);
    q2 = min(1/sqrt(3*N1*N2), 0.3232);

    % 校正因子 c(q) (对应论文式(6))
    ck = @(N, q) (pi * cot(pi*q/N) - pi*N*cot(pi*q)) / (N * cos(pi*q*(N-1)/N)^2);
    c1 = ck(N1, q1);
    c2 = ck(N2, q2);

    % --- 初始化存储变量 ---
    coarse_p = zeros(L, 2);
    coarse_delta = zeros(L, 2);
    coarse_amps_matrix = zeros(L, T, 'like', 1j);
    
    S_current = S_tensor;

    %% ========================================================================
    %  阶段一: 粗估计 (Coarse Estimation)
    % =========================================================================
    % 当前单流/单路径调用时 L=1；函数仍保留 L>1 的多频点扩展能力。
    for l = 1:L
        % 1. 寻找整数频率 (p1, p2)
        power_spectrum = sum(abs(fft2(S_current)).^2, 3);
        [~, max_idx] = max(power_spectrum(:));
        [p1_idx, p2_idx] = ind2sub(size(power_spectrum), max_idx);% p1_idx：第一维峰值下标，范围 1 到 N1
                                                                  % p2_idx：第二维峰值下标，范围 1 到 N2
        
        p1 = p1_idx - 1;% MATLAB 下标从 1 开始，而 DFT 频率索引通常从 0 开始，所以减 1
        if p1 > N1/2, p1 = p1 - N1; end% 高频部分转换成负频率
        p2 = p2_idx - 1;
        if p2 > N2/2, p2 = p2 - N2; end
        
        coarse_p(l, :) = [p1, p2];% 保存第 l 个二维分量的整数频率

        % 2. 估计小数频率 (delta1, delta2) - Sum of Ratios 方法
        % 维度1
        exp_term1_p = exp(-1j * 2 * pi * n1_vec * (p1 + q1) / N1);% N1×1，第一维 p1+q1 处的 DFT 核
        exp_term1_m = exp(-1j * 2 * pi * n1_vec * (p1 - q1) / N1);% 第一维 p1-q1 处的 DFT 核
        exp_term2   = exp(-1j * 2 * pi * n2_vec * p2 / N2);% 1×N2，固定第二维频率为整数频率 p2
        
        S_plus_q1_vec = squeeze(sum(sum(S_current .* (exp_term1_p * exp_term2), 1), 2));% 对二维阵列的所有采样点求和，因此得到每个快拍在二维频率点处的DFT系数
        S_minus_q1_vec = squeeze(sum(sum(S_current .* (exp_term1_m * exp_term2), 1), 2));

        ratios1 = (S_plus_q1_vec - S_minus_q1_vec) ./ (S_plus_q1_vec + S_minus_q1_vec);% 对每个快拍计算比率，是Q-shift型估计器的比率
        ratios1(isnan(ratios1) | isinf(ratios1)) = 0; % 处理分母为0的情况
        delta1 = (1/c1) * real(mean(ratios1));% 第一维的小数频率修正量
        
        % 维度2
        exp_term1   = exp(-1j * 2 * pi * n1_vec * p1 / N1);
        exp_term2_p = exp(-1j * 2 * pi * n2_vec * (p2 + q2) / N2);
        exp_term2_m = exp(-1j * 2 * pi * n2_vec * (p2 - q2) / N2);
        
        S_plus_q2_vec = squeeze(sum(sum(S_current .* (exp_term1 * exp_term2_p), 1), 2));
        S_minus_q2_vec = squeeze(sum(sum(S_current .* (exp_term1 * exp_term2_m), 1), 2));
        
        ratios2 = (S_plus_q2_vec - S_minus_q2_vec) ./ (S_plus_q2_vec + S_minus_q2_vec);
        ratios2(isnan(ratios2) | isinf(ratios2)) = 0;
        delta2 = (1/c2) * real(mean(ratios2));
        
        coarse_delta(l, :) = [delta1, delta2];% 保存当前二维分量的两个小数修正量

        % 3. 估计复振幅向量
        f1_est = (p1 + delta1) / N1;% 整数频率+小数修正转换成归一化二维频率
        f2_est = (p2 + delta2) / N2;
        exp_total = exp(-1j * 2 * pi * (n1_vec * f1_est + n2_vec * f2_est));
        A_l_vec = squeeze(sum(sum(S_current .* exp_total, 1), 2)) / (N1 * N2);% 估计当前二维分量在每个快拍中的复振幅
        % 对于第t个快拍，将当前信号投影到估计的二维正弦基底上。
        coarse_amps_matrix(l, :) = A_l_vec;

        % 4. 剥离分量
        reconstructed_signal = exp(1j * 2 * pi * (n1_vec * f1_est + n2_vec * f2_est)) .* reshape(A_l_vec, [1, 1, T]);
        S_current = S_current - reconstructed_signal;
    end

    %% ========================================================================
    %  阶段二: 精细估计 (Fine Estimation)
    % =========================================================================
    fine_delta = coarse_delta;
    fine_amps_matrix = coarse_amps_matrix;

    for i = 1:Qopt
        for l = 1:L
            % 1. 分离目标分量
            S_interferers = zeros(N1, N2, T, 'like', 1j);
            for m = 1:L
                if m ~= l
                    p1_m = coarse_p(m, 1);
                    p2_m = coarse_p(m, 2);
                    delta1_m = fine_delta(m, 1); % 使用最新的delta值
                    delta2_m = fine_delta(m, 2);
                    A_m_vec = fine_amps_matrix(m, :);

                    f1_m_est = (p1_m + delta1_m) / N1;
                    f2_m_est = (p2_m + delta2_m) / N2;
                    
                    reconstructed_m = exp(1j * 2 * pi * (n1_vec * f1_m_est + n2_vec * f2_m_est)) .* reshape(A_m_vec, [1, 1, T]);
                    S_interferers = S_interferers + reconstructed_m;
                end
            end
            S_isolated = S_tensor - S_interferers;

            % 2. 更新小数频率 (delta1, delta2) - Sum of Ratios
            p1_l = coarse_p(l, 1);
            p2_l = coarse_p(l, 2);
            delta1_l_prev = fine_delta(l, 1);
            delta2_l_prev = fine_delta(l, 2);

            % 维度1
            exp_term1_p = exp(-1j * 2 * pi * n1_vec * (p1_l + delta1_l_prev + q1) / N1);
            exp_term1_m = exp(-1j * 2 * pi * n1_vec * (p1_l + delta1_l_prev - q1) / N1);
            exp_term2   = exp(-1j * 2 * pi * n2_vec * (p2_l + delta2_l_prev) / N2);
            
            S_plus_q1_vec = squeeze(sum(sum(S_isolated .* (exp_term1_p * exp_term2), 1), 2));
            S_minus_q1_vec = squeeze(sum(sum(S_isolated .* (exp_term1_m * exp_term2), 1), 2));
            
            ratios1 = (S_plus_q1_vec - S_minus_q1_vec) ./ (S_plus_q1_vec + S_minus_q1_vec);
            ratios1(isnan(ratios1) | isinf(ratios1)) = 0;
            delta1_update = (1/c1) * real(mean(ratios1));
            fine_delta(l, 1) = delta1_l_prev + delta1_update;

            % 维度2
            exp_term1   = exp(-1j * 2 * pi * n1_vec * (p1_l + fine_delta(l, 1)) / N1); % 使用更新后的delta1
            exp_term2_p = exp(-1j * 2 * pi * n2_vec * (p2_l + delta2_l_prev + q2) / N2);
            exp_term2_m = exp(-1j * 2 * pi * n2_vec * (p2_l + delta2_l_prev - q2) / N2);

            S_plus_q2_vec = squeeze(sum(sum(S_isolated .* (exp_term1 * exp_term2_p), 1), 2));
            S_minus_q2_vec = squeeze(sum(sum(S_isolated .* (exp_term1 * exp_term2_m), 1), 2));

            ratios2 = (S_plus_q2_vec - S_minus_q2_vec) ./ (S_plus_q2_vec + S_minus_q2_vec);
            ratios2(isnan(ratios2) | isinf(ratios2)) = 0;
            delta2_update = (1/c2) * real(mean(ratios2));
            fine_delta(l, 2) = delta2_l_prev + delta2_update;
            
            % 3. 更新复振幅向量
            f1_est_new = (p1_l + fine_delta(l, 1)) / N1;
            f2_est_new = (p2_l + fine_delta(l, 2)) / N2;
            exp_total_new = exp(-1j * 2 * pi * (n1_vec * f1_est_new + n2_vec * f2_est_new));
            A_l_vec_new = squeeze(sum(sum(S_isolated .* exp_total_new, 1), 2)) / (N1 * N2);
            fine_amps_matrix(l, :) = A_l_vec_new;
        end
    end

    %% ========================================================================
    %  输出最终结果
    % =========================================================================
    estimated_freqs = zeros(L, 2);
    for l = 1:L
        estimated_freqs(l, 1) = (coarse_p(l, 1) + fine_delta(l, 1)) * Fs1 / N1;
        estimated_freqs(l, 2) = (coarse_p(l, 2) + fine_delta(l, 2)) * Fs2 / N2;
    end
end
