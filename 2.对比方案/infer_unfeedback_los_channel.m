% function [H_unfeedback_est,gain_hat,H_model_all,nu_unfeedback,index_dft_unfeedback] = ...
%     infer_unfeedback_los_channel(p_ue_est,index_unfeedback,index_known,H_feedback, ...
%                                  pb,Rb,R_ue,p_ue_lcs,p_lcs,k_c)
function [H_unfeedback_est,H_model_all,nu_unfeedback,index_dft_unfeedback] = ...
    infer_unfeedback_los_channel(p_ue_est,index_unfeedback,index_known,H_feedback, ...
                                 pb,Rb,R_ue,p_ue_lcs,p_lcs,k_c)
% 使用估计用户位置重构未反馈 UPA 的纯 LOS 信道。
%
% 输入：
%   p_ue_est         - 估计的用户中心位置，3×1
%   index_unfeedback - 未反馈 UPA 的索引
%   index_known      - 已反馈 UPA 的索引，用于估计公共复增益
%   H_feedback      - 整体已反馈/带噪信道，Nu×(B*Nb)；仅使用 index_known 块
%   pb,Rb            - UPA 中心位置和局部坐标到全局坐标的旋转矩阵
%   R_ue,p_ue_lcs    - UE 姿态和 UE 阵元相对中心的位置
%   p_lcs,k_c        - UPA 阵元局部坐标和波数
%
% 输出：
%   H_unfeedback_est - 重构的未反馈信道，Nu×Nb×numel(index_unfeedback)
%   gain_hat         - 由已反馈 AP 估计的公共复增益
%   H_model_all      - 未乘公共增益的所有 UPA 几何 LOS 信道
%   nu_unfeedback    - 未反馈 AP 的预测二维空间频率，行对应 AP，列为
%                      [nu_v, nu_h]，单位 cycles/element
%   index_dft_unfeedback - 按 O1/O2 过采样 DFT 码本得到的索引

    

    [~,B] = size(pb);
    Nb = size(p_lcs,2);
    Nu = size(p_ue_lcs,2);
    
    % p_ue_est = [25;25;0];

    H_model_all = complex(zeros(Nu,Nb,B));
    for b = 1:B
        d_bu = pb(:,b) - p_ue_est;       % UE -> AP
        r_bu = norm(d_bu);
       
        d_hat = d_bu / r_bu;
        d_hat_bLCS = Rb(:,:,b).' * d_hat;
        d_hat_uLCS = R_ue.' * d_hat;
        phase_center = exp(-1j*k_c*r_bu) / r_bu;

        for m = 1:Nu
            phase_UE = exp(1j*k_c * (d_hat_uLCS.' * p_ue_lcs(:,m)));
            for n = 1:Nb
                phase_AP = exp(-1j*k_c * (d_hat_bLCS.' * p_lcs(:,n)));
                H_model_all(m,n,b) = phase_center * phase_AP * phase_UE;
            end
        end
    end


    H_unfeedback_est = complex(zeros(Nu,Nb,numel(index_unfeedback)));
    nu_unfeedback = zeros(numel(index_unfeedback),2);
    index_dft_unfeedback = zeros(numel(index_unfeedback),1);
    for ii = 1:numel(index_unfeedback)

        b = index_unfeedback(ii);
        H_unfeedback_est(:,:,ii) =  H_model_all(:,:,b);

    end


end