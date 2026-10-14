% function [p_best,index_dft_infer,info] = solve_position_from_dft_index( ...
%     index_dft,index_known,pb,Rb,p_gcs,Nh,Nv,Oy,Oz,lambda_c,d,dft_hv)
function [p_best] = solve_position_from_dft_index( ...
    index_dft,index_known,pb,Rb,p_gcs,Nh,Nv,Oy,Oz,lambda_c,d,dft_hv)


% 取出已知基站的DFT索引
index_dft_knownb = index_dft(index_known);
B_known = length(index_known);
% 把一个总的DFT拆成y方向和z方向两个索引
index_dfty_knownb = mod(index_dft_knownb-1,Nh * Oy);
index_dftz_knownb = floor((index_dft_knownb-1)/(Nv * Oz));

%% fftshift processing the negative frequency ，负频移位，恢复真实角度
%% recover theta ，恢复仰角
indexz_negative = find(index_dftz_knownb/(Nv * Oz) > 0.5);
index_dftz_knownb(indexz_negative) = index_dftz_knownb(indexz_negative) - (Nv * Oz);
theta_est = acos(index_dftz_knownb/(Nv * Oz) * (lambda_c/d));

%% recover phi，恢复方位角
indexy_negative = find(index_dfty_knownb/(Nh * Oy) > 0.5);
index_dfty_knownb(indexy_negative) = index_dfty_knownb(indexy_negative) - (Nh * Oy);
phi_est = asin(index_dfty_knownb/(Nh * Oy) * (lambda_c/d) ./ sin(theta_est));

theta_est = reshape(theta_est,B_known,1);
phi_est = reshape(phi_est,B_known,1);

d_ue_b_est = [sin(theta_est).*cos(phi_est), ...
    sin(theta_est).*sin(phi_est), ...
    cos(theta_est)].';

% Convert each known AP direction from its local coordinate system to
% the global coordinate system and solve the multi-line least-squares
% intersection, exactly as in the original code.
R_inv = zeros(3,3);
res_inv = zeros(3,1);
for ii = 1:B_known
    b = index_known(ii);
    dg_ue_b_est = Rb(:,:,b) * d_ue_b_est(:,ii);
    proj_b = eye(3) - dg_ue_b_est * dg_ue_b_est.';
    R_inv = R_inv + proj_b;
    res_inv = res_inv + proj_b * pb(:,b);
end

p_best = R_inv \ res_inv;


end
