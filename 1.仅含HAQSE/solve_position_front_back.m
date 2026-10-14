%% 局部函数：枚举平面 UPA 的前后方向并做方向线交汇
function [p_best,sign_best,score_best,t_best, q_local_best,q_global_best] = ...
         solve_position_front_back(pb,Rb,index_known,q_y,q_z,q_x_abs)

    K = numel(index_known);% 参与定位的UPA数量
    n_comb = 2^K;% 每一个 UPA 都有两种 q_x 符号选择

    p_best = nan(3,1);% 初始化最优用户位置
    sign_best = ones(1,K);
    score_best = inf;
    t_best = nan(1,K);% 保存最优方向线上的距离参数：
    q_local_best = nan(3,K);% 保存最优组合下每个 UPA 的局部坐标方向
    q_global_best = nan(3,K);% 保存最优组合下每个 UPA 的全局坐标方向

    for mask = 0:(n_comb-1)
        g = zeros(3,K);              % AP -> UE 全局方向
        q_local_candidate = zeros(3,K);
        q_global_candidate = zeros(3,K);
        sign_candidate = ones(1,K);

        for ii = 1:K
            b = index_known(ii);
            sx = 1;
            if bitget(mask,ii)
                sx = -1;
            end
            sign_candidate(ii) = sx;

            q_local_candidate(:,ii) = [ sx*q_x_abs(b); q_y(b); q_z(b)];

            % q_local_candidate 是 UE -> AP；定位直线使用 AP -> UE
            q_global_candidate(:,ii) = Rb(:,:,b) * q_local_candidate(:,ii);
            g(:,ii) = -q_global_candidate(:,ii);
            g(:,ii) = g(:,ii) / norm(g(:,ii));
        end

        A = zeros(3,3);
        rhs = zeros(3,1);
        P_all = zeros(3,3,K);% 保存每个UPA的投影矩阵，后面计算定位残差时还会使用

        for ii = 1:K
            b = index_known(ii);
            P = eye(3)-g(:,ii)*g(:,ii).';% 垂直于方向 gb 的投影矩阵
            P_all(:,:,ii) = P;
            A = A + P;
            rhs = rhs + P*pb(:,b);
        end

        if rcond(A) < 1e-12
            continue;
        end

        p_candidate = A\rhs;
        t_candidate = zeros(1,K);
        residual = 0;

        for ii = 1:K
            b = index_known(ii);
            delta = p_candidate-pb(:,b);
            t_candidate(ii) = g(:,ii).' * delta;
            residual = residual + norm(P_all(:,:,ii)*delta)^2;
        end

        % 对位于方向线反向的候选组合施加惩罚
        negative_t = max(0,-t_candidate);
        score = residual + 1e6*sum(negative_t.^2);

        if score < score_best
            p_best = p_candidate;
            sign_best = sign_candidate;
            score_best = score;
            t_best = t_candidate;
            q_local_best = q_local_candidate;
            q_global_best = q_global_candidate;
        end
    end

end