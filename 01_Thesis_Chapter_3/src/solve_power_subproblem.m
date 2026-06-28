function P_opt = solve_power_subproblem(Q, par)
    % SOLVE_POWER_SUBPROBLEM 固定轨迹优化干扰功率
    % 基于论文 Eq. 19, 20, 21
    
    P_opt = zeros(1, par.N);
    
    % 1. 计算 Lambert W 函数相关的常数 K (Eq. 18)
    term = (par.p_prior - par.epsilon_det) / par.p_prior;
    if term <= 0
        error('检测约束参数设置错误，无法满足');
    end
    arg_W = (par.p_prior * log(term) / (1 - par.p_prior)) * term^(par.p_prior / (1 - par.p_prior));
    K_val = 1 + lambertw(arg_W) / log(term) - (1 - par.p_prior) / par.p_prior;
    
    for n = 1:par.N
        d_AU = norm(Q(:,n) - par.pos_Alice_Est); % d_Alice_UAV
        d_UL = norm(Q(:,n) - par.pos_LR);    % d_UAV_LR
        
        % --- 约束 A: 检测能力约束 (Eq. 17 / Eq. 21) ---
        % P_UAV <= (beta0 * P_Alice) / (rho * lambda * d_AU^alpha * K)
        % 注意: 论文 Eq 21 实际上就是检测约束取等号时的值，同时也是最大化中断概率的最优值(如果不管其他约束)
        P_det_limit = (par.beta0 * par.P_Alice) / ...
                      (par.rho * par.lambda_uu * (d_AU^par.alpha0) * K_val);
                      
        % --- 约束 B: 合法用户 QoS 约束 (Eq. 20) ---
        % R_LR >= R_min => SINR >= gamma_th
        gamma_th = 2^par.R_min_LR - 1;
        
        % 计算 LR 处的有用信号和来自 Alice 的干扰
        h_LT_LR_sq    = par.beta0 / norm(par.pos_LT - par.pos_LR)^par.alpha0;
        h_Alice_LR_sq = par.beta0 / norm(par.pos_Alice_Est - par.pos_LR)^par.alpha0;
        h_UAV_LR_sq   = par.beta0 / (d_UL^par.alpha0);
        
        % (P_LT * h_LT_LR) / (P_A * h_A_LR + P_U * h_U_LR + sigma) >= gamma
        % 推导 P_U 的上界:
        numerator = (par.P_LT * h_LT_LR_sq / gamma_th) - par.P_Alice * h_Alice_LR_sq - par.sigma2_noise;
        
        if numerator < 0
            P_qos_limit = 0; % 即使不干扰也无法满足 QoS
        else
            P_qos_limit = numerator / h_UAV_LR_sq;
        end
        
        % --- 综合求解 ---
        % 取所有限制的最小值
        P_opt(n) = min([par.P_UAV_max, P_det_limit, P_qos_limit]);
        P_opt(n) = max(P_opt(n), 0); % 保证非负
    end
end