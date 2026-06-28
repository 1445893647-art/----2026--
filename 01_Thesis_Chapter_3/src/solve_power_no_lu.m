function P_opt = solve_power_no_lu(Q, par)
    % SOLVE_POWER_NO_LU 固定轨迹优化干扰功率 (不考虑合法用户)
    
    P_opt = zeros(1, par.N);
    
    % 计算 Lambert W 相关常数 K
    term = (par.p_prior - par.epsilon_det) / par.p_prior;
    arg_W = (par.p_prior * log(term) / (1 - par.p_prior)) * term^(par.p_prior / (1 - par.p_prior));
    K_val = 1 + lambertw(arg_W) / log(term) - (1 - par.p_prior) / par.p_prior;
    
    for n = 1:par.N
        % 使用估计位置
        d_AU = norm(Q(:,n) - par.pos_Alice_Est); 
        
        % --- 约束 A: 检测能力约束 ---
        % 仅受限于检测错误概率
        P_det_limit = (par.beta0 * par.P_Alice) / ...
                      (par.rho * par.lambda_uu * (d_AU^par.alpha0) * K_val);
                      
        % --- 综合求解 (移除 QoS 约束) ---
        % P_opt = min(P_max, P_detection)
        P_opt(n) = min([par.P_UAV_max, P_det_limit]);
        P_opt(n) = max(P_opt(n), 0);
    end
end