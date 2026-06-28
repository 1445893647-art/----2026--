function [Q_new, status] = solve_traj_no_lu(P, Q_old, par)
    % SOLVE_TRAJ_NO_LU 固定功率优化轨迹 (不考虑合法用户)
    
    % 计算 Alice 距离限制圆半径
    R2_max_Alice = zeros(1, par.N);
    
    term = (par.p_prior - par.epsilon_det) / par.p_prior;
    arg_W = (par.p_prior * log(term) / (1 - par.p_prior)) * term^(par.p_prior / (1 - par.p_prior));
    K_val = 1 + lambertw(arg_W) / log(term) - (1 - par.p_prior) / par.p_prior;
    
    for n = 1:par.N
        if P(n) > 1e-9
            val = (par.beta0 * par.P_Alice) / (par.rho * P(n) * par.lambda_uu * K_val);
            R2_max_Alice(n) = val^(2/par.alpha0);
        else
            R2_max_Alice(n) = 1e10;
        end
    end
    
    % === CVX 求解 ===
    cvx_begin quiet
        variable q(3, par.N)
        
        % 目标函数: 最小化到 Bob (估计位置) 的距离
        expression dist_Bob_sq(par.N)
        for n = 1:par.N
            dist_Bob_sq(n) = pow_pos(norm(q(:,n) - par.pos_Bob_Est), 2);
        end
        minimize( sum(dist_Bob_sq) )
        
        subject to
            % 1. 初始位置
            q(:,1) == par.pos_UAV_start;
            
            % 2. 机动性约束
            q(3, :) >= par.H_min;
            for n = 2:par.N
                pow_pos(norm(q(:,n) - q(:,n-1)), 2) <= (par.V_max * par.delta_t)^2;
            end
            
            % 3. 检测能力约束 (Alice 圆内)
            for n = 1:par.N
                pow_pos(norm(q(:,n) - par.pos_Alice_Est), 2) <= R2_max_Alice(n);
            end
            
            % --- 移除: LR 保护约束 ---
            % 这里不再有关于 pos_LR 的约束
            
    cvx_end
    
    Q_new = q;
    status = cvx_status;
end