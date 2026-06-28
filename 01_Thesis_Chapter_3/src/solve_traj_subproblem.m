function [Q_new, status] = solve_traj_subproblem(P, Q_old, par)
    % SOLVE_TRAJ_SUBPROBLEM 固定功率优化轨迹 (SCA + CVX)
    % 对应论文 Problem P5
    
    % 预先计算约束所需的辅助变量 R_max_Alice (Eq. 23) 和 R_min_LR (Eq. 25)
    R2_max_Alice = zeros(1, par.N);
    R2_min_LR    = zeros(1, par.N);
    
    % 计算 Lambert K (同功率子问题)
    term = (par.p_prior - par.epsilon_det) / par.p_prior;
    arg_W = (par.p_prior * log(term) / (1 - par.p_prior)) * term^(par.p_prior / (1 - par.p_prior));
    K_val = 1 + lambertw(arg_W) / log(term) - (1 - par.p_prior) / par.p_prior;
    
    gamma_th = 2^par.R_min_LR - 1;
    h_LT_LR_sq    = par.beta0 / norm(par.pos_LT - par.pos_LR)^par.alpha0;
    h_Alice_LR_sq = par.beta0 / norm(par.pos_Alice_Est - par.pos_LR)^par.alpha0;
    
    for n = 1:par.N
        % 1. Alice 距离限制圆半径 (Eq. 23)
        if P(n) > 1e-9
            val = (par.beta0 * par.P_Alice) / (par.rho * P(n) * par.lambda_uu * K_val);
            R2_max_Alice(n) = val^(2/par.alpha0);
        else
            R2_max_Alice(n) = 1e10; % 无穷大
        end
        
        % 2. LR 保护区域半径 (Eq. 25)
        denom = (par.P_LT * h_LT_LR_sq / gamma_th) - par.P_Alice * h_Alice_LR_sq - par.sigma2_noise;
        if denom > 0
            % Eq 25: d^2 >= (P_UAV * beta / denom)^(2/alpha)
            val_lr = (par.rho * P(n) * par.beta0) / denom;
            R2_min_LR(n) = val_lr^(2/par.alpha0);
        else
            R2_min_LR(n) = 0;
        end
    end
    
    % === CVX 求解 ===
    cvx_begin quiet
        variable q(3, par.N)
        
        % 目标函数: 最小化到 Bob 的距离 (等效于最大化中断概率)
        expression dist_Bob_sq(par.N)
        for n = 1:par.N
            dist_Bob_sq(n) = pow_pos(norm(q(:,n) - par.pos_Bob_Est), 2);
        end
        minimize( sum(dist_Bob_sq) )
        
        subject to
            % 1. 初始位置
            q(:,1) == par.pos_UAV_start;
            
            % 2. 机动性约束 (Eq. 1, 2)
            q(3, :) >= par.H_min; % 高度
            for n = 2:par.N
                pow_pos(norm(q(:,n) - q(:,n-1)), 2) <= (par.V_max * par.delta_t)^2;
            end
            
            % 3. 检测能力约束 (Alice 圆内) (Eq. 22)
            for n = 1:par.N
                pow_pos(norm(q(:,n) - par.pos_Alice_Est), 2) <= R2_max_Alice(n);
            end
            
            % 4. LR 保护约束 (SCA 线性化 - Eq. 26)
            % 原始: ||q - q_LR||^2 >= R^2 (非凸)
            % SCA: ||q_old - q_LR||^2 + 2(q_old - q_LR)'(q - q_old) >= R^2
            for n = 1:par.N
                dist_old_sq = norm(Q_old(:,n) - par.pos_LR)^2;
                term_linear = 2 * (Q_old(:,n) - par.pos_LR)' * (q(:,n) - Q_old(:,n));
                
                dist_old_sq + term_linear >= R2_min_LR(n);
            end
    cvx_end
    
    Q_new = q;
    status = cvx_status;
end