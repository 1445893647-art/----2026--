function Gamma_req = get_gamma_req(P)
    % GET_GAMMA_REQ 根据目标检测错误概率 xi_req 反解所需的最小信干比 Gamma_req
    % 基于论文公式: 最优检测错误概率 xi* 是 gamma_w 的单调递减函数
    % 采用 Lambert W 函数计算最优阈值
    
    target_xi = P.xi_req;
    
    % 二分法或 fzero 求解 Gamma_req
    % 范围: 1e-4 (极小SIR) 到 1e4 (极大SIR)
    % xi(gamma) 是单调减函数
    
    try
        % 寻找 xi_exact(gamma) - target_xi = 0 的根
        Gamma_req = fzero(@(g) calc_xi_exact(g, P) - target_xi, [1e-4, 1e5]);
    catch
        % 如果解不出 (例如要求太严苛或太宽松)，返回边界值
        val_low = calc_xi_exact(1e-4, P);
        if val_low < target_xi
            Gamma_req = 1e-4; % 即使很小的 SIR 也能满足 (几乎不可能)
        else
            Gamma_req = 1e5; % 需要极高的 SIR
        end
        warning('get_gamma_req: fzero failed, using boundary value %.2e', Gamma_req);
    end
end

function xi = calc_xi_exact(gamma_w, P)
    % 计算给定 SIR (gamma_w) 下的精确最小检测错误概率
    % 假设归一化模型: eta_I = 1, eta_S = gamma_w
    
    Nr = P.N_r;
    pi0 = 0.5; pi1 = 0.5; % 假设等先验
    
    % 1. 计算分布参数 (矩匹配)
    eta_I = 1;
    eta_S = gamma_w;
    
    sum_eta = eta_S + eta_I;
    sq_sum  = eta_S^2 + eta_I^2;
    
    nu_Z  = (Nr * sum_eta^2) / (sq_sum + 1e-30);
    eta_Z = sq_sum / (sum_eta + 1e-30);
    
    % H0: Gamma(Nr, eta_I)
    % H1: Gamma(nu_Z, eta_Z)
    
    % 2. 求解最优阈值 u* (归一化阈值)
    % 方程: A*u - B*ln(u) = C
    % A = 1/eta_I - 1/eta_Z
    % B = Nr - nu_Z
    % C = ln(pi0/pi1) + ln(G(nu_Z)/G(Nr)) + nu_Z*ln(eta_Z) - Nr*ln(eta_I)
    
    A = (1/eta_I) - (1/eta_Z);
    B = Nr - nu_Z;
    
    % 使用 gammaln 防止溢出
    log_term = log(pi0/pi1) + gammaln(nu_Z) - gammaln(Nr) ...
               + nu_Z*log(eta_Z) - Nr*log(eta_I);
    C = log_term;
    
    % 求解 u:
    % A*u - B*ln(u) - C = 0
    % u - (B/A)ln(u) = C/A
    % let t = -A/B * u  => u = -B/A * t
    % -B/A * t * A - B * ln(-B/A * t) = C
    % -B*t - B*ln(t) - B*ln(-B/A) = C
    % t + ln(t) = -C/B - ln(-B/A)
    % exp(t) * t = exp( -C/B - ln(-B/A) )
    % t = LambertW( ... )
    
    u_opt = nan;
    
    if abs(A) < 1e-10 % eta_I approx eta_Z (gamma -> 0)
        u_opt = exp(-C/B); % 退化解
    elseif abs(B) < 1e-6 % Nr approx nu_Z
        u_opt = C/A; 
    else
        % Lambert W 参数计算
        % z = exp( -C/B - log(abs(B/A)) ) * sign(-A/B) ??
        % 仔细推导:
        % (A/B)u - ln(u) = C/B
        % (A/B)u + ln(1/u) = C/B
        % let x = (A/B)u  => u = (B/A)x
        % x - ln(B/A x) = C/B
        % x - ln(x) = C/B + ln(B/A)
        % ln(x) - x = - (C/B + ln(B/A))
        % x * exp(-x) = exp( - (C/B + ln(B/A)) )
        % -x * exp(-x) = -exp(...)
        % -x = W( -exp(...) )
        % u = (B/A) * (-W(...)) = (-B/A) * W( - (A/B) * exp(-C/B) )
        
        val = - (A/B) * exp(-C/B);
        
        % Lambert W 主分支 (通常对应最优解)
        try
            w_val = lambertw(0, real(val));
            u_opt = -(B/A) * w_val;
        catch
            u_opt = nan;
        end
    end
    
    % 如果闭式解失效（数值问题），回退到 fminbnd
    if isnan(u_opt) || u_opt < 0 || isinf(u_opt)
        % 搜索范围：均值之间
        mu0 = Nr * eta_I;
        mu1 = nu_Z * eta_Z;
        calc_err = @(u) pi0*gammainc(u/eta_I, Nr, 'upper') + pi1*gammainc(u/eta_Z, nu_Z);
        [u_opt, xi_val] = fminbnd(calc_err, 0, mu1*2);
        xi = xi_val;
    else
        % 3. 计算检测错误概率
        % P_FA = gammainc(u/eta_I, Nr, 'upper')
        % P_MD = gammainc(u/eta_Z, nu_Z) (lower)
        
        P_FA = gammainc(u_opt/eta_I, Nr, 'upper');
        P_MD = gammainc(u_opt/eta_Z, nu_Z);
        
        xi = pi0 * P_FA + pi1 * P_MD;
    end
end