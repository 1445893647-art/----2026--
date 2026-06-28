function P_out_avg = calc_outage_prob(Q, P, par)
    % CALC_OUTAGE_PROB 计算平均隐蔽通信中断概率 (Eq. 7)
    % Q: 3xN 轨迹矩阵
    % P: 1xN 功率向量
    
    P_out_sum = 0;
    
    for n = 1:par.N
        % 计算各链路距离
        % 真实距离
        d_AB = norm(par.pos_Alice_True - par.pos_Bob_True);
        d_UB = norm(Q(:,n) - par.pos_Bob_True);
        d_LB = norm(par.pos_LT - par.pos_Bob_True);
        
        % 根据公式 (7) 计算 k1, k2, k3
        % 注意: 信道增益 |h|^2 均值为 beta0 / d^alpha
        % k = (beta0 * Power * sigma^2_ch) / (Constant * d^alpha)
        % 此处 sigma^2_ch 归一化为 1
        
        k1 = (par.beta0 * par.P_Alice) / ((2^par.r_target - 1) * d_AB^par.alpha0);
        k2 = (par.beta0 * P(n)) / (d_UB^par.alpha0);
        k3 = (par.beta0 * par.P_LT) / (d_LB^par.alpha0);
        
        % Eq. 7: Closed-form Outage Probability
        term = (k1^2) / ((k1 + k2) * (k1 + k3));
        P_out_n = 1 - term * exp(-par.sigma2_noise / (2*k1));
        
        P_out_sum = P_out_sum + P_out_n;
    end
    
    P_out_avg = P_out_sum / par.N;
end