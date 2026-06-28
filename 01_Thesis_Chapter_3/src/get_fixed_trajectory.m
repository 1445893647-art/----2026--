function Q_fixed = get_fixed_trajectory(par)
    % GET_FIXED_TRAJECTORY 生成固定轨迹 (P&NT 方案)
    % 策略: 以匀速直线从起点飞向 Bob 的正上方 (高度 H_min)
    
    Q_fixed = zeros(3, par.N);
    
    % 设定终点: Bob 的 XY 坐标，高度为无人机最低飞行高度
    target_pos = par.pos_Bob_Est;
    target_pos(3) = par.H_min; 
    
    % 计算总位移向量
    vec_total = target_pos - par.pos_UAV_start;
    dist_total = norm(vec_total);
    
    % 计算匀速飞行的速度向量
    % 检查是否需要在 T 内到达，或者以 V_max 飞行
    % 论文中通常假设 T 时间内刚到达或飞越，这里采用 T 时间内匀速到达
    vel_req = dist_total / par.T;
    
    % 检查速度约束 (如果需要的速度超过最大速度，则限制为 V_max)
    if vel_req > par.V_max
        warning('固定轨迹所需速度 (%.2f m/s) 超过最大速度，将按最大速度飞行', vel_req);
        vel_vec = (vec_total / dist_total) * par.V_max;
    else
        vel_vec = vec_total / par.T;
    end
    
    % 生成轨迹点
    for n = 1:par.N
        curr_t = n * par.delta_t;
        % P[n] = P_start + v * t
        pos_n = par.pos_UAV_start + vel_vec * curr_t;
        
        % 确保高度不低于 H_min (虽然终点已设为 H_min，防止数值误差)
        pos_n(3) = max(pos_n(3), par.H_min);
        
        Q_fixed(:, n) = pos_n;
    end
end