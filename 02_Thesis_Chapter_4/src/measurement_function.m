%% --- 辅助函数 1: 测量函数 h(x) [cite: 79-83] ---
function [y, params] = measurement_function(x_b, q_w, v_w, w_wb, P)
    % 提取状态
    p_b = x_b(1:3); v_b = x_b(4:6);
    
    % 几何关系 [cite: 62-64]
    diff_p = p_b - q_w;
    dist = norm(diff_p);
    theta = acos((p_b(3) - q_w(3)) / dist); % 俯仰角
    phi = atan2(diff_p(2), diff_p(1));      % 方位角
    
    % 相对速度
    diff_v = v_b - v_w;
    v_r = dot(diff_v, diff_p) / dist;
    
    % 导向矢量
    a_W = get_steering_vec(theta, phi, P.N_tx, P.N_ty);
    b_W = get_steering_vec(theta, phi, P.N_rx, P.N_ry);
    
    % 1. 回波向量 r_tilde [cite: 80]
    beta_n = P.zeta_0 / (dist^2);
    r_tilde = beta_n * b_W * (a_W' * w_wb);
    
    % 2. 时延 tau [cite: 81]
    tau = 2 * dist / P.c;
    
    % 3. 多普勒 mu [cite: 82]
    mu = 2 * v_r * P.fc / P.c;
    
    y = [r_tilde; tau; mu];
    
    % 返回中间参数供雅可比计算使用
    params.dist = dist; params.theta = theta; params.phi = phi;
    params.diff_p = diff_p; params.diff_v = diff_v; 
    params.v_r = v_r; params.a_W = a_W; params.b_W = b_W;
end


