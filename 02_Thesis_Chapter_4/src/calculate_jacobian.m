%% --- 辅助函数 2: 计算雅可比矩阵 H [cite: 103, 321-355] ---
function H = calculate_jacobian(x_pred, q_w, v_w, w_wb, P)
    % 获取当前预测值的测量参数
    [~, p] = measurement_function(x_pred, q_w, v_w, w_wb, P);
    
    % 提取参数
    d = p.dist; dx = p.diff_p(1); dy = p.diff_p(2); dz = p.diff_p(3);
    rho_sq = dx^2 + dy^2; rho = sqrt(rho_sq);
    rho = max(rho, 1e-8);
    rho_sq = max(rho_sq, 1e-12);
    
    %% 1. 对位置 q 的偏导 [cite: 325-337]
    % d_tau / d_q [cite: 326]
    dtau_dq = (2/ (P.c * d)) * [dx, dy, dz];
    
    % d_mu / d_q [cite: 329]
    term1 = 2 * P.fc / (P.c * d);
    dvx = x_pred(4) - v_w(1); dvy = x_pred(5) - v_w(2); dvz = x_pred(6) - v_w(3);
    dmu_dq_x = term1 * (dvx - p.v_r * dx / d); % 需注意公式329中的 Delta v 是相对速度
    dmu_dq_y = term1 * (dvy - p.v_r * dy / d);
    dmu_dq_z = term1 * (dvz - p.v_r * dz / d);
    dmu_dq = [dmu_dq_x, dmu_dq_y, dmu_dq_z];
    
    % d_r_tilde / d_q [cite: 331]
    % 需要计算 d_r / d_d, d_r / d_theta, d_r / d_phi
    
    % J_d, J_theta, J_phi [cite: 332-334]
    J_d = (1/d) * [dx, dy, dz];
    J_theta = (1/(rho * d^2)) * [dx*dz, dy*dz, -rho_sq]; % 注意附录公式333可能有笔误，分母通常是 d^2
    J_phi = (1/rho_sq) * [-dy, dx, 0];
    
    % d_r / d_d [cite: 335]
    beta_n = P.zeta_0 / (d^2);
    u = p.b_W * (p.a_W' * w_wb); % r = beta * u
    dr_dd = -2/d * (beta_n * u);
    
    % d_r / d_psi (psi = theta, phi) [cite: 336-337]
    [daW_dtheta, dbW_dtheta] = get_steering_deriv(p.theta, p.phi, P, 'theta');
    [daW_dphi, dbW_dphi] = get_steering_deriv(p.theta, p.phi, P, 'phi');
    
    dr_dtheta = beta_n * (dbW_dtheta * (p.a_W' * w_wb) + p.b_W * (daW_dtheta' * w_wb));
    dr_dphi   = beta_n * (dbW_dphi   * (p.a_W' * w_wb) + p.b_W * (daW_dphi'   * w_wb));
    
    dr_dq = dr_dd * J_d + dr_dtheta * J_theta + dr_dphi * J_phi;
    
    %% 2. 对速度 v 的偏导 [cite: 351]
    dtau_dv = [0, 0, 0];
    dmu_dv = (2 * P.fc / (P.c * d)) * [dx, dy, dz]; % [cite: 354]
    dr_dv = zeros(P.N_r, 3); % 回波向量与速度无关
    
    %% 3. 组装雅可比矩阵 H [cite: 103]
    % H 维度: (Nr + 2) x 6
    H_q = [dr_dq; dtau_dq; dmu_dq]; % (Nr+2) x 3
    H_v = [dr_dv; dtau_dv; dmu_dv]; % (Nr+2) x 3
    H = [H_q, H_v];
end