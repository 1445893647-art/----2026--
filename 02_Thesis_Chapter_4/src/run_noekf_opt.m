function out = run_noekf_opt(P, x_true, F_mat, Qx, mode)
% mode = "oracle" or "lag1"
% This baseline does NOT run EKF and does NOT generate noisy measurements.
% It still propagates covariance deterministically (needed if your solvers use tracking/PCRB constraints).

N  = P.N_slots;
Nr = P.N_r;

% ---- init Willie state (same as main) ----
q_w_curr = P.q_w_0;
v_w_curr = [0;0;0];
P_J_curr = P.P_J;

w_wb_curr = get_steering_vec(0, 0, P.N_tx, P.N_ty);
w_wb_curr = w_wb_curr / norm(w_wb_curr);

% dummy Alice beam for metrics (same as your main)
w_ab_dummy = ones(P.N_a, 1);

% ---- init covariance (only for constraints) ----
P_cov = P.P_cov;

% ---- logs ----
history_P_out   = zeros(1, N);
history_xi      = zeros(1, N);
history_dist_wb = zeros(1, N);
history_P_J     = zeros(1, N);
q_w_hist        = zeros(3, N);
q_w_hist(:,1)   = q_w_curr;

% BCD params (keep same as main for fair compare)
max_bcd_iter = 10;
tol_bcd = 1e-3;

for n = 2:N
    % choose optimization state
    if mode == "oracle"
        x_opt = x_true(:,n);
    else % "lag1"
        x_opt = x_true(:,max(n-1,1));
    end

    % covariance prediction (for constraint use)
    P_pred = F_mat * P_cov * F_mat' + Qx;

    % ---- BCD optimize for this slot ----
    q_w_prev_slot = q_w_curr;
    q_w_opt = q_w_curr;
    v_w_opt = v_w_curr;
    P_J_opt = P_J_curr;
    w_wb_opt = w_wb_curr;

    for iter = 1:max_bcd_iter
        q_w_old = q_w_opt;
        P_J_old = P_J_opt;
        w_wb_old = w_wb_opt;
        % 2.3 trajectory (use x_opt)
        q_w_opt = solve_trajectory_opt(x_opt, P_pred, q_w_prev_slot, v_w_opt, P_J_opt, w_wb_opt, P);

        % 2.1 power
        P_J_opt = solve_power_opt(q_w_opt, P);

        % 2.2 beamforming (IMPORTANT: use x_opt, not x_true and not EKF)
        try
            w_wb_opt = solve_beamforming_opt(x_opt, P_pred, v_w_opt, q_w_opt, P_J_opt, w_wb_opt, P);
        catch ME
            % 保持上一次的波束
            warning('Beamforming failed: %s', ME.message);
            disp(getReport(ME,'extended'));   % <-- 这行会把真实栈信息打印出来
            q_b_tmp = x_opt(1:3);
            [th_tmp, ph_tmp, ~] = get_angle_dist(q_w_opt, q_b_tmp);
            w_wb_opt = get_steering_vec(th_tmp, ph_tmp, P.N_tx, P.N_ty);
            w_wb_opt = w_wb_opt / norm(w_wb_opt);
        end

        
        % update velocity
        v_w_opt = (q_w_opt - q_w_prev_slot) / P.dt;

        % convergence
        diff_q = norm(q_w_opt - q_w_old);
        diff_P = abs(P_J_opt - P_J_old);
        diff_w = norm(w_wb_opt - w_wb_old);            
        if diff_q < tol_bcd && diff_P < tol_bcd && diff_w < tol_bcd
            fprintf('Slot %d iter: %d, diff_q: %.6g, diff_P: %.6g, diff_w: %.6g\n',n, iter,diff_q,diff_P,diff_w);
            break; % 已收敛
        end
    end

    % execute
    q_w_curr  = q_w_opt;
    v_w_curr  = v_w_opt;
    P_J_curr  = P_J_opt;
    w_wb_curr = w_wb_opt;

    q_w_hist(:,n) = q_w_curr;
    history_P_J(n)= P_J_curr;

    % ---- deterministic covariance update (no innovation, no sampling) ----
    % Linearize at x_opt and current action, build H_r and Qy_r same way as main.
    H_jac = calculate_jacobian(x_opt, q_w_curr, v_w_curr, w_wb_curr, P);
    [~, params_opt] = measurement_function(x_opt, q_w_curr, v_w_curr, w_wb_curr, P);

    beta_r = P.zeta_0 / (params_opt.dist^2);
    delta_n = params_opt.a_W' * w_wb_curr;
    denom = P.G * abs(beta_r)^2 * abs(delta_n)^2 * P_J_curr;

    var_1 = (P.C1^2 * P.sigma_r_sq) / (P.G * P_J_curr);
    var_2 = (P.C2^2 * P.sigma_r_sq) / denom;
    var_3 = (P.C3^2 * P.sigma_r_sq) / denom;

    H_r = [real(H_jac(1:Nr,:));
           imag(H_jac(1:Nr,:));
           real(H_jac(Nr+1:end,:))];

    Qy_r = blkdiag((var_1/2)*eye(Nr), (var_1/2)*eye(Nr), var_2, var_3);

    S = H_r * P_pred * H_r' + Qy_r;
    K = P_pred * H_r' / S;

    % Joseph form covariance update (stable, no need for y)
    I6 = eye(6);
    P_cov = (I6 - K*H_r) * P_pred * (I6 - K*H_r)' + K*Qy_r*K';
    P_cov = (P_cov + P_cov')/2;

    % ---- metrics: ALWAYS evaluate using TRUE Bob position at slot n ----
    x_b_true = x_true(:,n);
    [p_out_val, xi_val] = calc_metrics(q_w_curr, x_b_true(1:3), P_J_curr, w_wb_curr, w_ab_dummy, P);

    history_P_out(n)   = p_out_val;
    history_xi(n)      = xi_val;
    history_dist_wb(n) = norm(q_w_curr - x_b_true(1:3));
    % 进度
    fprintf('[%s] Slot %d | Dist_wb: %.1fm | P_J: %.2fW | P_out: %.2f\n', ...
        char(mode), n, history_dist_wb(n), P_J_curr, p_out_val);
end

out.P_out   = history_P_out;
out.xi      = history_xi;
out.dist_wb = history_dist_wb;
out.P_J     = history_P_J;
out.q_w_hist= q_w_hist;
out.mode    = mode;
end
