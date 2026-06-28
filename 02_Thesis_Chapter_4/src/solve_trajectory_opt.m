function q_w_opt = solve_trajectory_opt(x_pred, P_cov_pred, q_w_prev, v_w, P_J, w_wb, P)
%SOLVE_TRAJECTORY_OPT Trajectory subproblem (P6) - CVX with soft constraints.
%
% This implementation follows the paper's idea of turning the Willie
% detection constraint into an explicit *distance-to-Alice* constraint:
%   ||q_a - q_w|| <= R_detect
% where
%   R_detect = (P_a*beta0 / (rho*P_J*Gamma_req))^(1/alpha_aw).
% (Gamma_req is solved from the Bhattacharyya bound once, e.g. in config.)
%
% The PCRB/tracking constraint here is implemented in the same spirit as your
% existing code: we derive a conservative maximum allowable distance to Bob
% (d_pcrb_max) from the current linearization point q_w_prev.

    % --- 1) Geometry / linearization point ---
    q_b = x_pred(1:3);

    % Use the same measurement model to get dist and beam-gain at q_w_prev
    [~, meas_p] = measurement_function(x_pred, q_w_prev, v_w, w_wb, P);
    dist_prev = max(meas_p.dist, 1e-3);

    % Beam matching gain |a^H w|^2 (appears in tau/mu noise)
    delta_sq = abs(meas_p.a_W' * w_wb)^2;
    delta_sq = max(delta_sq, 1e-12);

    % Jacobian (complex)
    H_jac = calculate_jacobian(x_pred, q_w_prev, v_w, w_wb, P);

    % --- 2) Build D and C matrices for PCRB constraint ---
    % Q1^{-1}: echo-vector noise (depends on P_J)
    var_1 = (P.C1 * P.sigma_r_sq) / (P.G * P_J);
    Q1_inv = diag([ones(P.N_r, 1) / var_1; 0; 0]);

    % Q2^{-1}: delay/doppler noise (depends on |beta|^2, |delta|^2, P_J)
    beta_sq = (P.zeta_0 / (dist_prev^2))^2; % |beta|^2
    denom_common = P.G * beta_sq * delta_sq * P_J;
    var_2 = (P.C2^2 * P.sigma_r_sq) / denom_common;
    var_3 = (P.C3^2 * P.sigma_r_sq) / denom_common;
    Q2_inv = diag([zeros(P.N_r, 1); 1/var_2; 1/var_3]);

    % D = inv(M_pred) + H' Q1^{-1} H
    M_inv = inv(P_cov_pred + 1e-9 * eye(6));
    D_mat = M_inv + H_jac' * Q1_inv * H_jac;
    D_mat = (D_mat + D_mat') / 2;

    % C = H' Q2^{-1} H
    C_mat = H_jac' * Q2_inv * H_jac;
    C_mat = (C_mat + C_mat') / 2;

    % --- 3) Solve conservative d_pcrb_max from linearization point ---
    % Generalized eigen-decomposition: C v = lambda D v
    [U_gen, Lambda_gen] = eig(C_mat, D_mat);
    lambda_vec = real(diag(Lambda_gen));

    % Normalize eigenvectors so that v' D v = 1
    l_vec = zeros(6, 1);
    for k = 1:6
        v_k = U_gen(:, k);
        denom = real(v_k' * D_mat * v_k);
        denom = max(denom, 1e-12);
        v_norm = v_k / sqrt(denom);
        l_vec(k) = norm(v_norm)^2;
        U_gen(:, k) = v_norm;
    end

    % Here we keep your original scalar form:
    %   sum( l_k / (1 + mu_k / d^4) ) <= eta_th
    % with mu_k calibrated from the current point.
    mu_vec = lambda_vec * (dist_prev^4);

    pcrb_func = @(d) sum(l_vec ./ (1 + mu_vec ./ (d^4 + 1e-12))) - P.eta_th;

    try
        d_pcrb_max = fzero(pcrb_func, [0.1, 2000]);
    catch
        if pcrb_func(0.1) > 0
            d_pcrb_max = 0.1;
        else
            d_pcrb_max = 2000;
        end
    end

    % --- 4) Detection constraint radius (paper closed-form) ---
    if isfield(P, 'Gamma_req') && ~isempty(P.Gamma_req)
        Gamma_req = P.Gamma_req;
    else
        Gamma_req = get_gamma_req(P);
    end
    R_detect = (P.P_a * P.beta_0 / (P.rho * P_J * Gamma_req))^(1 / P.alpha_aw);

    % Numerical safety
    R_detect = max(R_detect, 0);

    % --- 5) CVX step: move closer to Bob under mobility/height/detection/PCRB ---
    cvx_begin quiet
        variable q_next(3, 1)
        variable slack_det nonnegative
        variable slack_pcrb nonnegative

        % Objective: get close to Bob, but strongly penalize infeasibility
        obj_dist = norm(q_next - q_b);
        penalty = 1e3 * (slack_det + slack_pcrb);
        minimize( obj_dist + penalty )

        subject to
            norm(q_next - q_w_prev) <= P.V_max * P.dt;
            q_next(3) >= P.H_min;

            % Detection constraint (soft)
            norm(q_next - P.q_a) <= R_detect + slack_det;

            % Tracking/PCRB constraint (soft)
            norm(q_next - q_b) <= d_pcrb_max + slack_pcrb;
    cvx_end

    % --- 6) Output ---
    if strcmp(cvx_status, 'Solved') || strcmp(cvx_status, 'Inaccurate/Solved')
        q_w_opt = q_next;
    else
        warning('Trajectory optimization failed. Coasting.');
        q_w_opt = q_w_prev;
        q_w_opt(3) = max(q_w_opt(3), P.H_min);
    end
end
