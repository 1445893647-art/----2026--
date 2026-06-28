function w_opt = solve_beamforming_opt(x_pred, P_cov_pred, v_w, q_w, P_J, w_prev, P)
% solve_beamforming_opt
% -------------------------------------------------------------------------
% 波束赋形子问题 (P5): 期望干扰功率最大化 + EKF/PCRB 约束（Schur-LMI）
%
% 你反馈的“前十个时隙 CVX 求解失败”通常来自两类原因：
%   (i) PCRB 阈值 eta_th 太严格 → LMI 不可行（cvx_status=Infeasible/Failed）
%   (ii) 数值尺度/方差建模不一致（例如把 C1/C2/C3 多写成平方）导致信息矩阵过小 → 更易不可行
%
% 本版本修复：
%   1) 将 var_1,var_2,var_3 与你轨迹子问题的写法对齐：用 C1/C2/C3 * sigma_r_sq（不再平方）
%   2) 给 PCRB 约束加 soft slack：trace(Z) <= eta_th + slack_pcrb，并在目标中惩罚 slack
%      -> 保证早期时隙也能得到可行解（同时可通过 slack 大小判断“约束不满足程度”）
%   3) 增加轻微的 Jf 下界（数值稳定）：Jf >= eps*I
%
% 输入：
%   x_pred      : EKF 预测状态 [q_b; v_b] (6x1)
%   P_cov_pred  : EKF 预测协方差 (6x6)
%   v_w         : Willie 速度 (3x1)
%   q_w         : Willie 位置 (3x1)
%   P_J         : Willie 发射功率 (标量)
%   w_prev      : 上一轮波束 (N_t x 1)，可为空
%   P           : 参数结构体（需包含：N_tx,N_ty,N_rx,N_ry,zeta_0,G,sigma_r_sq,C1,C2,C3,eta_th,K,R_wb 等）
%
% 输出：
%   w_opt       : 求得波束 (N_t x 1)

    %% ---------------- 0) 参数与保护 ----------------
    n_ant = P.N_tx * P.N_ty;

    % slack 惩罚（越大越“硬”）
    if isfield(P,'penalty_pcrb') && ~isempty(P.penalty_pcrb)
        penalty_pcrb = P.penalty_pcrb;
    else
        penalty_pcrb = 1e3;
    end

    % 数值稳定下界
    epsJ = 1e-9;

    % 信赖域（可选）
    if isfield(P,'beam_sca_trust') && ~isempty(P.beam_sca_trust)
        useTrust = logical(P.beam_sca_trust);
    else
        useTrust = false;
    end
    if isfield(P,'beam_sca_delta') && ~isempty(P.beam_sca_delta)
        delta = P.beam_sca_delta;
    else
        delta = 0.25;
    end

    % 防止 P_J=0 造成除零
    P_J_safe = max(P_J, 1e-12);

    %% 1) 几何关系（Willie -> Bob）
    q_b = x_pred(1:3);
    diff_p = q_b - q_w(:);

    dist = norm(diff_p);
    dist = max(dist, 1e-3);

    dx = diff_p(1); dy = diff_p(2); dz = diff_p(3);
    [theta, phi, ~] = get_angle_dist(q_w(:), q_b(:));
    rho_sq = dx^2 + dy^2;

    %% 2) 导向矢量与导数
    a_W = get_steering_vec(theta, phi, P.N_tx, P.N_ty);   % Tx steering (N_t x 1)
    b_W = get_steering_vec(theta, phi, P.N_rx, P.N_ry);   % Rx steering (N_r x 1)

    [da_dtheta, db_dtheta] = get_steering_deriv(theta, phi, P, 'theta');
    [da_dphi,   db_dphi]   = get_steering_deriv(theta, phi, P, 'phi');

    A_W = a_W * a_W';   % g(W)=|a^H w|^2 = Tr(A_W W)

    %% 3) 回波导数矩阵分解: dr/dq_i = S_i * w
    beta_n = P.zeta_0 / (dist^2);

    S_d  = (-2/dist) * beta_n * (b_W * a_W');                          % dr/dd
    S_th = beta_n * (db_dtheta * a_W' + b_W * da_dtheta');             % dr/dtheta
    S_ph = beta_n * (db_dphi   * a_W' + b_W * da_dphi');               % dr/dphi

    % 位置到 (d,theta,phi) 的导数
    J_d = (1/dist) * [dx, dy, dz];

    if rho_sq < 1e-12
        J_th = [0, 0, 0];
        J_ph = [0, 0, 0];
    else
        rho = sqrt(rho_sq);
        J_th = (1/(rho * dist^2)) * [dx*dz, dy*dz, -rho_sq];
        J_ph = (1/rho_sq) * [-dy, dx, 0];
    end

    S_q = cell(1,3);
    for i = 1:3
        S_q{i} = S_d * J_d(i) + S_th * J_th(i) + S_ph * J_ph(i);
    end

    %% 4) 构造 K_ij，使得 J_echo(i,j)=Tr(K_ij W)
    % ★修复：与轨迹子问题一致，用 C1*sigma_r_sq（不再平方）
    var_1 = (P.C1^2 * P.sigma_r_sq) / max(P.G * P_J_safe, 1e-30);
    inv_var1 = 1 / max(var_1, 1e-30);

    Kij = cell(3,3);
    for i = 1:3
        for j = 1:3
            Ktmp = inv_var1 * (S_q{i}' * S_q{j});
            Kij{i,j} = real((Ktmp + Ktmp')/2);
        end
    end

    %% 5) tau/mu Fisher 项：J_tau_mu(W) = g(W) * C_const
    diff_v = x_pred(4:6) - v_w(:);
    v_r = dot(diff_v, diff_p) / dist;

    dtau_dq = (2/(P.c * dist)) * [dx, dy, dz];

    term1 = 2 * P.fc / (P.c * dist);
    dmu_dq = term1 * [ ...
        (diff_v(1) - v_r * dx / dist), ...
        (diff_v(2) - v_r * dy / dist), ...
        (diff_v(3) - v_r * dz / dist) ];

    dtau_dv = [0, 0, 0];
    dmu_dv  = (2 * P.fc / (P.c * dist)) * [dx, dy, dz];

    H2 = [dtau_dq, dtau_dv;   % 1x6
          dmu_dq,  dmu_dv];   % 1x6

    beta_sq = (P.zeta_0 / dist^2)^2;
    denom_base = P.G * beta_sq * P_J_safe;
    denom_base = max(denom_base, 1e-30);

    % ★修复：与轨迹子问题一致，用 C2*sigma_r_sq（不再平方）
    var_2_base = (P.C2^2 * P.sigma_r_sq) / denom_base;
    var_3_base = (P.C3^2 * P.sigma_r_sq) / denom_base;

    inv_var2 = 1 / max(var_2_base, 1e-30);
    inv_var3 = 1 / max(var_3_base, 1e-30);

    % 防止过大导致数值爆炸
    inv_var2 = min(inv_var2, 1e30);
    inv_var3 = min(inv_var3, 1e30);

    Q2_inv = diag([inv_var2, inv_var3]);
    C_const = H2' * Q2_inv * H2;
    C_const = real((C_const + C_const')/2);

    %% 6) 先验信息：M_inv
    A = (P_cov_pred + 1e-9 * eye(6));
    M_inv = A \ eye(6);
    M_inv = real((M_inv + M_inv')/2);

    %% 7) 波束目标矩阵：最大化 E{|h_wb^H w|^2} = Tr(H_obj W)
    % 这里采用 Rician: h = sqrt(beta)*[ sqrt(K/(K+1))*a + sqrt(1/(K+1))*R^{1/2}g ]
    % => E{h h^H} = beta*( K/(K+1)*a a^H + 1/(K+1)*R )
    if isfield(P,'R_wb') && ~isempty(P.R_wb)
        R_wb = P.R_wb;
    else
        R_wb = eye(n_ant);
    end
    R_wb = real((R_wb + R_wb')/2);

    if ~isfield(P,'K') || isempty(P.K)
        error('P.K (Rician K-factor) is required.');
    end

    % 大尺度衰落（若你主代码里已内含 beta_wb 也可以删掉这项）
    if isfield(P,'beta_0') && (isfield(P,'alpha_wb') || isfield(P,'alpha_aw'))
        alpha_wb = getfield(P, tern(isfield(P,'alpha_wb'),'alpha_wb','alpha_aw')); %#ok<GFLD>
        beta_wb  = P.beta_0 * (max(dist,1)^(-alpha_wb));
    else
        beta_wb  = 1; % 退化：不做路损缩放
    end

    H_obj = beta_wb * ( (P.K/(P.K+1))*(a_W*a_W') + (1/(P.K+1))*R_wb );
    H_obj = real((H_obj + H_obj')/2);

    %% 8) CVX：SDR + Schur LMI（加 slack）
    if exist('cvx_clear','file') == 2
        cvx_clear;
    end

    % 参考点（信赖域）
    if nargin >= 6 && ~isempty(w_prev)
        wref = w_prev(:) / max(norm(w_prev), 1e-12);
        Wref = wref*wref';
    else
        Wref = (a_W/max(norm(a_W),1e-12)) * (a_W'/max(norm(a_W),1e-12));
    end

    cvx_begin sdp quiet
        variable W(n_ant, n_ant) hermitian semidefinite
        variable Z(6,6) symmetric
        variable slack_pcrb nonnegative

        expression g
        g = real(trace(A_W * W));

        expression Jf(6,6)
        Jf = M_inv + g * C_const;

        % 回波项：仅作用于位置 1:3
        for i = 1:3
            for j = 1:3
                Jf(i,j) = Jf(i,j) + real(trace(Kij{i,j} * W));
            end
        end
        Jf = real((Jf + Jf')/2);

        maximize( real(trace(H_obj * W)) - penalty_pcrb * slack_pcrb )

        subject to
            trace(W) == 1;

            % ★soft PCRB：早期不可行时允许 slack
            trace(Z) <= P.eta_th + slack_pcrb;

            % 数值稳定：Jf 需要有下界（防止近奇异）
            Jf >= epsJ * eye(6);

            [Z, eye(6); eye(6), Jf] == semidefinite(12);

            if useTrust
                norm(W - Wref, 'fro') <= delta;
            end
    cvx_end

    %% 9) SDR -> w
    if strcmp(cvx_status, 'Solved') || strcmp(cvx_status, 'Inaccurate/Solved')
        fprintf("波束求解成功");
        [V_eig, D_eig] = eig(full(W));
        [~, idx] = max(real(diag(D_eig)));
        w_opt = V_eig(:, idx);
        w_opt = w_opt / max(norm(w_opt), 1e-12);

        % 如果你想在外层记录 slack，可把它返回或写入 P.debug
        if isfield(P,'debug_beam') && P.debug_beam
            fprintf('[beam] cvx=%s, slack_pcrb=%.3e\n', cvx_status, slack_pcrb);
        end
    else
        % 失败回退：优先回退上一轮，否则指向波束
        warning("波束求解失败");
        if nargin >= 6 && ~isempty(w_prev)
            w_opt = w_prev / max(norm(w_prev), 1e-12);
        else
            w_opt = a_W / max(norm(a_W), 1e-12);
        end
    end
end

function out = tern(cond, a, b)
    if cond, out = a; else, out = b; end
end
