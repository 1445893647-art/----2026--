function [q_w_opt, v_w_opt, P_J_opt, w_wb_opt, info] = solve_joint_sca_opt(x_pred, P_pred, q_w_prev, P_J_init, w_wb_init, P)
% ---- CVX 保护：避免上一次 CVX 异常退出导致的“已有 cvx problem”警告 ----
try
    cvx_clear;  % 清空可能残留的 CVX problem（warning: A non-empty cvx problem...）
catch
end

%======================================================================
% solve_joint_sca_opt.m  (v4: 检测约束用“乘积线性化”替代 1/d^alpha 线性化)
%
% 说明（针对“rho 变大后 SCA 原地摇摆”）：
%  - 你论文的检测约束来自：gamma_w >= Gamma_req（由检测错误概率约束等价得到）
%    且自干扰项服从 Gamma，均值 eta_I = rho * P_J（rho 为残留自干扰系数）。
%  - rho 变大 => 检测约束更紧：必须减小 P_J 或者靠近 Alice（增大 beta_aw）。
%  - 若在 SCA 中用 PJ <= c_det * ||qa-q||^{-alpha} 的倒数形式做一阶“内逼近”，会非常保守，
%    常导致 q 不敢远离 Alice（看起来“原地摇摆”）。
%  - 更稳的 SCA：改用等价形式 PJ * ||qa-q||^{alpha} <= const_det，并对“乘积”做线性化。
%
% 等价变换（与你初稿 Lemma 1 的距离约束一致）：
%   ||qa-q|| <= (Pa*beta0/(rho*PJ*Gamma_req))^(1/alpha)
%   <=> PJ * ||qa-q||^{alpha} <= Pa*beta0/(rho*Gamma_req) = const_det
%======================================================================

info = struct();
info.cvx_status = "";
info.iters = 0;
info.backtrack_trials = 0;
info.slack_ekf = NaN;
info.det_ub = NaN;

q_b_pred = x_pred(1:3);
Nt = P.N_t;

% SCA 参数
if ~isfield(P,'sca'), P.sca = struct(); end
maxIter = getfield_def(P.sca,'max_iter', 10);
tol     = getfield_def(P.sca,'tol', 1e-3);
eps_q   = getfield_def(P.sca,'eps_q', 5e-2);
trust_q0= getfield_def(P.sca,'trust_q', 15);
trust_P0= getfield_def(P.sca,'trust_PJ', 5.0);   % 建议更大：允许 PJ 随 q 明显下降
trust_W0= getfield_def(P.sca,'trust_Wfro', 1.5);
eps_pd  = getfield_def(P.sca,'eps_pd', 1e-9);

lambda_track = getfield_def(P.sca,'lambda_track', 0.0); % baseline 建议 0
penalty_pcrb = getfield_def(P,'penalty_pcrb', 1e4);

% 初值
q0  = q_w_prev;
PJ0 = min(max(P_J_init,0), P.P_max);

if nargin >= 5 && ~isempty(w_wb_init)
    w0 = w_wb_init(:); w0 = w0/(norm(w0)+1e-12);
else
    [th0, ph0, ~] = get_angle_dist(q0, q_b_pred);
    w0 = get_steering_vec(th0, ph0, P.N_tx, P.N_ty);
    w0 = w0/(norm(w0)+1e-12);
end
W0 = w0*w0';

% 先验信息矩阵
B0 = inv(P_pred + 1e-9*eye(6));

% 检测约束常数（等价形式）
qa = P.q_a(:);
const_det = (P.P_a * P.beta_0) / (P.rho * P.Gamma_req + 1e-12);  % Pa*beta0/(rho*Gamma_req)

% PJ0 截断到“真实上界”（以免一上来 infeasible）
r0 = norm(qa - q0) + 1e-12;
u0 = r0^(P.alpha_aw);
PJ_ub0 = min(P.P_max, const_det / u0);
PJ0 = min(PJ0, 0.95*PJ_ub0);
info.det_ub = PJ_ub0;

% 回溯信赖域
trust_q = trust_q0; trust_P = trust_P0; trust_W = trust_W0;

for it = 1:maxIter
    info.iters = it;

    v0 = (q0 - q_w_prev)/P.dt;

    % ===== 冻结雅可比，构造 K0, D0（保持“第一块只含 PJ”结构）=====
    H_jac = calculate_jacobian(x_pred, q0, v0, w0, P);
    Nr = P.N_r;
    H_r = [real(H_jac(1:Nr,:));
           imag(H_jac(1:Nr,:));
           real(H_jac(Nr+1:end,:))];

    H1 = H_r(1:2*Nr,:);
    h2 = H_r(2*Nr+1,:);
    h3 = H_r(2*Nr+2,:);

    K0 = (2*P.G)/(P.C1^2 * P.sigma_r_sq) * (H1' * H1);
    D0 = (P.G)/(P.sigma_r_sq) * ( (1/(P.C2^2))*(h2'*h2) + (1/(P.C3^2))*(h3'*h3) );
    K0 = (K0 + K0')/2;
    D0 = (D0 + D0')/2;

    % ===== g(q,W)=|beta_r|^2*|delta|^2 线性化 =====
    [g0, c0, A0] = g_c_A(q0, q_b_pred, W0, P);
    gradg = grad_numeric(@(q) g_only(q, q_b_pred, W0, P), q0, eps_q);

    % ===== 目标代理：PJ*s(q,W) =====
    [s0, beta_wb0, Aobj0] = s_interf(q0, q_b_pred, W0, P);
    grads = grad_numeric(@(q) s_only(q, q_b_pred, W0, P), q0, eps_q);

    solved = false;
    for trial = 1:6
        info.backtrack_trials = trial;

        % 线性化点处 u0（u=||qa-q||^alpha）
        r0 = norm(qa - q0) + 1e-12;
        u0 = r0^(P.alpha_aw);

        s_ekf_val = NaN; % 若 CVX 未成功运行/报错，用 NaN 占位
            cvx_clear;  % 防止残留 CVX 模型

        try

            cvx_begin sdp quiet
            variable W(Nt,Nt) hermitian semidefinite
            variables q(3,1) PJ
            variable s_ekf
            variable r nonnegative
            variable u nonnegative

            s_ekf >= 0;

            % 信赖域
            norm(q - q0, 2) <= trust_q;
            abs(PJ - PJ0) <= trust_P;
            norm(W - W0, 'fro') <= trust_W;

            trace(W) == 1;
            0 <= PJ <= P.P_max;

            norm(q - q_w_prev, 2) <= P.V_max * P.dt;
            q(3) >= P.H_min;

            % ===== C3 检测约束：PJ * ||qa-q||^alpha <= const_det =====
            r >= norm(q - qa, 2);
            u >= pow_pos(r, P.alpha_aw);

            % 乘积线性化：PJ*u ≈ PJ0*u + u0*PJ - PJ0*u0  <= const_det
            PJ0*u + u0*PJ - PJ0*u0 <= const_det;

            % ===== EKF 约束（允许 slack）=====
            g_aff = c0 * real(trace(A0 * W)) + gradg'*(q - q0);
            g_aff >= 0;

            t_hat = PJ0 * g_aff + g0 * PJ - PJ0 * g0;
            t_hat >= 0;

            J = B0 + PJ*K0 + t_hat*D0;
            J >= eps_pd * eye(6);

            trace_inv(J) <= P.eta_th + s_ekf;

            % ===== 目标 =====
            s_aff = beta_wb0 * ( (P.K/(P.K+1))*real(trace(Aobj0*W)) + (1/(P.K+1))*real(trace(P.R_wb*W)) ) ...
                    + grads'*(q - q0);
            obj_hat = PJ0 * s_aff + s0 * PJ - PJ0 * s0;

            track_term = -lambda_track * norm(q - q_b_pred, 2);

            maximize( obj_hat + track_term - penalty_pcrb*s_ekf )
            cvx_end
            % 取出 slack 数值（cvx 变量在 cvx_end 后可读取）
            try s_ekf_val = s_ekf; catch, s_ekf_val = NaN; end
        catch ME_cvx
            % 若 CVX 内部报错（如 pow_pos/Disciplined Convex Programming 不满足），清空并进入回溯
            cvx_clear;
            cvx_status = "Error";
        end
        info.cvx_status = cvx_status;
        info.slack_ekf  = s_ekf_val;

        if contains(cvx_status,'Solved')
            solved = true;
            break;
        end

        % 回溯：缩小信赖域 & 重置 PJ0 到检测上界
        trust_q = max(1.0, trust_q * 0.5);
        trust_P = max(0.1, trust_P * 0.5);
        trust_W = max(0.2, trust_W * 0.7);

        r0 = norm(qa - q0) + 1e-12;
        u0 = r0^(P.alpha_aw);
        PJ_ub0 = min(P.P_max, const_det / u0);
        PJ0 = min(PJ0, 0.95*PJ_ub0);
    end

    if ~solved
        warning('SCA-CVX failed. status=%s', cvx_status);
        break;
    end

    % 从 W 恢复 w
    w1 = principal_eigvec(W);
    w1 = w1/(norm(w1)+1e-12);

    q1 = q;
    PJ1 = PJ;

    % 解后真实 C3 投影
    r1 = norm(qa - q1) + 1e-12;
    u1 = r1^(P.alpha_aw);
    PJ_ub1 = min(P.P_max, const_det / u1);
    PJ1 = min(PJ1, 0.999*PJ_ub1);

    % 收敛
    rel_q  = norm(q1 - q0)/(norm(q0)+1e-9);
    rel_PJ = abs(PJ1 - PJ0)/(abs(PJ0)+1e-9);

    q0 = q1; PJ0 = PJ1; w0 = w1; W0 = W;

    if max(rel_q, rel_PJ) < tol
        break;
    end
end

q_w_opt = q0;
v_w_opt = (q_w_opt - q_w_prev)/P.dt;
P_J_opt = PJ0;
w_wb_opt = w0;

end

%========================== Helper functions ==========================

function v = getfield_def(s, name, def)
if isstruct(s) && isfield(s,name), v = s.(name); else, v = def; end
end

function w = principal_eigvec(W)
W = (W+W')/2;
[V,D] = eig(W);
[~,idx] = max(real(diag(D)));
w = V(:,idx);
end

function g = grad_numeric(fun, q0, eps_q)
g = zeros(3,1);
for i = 1:3
    dq = zeros(3,1); dq(i)=eps_q;
    g(i) = (fun(q0 + dq) - fun(q0 - dq)) / (2*eps_q);
end
end

function val = g_only(q, qb, W0, P)
[val,~,~] = g_c_A(q, qb, W0, P);
end

function [g0, c0, A0] = g_c_A(q, qb, W0, P)
[th, ph, d] = get_angle_dist(q, qb);
a = get_steering_vec(th, ph, P.N_tx, P.N_ty);
a = a/(norm(a)+1e-12);
A0 = a*a';
beta_r = P.zeta_0 / (d^2 + 1e-12);
c0 = abs(beta_r)^2;
s0 = real(trace(A0 * W0));
g0 = c0 * s0;
end

function val = s_only(q, qb, W0, P)
[val,~,~] = s_interf(q, qb, W0, P);
end

function [s0, beta_wb0, A0] = s_interf(q, qb, W0, P)
[th, ph, d] = get_angle_dist(q, qb);
a = get_steering_vec(th, ph, P.N_tx, P.N_ty);
a = a/(norm(a)+1e-12);
A0 = a*a';
beta_wb0 = P.beta_0 * (d + 1e-12)^(-P.alpha_wb);
s_gain = (P.K/(P.K+1))*real(trace(A0*W0)) + (1/(P.K+1))*real(trace(P.R_wb*W0));
s0 = beta_wb0 * s_gain;
end