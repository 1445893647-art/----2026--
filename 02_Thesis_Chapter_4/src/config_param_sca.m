% config_param_sca.m
% 在原 config_param.m 基础上补齐 SCA 参数（硬检测约束 + 回溯版）

config_param;  % 先加载原参数（会创建 global P）

if ~isfield(P,'sca'), P.sca = struct(); end

% 迭代与数值差分
P.sca.max_iter = 10;
P.sca.tol      = 1e-3;
P.sca.eps_q    = 5e-2;

% 信赖域（若仍出现 infeasible，可增大 trust_PJ 或减少 max_iter）
P.sca.trust_q    = 15;     % (m)
P.sca.trust_PJ   = 1.0;    % (W)
P.sca.trust_Wfro = 1.5;

% 信息矩阵正定裕量
P.sca.eps_pd = 1e-9;

% EKF 更新数值保护（用于 main_ekf_sim_SCA）
P.sca.var_floor = 1e-10;

% 可选：轻微追踪 Bob（0 表示不启用；建议 0.01~0.1）
P.sca.lambda_track = 0.02;

% 若未定义 Willie-Bob 路损指数，给默认
if ~isfield(P,'alpha_wb')
    P.alpha_wb = 2.2;
end

% EKF slack 惩罚系数（越大越逼近满足 EKF）
if ~isfield(P,'penalty_pcrb')
    P.penalty_pcrb = 1e4;
end
