% =========================================================================
%  Result Visualization & Metric Calculation (Updated with R_LR)
%  (Run this after main.m)
% =========================================================================

% 检查必要变量
if ~exist('par', 'var') || ~exist('Q_curr', 'var') || ~exist('P_curr', 'var')
    error('请先运行 main.m 以获取优化结果 (par, Q_curr, P_curr)');
end

%% 1. 计算逐时隙的性能指标 (Per-slot Metrics)
P_out_vec = zeros(1, par.N); % 中断概率
Xi_vec    = zeros(1, par.N); % 检测错误概率
R_LR_vec  = zeros(1, par.N); % 合法用户速率 [新增]

% --- 预计算 R_LR 相关的常量 (信号与固定干扰) ---
% 信号: LT -> LR
d_LT_LR = norm(par.pos_LT - par.pos_LR);
h_LT_LR_sq = par.beta0 / d_LT_LR^par.alpha0;
Signal_LR = par.P_LT * h_LT_LR_sq;

% 干扰: Alice -> LR
d_Alice_LR = norm(par.pos_Alice_True - par.pos_LR);
h_Alice_LR_sq = par.beta0 / d_Alice_LR^par.alpha0;
Interference_Alice = par.P_Alice * h_Alice_LR_sq;

% 检测概率计算所需常量
term_inner = (par.p_prior - par.epsilon_det)/par.p_prior;

for n = 1:par.N
    % --- (1) 计算中断概率 P_out[n] ---
    d_AB = norm(par.pos_Alice_True - par.pos_Bob_True);
    d_UB = norm(Q_curr(:,n) - par.pos_Bob_True);
    d_LB = norm(par.pos_LT - par.pos_Bob_True);
    
    k1 = (par.beta0 * par.P_Alice) / ((2^par.r_target - 1) * d_AB^par.alpha0);
    k2 = (par.beta0 * P_curr(n)) / (d_UB^par.alpha0);
    k3 = (par.beta0 * par.P_LT) / (d_LB^par.alpha0);
    
    term_out = (k1^2) / ((k1 + k2) * (k1 + k3));
    P_out_vec(n) = 1 - term_out * exp(-par.sigma2_noise / (2*k1));
    
    % --- (2) 计算检测错误概率 Xi[n] ---
    d_AU = norm(Q_curr(:,n) - par.pos_Alice_True);
    d_LU = norm(par.pos_LT - Q_curr(:,n));
    
    zeta_au = (par.beta0 * par.P_Alice) / (d_AU^par.alpha0);
    zeta_uu = par.rho * P_curr(n) * par.lambda_uu;
    A_n     = (par.P_LT * par.beta0) / (d_LU^par.alpha0);
    
    if zeta_uu < 1e-12
        Xi_vec(n) = 1 - par.p_prior; 
    else
        ln_term = log( (1/par.p_prior + 1)*(zeta_au/zeta_uu) - 1/par.p_prior + 2 );
        coef = (zeta_au * zeta_uu) / (zeta_au + zeta_uu);
        tau_star = coef * ln_term + par.sigma2_noise + A_n;
        
        denom_w = zeta_uu - zeta_au;
        if abs(denom_w) < 1e-9, denom_w = 1e-9; end
        
        omega1 = 1 - par.p_prior * (1 + zeta_uu/denom_w);
        omega2 = par.p_prior * (zeta_au/denom_w);
        
        term_exp1 = exp( -(tau_star - par.sigma2_noise - A_n) / (2*zeta_uu) );
        term_exp2 = exp( -(tau_star - par.sigma2_noise - A_n) / (2*zeta_au) );
        Xi_vec(n) = real(par.p_prior + omega1 * term_exp1 + omega2 * term_exp2);
    end
    
    % --- (3) 计算合法用户速率 R_LR[n] [新增] ---
    % 干扰: UAV -> LR (随时间变化)
    d_UAV_LR = norm(Q_curr(:,n) - par.pos_LR);
    h_UAV_LR_sq = par.beta0 / d_UAV_LR^par.alpha0;
    Interference_UAV = P_curr(n) * h_UAV_LR_sq;
    
    % SINR = Signal / (I_Alice + I_UAV + Noise)
    SINR_LR = Signal_LR / (Interference_Alice + Interference_UAV + par.sigma2_noise);
    R_LR_vec(n) = log2(1 + SINR_LR);
end

% 时间轴
time_axis = (0:par.N-1) * par.delta_t;
set(0, 'DefaultAxesFontSize', 12);
set(0, 'DefaultLineLineWidth', 1.5);

%% 2. 绘图 (Plotting)

% 图 1: 收敛图
figure('Name', 'Convergence');
yyaxis left; plot(1:length(history.obj), history.obj, '-s');
xlabel('Iter'); ylabel('Outage Prob');
yyaxis right; semilogy(1:length(history.diff_Q), history.diff_Q, '-o');
ylabel('Traj Diff'); title('Convergence');

% 图 2: 3D 轨迹
figure('Name', '3D Trajectory');
plot3(par.pos_Alice_True(1), par.pos_Alice_True(2), par.pos_Alice_True(3), 'rp', 'MarkerSize',12,'MarkerFaceColor','r'); hold on; grid on;
plot3(par.pos_Bob_True(1), par.pos_Bob_True(2), par.pos_Bob_True(3), 'bs', 'MarkerSize',12,'MarkerFaceColor','b');
plot3(par.pos_LT(1), par.pos_LT(2), par.pos_LT(3), 'gd', 'MarkerSize',12,'MarkerFaceColor','g');
plot3(par.pos_LR(1), par.pos_LR(2), par.pos_LR(3), 'mo', 'MarkerSize',12,'MarkerFaceColor','m');
plot3(Q_curr(1,:), Q_curr(2,:), Q_curr(3,:), 'k.-');
plot3(Q_curr(1,1), Q_curr(2,1), Q_curr(3,1), 'k^', 'MarkerSize',10,'MarkerFaceColor','k');
legend('Alice','Bob','LT','LR','Traj','Start'); view(-30, 30);
title('UAV 3D Trajectory');

% 图 3: 中断概率
figure('Name', 'Outage Probability');
plot(time_axis, P_out_vec, 'b-o'); grid on;
xlabel('Time (s)'); ylabel('P_{out}'); title('Outage Probability vs Time');

% 图 4: 检测错误概率
figure('Name', 'Detection Error Probability');
plot(time_axis, Xi_vec, 'r-d'); hold on;
yline(par.epsilon_det, 'k--', 'Constraint \epsilon'); grid on;
xlabel('Time (s)'); ylabel('\xi^*'); title('Detection Error Prob vs Time');

% === 图 5: 合法用户通信速率 (新增图表) ===
figure('Name', 'Legitimate User Rate');
plot(time_axis, R_LR_vec, 'g-s', 'LineWidth', 1.5, 'MarkerSize', 4); hold on;

% 绘制约束阈值线 R_min
yline(par.R_min_LR, 'r--', 'LineWidth', 2, 'DisplayName', 'Min Rate R_{min}');

grid on;
xlabel('时间 (Time / s)');
ylabel('合法用户速率 R_{LR} (bit/s/Hz)');
title('合法用户通信速率随时间变化');
legend('Rate R_{LR}', 'Constraint R_{min}');
% 设置 Y 轴下限，确保能看到约束线
ylim([min(min(R_LR_vec), par.R_min_LR)*0.9, max(R_LR_vec)*1.1]);