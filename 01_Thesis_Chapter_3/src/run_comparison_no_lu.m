% =========================================================================
%  Run Scheme WITHOUT Legitimate User Constraints
% =========================================================================

% 确保 par 已存在 (如果未运行 main.m，请先加载配置)
if ~exist('par', 'var'), par = config_params(); end

% 初始化
Q_no_lu = zeros(3, par.N);
for n = 1:par.N
    direction = (par.pos_Bob_Est - par.pos_UAV_start); direction(3) = 0; 
    step_vec = (direction / norm(direction)) * (par.V_max * 0.5 * par.delta_t);
    Q_no_lu(:, n) = par.pos_UAV_start + step_vec * (n-1);
    Q_no_lu(3, n) = max(Q_no_lu(3, n), par.H_min);
end
P_no_lu = par.P_UAV_max * ones(1, par.N);

fprintf('开始运行无合法用户约束方案...\n');

% --- AO 循环 ---
for iter = 1:par.MaxIter
    Q_prev = Q_no_lu;
    P_prev = P_no_lu;
    
    % 1. 优化功率 (无 LU)
    P_no_lu = solve_power_no_lu(Q_no_lu, par);
    
    % 2. 优化轨迹 (无 LU)
    [Q_no_lu, ~] = solve_traj_no_lu(P_no_lu, Q_prev, par);
    
    % 收敛判断
    if norm(P_no_lu - P_prev)/norm(P_prev) < par.tol_var && ...
       norm(Q_no_lu - Q_prev, 'fro')/norm(Q_prev, 'fro') < par.tol_var
        break;
    end
end

% --- 计算该方案下的性能 ---
P_out_no_lu_vec = zeros(1, par.N);
R_LR_no_lu_vec  = zeros(1, par.N);

% 计算常数
d_LT_LR = norm(par.pos_LT - par.pos_LR);
Signal_LR = par.P_LT * (par.beta0 / d_LT_LR^par.alpha0);

for n = 1:par.N
    % 1. 中断概率 (使用真实位置)
    d_AB = norm(par.pos_Alice_True - par.pos_Bob_True);
    d_UB = norm(Q_no_lu(:,n) - par.pos_Bob_True);
    d_LB = norm(par.pos_LT - par.pos_Bob_True);
    
    k1 = (par.beta0 * par.P_Alice) / ((2^par.r_target - 1) * d_AB^par.alpha0);
    k2 = (par.beta0 * P_no_lu(n)) / (d_UB^par.alpha0);
    k3 = (par.beta0 * par.P_LT) / (d_LB^par.alpha0);
    
    term = (k1^2) / ((k1 + k2) * (k1 + k3));
    P_out_no_lu_vec(n) = 1 - term * exp(-par.sigma2_noise / (2*k1));
    
    % 2. 合法用户速率 (看看它到底下降了多少)
    d_Alice_LR = norm(par.pos_Alice_True - par.pos_LR); % 用 True 比较好
    I_Alice = par.P_Alice * (par.beta0 / d_Alice_LR^par.alpha0);
    
    d_UAV_LR = norm(Q_no_lu(:,n) - par.pos_LR);
    I_UAV = P_no_lu(n) * (par.beta0 / d_UAV_LR^par.alpha0);
    
    SINR = Signal_LR / (I_Alice + I_UAV + par.sigma2_noise);
    R_LR_no_lu_vec(n) = log2(1 + SINR);
end

%% === 绘图对比 ===
time_axis = (0:par.N-1) * par.delta_t;

% 图 1: 中断概率对比
figure('Name', 'Impact of LU Constraint on Outage');
plot(time_axis, P_out_vec, 'r-o', 'LineWidth', 1.5, 'DisplayName', 'Proposed (With LU Constraint)'); hold on;
plot(time_axis, P_out_no_lu_vec, 'k--x', 'LineWidth', 1.5, 'DisplayName', 'Benchmark (No LU Constraint)');
xlabel('Time (s)'); ylabel('Outage Probability');
legend; title('中断概率对比：是否考虑合法用户');
grid on;

% 图 2: 合法用户速率对比 (关键验证)
figure('Name', 'Impact on Legitimate User Rate');
plot(time_axis, R_LR_vec, 'g-s', 'LineWidth', 1.5, 'DisplayName', 'Proposed R_{LR}'); hold on;
plot(time_axis, R_LR_no_lu_vec, 'm--^', 'LineWidth', 1.5, 'DisplayName', 'No Constraint R_{LR}');
yline(par.R_min_LR, 'r--', 'LineWidth', 2, 'DisplayName', 'Min Rate R_{min}');
xlabel('Time (s)'); ylabel('Rate (bit/s/Hz)');
legend; title('合法用户速率对比');
grid on;

% 图 3: 轨迹对比
figure('Name', 'Trajectory Comparison');
plot(par.pos_LR(1), par.pos_LR(2), 'mo', 'MarkerSize', 10, 'MarkerFaceColor', 'm', 'DisplayName', 'Legitimate Receiver'); hold on;
plot(Q_curr(1,:), Q_curr(2,:), 'r-', 'LineWidth', 1.5, 'DisplayName', 'Proposed Traj');
plot(Q_no_lu(1,:), Q_no_lu(2,:), 'k--', 'LineWidth', 1.5, 'DisplayName', 'No LU Constraint Traj');
legend; title('飞行轨迹对比 (2D投影)');
xlabel('X (m)'); ylabel('Y (m)');
grid on; axis equal;