% =========================================================================
%  Joint Trajectory and Power Optimization (Variable Convergence Version)
% =========================================================================

clear; clc; close all;

%% 1. 加载参数
par = config_params();

%% 1.5 生成真实位置 (True Positions with Uncertainty)
rng('shuffle'); % 确保随机性
% 生成高斯噪声 N(0, epsilon)
% 仅在 X, Y 平面添加误差 (Z轴通常为0)
std_dev = sqrt(par.epsilon_pos_var);

noise_Alice = std_dev * randn(3, 1); noise_Alice(3) = 0;
noise_Bob   = std_dev * randn(3, 1); noise_Bob(3)   = 0;

% 真实位置 = 估计位置 + 噪声
par.pos_Alice_Est = par.pos_Alice_True + noise_Alice;
par.pos_Bob_Est   = par.pos_Bob_True   + noise_Bob;

fprintf('位置误差生成完毕 (方差 %.1f):\n', par.epsilon_pos_var);
fprintf('Alice 偏移: [%.2f, %.2f]\n', noise_Alice(1), noise_Alice(2));
fprintf('Bob   偏移: [%.2f, %.2f]\n', noise_Bob(1), noise_Bob(2));

%% 2. 初始化变量
% 轨迹初始化：生成一个稍微有些偏移的直线轨迹，避免 SCA 初始点问题
Q_curr = zeros(3, par.N);
for n = 1:par.N
    % 简单直线飞向 Bob 方向，但保持一定高度
    direction = (par.pos_Bob_Est - par.pos_UAV_start);
    direction(3) = 0; 
    step_vec = (direction / norm(direction)) * (par.V_max * 0.5 * par.delta_t);
    
    Q_curr(:, n) = par.pos_UAV_start + step_vec * (n-1);
    Q_curr(3, n) = max(Q_curr(3, n), par.H_min);
end

% 功率初始化
P_curr = par.P_UAV_max * ones(1, par.N);

% 记录器
history.obj = [];
history.diff_Q = [];
history.diff_P = [];

fprintf('开始优化...\n');
fprintf('%-5s | %-12s | %-12s | %-12s\n', 'Iter', 'Outage Prob', 'Diff Traj', 'Diff Power');
fprintf('----------------------------------------------------------\n');

%% 3. AO 主循环
for iter = 1:par.MaxIter
    
    Q_prev = Q_curr;
    P_prev = P_curr;
    
    % --- Step 1: 优化功率 ---
    P_curr = solve_power_subproblem(Q_curr, par);
    
    % --- Step 2: 优化轨迹 ---
    [Q_curr, status] = solve_traj_subproblem(P_curr, Q_prev, par);
    
    if ~strcmp(status, 'Solved') && ~strcmp(status, 'Inaccurate/Solved')
        warning('轨迹子问题未最优求解: %s', status);
    end
    
    % --- Step 3: 计算指标与收敛性 ---
    % 计算当前目标函数值
    obj_val = calc_outage_prob(Q_curr, P_curr, par);
    history.obj = [history.obj, obj_val];
    
    % 计算变量变化量 (归一化差异)
    diff_Q = norm(Q_curr - Q_prev, 'fro') / norm(Q_prev, 'fro');
    diff_P = norm(P_curr - P_prev) / norm(P_prev);
    
    history.diff_Q = [history.diff_P, diff_Q];
    history.diff_P = [history.diff_P, diff_P];
    
    fprintf('%-5d | %-12.4f | %-12.4e | %-12.4e\n', iter, obj_val, diff_Q, diff_P);
    
    % --- Step 4: 判断收敛 (基于变量) ---
    if iter > 1 && diff_Q < par.tol_var && diff_P < par.tol_var
        fprintf('算法在第 %d 次迭代收敛 (变量准则)。\n', iter);
        break;
    end
end

%% 4. 绘图结果
figure('Position', [100, 100, 1000, 400]);

% 子图1: 3D 轨迹
subplot(1, 2, 1);
plot3(par.pos_Alice_True(1), par.pos_Alice_True(2), par.pos_Alice_True(3), 'rx', 'LineWidth', 2, 'MarkerSize', 10); hold on;
plot3(par.pos_Bob_True(1), par.pos_Bob_True(2), par.pos_Bob_True(3), 'bx', 'LineWidth', 2, 'MarkerSize', 10);
plot3(par.pos_LT(1), par.pos_LT(2), par.pos_LT(3), 'gx', 'LineWidth', 2, 'MarkerSize', 10);
plot3(par.pos_LR(1), par.pos_LR(2), par.pos_LR(3), 'mx', 'LineWidth', 2, 'MarkerSize', 10);
plot3(Q_curr(1,:), Q_curr(2,:), Q_curr(3,:), 'k-o', 'LineWidth', 1.5, 'MarkerSize', 4);
grid on;
legend('Alice', 'Bob', 'LT', 'LR', 'UAV Trajectory');
title('无人机 3D 飞行轨迹');
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
view(2); % 俯视图查看更清晰

% 子图2: 功率分配
subplot(1, 2, 2);
t_axis = (0:par.N-1) * par.delta_t;
plot(t_axis, P_curr, 'r-o', 'LineWidth', 1.5);
grid on;
title('无人机干扰功率分配');
xlabel('时间 (s)'); ylabel('功率 (W)');

% 额外图: 收敛曲线
figure;
yyaxis left
plot(history.obj, '-s', 'LineWidth', 2);
ylabel('平均中断概率');
yyaxis right
plot(history.diff_Q, '--', 'LineWidth', 1.5); hold on;
% plot(history.diff_P, ':', 'LineWidth', 1.5);
set(gca, 'YScale', 'log');
ylabel('变量相对变化量 (Log Scale)');
legend('Objective', 'Traj Diff');
title('算法收敛性能');