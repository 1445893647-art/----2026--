% =========================================================================
%  Analysis: Average Outage Probability vs Detection Constraint epsilon
%  (Runs the full Joint Optimization for each epsilon value)
% =========================================================================

clear; clc; close all;

%% 1. 准备工作
% 加载默认参数
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

% 定义 epsilon 的扫描范围 (0.05 到 0.30，步长 0.05)
epsilon_vec = 0.05 : 0.05 : 0.30;
num_points = length(epsilon_vec);

% 预分配结果数组
P_out_results = zeros(1, num_points);

% 算法设置
par.MaxIter = 15;   % 每个 epsilon 点的内部迭代次数 (通常10-15次即可收敛)
par.tol_var = 1e-4; % 收敛阈值

fprintf('开始仿真 P_out 随 epsilon 的变化...\n');
fprintf('%-10s | %-15s | %-10s\n', 'Epsilon', 'Avg P_out', 'Status');
fprintf('--------------------------------------------\n');

%% 2. 循环仿真 (Loop over epsilon)
for k = 1:num_points
    
    % --- 更新当前 epsilon ---
    current_eps = epsilon_vec(k);
    par.epsilon_det = current_eps;
    
    % --- 必须检查 epsilon 是否合法 ---
    % 朗伯 W 函数要求 (p - epsilon)/p > 0，即 epsilon < p
    if current_eps >= par.p_prior
        warning('Epsilon %.2f >= Prior p %.2f, 跳过此点', current_eps, par.p_prior);
        P_out_results(k) = NaN;
        continue;
    end
    
    % --- 初始化 (每次都重置，保证独立性) ---
    % 轨迹初始化：直线
    Q_curr = zeros(3, par.N);
    for n = 1:par.N
        direction = (par.pos_Bob_Est - par.pos_UAV_start);
        direction(3) = 0; 
        step_vec = (direction / norm(direction)) * (par.V_max * 0.5 * par.delta_t);
        Q_curr(:, n) = par.pos_UAV_start + step_vec * (n-1);
        Q_curr(3, n) = max(Q_curr(3, n), par.H_min);
    end
    
    % 功率初始化：满功率
    P_curr = par.P_UAV_max * ones(1, par.N);
    
    % --- AO 联合优化算法 ---
    for iter = 1:par.MaxIter
        Q_prev = Q_curr;
        P_prev = P_curr;
        
        % 1. 优化功率
        P_curr = solve_power_subproblem(Q_curr, par);
        
        % 2. 优化轨迹
        [Q_curr, status] = solve_traj_subproblem(P_curr, Q_prev, par);
        
        % 3. 简单收敛检查 (基于变量)
        diff_Q = norm(Q_curr - Q_prev, 'fro') / norm(Q_prev, 'fro');
        diff_P = norm(P_curr - P_prev) / norm(P_prev);
        
        if diff_Q < par.tol_var && diff_P < par.tol_var
            break; 
        end
    end
    
    % --- 记录最终结果 ---
    final_P_out = calc_outage_prob(Q_curr, P_curr, par);
    P_out_results(k) = final_P_out;
    
    fprintf('%-10.2f | %-15.4f | Converged in %d iters\n', current_eps, final_P_out, iter);
end

%% 3. 绘图 (Plotting)
figure('Name', 'P_out vs Epsilon');

% 绘制曲线
plot(epsilon_vec, P_out_results, 'r-s', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'r');

grid on;
xlabel('检测错误概率约束 \epsilon (Detection Error Probability Constraint)');
ylabel('平均隐蔽通信中断概率 \bar{P}_{out} (Average Outage Probability)');
title('中断概率随检测性能约束的变化');

% 设置坐标轴范围美化
xlim([0.04, 0.31]);
ylim([min(P_out_results)*0.95, 1.0]); % 动态调整Y轴，上限为1

% 添加图例
legend('Proposed Joint Optimization', 'Location', 'SouthEast');

% 标注数据点 (可选)
for i = 1:length(epsilon_vec)
    text(epsilon_vec(i), P_out_results(i), sprintf('  %.3f', P_out_results(i)), ...
         'VerticalAlignment', 'bottom', 'FontSize', 9);
end