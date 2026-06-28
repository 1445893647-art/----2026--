% main_ekf_sim.m
clc; clear; close all;
config_param; % 加载参数
seed = P.seed;          % 你随便选一个整数
rng(seed, 'twister');     % 固定 rand / randn / randi
%% 1. 轨迹生成 (Ground Truth)
% 状态转移矩阵 F (匀速模型) [cite: 94]
I3 = eye(3); O3 = zeros(3);
F_mat = [I3, P.dt*I3; O3, I3]; 

% 过程噪声协方差 Qx [cite: 58]
Qx = diag([P.sigma_qx^2, P.sigma_qy^2, P.sigma_qz^2, ...
           P.sigma_vx^2, P.sigma_vy^2, P.sigma_vz^2]);

% 初始化存储
x_true = zeros(6, P.N_slots);
x_true(:,1) = P.x_b_0;

% Willie 的轨迹存储 (作为优化变量，初始设为全 0，随仿真填充)
q_w_hist = zeros(3, P.N_slots);
q_w_hist(:,1) = P.q_w_0;

xi_th = 0.05:0.05:0.30;
average_P_out = zeros(1,length(xi_th));
% 初始化结果存储
avg_pout_ekf    = zeros(1, length(xi_th));
avg_pout_oracle = zeros(1, length(xi_th));
avg_pout_lag1   = zeros(1, length(xi_th));

% 生成真实轨迹
for n = 2:P.N_slots
    % Bob 运动 (直线 + 轻微扰动)
    x_true(:,n) = F_mat * x_true(:,n-1) + sqrt(Qx) * randn(6,1);
end

for i = 1:length(xi_th)
    rng(seed, 'twister');
    P.xi_req = xi_th(i);
    P.Gamma_req = get_gamma_req(P);
    sum_P_out = 0;
    %% 2. EKF 初始化
    x_est = zeros(6, P.N_slots);
    x_est(:,1) = x_true(:,1) + P.x_b_0_error; % 初始估计有一定的误差
    P_cov = P.P_cov;  % 初始误差协方差矩阵

    %% 3. 性能记录初始化
    history_P_out = zeros(1, P.N_slots);
    history_xi = zeros(1, P.N_slots);
    history_dist_wb = zeros(1, P.N_slots);
    history_P_J = zeros(1, P.N_slots);
    history_xi_exact = zeros(1, P.N_slots); % 新增存储真实 xi

    %% 3. EKF 循环 
    fprintf('开始 EKF 跟踪与联合优化仿真...\n');

    % 第1个时隙的初始值
    q_w_curr = P.q_w_0; 
    % 初始速度 (假设初始静止或给定初速度)
    v_w_curr = [0; 0; 0];

    P_J_curr = P.P_J;   
    w_wb_curr = get_steering_vec(0, 0, P.N_tx, P.N_ty); 
    w_wb_curr = w_wb_curr / norm(w_wb_curr);

    % Alice 波束 (简化假设：Alice 始终尝试对准 Bob 的预测/真实位置)
    % 在实际对抗中，Alice 也是动态的，这里假设 Alice 使用理想波束作为基准
    w_ab_dummy = ones(P.N_a, 1);


    for n = 2:P.N_slots    
    % =======================================================
        % 记录上一时隙的最终位置 (用于计算位移和速度)
        q_w_prev_slot = q_w_curr;
        v_w_prev_slot = v_w_curr;

        % Step 1: EKF 状态预测 (Predict)
        % =======================================================
        x_pred = F_mat * x_est(:,n-1);
        P_pred = F_mat * P_cov * F_mat' + Qx;

        % =======================================================
        % Step 2: BCD 交替优化 (Joint Optimization)
        % =======================================================
        % 目标：基于预测状态 x_pred，找到当前时隙最优的 {P_J, w_wb, q_w}

        % BCD 迭代参数
        max_bcd_iter = 10;       % 最大迭代次数 (实时性要求，通常 3-5 次即可)
        tol_bcd = 1e-3;         % 收敛容差

        % 初始化本时隙的优化变量 (利用上一时隙解作为初始点 - Hot Start)
        % q_w_opt 初始值设为上一位置 (Willie 还没动)
        q_w_opt = q_w_curr; 
        v_w_opt = v_w_curr; % 初始速度猜测
        P_J_opt = P_J_curr;
        w_wb_opt = w_wb_curr;
        for iter = 1:max_bcd_iter
            % 记录旧变量用于判断收敛
            q_w_old = q_w_opt;
            P_J_old = P_J_opt;
            w_wb_old = w_wb_opt;
            % --- 2.3 轨迹优化子问题 ---
            % 输入: 最新的 P_J_opt 和 w_wb_opt
            % 注意：这里输入的是 q_w_curr (上一时隙位置) 作为起始点，优化出 q_w_opt (本时隙位置)
            q_w_new = solve_trajectory_opt(x_pred, P_pred, q_w_prev_slot, v_w_opt, P_J_opt, w_wb_opt, P);
            q_w_opt = q_w_new;
            % --- 2.1 功率优化子问题 ---
            % 输入: 当前的轨迹 q_w_opt (此时 w_wb 不影响功率闭式解)
            P_J_opt = solve_power_opt(q_w_opt, P);

            % --- 2.2 波束赋形优化子问题 ---
            % 输入: 最新的 P_J_opt 和 q_w_opt
            try
                % w_prev 传入上一轮的 w_wb_opt 用于线性化 H 矩阵
                w_wb_new = solve_beamforming_opt(x_pred, P_pred, v_w_opt, q_w_opt, P_J_opt, w_wb_opt, P);
                w_wb_opt = w_wb_new;
            catch
                % 保持上一次的波束
                warning('CVX 未运行或求解失败，保持上一次的波束。');
                q_b_tmp = x_pred(1:3);
                [th_tmp, ph_tmp, ~] = get_angle_dist(q_w_opt, q_b_tmp);
                w_wb_opt = get_steering_vec(th_tmp, ph_tmp, P.N_tx, P.N_ty); 
            end

            

            % 更新速度
            % 速度 = (当前优化出的位置 - 上一时隙位置) / delta_T
            v_w_opt = (q_w_opt - q_w_prev_slot) / P.dt;

            % --- 2.4 收敛性检查 ---
            diff_q = norm(q_w_opt - q_w_old);
            diff_P = abs(P_J_opt - P_J_old);
            diff_w = norm(w_wb_opt - w_wb_old);            
            if diff_q < tol_bcd && diff_P < tol_bcd && diff_w < tol_bcd
                fprintf('Slot %d iter: %d, diff_q: %.6g, diff_P: %.6g, diff_w: %.6g\n',n, iter,diff_q,diff_P,diff_w);
                break; % 已收敛
            end
        end
        % =======================================================
        % Step 3: 执行动作与环境交互 (Execute & Measure)
        % =======================================================
        % 更新系统状态为 BCD 优化后的结果
        q_w_curr = q_w_opt;
        v_w_curr = v_w_opt; % 更新由于轨迹变化产生的真实速度
        P_J_curr = P_J_opt;
        w_wb_curr = w_wb_opt;    

        % 记录轨迹
        q_w_hist(:, n) = q_w_curr;
        history_P_J(n) = P_J_curr;

        % =======================================================

        % --- B. 模拟测量 (Measurement Generation) ---

        % 基于 Bob 真实位置计算无噪声测量值
        x_b_true = x_true(:,n);
        %接收 params_true 以计算真实信噪比
        [y_clean, params_true] = measurement_function(x_b_true, q_w_curr, v_w_curr, w_wb_curr, P);

        % 3.动态计算测量噪声 Qy (基于真实状态生成仿真噪声) 
        % 计算反射系数 beta_r = zeta_0 / d^2
        beta_r_true = P.zeta_0 / (params_true.dist^2);

        % 计算波束成形增益因子 delta_n = a_W' * w_wb
        delta_n_true = params_true.a_W' * w_wb_curr;

        % 计算公共项分母: G * |beta|^2 * |delta|^2 * P_J
        denom_true = P.G * abs(beta_r_true)^2 * abs(delta_n_true)^2 * P_J_curr;

        % 根据公式 85 计算方差
        % 生成测量噪声 Qy [cite: 85]
        % 注意：噪声方差与回波强度有关，这里简化处理，假设信噪比足够高
        % 实际应根据公式 85 动态计算 sigma
        var_1_true = (P.C1^2 * P.sigma_r_sq) / (P.G * P_J_curr); % sigma_1^2
        var_2_true = (P.C2^2 * P.sigma_r_sq) / denom_true;    % sigma_tau^2
        var_3_true = (P.C3^2 * P.sigma_r_sq) / denom_true;    % sigma_mu^2

        % 4. 添加噪声得到观测值 y [cite: 91]
        noise_vec = [sqrt(var_1_true/2)*(randn(P.N_r,1) + 1j*randn(P.N_r,1)); ...
                     sqrt(var_2_true)*randn; ...
                     sqrt(var_3_true)*randn];
        y_meas = y_clean + noise_vec;

        % --- C. 雅可比矩阵计算 (Linearization) [cite: 111] ---
        H_jac = calculate_jacobian(x_pred, q_w_curr, v_w_curr, w_wb_curr, P);

        % --- D. 卡尔曼增益与更新 [cite: 115-119] ---
        % 注意：由于测量值包含复数 (回波)，MATLAB 处理复数求逆通常没问题
        % 但为了严谨，通常将实部虚部展开。这里直接使用复数运算。
         % 滤波器不知道真实位置，只能用预测位置 x_pred 来估算当前的噪声水平
        [~, params_pred] = measurement_function(x_pred, q_w_curr, v_w_curr, w_wb_curr, P);

        beta_r_pred = P.zeta_0 / (params_pred.dist^2);
        delta_n_pred = params_pred.a_W' * w_wb_curr;
        denom_pred = P.G * abs(beta_r_pred)^2 * abs(delta_n_pred)^2 * P_J_curr;

        var_1_p = (P.C1^2 * P.sigma_r_sq) / (P.G * P_J_curr);
        var_2_p = (P.C2^2 * P.sigma_r_sq) / denom_pred;
        var_3_p = (P.C3^2 * P.sigma_r_sq) / denom_pred;

        % --- 复数测量处理：将回波向量拆成实部/虚部，构造纯实EKF ---
        Nr = P.N_r;

        % 预测测量
        [y_pred, ~] = measurement_function(x_pred, q_w_curr, v_w_curr, w_wb_curr, P);

        % 观测/预测：复 -> 实拼接 (2Nr+2 维)
        y_meas_r = [real(y_meas(1:Nr)); imag(y_meas(1:Nr)); real(y_meas(Nr+1:end))];
        y_pred_r = [real(y_pred(1:Nr)); imag(y_pred(1:Nr)); real(y_pred(Nr+1:end))];
        innovation_r = y_meas_r - y_pred_r;

        % 雅可比：复 -> 实拼接
        H_r = [real(H_jac(1:Nr,:)); imag(H_jac(1:Nr,:)); real(H_jac(Nr+1:end,:))];

        % 预测测量噪声协方差：回波实/虚各占一半方差
        Qy_r = blkdiag((var_1_p/2)*eye(Nr), (var_1_p/2)*eye(Nr), var_2_p, var_3_p);

        % 卡尔曼增益
        S = H_r * P_pred * H_r' + Qy_r;
        K = P_pred * H_r' / S;

        % 更新状态与协方差
        x_est(:,n) = x_pred + K * innovation_r;
        P_cov = (eye(6) - K * H_r) * P_pred;

        % 强制 P_cov 对称正定 (防止数值问题)
        P_cov = (P_cov + P_cov') / 2;

        % Step 5: 计算性能指标 (Metrics)
        % =======================================================
        % 使用更新后的状态和参数计算
        [p_out_val, xi_val] = calc_metrics(q_w_curr, x_b_true(1:3), P_J_curr, w_wb_curr, w_ab_dummy, P);

        history_P_out(n) = p_out_val;
        history_xi(n) = xi_val;
        history_dist_wb(n) = norm(q_w_curr - x_b_true(1:3));
        sum_P_out = sum_P_out+p_out_val;
        % 进度
        fprintf('xi_th: %.2f | Slot %d/%d | Dist_wb: %.1fm | P_J: %.2fW | P_out: %.2f\n', ...
                P.xi_req, n, P.N_slots, history_dist_wb(n), P_J_curr, p_out_val);
    end
    average_P_out(i) = sum_P_out/P.N_slots;
    avg_pout_ekf(i) = sum_P_out/P.N_slots;
    %res_oracle = run_noekf_opt(P, x_true, F_mat, Qx, "oracle");
    %avg_pout_oracle(i) = mean(res_oracle.P_out);
    %res_lag1 = run_noekf_opt(P, x_true, F_mat, Qx, "lag1");
    %avg_pout_lag1(i) = mean(res_lag1.P_out);
end



%% 4. 绘图结果
set(0, 'DefaultAxesFontSize', 12);
set(0, 'DefaultLineLineWidth', 1.5);
set(0, 'DefaultTextInterpreter', 'latex');
set(0, 'DefaultLegendInterpreter', 'latex');

figure('Position', [300, 300, 700, 500], 'Name', 'Average P_out vs Delta');

%plot(xi_th, avg_pout_oracle, 'b--s', 'MarkerSize', 8, 'LineWidth', 1.5, 'DisplayName', 'Oracle (Upper Bound)'); hold on;
plot(xi_th, avg_pout_ekf,    'r-o',  'MarkerSize', 8, 'LineWidth', 2,   'DisplayName', 'Proposed EKF');
%plot(xi_th, avg_pout_lag1,   'm:^',  'MarkerSize', 8, 'LineWidth', 1.5, 'DisplayName', 'Lag-1 (Benchmark)');

grid on;
xlabel('Detection Error Threshold $\delta$', 'Interpreter', 'latex');
ylabel('Average Outage Probability $\bar{P}_{\text{out}}$', 'Interpreter', 'latex');
title('\textbf{Covertness vs. Detection Requirement}', 'Interpreter', 'latex');
legend('Location', 'best', 'Interpreter', 'latex');
xlim([min(xi_th), max(xi_th)]);
ylim([0, 1.0]);

% 保存数据
% 定义保存路径
save_dir = 'Sweep_Data'; % 文件夹名称
if ~exist(save_dir, 'dir')
    mkdir(save_dir); % 如果文件夹不存在，自动创建
end

% 生成带时间戳的文件名 (格式: sweep_results_xi_年月日_时分秒.csv)
timestamp = datestr(now, 'yyyymmdd_HHMMSS');
filename = sprintf('sweep_results_xi_%s.csv', timestamp);
full_path = fullfile(save_dir, filename);

% 构建表格并保存
results_sweep = table(xi_th', avg_pout_ekf', avg_pout_oracle', avg_pout_lag1', ...
    'VariableNames', {'Delta', 'Avg_Pout_EKF', 'Avg_Pout_Oracle', 'Avg_Pout_Lag1'});

writetable(results_sweep, full_path);
fprintf('结果文件已保存至: %s\n', full_path);
