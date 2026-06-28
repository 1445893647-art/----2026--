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

% 生成真实轨迹
for n = 2:P.N_slots
    % Bob 运动 (直线 + 轻微扰动)
    x_true(:,n) = F_mat * x_true(:,n-1) + sqrt(Qx) * randn(6,1);
end

%% 2. EKF 初始化
x_est = zeros(6, P.N_slots);
x_est(:,1) = x_true(:,1) + P.x_b_0_error; % 初始估计有一定的误差
P_cov = P.P_cov; % 初始误差协方差矩阵

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

% BCD 迭代参数
max_bcd_iter = 10;       % 最大迭代次数 (实时性要求，通常 3-5 次即可)
tol_bcd = 1e-6;         % 收敛容差
pout_bcd_hist = zeros(P.N_slots, max_bcd_iter); % 每个时隙、每次迭代的 P_out（用于画平均收敛曲线）
    
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
        catch ME
            % 保持上一次的波束
            warning('Beamforming failed: %s', ME.message);
            disp(getReport(ME,'extended'));   % <-- 这行会把真实栈信息打印出来
            q_b_tmp = x_pred(1:3);
            [th_tmp, ph_tmp, ~] = get_angle_dist(q_w_opt, q_b_tmp);
            w_wb_opt = get_steering_vec(th_tmp, ph_tmp, P.N_tx, P.N_ty); 
        end
        
       
        
        % 更新速度
        % 速度 = (当前优化出的位置 - 上一时隙位置) / delta_T
        v_w_opt = (q_w_opt - q_w_prev_slot) / P.dt;
        
        [pout_tmp, ~] = calc_metrics(q_w_opt, x_true(1:3,n), P_J_opt, w_wb_opt, w_ab_dummy, P);
        pout_bcd_hist(n, iter) = pout_tmp;
        
        % --- 2.4 收敛性检查 ---
        diff_q = norm(q_w_opt - q_w_old);
        diff_P = abs(P_J_opt - P_J_old);
        diff_w = norm(w_wb_opt - w_wb_old);
        
        if diff_q < tol_bcd && diff_P < tol_bcd && diff_w < tol_bcd
            fprintf('Slot %d iter: %d, diff_q: %.6g, diff_P: %.6g, diff_w: %.6g\n',n, iter,diff_q,diff_P,diff_w);
            
            break; % 已收敛
        end
    end
    if iter < max_bcd_iter && ~isnan(pout_bcd_hist(n, iter))
        pout_bcd_hist(n, iter+1:end) = pout_bcd_hist(n, iter);
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
    
    % 进度
    fprintf('Slot %d/%d | Dist_wb: %.1fm | P_J: %.2fW | P_out: %.2f\n', ...
            n, P.N_slots, history_dist_wb(n), P_J_curr, p_out_val);
end


%% 3.5 交替优化收敛曲线：平均中断概率 vs 迭代次数
% 说明：每个时隙内部有一段 BCD 迭代；这里把“第 i 次迭代后的 P_out”对所有时隙取平均，观察收敛趋势。
avg_pout_bcd = mean(pout_bcd_hist(1:end,:), 1, 'omitnan');
figure('Name','BCD Convergence: Avg Pout vs Iter','NumberTitle','off');
plot(1:max_bcd_iter, avg_pout_bcd, 'o-', 'LineWidth', 2);
ylim([0, 1.05]);
grid on; xlabel('迭代次数'); ylabel('平均中断概率');
title('交替优化收敛曲线：平均中断概率 vs 迭代次数');

%% 4. 绘图结果（完全分开：每个指标一张图，轨迹除外）
%% EKF 跟踪对比图：角度/距离/速度 + x/y/z

t = (0:P.N_slots-1) * P.dt;   % 时间轴 (s)

theta_true = zeros(1, P.N_slots);
theta_est  = zeros(1, P.N_slots);
phi_true = zeros(1, P.N_slots);
phi_est  = zeros(1, P.N_slots);
dist_true  = zeros(1, P.N_slots);
dist_est   = zeros(1, P.N_slots);
speed_true = zeros(1, P.N_slots);
speed_est  = zeros(1, P.N_slots);

for k = 1:P.N_slots
    q_wk = q_w_hist(:,k);

    % 真实角度/距离（用你已有的 get_angle_dist 定义：theta=俯仰, phi=方位）
    [th_t, ph_t, d_t] = get_angle_dist(q_wk, x_true(1:3,k));
    theta_true(k) = th_t;
    phi_true(k)=ph_t;
    dist_true(k)  = d_t;

    % 估计角度/距离
    [th_e, ph_e, d_e] = get_angle_dist(q_wk, x_est(1:3,k));
    theta_est(k) = th_e;
    phi_est(k) = ph_e;
    dist_est(k)  = d_e;

    % 速度（这里用速度向量模长；如果你想画“径向速度”，我也可以给对应代码）
    speed_true(k) = norm(x_true(4:6,k));
    speed_est(k)  = norm(x_est(4:6,k));
end

% ===== Figure A: 角度/距离/速度（三子图，真实 vs 估计）=====
figure('Name','EKF Tracking: Angle/Distance/Speed','NumberTitle','off');

subplot(4,1,1);
plot(t, theta_true, 'k-', 'LineWidth', 1.5); hold on;
plot(t, theta_est,  'k--','LineWidth', 1.5);
grid on; ylabel('仰角 (rad)');
legend('实际仰角','估计仰角','Location','best');

subplot(4,1,2);
plot(t, phi_true, 'k-', 'LineWidth', 1.5); hold on;
plot(t, phi_est,  'k--','LineWidth', 1.5);
grid on; ylabel('俯角 (rad)');
legend('实际俯角','估计俯角','Location','best');

subplot(4,1,3);
plot(t, dist_true, 'k-', 'LineWidth', 1.5); hold on;
plot(t, dist_est,  'k--','LineWidth', 1.5);
grid on; ylabel('距离 (m)');
legend('实际距离','估计距离','Location','best');

subplot(4,1,4);
plot(t, speed_true, 'k-', 'LineWidth', 1.5); hold on;
plot(t, speed_est,  'k--','LineWidth', 1.5);
grid on; ylabel('速度 (m/s)'); xlabel('时间 (s)');
legend('实际速度','估计速度','Location','best');


% ===== Figure B: x/y/z（三子图，真实 vs 估计）=====
figure('Name','EKF Tracking: X/Y/Z','NumberTitle','off');

subplot(3,1,1);
plot(t, x_true(1,:), 'k-', 'LineWidth', 1.5); hold on;
plot(t, x_est(1,:),  'k--','LineWidth', 1.5);
grid on; ylabel('X (m)');
legend('实际 X','估计 X','Location','best');

subplot(3,1,2);
plot(t, x_true(2,:), 'k-', 'LineWidth', 1.5); hold on;
plot(t, x_est(2,:),  'k--','LineWidth', 1.5);
grid on; ylabel('Y (m)');
legend('实际 Y','估计 Y','Location','best');

subplot(3,1,3);
plot(t, x_true(3,:), 'k-', 'LineWidth', 1.5); hold on;
plot(t, x_est(3,:),  'k--','LineWidth', 1.5);
grid on; ylabel('Z (m)'); xlabel('时间 (s)');
legend('实际 Z','估计 Z','Location','best');

% Figure 1: 3D 轨迹（保留）
figure('Name','Trajectory 3D','NumberTitle','off','Position',[100, 100, 850, 650]);
plot3(x_true(1,:), x_true(2,:), x_true(3,:), 'b-', 'LineWidth', 1.5); hold on;
plot3(q_w_hist(1,:), q_w_hist(2,:), q_w_hist(3,:), 'r-o', 'LineWidth', 1.5, ...
      'MarkerIndices', 1:5:P.N_slots);
plot3(P.q_a(1), P.q_a(2), P.q_a(3), 'kp', 'MarkerSize', 12, 'MarkerFaceColor', 'g');
legend('Bob 真实轨迹', 'Willie 优化轨迹', 'Alice 位置', 'Location', 'best');
grid on; xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
title('无人机三维动态对抗轨迹');
axis equal; view(3);

% Figure 2: P_out 单独一张
figure('Name','Pout','NumberTitle','off');
plot(1:P.N_slots, history_P_out, 'r-', 'LineWidth', 2);
ylim([0, 1.1]);
grid on; xlabel('时隙 n'); ylabel('P_{out}');
title('隐蔽中断概率 P_{out}');

% Figure 3: xi 单独一张（带阈值线）
figure('Name', 'Detection Error Prob Comparison', 'NumberTitle', 'off');
hold on; box on; grid on;
plot(1:P.N_slots, history_xi, 'b--', 'LineWidth', 1.5, 'DisplayName', '检测错误概率 (\xi^*)');

% 3. 绘制目标阈值线
yline(P.xi_req, 'k:', 'LineWidth', 1.5, 'DisplayName', '检测阈值 (\delta)');

xlabel('时隙 n');
ylabel('检测错误概率');
legend('Location', 'best');
title('检测性能');
ylim([0, 0.3]); % 错误概率最高 0.5 (瞎猜)

% Figure 4: P_J 单独一张
figure('Name','P_J','NumberTitle','off');
plot(1:P.N_slots, history_P_J, 'm-', 'LineWidth', 2);
grid on; xlabel('时隙 n'); ylabel('P_J (W)');
title('干扰功率 P_J');

% Figure 5: Willie-Bob 距离 单独一张
figure('Name','Distance W-B','NumberTitle','off');
plot(1:P.N_slots, history_dist_wb, 'k-.', 'LineWidth', 1.5);
grid on; xlabel('时隙 n'); ylabel('Distance (m)');
title('Willie-Bob 相对距离');


%% 保存数据（仿真结束后）
if isfield(P,'save_data') && P.save_data
    run_id = datestr(now,'yyyymmdd_HHMMSS');
    outdir = fullfile(P.output_root, ['run_' run_id]);
    if ~exist(outdir,'dir'), mkdir(outdir); end

    % 时间轴
    t = (0:P.N_slots-1)' * P.dt;   % Nx1

    % 组织结果（mat 里建议把能用的都存进去）
    results = struct();
    results.P = P;
    results.run_id = run_id;
    results.t = t;

    % 状态/轨迹
    results.x_true = x_true;       % 6xN
    results.x_est  = x_est;        % 6xN
    results.q_w_hist = q_w_hist;   % 3xN

    % 关键指标
    results.history_P_out   = history_P_out;
    results.history_xi      = history_xi;
    results.history_P_J     = history_P_J;
    results.history_dist_wb = history_dist_wb;
    % ===== 对比方案（Oracle / Lag-1）======
    % 建议只存关键曲线，避免 oracle/lag1 里重复存一堆 P/x_true 导致 results.mat 变很大
    if exist('oracle','var') && isstruct(oracle)
        results.compare.oracle = struct();
        if isfield(oracle,'P_out'),   results.compare.oracle.P_out   = oracle.P_out(:); end
        if isfield(oracle,'xi'),      results.compare.oracle.xi      = oracle.xi(:); end
        if isfield(oracle,'P_J'),     results.compare.oracle.P_J     = oracle.P_J(:); end
        if isfield(oracle,'dist_wb'), results.compare.oracle.dist_wb = oracle.dist_wb(:); end
        if isfield(oracle,'q_w_hist'),results.compare.oracle.q_w_hist= oracle.q_w_hist; end
    end

    if exist('lag1','var') && isstruct(lag1)
        results.compare.lag1 = struct();
        if isfield(lag1,'P_out'),   results.compare.lag1.P_out   = lag1.P_out(:); end
        if isfield(lag1,'xi'),      results.compare.lag1.xi      = lag1.xi(:); end
        if isfield(lag1,'P_J'),     results.compare.lag1.P_J     = lag1.P_J(:); end
        if isfield(lag1,'dist_wb'), results.compare.lag1.dist_wb = lag1.dist_wb(:); end
        if isfield(lag1,'q_w_hist'),results.compare.lag1.q_w_hist= lag1.q_w_hist; end
    end
    
    % 若你已加 EKF 协方差记录（有就存，没有就跳过）
    if exist('history_traceP','var'),    results.history_traceP = history_traceP; end
    if exist('history_sigma_pos','var'), results.history_sigma_pos = history_sigma_pos; end
    if exist('history_sigma_vel','var'), results.history_sigma_vel = history_sigma_vel; end

    % 波束一般是复数/维度多：建议只存 mat
    if exist('w_wb_hist','var'), results.w_wb_hist = w_wb_hist; end

    save(fullfile(outdir,'results.mat'),'results','-v7.3');

    % ===== 导出 CSV（便于你直接发我）=====
    if isfield(P,'save_csv') && P.save_csv
        % 1) 指标汇总：n, t, P_out, xi, P_J, distWB
        metrics = [ (1:P.N_slots)', t, ...
                    history_P_out(:), history_xi(:), history_P_J(:), history_dist_wb(:) ];
        write_csv_safe(fullfile(outdir,'history_metrics.csv'), metrics);

        % 2) 真实/估计状态：n, t, x,y,z,vx,vy,vz
        st_true = [ (1:P.N_slots)', t, x_true.' ];  % Nx(2+6)
        st_est  = [ (1:P.N_slots)', t, x_est.'  ];
        write_csv_safe(fullfile(outdir,'state_true.csv'), st_true);
        write_csv_safe(fullfile(outdir,'state_est.csv'),  st_est);

        % 3) Willie 轨迹：n, t, qx,qy,qz
        traj_w = [ (1:P.N_slots)', t, q_w_hist.' ]; % Nx(2+3)
        write_csv_safe(fullfile(outdir,'traj_willie.csv'), traj_w);

        % 4) （可选）EKF 误差：n,t,ex,ey,ez,evx,evy,evz
        if ~isempty(x_est) && ~isempty(x_true)
            err = (x_est - x_true).';
            err_out = [ (1:P.N_slots)', t, err ];
            write_csv_safe(fullfile(outdir,'ekf_error.csv'), err_out);
        end
        % 5)（新增）对比方案：n, t, Pout_EKF, Pout_oracle, Pout_lag1, (可选) xi/PJ/dist
        if exist('oracle','var') && isstruct(oracle) && exist('lag1','var') && isstruct(lag1) ...
                && isfield(oracle,'P_out') && isfield(lag1,'P_out')

            Pout_oracle = oracle.P_out(:);
            Pout_lag1   = lag1.P_out(:);

            % 保证长度一致
            L = min([P.N_slots, numel(history_P_out), numel(oracle.P_out), numel(lag1.P_out)]);

            nvec = (1:L)';           % Lx1
            tvec = t(1:L);  tvec = tvec(:);

            p_ekf    = history_P_out(1:L); p_ekf    = p_ekf(:);
            p_oracle = oracle.P_out(1:L);  p_oracle = p_oracle(:);
            p_lag1   = lag1.P_out(1:L);    p_lag1   = p_lag1(:);

            comp = [nvec, tvec, p_ekf, p_oracle, p_lag1];

            % 可选追加 xi / PJ / dist（都统一成列）
            if isfield(oracle,'xi') && isfield(lag1,'xi') && numel(history_xi) >= L
                xi_ekf = history_xi(1:L); xi_ekf = xi_ekf(:);
                xi_o   = oracle.xi(1:L);  xi_o   = xi_o(:);
                xi_l   = lag1.xi(1:L);    xi_l   = xi_l(:);
                comp = [comp, xi_ekf, xi_o, xi_l];
            end

            if isfield(oracle,'P_J') && isfield(lag1,'P_J') && numel(history_P_J) >= L
                pj_ekf = history_P_J(1:L); pj_ekf = pj_ekf(:);
                pj_o   = oracle.P_J(1:L);  pj_o   = pj_o(:);
                pj_l   = lag1.P_J(1:L);    pj_l   = pj_l(:);
                comp = [comp, pj_ekf, pj_o, pj_l];
            end

            if isfield(oracle,'dist_wb') && isfield(lag1,'dist_wb') && numel(history_dist_wb) >= L
                d_ekf = history_dist_wb(1:L); d_ekf = d_ekf(:);
                d_o   = oracle.dist_wb(1:L);  d_o   = d_o(:);
                d_l   = lag1.dist_wb(1:L);    d_l   = d_l(:);
                comp = [comp, d_ekf, d_o, d_l];
            end

            write_csv_safe(fullfile(outdir,'history_compare.csv'), comp);
        end
    end

    fprintf('数据已保存到：%s\n', outdir);
end


%% 本地函数：兼容不同 MATLAB 版本写 CSV
function write_csv_safe(fname, M)
    if exist('writematrix','file') == 2
        writematrix(M, fname);
    else
        % 老版本没有 writematrix
        csvwrite(fname, M);
    end
end
