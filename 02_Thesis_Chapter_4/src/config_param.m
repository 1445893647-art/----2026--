% config_param.m
% 系统参数配置

global P
P.seed = 33455;             %随机种子13455、23455、

%% 1. 基础物理常数
P.c = 3e8;                  % 光速 (m/s)
P.fc = 28e9;                % 载波频率 (Hz) (假设毫米波)
P.lambda = P.c / P.fc;      % 波长 (m)
P.d = P.lambda / 2;         % 天线阵元间距 [cite: 60]

%% 2. 天线阵列参数 (Willie) [cite: 60]
% UPA (Uniform Planar Array)
P.N_tx = 4; P.N_ty = 4;     % 发射天线行列数
P.N_t = P.N_tx * P.N_ty;    % 发射天线总数
P.N_rx = 4; P.N_ry = 4;     % 接收天线行列数
P.N_r = P.N_rx * P.N_ry;    % 接收天线总数

%% 3. 时间与仿真参数 [cite: 8]
P.T = 10;                   % 总飞行时间 (s)
P.N_slots = 20;            % 总时隙数
P.dt = P.T / P.N_slots;     % 时隙间隔 Delta T

%% 4. 功率与信道参数
P.P_J = 1;                  % Willie 发射功率 (W) (本阶段固定)
P.sigma_r_sq = 1;       % 接收噪声功率 sigma_r^2 [cite: 76]
P.beta_0 = 1e-3;            % 参考距离处的路径损耗 [cite: 16]
P.zeta_0 = 50;            % 雷达反射系数 (简化) [cite: 76]
P.G = 10;                  % 匹配滤波增益 [cite: 85]

%% 5. EKF 噪声协方差参数
% 过程噪声 (Process Noise) Qx [cite: 58]
P.sigma_qx = 0.1; P.sigma_qy = 0.1; P.sigma_qz = 0.1;
P.sigma_vx = 0.01; P.sigma_vy = 0.01; P.sigma_vz = 0.01;

% 初始偏差要小：
P.x_b_0_error = [1.8; 1.8; 1; 0.5; 0.5; 0.1]; % 初始位置差1米，速度差0.5米/秒


% 测量噪声系数 (Measurement Noise Coeffs) [cite: 85]
P.C1 = 1; P.C2 = 6.7*1e-7; P.C3 = 2*1e4; 

%% 6. 初始状态设置 (用于生成轨迹)
% Bob 初始状态 [x, y, z, vx, vy, vz]
P.x_b_0 = [15, 15, 110, 12.4, 12.4, 0]'; 

% Willie 初始状态 (假设直线运动)
%P.q_w_0 = [50, 50, 120]';
P.q_w_0 =[0; 0; 80];

P.v_w = [5, 5, 0]'; % Willie 的恒定速度

%% 7. 隐蔽通信与优化参数
% Alice 参数
P.N_a = 32;                  % Alice 天线数
P.q_a = [80, 30, 0]';         % Alice 位置 
P.P_a_dB = 0; 
P.P_a = 10^(P.P_a_dB/10);                  % Alice 发射功率 (W)
P.P_max = 1;                % Willie 最大干扰功率 (W)

% 阈值参数
P.r_th = 1;                 % 目标隐蔽速率 (bits/s/Hz)
P.xi_req = 0.10;             % 检测错误概率上限 delta
P.eta_th = 1;            % 最大容忍跟踪 MSE (PCRB 阈值)     % 并不是硬约束，而是目标函数的方向

% 预计算 Gamma_req (用于检测约束 xi^u <= delta)
% 由 Bhattacharyya 系数单调性可知，B(gamma)=2*delta (pi0=pi1=0.5) 有唯一解。
% 这里直接数值反解得到 Gamma_req。
P.Gamma_req = get_gamma_req(P);

% 路径损耗指数
P.alpha_ab = 2.2;           % Alice-Bob
P.alpha_aw = 2.2;           % Alice-Willie
%% 8. 信道衰落参数 (莱斯信道)
P.K_dB = 10;                % 莱斯因子 (dB)
P.K = 10^(P.K_dB/10);       % 莱斯因子 (线性值)

P.V_max = 20;       % 最大速度 (m/s)
P.H_min = 50;       % 最小高度 (m)
P.alpha_aw = 2.2;   % 路径损耗指数

% Bob 接收机噪声标准差 (用于计算 P_out)
P.sigma_b = sqrt(1e-13);  % 假设噪声功率为 -100dBm 左右

% Willie 接收机噪声标准差 (用于计算检测错误概率 xi)
P.sigma_w = sqrt(1e-13);

%自干扰系数
P.rho = 2.6e-7;
P.rho_dB = -66;                % 莱斯因子 (dB)
P.rho = 10^(P.rho_dB/10);       % 莱斯因子 (线性值)

rho_x = 0.7;     % x 方向相关系数（0~1）
rho_y = 0.7;     % y 方向相关系数（0~1）

Rx = toeplitz(rho_x.^(0:P.N_tx-1));
Ry = toeplitz(rho_y.^(0:P.N_ty-1));

P.R_wb = kron(Rx, Ry);      % (Nt x Nt), Hermitian PSD

% 可选：按 trace 归一化（让 tr(R)=Nt）
P.R_wb = P.R_wb * (P.N_t / trace(P.R_wb));

% 可选：数值对称化
%P.R_wb = (P.R_wb + P.R_wb')/2;

P.P_cov = diag([2, 2, 2, 0.5, 0.5, 0.1].^2);
% ===== 保存数据设置 =====
P.save_data   = true;          % 是否保存数据
P.save_csv    = true;          % 是否导出 csv
P.output_root = 'sim_outputs'; % 输出根目录
P.beamforming_fix = 0.17;
P.penalty_pcrb = 1e3;   % 和轨迹子问题同量级即可