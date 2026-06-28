function par = config_params()
    % CONFIG_PARAMS 设置仿真场景和算法参数
    
    %% === 物理场景参数 ===
    % 坐标设置 (x, y, z) [米]
    par.pos_Alice_True = [500, 300, 0]';
    par.pos_Bob_True   = [200, 300, 0]';
    par.pos_LT    = [500, 100, 0]';
    par.pos_LR    = [300, 100, 0]';
    par.pos_UAV_start = [400, 0, 50]';
    
    % 无人机运动参数
    par.T = 15;             % 飞行周期 (s) [cite: 16]
    par.H_min = 20;         % 最低飞行高度 (m) [cite: 16]
    par.V_max = 30;         % 最大速度 (m/s) [cite: 16]
    par.delta_t = 0.5;      % 时隙长度 (s)
    par.N = par.T / par.delta_t; % 时隙数量
    
    %% === 通信与信道参数 ===
    % 功率转线性值 function
    dBm2Watt = @(x) 10^((x-30)/10);
    
    par.sigma2_noise = dBm2Watt(-50);  % 噪声功率 [cite: 16]
    par.P_Alice      = dBm2Watt(7);   % Alice 发射功率 [cite: 16]
    par.P_LT         = dBm2Watt(10);   % LT 发射功率 [cite: 16]
    par.P_UAV_max    = dBm2Watt(30);   % 无人机最大干扰功率 [cite: 16]
    par.lambda_uu    = dBm2Watt(-20);  % 自干扰系数 [cite: 16]
    
    par.alpha0 = 2.2;       % 路径损耗指数 [cite: 16]
    par.beta0  = 0.633;     % 参考信道增益 [cite: 16]
    par.rho    = 1.0;       % 自干扰消除系数 
    
    %% === 约束阈值 ===
    par.epsilon_det = 0.10; % 检测错误概率约束 [cite: 16]
    par.p_prior     = 0.9;  % Alice 发送概率 [cite: 16]
    par.R_min_LR    = 1;    % LR 最小速率 (bit/s/Hz) [cite: 16]
    par.r_target    = 1;    % Bob 隐蔽目标速率 [cite: 16]
    
    %% === 算法控制参数 ===
    par.MaxIter = 30;       % 最大迭代次数
    par.tol_var = 1e-4;     % 变量收敛阈值 (Norm difference)
    
    %% === 误差参数 ===
    par.epsilon_pos_var = 5.0;
end