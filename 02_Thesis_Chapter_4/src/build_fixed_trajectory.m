function q_fix = build_fixed_trajectory(P, q_goal1, q_goal2)
%BUILD_FIXED_TRAJECTORY 生成两阶段固定轨迹 (自动时间分配版)
%
% 用法:
%   q_fix = build_fixed_trajectory(P, q_goal1, q_goal2)
%
% 逻辑:
%   1. 初始时刻 Willie 从 P.q_w_0 出发。
%   2. 全速飞向 q_goal1 (第一阶段)。
%   3. 一旦到达 q_goal1 (距离 < 阈值)，立即切换目标，全速飞向 q_goal2 (第二阶段)。
%   4. 如果到达 q_goal2，则在原地悬停。
%   5. 整个过程受 P.V_max (最大速度) 和 P.H_min (最小高度) 约束。

    % --- 1. 参数检查与初始化 ---
    if nargin < 2 || isempty(q_goal1)
        % 默认目标1: Bob 初始位置上方
        q_goal1 = [P.x_b_0(1); P.x_b_0(2); max(P.H_min, P.q_w_0(3))];
    end
    
    if nargin < 3
        q_goal2 = []; % 如果未提供第二个点，则仅执行单阶段任务
    end
    
    N = P.N_slots;
    q_fix = zeros(3, N);
    q_fix(:,1) = P.q_w_0(:); % 初始位置
    
    % 每个时隙的最大移动距离
    % 注意: 请确保 P 结构体中有 P.V_max (例如 10m/s 或 20m/s)
    if isfield(P, 'V_max')
        max_step = P.V_max * P.dt;
    else
        % 如果未定义 V_max，使用一个合理的默认值 (例如根据 P.v_w 推断)
        max_step = 10 * P.dt; 
    end

    % 状态标志: 是否已经到达过目标 1
    has_reached_goal1 = false;
    
    % --- 2. 轨迹生成循环 ---
    for n = 2:N
        q_prev = q_fix(:,n-1);
        
        % --- 状态机逻辑 ---
        if ~has_reached_goal1
            % [阶段 1] 尚未到达目标 1 -> 飞向 goal1
            current_target = q_goal1;
            
            % 检查距离
            d_vec = current_target - q_prev;
            dist = norm(d_vec);
            
            if dist <= 1e-3 % 认为已到达 (容差 1mm)
                has_reached_goal1 = true;
                % 到达的瞬间，本时隙是否还要动？
                % 如果已经重合，则直接看是否有 goal2
                if ~isempty(q_goal2)
                    current_target = q_goal2; % 立即切换
                else
                    current_target = q_prev;  % 没目标了，悬停
                end
            end
        else
            % [阶段 2] 已经到达过目标 1 -> 飞向 goal2
            if ~isempty(q_goal2)
                current_target = q_goal2;
            else
                current_target = q_prev; % 悬停
            end
        end
        
        % --- 运动逻辑 ---
        % 计算指向当前目标的向量
        d_vec = current_target - q_prev;
        dist = norm(d_vec);
        
        if dist < 1e-9
            % 已经完全重合 -> 悬停
            q_fix(:,n) = q_prev;
        else
            % 未重合 -> 移动
            step_len = min(max_step, dist); % 不能超速，也不能冲过头
            step_vec = step_len * (d_vec / dist);
            q_fix(:,n) = q_prev + step_vec;
        end
        
        % 再次检查：如果这一步走完刚好到达目标1，更新标志位
        % 这样下一个时隙就会自动切换到 goal2
        if ~has_reached_goal1 && norm(q_fix(:,n) - q_goal1) <= 1e-3
            has_reached_goal1 = true;
        end
        
        % --- 约束限制 ---
        % 强制高度约束 (例如不能低于 H_min)
        q_fix(3,n) = max(q_fix(3,n), P.H_min);
    end
end