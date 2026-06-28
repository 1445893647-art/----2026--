function [theta, phi, dist] = get_angle_dist(q_w, q_b)
% GET_ANGLE_DIST 计算两点间的距离和角度
% 输入:
%   q_w: Willie 的位置坐标 [xw; yw; zw]
%   q_b: Bob 的位置坐标 [xb; yb; zb]
% 输出:
%   theta: 俯仰角 (rad) [cite: 64]
%   phi:   方位角 (rad) [cite: 64]
%   dist:  距离 (m) [cite: 63]

    % 计算位置差向量
    diff_p = q_b - q_w;
    
    % 1. 计算距离
    dist = norm(diff_p);
    
    % 2. 计算方位角 phi (Azimuth)
    % 使用 atan2 可以自动处理所有四个象限的情况 (-pi 到 pi)
    phi = atan2(diff_p(2), diff_p(1)); 
    
    % 3. 计算俯仰角 theta (Elevation)
    % 根据论文公式: theta = arccos((Hb - Hw) / d)
    % 这里 diff_p(3) 即为 zb - zw
    if dist < 1e-9
        theta = 0; % 防止除以零
    else
        % 限制输入范围在 [-1, 1] 之间，防止数值误差导致复数结果
        val = diff_p(3) / dist;
        val = max(min(val, 1), -1);
        theta = acos(val);
    end
end