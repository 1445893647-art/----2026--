%% --- 辅助函数 3: UPA 导向矢量生成 [cite: 67-73] ---
function a = get_steering_vec(theta, phi, Nx, Ny)
    % 水平与垂直分量
    idx_x = (0:Nx-1).';
    idx_y = (0:Ny-1).';
    
    % 公式 [cite: 71-72]
    ax = exp(-1j * pi * idx_x * sin(theta) * cos(phi));
    ay = exp(-1j * pi * idx_y * sin(theta) * sin(phi));
    
    % Kronecker 积 [cite: 67]
    a = kron(ax, ay);
end