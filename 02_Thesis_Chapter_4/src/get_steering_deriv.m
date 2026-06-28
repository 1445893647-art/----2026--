%% --- 辅助函数 4: 导向矢量导数 [cite: 339-346] ---
function [da_tx, da_rx] = get_steering_deriv(theta, phi, P, type)
    % 此函数同时计算发射和接收导向矢量的导数
    % type: 'theta' or 'phi'
    
    % 基础导向矢量
    atx_x = exp(-1j * pi * (0:P.N_tx-1).' * sin(theta) * cos(phi));
    atx_y = exp(-1j * pi * (0:P.N_ty-1).' * sin(theta) * sin(phi));
    
    arx_x = exp(-1j * pi * (0:P.N_rx-1).' * sin(theta) * cos(phi));
    arx_y = exp(-1j * pi * (0:P.N_ry-1).' * sin(theta) * sin(phi));
    
    % 对角矩阵 D
    D_tx = diag(0:P.N_tx-1); D_ty = diag(0:P.N_ty-1);
    D_rx = diag(0:P.N_rx-1); D_ry = diag(0:P.N_ry-1);
    
    % 空间相位偏导因子 Phi [cite: 343-346]
    if strcmp(type, 'theta')
        Phi_x = cos(theta) * cos(phi);
        Phi_y = cos(theta) * sin(phi);
    else % phi
        Phi_x = -sin(theta) * sin(phi);
        Phi_y = sin(theta) * cos(phi);
    end
    
    % 分量导数 [cite: 341]
    datx_x = -1j * pi * Phi_x * D_tx * atx_x;
    datx_y = -1j * pi * Phi_y * D_ty * atx_y;
    
    darx_x = -1j * pi * Phi_x * D_rx * arx_x;
    darx_y = -1j * pi * Phi_y * D_ry * arx_y;
    
    % 链式法则 [cite: 339]
    % da = da_x kron a_y + a_x kron da_y
    da_tx = kron(datx_x, atx_y) + kron(atx_x, datx_y);
    da_rx = kron(darx_x, arx_y) + kron(arx_x, darx_y);
end