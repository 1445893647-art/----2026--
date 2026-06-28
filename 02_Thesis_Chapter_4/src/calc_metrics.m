function [P_out, xi_star, tau_opt] = calc_metrics_fixR(q_w, q_b, P_J, w_wb, w_ab, P)
%CALC_METRICS_FIXR  与“相关散射 R_wb”一致的指标计算（修复版）
% - 关键修复：Willie->Bob 的等效标量信道 g = h_wb^H w_wb 在相关散射下满足
%       g ~ CN(mu, sigma_g2),
%   其中 sigma_g2 = beta_wb*(1/(K+1))*(w^H R_wb w)，而不是 beta_wb*(1/(K+1))*||w||^2。
%   你之前的 calc_metrics.m 用的是 ||w||^2（等价于假设 R_wb = I），
%   这会导致：你在波束优化里引入 R_wb（各向异性散射）但评估指标仍按各向同性算，
%   从而出现“优化不起作用/看起来还不如固定指向波束”的现象。
%
% - 同时：强制 w_wb 单位范数，避免 baseline 因为未归一化而出现“中断概率=1饱和”。

    %#ok<NASGU> % w_ab kept for future extensions

    % 强制单位范数（防止固定波束方案功率被放大）
    %w_wb = w_wb(:);
    %w_wb = w_wb / max(norm(w_wb), 1e-12);

    % === 1) Outage probability at Bob ===
    d_aw = norm(P.q_a - q_w);
    d_ab = norm(P.q_a - q_b);
    d_wb = norm(q_w - q_b);

    [theta, phi, ~] = get_angle_dist(q_w, q_b);

    % Path loss
    d_ab = max(d_ab, 1);
    d_wb = max(d_wb, 1);
    beta_ab = P.beta_0 * (d_ab ^ (-P.alpha_ab));

    % alpha_wb 优先，否则沿用 alpha_aw（保持与你旧代码兼容）
    if isfield(P,'alpha_wb') && ~isempty(P.alpha_wb)
        alpha_wb = P.alpha_wb;
    else
        alpha_wb = P.alpha_aw;
    end
    beta_wb = P.beta_0 * (d_wb ^ (-alpha_wb));

    % Alice-to-Bob (Rayleigh) + MRT: Z = ||h_ab||^2 ~ Gamma(N_a, beta_ab)
    Pa = P.P_a;
    Omega_ab = Pa * beta_ab;

    gamma_th = 2^P.r_th - 1;

    lambda = gamma_th / (Omega_ab + 1e-30);
    c = lambda * P.sigma_b^2;
    b = lambda * P_J;

    % Willie->Bob 等效标量信道 g
    % mu = sqrt(beta)*sqrt(K/(K+1)) * a^H w
    % sigma_g2 = beta*(1/(K+1)) * w^H R w
    a = get_steering_vec(theta, phi, P.N_tx, P.N_ty);   % 不强制归一化：保留阵列增益
    gain_los = abs(a' * w_wb)^2;

    if isfield(P,'R_wb') && ~isempty(P.R_wb)
        R = P.R_wb;
        sigma_g2 = beta_wb * (1/(P.K + 1)) * real(w_wb' * R * w_wb);
        sigma_g2 = max(sigma_g2, 1e-30);
    else
        % 退化为各向同性散射：R=I
        sigma_g2 = beta_wb * (1/(P.K + 1)) * (w_wb' * w_wb); % ==1
    end

    mu2 = beta_wb * (P.K/(P.K + 1)) * gain_los;

    Phi   = b * sigma_g2;
    Theta = b * mu2;
    nu    = mu2 / (sigma_g2 + 1e-30);

    % Success probability under MRT (Gamma shape = N_a)
    Na = P.N_a;
    onePlusPhi = 1 + Phi;
    argL = -nu / onePlusPhi;
    sumTerm = 0;

    for m = 0:(Na-1)
        for k = 0:m
            mk = m - k;
            if mk == 0
                coeff1 = 1;
            elseif c == 0
                coeff1 = 0;
            else
                coeff1 = exp(mk*log(c) - gammaln(mk + 1));
            end

            if k == 0
                coeff2 = 1 / onePlusPhi;
            elseif Phi == 0
                coeff2 = 0;
            else
                coeff2 = exp(k*log(Phi) - (k+1)*log(onePlusPhi));
            end

            Lk = laguerreL_int(k, argL);
            sumTerm = sumTerm + coeff1 * coeff2 * Lk;
        end
    end

    P_s = exp(-c) * exp( -Theta / onePlusPhi ) * sumTerm;
    P_out = 1 - P_s;
    P_out = min(max(real(P_out), 0), 1);

    % === 2) Willie detection error xi*（沿用你旧版本） ===
    d_aw = max(d_aw, 1);
    beta_aw = P.beta_0 * (d_aw ^ (-P.alpha_aw));

    eta_S_phys = P.P_a * beta_aw;
    eta_I_phys = P.rho * P_J + 1e-30;
    gamma_w = eta_S_phys / eta_I_phys;

    xi = calc_xi_exact_internal(gamma_w, P.N_r);
    xi_star = xi;
    tau_opt = NaN;
end

function xi = calc_xi_exact_internal(gamma_w, Nr)
    pi0 = 0.5; pi1 = 0.5;
    eta_I = 1;
    eta_S = gamma_w;

    sum_eta = eta_S + eta_I;
    sq_sum  = eta_S^2 + eta_I^2;
    nu_Z    = (Nr * sum_eta^2) / (sq_sum + 1e-30);
    eta_Z   = sq_sum / (sum_eta + 1e-30);

    A = (1/eta_I) - (1/eta_Z);
    B = Nr - nu_Z;
    C = log(pi0/pi1) + gammaln(nu_Z) - gammaln(Nr) + nu_Z*log(eta_Z) - Nr*log(eta_I);

    u_opt = nan;
    if abs(B) < 1e-6
        u_opt = C/A;
    else
        val = -(A/B)*exp(-C/B);
        try
            w_val = lambertw(0, real(val));
            u_opt = -(B/A)*w_val;
        catch
            u_opt = nan;
        end
    end

    if isnan(u_opt) || u_opt < 0
        calc_err = @(u) pi0*gammainc(u/eta_I, Nr, 'upper') + pi1*gammainc(u/eta_Z, nu_Z);
        [~, xi] = fminbnd(calc_err, 0, nu_Z*eta_Z*2);
    else
        P_FA = gammainc(u_opt/eta_I, Nr, 'upper');
        P_MD = gammainc(u_opt/eta_Z, nu_Z);
        xi = pi0 * P_FA + pi1 * P_MD;
    end
end

function Lk = laguerreL_int(k, x)
    if k == 0
        Lk = 1;
        return;
    end
    Lk = 0;
    for i = 0:k
        comb = exp(gammaln(k+1) - gammaln(i+1) - gammaln(k-i+1));
        Lk = Lk + comb * ((-x)^i) / exp(gammaln(i+1));
    end
end
