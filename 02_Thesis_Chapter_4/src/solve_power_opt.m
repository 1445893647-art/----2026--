function P_J_opt = solve_power_opt(q_w, P)
%SOLVE_POWER_OPT Power optimization subproblem (paper closed-form).
%
% The paper uses the Bhattacharyya upper bound:
%   xi^u = sqrt(pi0*pi1) * B(gamma_w) <= delta
% with gamma_w = eta_S / eta_I = (P_a*beta_aw) / (rho*P_J).
% Since B(gamma_w) is strictly decreasing, the constraint is equivalent to
%   gamma_w >= Gamma_req  (Gamma_req solved from B(Gamma_req)=2*delta for pi0=pi1=0.5)
% which yields an upper bound on P_J:
%   P_J <= (P_a * beta_aw) / (rho * Gamma_req).
% To maximize jamming effectiveness (and thus Bob outage), we take equality
% and clip by power budget P_max.

    % Geometry
    d_aw = norm(P.q_a - q_w);
    d_aw = max(d_aw, 1); % avoid singular path-loss

    % Large-scale path loss (Alice -> Willie)
    beta_aw = P.beta_0 * (d_aw ^ (-P.alpha_aw));

    % Required SINR Gamma_req for the Bhattacharyya constraint
    if isfield(P, 'Gamma_req') && ~isempty(P.Gamma_req)
        Gamma_req = P.Gamma_req;
    else
        Gamma_req = get_gamma_req(P);
    end

    % Closed-form power upper bound
    P_J_upper = (P.P_a * beta_aw) / (P.rho * Gamma_req);

    % Optimal power: saturate at the upper bound and max budget
    P_J_opt = min(P.P_max, P_J_upper);

    % Numerical safety
    P_J_opt = max(P_J_opt, 1e-6);
end
