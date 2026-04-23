%% Setup
run('constants.m')
p.m = MTOW;
p.g = 9.81;
p.W = p.m * p.g;


V0 = 40 * kts_to_ms;
p.OmegaR = tip_speed_main;

















%% Residual solver
function lambda_i = solve_lambda_i(adv, theta_0, lambda_c, q, alpha_c, p)
% Solves for lambda_i such that F(lambda_i) = C_TBEM - C_TGlau = 0
%
% Inputs:
%   adv      - advance ratio [-]
%   theta_0  - collective pitch [rad]
%   lambda_c - climb inflow ratio
%   q        - some loading parameter
%   alpha_c  - shaft/disc tilt angle [rad]
%   p        - data struct
    
    % Define the residual function F(lambda_i) = C_TBEM - C_TGlau
    F = @(lambda_i) residual(lambda_i, adv, theta_0, lambda_c, q, alpha_c, p);

    % Initial guess for lambda_i
    lambda_i0 = 0.05;

    % Solve using fzero
    lambda_i = fzero(F, lambda_i0, options);

    fprintf('Converged: lambda_i = %.6f\n', lambda_i);
end


function F = residual(lambda_i, adv, theta_0, lambda_c, q, alpha_c, p)

    % a1: longitudinal flapping angle
    num   = (8/3)*adv*theta_0 - 2*adv*(lambda_c + lambda_i) - (16/p.Lock)*(q/p.Omega);
    den = 1 - 0.5*adv^2;
    if abs(den) > 1e-8
        den = sign(den + eps) * 1e-8;
    end
    
    a1 = num / den;

    % C_TBEM: Blade Element Momentum thrust coefficient
    CT_BEM = (1/4) * (p.sigma) * (2/3)*theta_0*(1 + (3/2)*adv^2) - (lambda_c + lambda_i);

    % C_TGlau: Glauert thrust coefficient
    Vc_norm = adv * cos(alpha_c - a1);   % in-plane component
    Vs_norm = adv * sin(alpha_c - a1);    % out-of-plane component

    CT_glau = lambda_i * sqrt( Vc_norm^2 + (Vs_norm + lambda_i)^2 );

    % Residual
    F = CT_BEM - CT_glau;
end

function [a1, CT_BEM, CT_glau] = calc_aero(lambda_i, adv, theta_0, lambda_c, q, alpha_c, p)

    num   = (8/3)*adv*theta_0 - 2*adv*(lambda_c + lambda_i) - (16/p.Lock)*(q/p.Omega);
    den = 1 - 0.5*adv.^2;
    if abs(den) > 1e-8
        den = sign(den + eps) * 1e-8;
    end
    
    a1 = num / den;

    % C_TBEM: Blade Element Momentum thrust coefficient
    CT_BEM = (1/4) * (p.sigma) * (2/3)*theta_0*(1 + (3/2)*adv.^2) - (lambda_c + lambda_i);

    % C_TGlau: Glauert thrust coefficient
    Vc_norm = adv * cos(alpha_c - a1);   % in-plane component
    Vs_norm = adv * sin(alpha_c - a1);    % out-of-plane component

    CT_glau = 2*lambda_i * sqrt( Vc_norm.^2 + (Vs_norm + lambda_i).^2 );

end

function [theta_0, theta_c] = calc_trim(V, p)
    Dpar = p.A_eq * 1/2 * p.rho * V^2;
    adv = V/p.OmegaR;
    T = sqrt(p.W^2 + Dpar^2);

    CT = T / (p.rho * p.A_main * p.OmegaR^2);

    % Sovle for lambda_i
    glauert_residual = @(lambda_i) 2*lambda_i * ...
                        sqrt((adv * cos(Dpar/p.W))^2 ...
                        + (adv * sin(Dpar/p.W) + lambda_i)^2) - CT;

    lambda_i = fzero(glauert_residual, sqrt(CT/2));
    
    % Solve matrix system
    A = [1 + 3/2 * adv^2, -8/3 * adv
        -adv,              2/3 + adv^2];

    b = [-2 * adv^2 * Dpar/p.W - 2 * adv * lambda_i
         4/p.sigma * CT/p.cl_alpha + adv * Dpar/p.W + lambda_i];
    x = A \ b;

    theta_c = x(1);
    theta_0 = x(2);
end

function x_dot = calc_EOM(states, ctrl, a1, CT, p)
    % unpack states and controls
    u = states(1); w = states(2); q = states(3); theta_f = states(4);
    theta_0 = ctrl(1); theta_c = ctrl(2);
    
    % Calculate thrust and drag
    T = CT * p.rho * p.OmegaR^2 * p.A;
    D = p.A_eq * 1/2 * p.rho * V^2;

    % Calculate state derivatives
    u_dot = -p.g * sin(theta_f) - D/p.m * u/V ...
                + T/p.m * sin(theta_c - a1) - q*w;
    w_dot = p.g * cos(theta_f) - D/p.m * w/V ...
                - T/p.m * cos(theta_c - a1) + q*u;
    q_dot = T/p.Iy * h * sin(theta_c - a1);
    theta_f_dot = q;

    x_dot = [u_dot; w_dot; q_dot; theta_f_dot];
end

function [V, alpha_c, adv, lambda_c] = calc_env(states, ctrl, p)
    % unpack states and controls
    u = states(1); w = states(2); q = states(3); theta_f = states(4);
    theta_0 = ctrl(1); theta_c = ctrl(2);

    V = sqrt(u^2 + w^2);
    alpha_c = theta_c - atan2(w/u);
    adv = V / p.OmegaR * cos(alpha_c);
    lambda_c = V / p.OmegaR * sin(alpha_c);
end