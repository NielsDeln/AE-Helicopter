clear; clc; close all;
run('constants.m');

%% Parameters
p.rho = rho;
p.g = g;
p.m = MTOW;
p.W = W;

p.R = R_main;                    % [m]
p.Omega_tip = tip_speed_main;    % [m/s] = Omega*R
p.Omega = p.Omega_tip / p.R;     % [rad/s]

p.sigma = solidity_main;
p.A = A_main;
p.A_eq = A_eq;

p.cl_alpha = 5.73;
p.Iy = 22713;
p.I_bl = 524.4986;
p.Lock = p.rho * p.cl_alpha * c_main * R_main^4 / p.I_bl;

p.h = 0.9;
p.hinge_offset = hinge_offset;

%% Hover trim
trim = hover_trim(p);

fprintf('Hover trim:\n');
fprintf('theta_0 = %.4f deg\n', rad2deg(trim.theta_0));
fprintf('theta_c = %.4f deg\n', rad2deg(trim.theta_c));
fprintf('lambda_i = %.6f\n', trim.lambda_i);
fprintf('CT = %.6f\n', trim.CT);
fprintf('T = %.2f N\n\n', trim.T);

%% Hover equilibrium state
states0.u = 0;
states0.w = 0;
states0.q = 0;
states0.theta_f = 0;
states0.V_z = 0;

ctrl.theta_0 = trim.theta_0;
ctrl.theta_c = trim.theta_c;

%% Perturbation sizes
du = 0.1;               % [m/s]
dq = deg2rad(0.1);      % [rad/s]

%% Baseline dynamics at hover
out0 = evaluate_dynamics(ctrl, states0, p);

fprintf('Baseline hover dynamics:\n');
fprintf('u_dot = %.8f m/s^2\n', out0.xdot.u_dot);
fprintf('q_dot = %.8f rad/s^2\n\n', out0.xdot.q_dot);

%% Xu and Mu from perturbation in u
states_p = states0;
states_m = states0;

states_p.u = +du;
states_m.u = -du;

out_pu = evaluate_dynamics(ctrl, states_p, p);
out_mu = evaluate_dynamics(ctrl, states_m, p);

Xu = (out_pu.X - out_mu.X) / (2*du) / p.m;
Mu = (out_pu.M - out_mu.M) / (2*du) / p.Iy;

%% Mq from perturbation in q
states_p = states0;
states_m = states0;

states_p.q = +dq;
states_m.q = -dq;

out_pq = evaluate_dynamics(ctrl, states_p, p);
out_mq = evaluate_dynamics(ctrl, states_m, p);

Mq = (out_pq.M - out_mq.M) / (2*dq) / p.Iy;

%% Phugoid approximation
if abs(Mq) < 1e-12
    error('Mq evaluated to zero. Phugoid approximation cannot be computed.');
end

omega_p_sq = -p.g * Mu / Mq;

if omega_p_sq <= 0
    warning('omega_p^2 <= 0. Check signs/derivatives.');
end

omega_p = sqrt(abs(omega_p_sq));
zeta_p = -(Xu + p.g * Mu / Mq^2) / (2 * omega_p);

A = Xu + p.g * Mu / Mq^2;
B = p.g * Mu / Mq;

disc = A^2 + 4*B;

lambda1 = 0.5 * (A + sqrt(disc));
lambda2 = 0.5 * (A - sqrt(disc));

T_period = 2*pi / omega_p;

%% Print results
fprintf('Stability derivatives:\n');
fprintf('Xu = %.8f 1/s\n', Xu);
fprintf('Mu = %.8f 1/(m s)\n', Mu);
fprintf('Mq = %.8f 1/s\n\n', Mq);

fprintf('Phugoid characteristics:\n');
fprintf('omega_p = %.8f rad/s\n', omega_p);
fprintf('zeta_p  = %.8f\n', zeta_p);
fprintf('lambda1 = %.8f %+ .8fi 1/s\n', real(lambda1), imag(lambda1));
fprintf('lambda2 = %.8f %+ .8fi 1/s\n', real(lambda2), imag(lambda2));
fprintf('Period  = %.4f s\n', T_period);

%% =============================
% Local functions
% =============================

function trim = hover_trim(p)
    T = p.W;
    CT = T / (p.rho * p.A * p.Omega_tip^2);
    lambda_i = sqrt(CT/2);

    theta_0 = (3/2) * (lambda_i + 4*CT/(p.cl_alpha*p.sigma));
    theta_c = 0;

    trim.theta_0 = theta_0;
    trim.theta_c = theta_c;
    trim.lambda_i = lambda_i;
    trim.CT = CT;
    trim.T = T;
end

function out = evaluate_dynamics(ctrl, states, p)
    [adv, lambda_c, alpha_c] = calc_env(ctrl, states, p);

    env.adv = adv;
    env.lambda_c = lambda_c;
    env.alpha_c = alpha_c;

    lambda_i = solve_lambda_i(ctrl, env, states, p);
    [a1, CT_glau, CT_BEM] = calc_aero(ctrl, env, lambda_i, states, p);

    [xdot, T, D] = calc_EOM(ctrl, a1, CT_glau, states, p);

    X = p.m * xdot.u_dot;
    M = p.Iy * xdot.q_dot;

    out.env = env;
    out.lambda_i = lambda_i;
    out.a1 = a1;
    out.CT_glau = CT_glau;
    out.CT_BEM = CT_BEM;
    out.xdot = xdot;
    out.T = T;
    out.D = D;
    out.X = X;
    out.M = M;
end

function [adv, lambda_c, alpha_c] = calc_env(ctrl, states, p)
    u = states.u;
    w = states.w;
    theta_c = ctrl.theta_c;

    V = sqrt(u^2 + w^2);

    if V < 1e-8
        adv = 0;
        lambda_c = 0;
        alpha_c = theta_c;
    else
        alpha_c = theta_c - atan2(w, u);
        adv = V / p.Omega_tip * cos(alpha_c);
        lambda_c = V / p.Omega_tip * sin(alpha_c);
    end
end

function lambda_i = solve_lambda_i(ctrl, env, states, p)
    theta_0 = ctrl.theta_0;
    q = states.q;
    adv = env.adv;
    lambda_c = env.lambda_c;
    alpha_c = env.alpha_c;

    F = @(lambda_i) residual(lambda_i, adv, theta_0, lambda_c, q, alpha_c, p);

    lambda_i0 = 0.05;
    lambda_i = fzero(F, lambda_i0);
end

function F = residual(lambda_i, adv, theta_0, lambda_c, q, alpha_c, p)
    num = (8/3)*adv*theta_0 ...
        - 2*adv*(lambda_c + lambda_i) ...
        - (16/p.Lock)*(q/p.Omega);

    den = 1 - 0.5*adv^2;
    if abs(den) < 1e-8
        den = sign(den + eps) * 1e-8;
    end

    a1 = num / den;

    CT_BEM = p.cl_alpha * 0.25 * p.sigma * ...
             ((2/3)*theta_0*(1 + 1.5*adv^2) - (lambda_c + lambda_i));

    Vc_norm = adv * cos(alpha_c - a1);
    Vs_norm = adv * sin(alpha_c - a1);

    CT_glau = 2 * lambda_i * sqrt(Vc_norm^2 + (Vs_norm + lambda_i)^2);

    F = CT_BEM - CT_glau;
end

function [a1, CT_glau, CT_BEM] = calc_aero(ctrl, env, lambda_i, states, p)
    theta_0 = ctrl.theta_0;
    q = states.q;
    adv = env.adv;
    lambda_c = env.lambda_c;
    alpha_c = env.alpha_c;

    num = (8/3)*adv*theta_0 ...
        - 2*adv*(lambda_c + lambda_i) ...
        - (16/p.Lock)*(q/p.Omega);

    den = 1 - 0.5*adv^2;
    if abs(den) < 1e-8
        den = sign(den + eps) * 1e-8;
    end

    a1 = num / den;

    CT_BEM = p.cl_alpha * 0.25 * p.sigma * ...
             ((2/3)*theta_0*(1 + 1.5*adv^2) - (lambda_c + lambda_i));

    Vc_norm = adv * cos(alpha_c - a1);
    Vs_norm = adv * sin(alpha_c - a1);

    CT_glau = 2 * lambda_i * sqrt(Vc_norm^2 + (Vs_norm + lambda_i)^2);
end

function [xdot, T, D] = calc_EOM(ctrl, a1, CT, states, p)
    u = states.u;
    w = states.w;
    q = states.q;
    theta_f = states.theta_f;
    V_z = states.V_z;
    theta_c = ctrl.theta_c;

    V = sqrt(u^2 + w^2);
    if V < 1e-8
        V = 1e-8;
    end

    T = CT * p.rho * p.Omega_tip^2 * p.A;
    D = p.A_eq * 0.5 * p.rho * V^2;

    if V < 0.1
        xdot.u_dot = -p.g * sin(theta_f) + T/p.m * sin(theta_c - a1);
        xdot.w_dot =  p.g * cos(theta_f) - T/p.m * cos(theta_c - a1);
    else
        xdot.u_dot = -p.g * sin(theta_f) ...
                   - D/p.m * u/V ...
                   + T/p.m * sin(theta_c - a1) ...
                   - q*w;

        xdot.w_dot =  p.g * cos(theta_f) ...
                   - D/p.m * w/V ...
                   - T/p.m * cos(theta_c - a1) ...
                   + q*u;
    end

    h_eff = p.h + p.hinge_offset * p.R;
    xdot.q_dot = -T/p.Iy * h_eff * sin(theta_c - a1);

    xdot.theta_f_dot = q;

    xdot.V_z = -xdot.w_dot * cos(theta_f) ...
               + w * sin(theta_f) * q ...
               + xdot.u_dot * sin(theta_f) ...
               + u * cos(theta_f) * q;

    xdot.pos = u * cos(theta_f) - w * sin(theta_f);
    xdot.h = V_z;
end