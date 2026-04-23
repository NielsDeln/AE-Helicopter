clear;
clc; 
close all;

%% Setup
run('constants.m')

%% Gains
Kp_0 = -0.01;
Ki_0 = 0;
Kd_0 = 0;
N_0 = 0;

Kp_c = -0.01;
Ki_c = 0;
Kd_c = 0;
N_c = 0;

%% Signal bus setup
% States definition
states.u = 0;
states.w = 0;
states.q = 0;
states.theta_f = 0;
states.V_z = 0;

Simulink.Bus.createObject(states);

% Controls definition
ctrl.theta_0 = 0;
ctrl.theta_c = 0;

Simulink.Bus.createObject(ctrl);

% Derivative definition
xdot.u_dot = 0;
xdot.w_dot = 0;
xdot.q_dot = 0;
xdot.theta_f_dot = 0;
xdot.V_z = 0;
xdot.pos = 0;
xdot.h = 0;

Simulink.Bus.createObject(xdot);

% Environment definition
env.adv = 0;
env.lambda_c = 0;
env.alpha_c = 0;

Simulink.Bus.createObject(env);

%% additional param
V0 = 40 * kts_to_ms;
p.rho = 1.225; % [kg/m^3]
p.g = 9.81; % [m/s^2]
p.m = MTOW; % [kg]
p.W = p.m * p.g; % [N]


p.OmegaR = tip_speed_main;
p.sigma = solidity_main;
p.A = A_main;
p.A_eq = A_eq;


p.cl_alpha = 5.73;
p.Iy = 22713;
p.I_bl = 524.4986;
p.Lock = p.rho * p.cl_alpha * c_main * R_main^4 / p.I_bl; % lock number
p.h = 0.9;

%% Initial conditions
D0 = p.A_eq * 1/2 * p.rho * V0^2;
theta_f_0 = atan2(-D0,p.W);
u0 = V0 * cos(theta_f_0);
w0 = V0 * sin(theta_f_0);
q0 = 0;

%% Trim calculation

[theta_0_trim, theta_c_trim] = calc_trim(V0, p);

function [theta_0, theta_c] = calc_trim(V, p)
    Dpar = p.A_eq * 1/2 * p.rho * V^2;
    adv = V/p.OmegaR;
    T = sqrt(p.W^2 + Dpar^2);

    CT = T / (p.rho * p.A * p.OmegaR^2);

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

%%
disp('Setup successful')