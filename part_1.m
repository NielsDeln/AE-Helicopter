clear;
clc; 
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Induced Velocity: questio 1.2 %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 0 Load constants and Design Params.
run('constants.m');

% atmospheric / flight assumptions
rho = 1.225;         % Air density at sea level [kg/m^3]
g   = 9.80665;          % Gravity [m/s^2]

% Choose helicopter mass for the calculation
m = MTOW;            % [kg]
W = m * g;           % Weight [N]

% Power calculations
sigma_m = solidity_main;
OmegaR = tip_speed_main;
C_T = W / (rho * A_main * OmegaR^2); % [-] Approxed as weight coeff.
C_Lbar = 6.6 * C_T / sigma_m; % [-] Medium lift coeff.
C_D_p = 0.015; % PROPER VALUE NEEDS TO BE FOUND
A_eq = 12.6 * ft_to_m^2; % Equivalent plate area Found in https://doi.org/10.2514/6.2024-1117

%% 1 Hover induced velocity
vi_hover = sqrt(W / (2 * rho * A_main));   % [m/s]

fprintf('Hover induced velocity: %.3f m/s\n', vi_hover);

%% 2 forward speed range and preallocate
V = linspace(0, 79.7, 300);   % Forward speed [m/s] (79.7 is max speed)

% preallocate vector vi
vi_glauert   = zeros(size(V));   % Numerical Glauert solution
vi_lowspeed  = zeros(size(V));   % Lowspeed closed-form solution

%% 3 Disc angle of attack calcualtion
Dpar = A_eq * 1/2 * rho * V.^2;
alpha_d = asin(Dpar / W);   % [rad] (flight path = 0 assumed)

%% 4 solve induced velocity in forward flight
for i = 1:length(V)
    Vi = V(i);
    alpha_di = alpha_d(i);

    if Vi == 0
        % At hover
        vi_glauert(i)   = vi_hover;
        vi_lowspeed(i)  = vi_hover;
    else
        % --- General Glauert implicit equation --
        % W = 2*rho*A*vi*sqrt((V*cos(alpha_d))^2 + (V*sin(alpha_d)+vi)^2)
        f = @(vi) 2 * rho * A_main * vi .* ...
            sqrt((Vi*cos(alpha_di)).^2 + (Vi*sin(alpha_di) + vi).^2) - W;

        % Use hover induced velocity as initial guess
        vi_guess = max(vi_hover/2, 0.1);

        vi_glauert(i) = fzero(f, vi_guess);

        % --- Low-speed analytical approximation from lecture slides --
        Vbar = Vi / vi_hover;
        vibar_low = sqrt(-Vbar^2/2 + sqrt(1 + Vbar^4/4));
        vi_lowspeed(i) = vibar_low * vi_hover;

    end
end

%% Plot
figure;
plot(V, vi_glauert, 'LineWidth', 1.8); hold on;
plot(V, vi_lowspeed, '--', 'LineWidth', 1.5);
grid on;
xlabel('Forward speed V [m/s]', 'FontSize', 14);
ylabel('Induced velocity v_i [m/s]', 'FontSize', 14);
title('Main rotor induced velocity versus forward speed', 'FontSize', 16);
legend('Glauert numerical', 'Low-speed approximation', ...
       'Location', 'northeast', 'FontSize', 12);
set(gca, 'FontSize', 12);

%% dimensionless plot
Vbar = V / vi_hover;
vibar_glauert = vi_glauert / vi_hover;
vibar_low     = vi_lowspeed / vi_hover;

figure;
plot(Vbar, vibar_glauert, 'LineWidth', 1.8); hold on;
plot(Vbar, vibar_low, '--', 'LineWidth', 1.5);
grid on;
xlabel('$\bar{V} = V / v_{i,h}$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$\bar{v}_i = v_i / v_{i,h}$', 'Interpreter', 'latex', 'FontSize', 14);
title('Non-dimensional induced velocity versus forward speed', 'FontSize', 16);
legend('Glauert numerical', 'Low-speed approximation', ...
       'Location', 'northeast', 'FontSize', 12);
set(gca, 'FontSize', 12);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Induced Velocity: questio 1.3 %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 1 Hover Calculations
% Actuator disc theory
Pideal = W * vi_hover; % Ideal power in hover for MTOW in W
FOM_ACT = 0.74; % Assumed Hover Figure of Merit (found in 1977 paper)
Phov_ACT = Pideal / FOM_ACT; % True hover power in ACT assuming an FOM

% Blade element theory
k_main = 1.15; % Assumption
Pi_main = k_main * Pideal; %
Pp_main = sigma_m * C_D_p / 8 * rho * OmegaR^3 * A_main;
Phov_BEM = Pi_main + Pp_main; % Power is induced plus profile drag
FOM_BEM = Pideal / Phov_BEM; % FOM according to BEM theory

%% 2 Forward Flight Calculations
% Main rotor
adv_ratio = V / OmegaR; % Array of advance ratios

Ppar_fw_main = Dpar .* V; % Main rotor parasitic power
Pi_fw_main = k_main * W * vi_glauert; % glauert induced velocity is taken
P_benett = Pp_main * (1 + 4.65 * adv_ratio.^2); % Using Benett approx.

P_main = P_benett + Pi_fw_main + Ppar_fw_main;

% Tail rotor
k_tr = 1.4;
l_tail = l_LOA - R_main - R_tail;
T_tail = P_main / (OmegaR/R_main * l_tail);
vi_tail = sqrt(T_tail / (2 * rho * A_tail));
Pi_fw_tail = 1.1 * k_tr * T_tail .* vi_tail;

P_total = P_main + Pi_fw_tail;

[P_min, idx_min] = min(P_total);
V_Pmin = V(idx_min);

fprintf('Minimum Power required: %.3f KW\n', P_min/1e3);
fprintf('Forward Velocity for minimum power required: %.3f m/s\n', V_Pmin);


figure;
plot(V, P_total, 'LineWidth', 1.8); hold on;
grid on;
xlabel('$V$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$P_{total}$', 'Interpreter', 'latex', 'FontSize', 14);
title('Total power versus forward speed', 'FontSize', 16);
legend('Total Power', 'Location', 'northeast', 'FontSize', 12);
ylim([0, 1e6]);
set(gca, 'FontSize', 12);