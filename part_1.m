clear;
clc; 
close all;

plots = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Induced Velocity: questio 1.2 %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 0 Load constants and Design Params.
run('constants.m');

% Choose helicopter mass for the calculation
m = MTOW;            % [kg]

% Power calculations
sigma_m = solidity_main;
OmegaR_m = tip_speed_main;
sigma_t = solidity_tail;
OmegaR_t = tip_speed_tail;


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
if plots
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
end

%% dimensionless plot
Vbar = V / vi_hover;
vibar_glauert = vi_glauert / vi_hover;
vibar_low     = vi_lowspeed / vi_hover;


if plots
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
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Induced Velocity: questio 1.3 %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 1 Hover Calculations
% Coefficient calculations
C_T_m = W / (rho * A_main * OmegaR_m^2); % [-] Approxed as weight coeff.
C_Lbar_m = 6.6 * C_T_m / sigma_m; % [-] Medium lift coeff.
C_Dp_m = 0.00568; % PROPER VALUE NEEDS TO BE FOUND

% Actuator disc theory
Pideal = W * vi_hover; % Ideal power in hover for MTOW [W]
FOM_ACT = 0.67; % Assumed Hover Figure of Merit (found in 1977 paper)
Phov_ACT = Pideal / FOM_ACT; % True hover power in ACT assuming an FOM


% Blade element theory
k_main = 1.15; % Assumption
Pi_main = k_main * Pideal; %
Pp_main = sigma_m * C_Dp_m / 8 * rho * OmegaR_m^3 * A_main;
Phov_BEM = Pi_main + Pp_main; % Power is induced plus profile drag
FOM_BEM = Pideal / Phov_BEM; % FOM according to BEM theory

%% 2 Forward Flight Calculations
% Main rotor
adv_ratio_m = V / OmegaR_m; % Array of advance ratios

Ppar_fw_m = Dpar .* V; % Main rotor parasitic power
Pi_fw_m = k_main * W * vi_glauert; % glauert induced velocity is taken
P_benett_m = Pp_main * (1 + 4.65 * adv_ratio_m.^2); % Using Benett approx.

P_main = P_benett_m + Pi_fw_m + Ppar_fw_m; % Main rotor power

% Tail rotor
adv_ratio_t = V / OmegaR_t;

% Tail rotor thrust
k_tr = 1.4;
l_tail = l_LOA - R_main - R_tail;
T_tail = P_main / (OmegaR_m/R_main * l_tail);

% coefficient calculations
C_T_t = T_tail / (rho * A_tail * OmegaR_t^2); % [-] Approxed as weight coeff.
C_Lbar_t = 6.6 * C_T_t / sigma_t; % [-] Medium lift coeff.
C_Dp_t = 0.0075; % value taken from graph in lecture notes

vi_tail = sqrt(T_tail / (2 * rho * A_tail));
Pi_fw_t = 1.1 * k_tr * T_tail .* vi_tail;

Pp_tail = sigma_t * C_Dp_t / 8 * rho * OmegaR_t^3 * A_tail;
P_benett_t = Pp_tail * (1 + 4.65 * adv_ratio_t.^2);

P_tail = Pi_fw_t + P_benett_t + Pp_tail; % Tail rotor power

% Total power 
% Total power
P_total = P_main + P_tail;

% Minimum power and corresponding speed
[P_min, idx_min] = min(P_total);
V_Pmin = V(idx_min);

fprintf('Minimum Power required: %.3f KW\n', P_min/1e3);
fprintf('Forward Velocity for minimum power required: %.3f m/s\n', V_Pmin);

% Compute power per speed (slope from origin) and find its minimum.
slope = P_total ./ V;
% Avoid division by zero at V=0 by setting slope(1) to large value
if V(1) == 0
    slope(1) = inf;
end
[slope_min, idx_range] = min(slope);
V_range = V(idx_range);
P_range = P_total(idx_range);

fprintf('Power for maximum range (tangent slope): %.3f KW\n', P_range/1e3);
fprintf('Forward Velocity for maximum range: %.3f m/s\n', V_range);

if plots
    figure;
    plot(V, P_total, 'LineWidth', 1.8); hold on;
    plot(V, P_main, 'LineWidth', 1.8);
    plot(V, P_tail, 'LineWidth', 1.8);
    plot(V, Ppar_fw_m, 'LineWidth', 1.8)
    plot(V, Pi_fw_m, 'LineWidth', 1.8)
    grid on;
    xlabel('$V$', 'Interpreter', 'latex', 'FontSize', 14);
    ylabel('$P_{total}$', 'Interpreter', 'latex', 'FontSize', 14);
    title('Total power versus forward speed', 'FontSize', 16);
    legend('Total Power', 'Location', 'northeast', 'FontSize', 12);
    ylim([0, 1e6]);
    set(gca, 'FontSize', 12);
end