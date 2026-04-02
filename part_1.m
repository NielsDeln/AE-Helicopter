clear;
clc; 
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Induced Velocity: questio 1.2 %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 0 Load constants
run('constants.m');

% atmospheric / flight assumptions
rho = 1.225;         % Air density at sea level [kg/m^3]
g   = 9.80665;          % Gravity [m/s^2]

% Choose helicopter mass for the calculation
m = MTOW;            % [kg]
W = m * g;           % Weight [N]

%% 1 Hover induced velocity
vi_hover = sqrt(W / (2 * rho * A_main));   % [m/s]

fprintf('Hover induced velocity: %.3f m/s\n', vi_hover);

%% 2 forward speed range and preallocate
V = linspace(0, 100, 300);   % Forward speed [m/s]

% preallocate vector vi
vi_glauert   = zeros(size(V));   % Numerical Glauert solution
vi_lowspeed  = zeros(size(V));   % Lowspeed closed-form solution

%% 3 Disc angle of attack assumption
alpha_d = 0;   % [rad]

%% 4 solve induced velocity in forward flight
for i = 1:length(V)
    Vi = V(i);

    if Vi == 0
        % At hover
        vi_glauert(i)   = vi_hover;
        vi_lowspeed(i)  = vi_hover;
    else
        % --- General Glauert implicit equation --
        % W = 2*rho*A*vi*sqrt((V*cos(alpha_d))^2 + (V*sin(alpha_d)+vi)^2)
        f = @(vi) 2 * rho * A_main * vi .* ...
            sqrt((Vi*cos(alpha_d))^2 + (Vi*sin(alpha_d) + vi).^2) - W;

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
sigma_m = solidity_main;
OmegaR = tip_speed_main;
C_D_p = input('Blade medium drag coefficient determined graphically:'); % DETERMINE GRAPHICALLY AS FUNCTION OF MED LIFT 

% Actuator disc theory
Pideal = W * vi_hover; % Ideal power in hover for MTOW in W
FOM_ACT = 0.6; % Assumed Figure of merit
Phov_ACT = Pideal / FOM_ACT; % True hover power in ACT assuming an FOM

% Blade element theory
k_main = 1.1; % Assumption
Pi_main = k_main * Pideal; %
Pp_main = sigma_m * C_D_p / 8 * rho * OmegaR^3 * A_main;
Phov_BEM = Pi_main + Pp_main; % Power is induced plus profile drag
FOM_BEM = Pideal / Phov_BEM; % FOM according to BEM theory

%% 2 Forward Flight Calculations
% Main rotor
adv_ratio = V / OmegaR; % Array of advance ratios
Pp_fw_main = Pp_main * (1 + adv_ratio^2); % Profile power of main rotors for forward velocities
Pd_fw_main = Pp_main * 2 * adv_ratio^2; % Main rotor drag power for forward flight
Pi_fw_main = 0; % ADD CALCULATION OF INDUCED POWER AS FUNC OF FWRD VEL.

P_benett = Pp_main * (1 + 4.65 * adv_ratio^2); % Using Benett approx.

