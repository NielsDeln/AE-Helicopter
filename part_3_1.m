clear;
clc; 
close all;

%% 0 Load constants and Design Params.
run('constants.m');
plots = true;

cl_alpha = 2*pi; % ASSUMPTION Lift slope curve
sigma = solidity_main;

V = linspace(0, 79.7, 300);   % Forward speed [m/s] (79.7 is max speed)
OmegaR = tip_speed_main;
adv = V / OmegaR;

Dpar = A_eq * 1/2 * rho * V.^2;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Trim Calculation: questio 3.1 %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Iterative Loops
% Caluclate Glauert thrust coefficient
T_glau = sqrt(W^2 + Dpar.^2);
CT_glau = T_glau / (rho * A_main * OmegaR^2);

% Preallocate
lambda_i = zeros(1, length(V));
theta_c  = zeros(1, length(V));
theta_0  = zeros(1, length(V));
lambda_c = zeros(1, length(V));
CT_ver   = zeros(1, length(V));

for k = 1:length(V)
    % Define iterated values
    CT_k = CT_glau(k);
    Dpar_k = Dpar(k);
    V_k = V(k);
    adv_k = adv(k);

    % Sovle for lambda_i
    glauert_residual = @(lambda_i) 2*lambda_i * ...
                        sqrt((adv_k * cos(Dpar_k/W))^2 ...
                        + (adv_k * sin(Dpar_k/W) + lambda_i)^2) - CT_k;

    lambda_i_k = fzero(glauert_residual, sqrt(CT_k/2));
    
    % Solve matrix system
    A = [1 + 3/2 * adv_k^2, -8/3 * adv_k
        -adv_k,              2/3 + adv_k^2];

    b = [-2 * adv_k^2 * Dpar_k/W - 2 * adv_k * lambda_i_k
         4/sigma * CT_k/cl_alpha + adv_k * Dpar_k/W + lambda_i_k];
    x = A \ b;

    theta_c_k = x(1);
    theta_0_k = x(2);

    lambda_c_k = adv_k * theta_c_k + adv_k * Dpar_k/W;

    CT_ver_k = cl_alpha * sigma/4 * (2/3 * theta_0_k * (1 + 3/2 *adv_k^2) ...
                - (lambda_c_k + lambda_i_k));

    % Store
    lambda_i(k) = lambda_i_k;
    theta_c(k)  = rad2deg(theta_c_k); % converted to deg
    theta_0(k)  = rad2deg(theta_0_k); % converted to deg
    lambda_c(k) = lambda_c_k;
    CT_ver(k)   = CT_ver_k;
end


%% Results and graphs
% Check result with other equation
max_CT_err = max(abs(CT_ver - CT_glau));
fprintf('Max blade-element C_T residual: %.2e  (should be ~0)\n', max_CT_err)

% Plot a1 and theta0 versus flight speed V
if plots
    figure;
    plot(V, theta_c, 'LineWidth', 1.5); hold on;
    plot(V, theta_0, 'LineWidth', 1.5);
    ylabel('Control deflection [deg]', 'Interpreter', 'latex', 'FontSize', 14);
    xlabel('$V$ [m/s]', 'Interpreter', 'latex', 'FontSize', 14);
    grid on;
    legend('Long. cyclic, \theta_c','Collective, \theta_0','Location','best');
    set(gca, 'FontSize', 12);
end

