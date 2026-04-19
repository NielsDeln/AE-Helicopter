clear;
clc; 
close all;

plots = false;

%% 0 Load constants and Design Params.
run('part_1.m');
cl_alpha = 2*pi; % ASSUMPTION Lift slope curve


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Trim Calculation: questio 3.1 %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Iterative Loops
% Caluclate Glauert thrust coefficient
T_glau = sqrt(W^2 + Dpar.^2);
C_T_glau = T_glau / (rho * A_main * OmegaR_m^2);

% Preallocate
lambda_i = zeros(1, length(V));
a1       = zeros(1, length(V));
theta0   = zeros(1, length(V));
lambda_c = zeros(1, length(V));

for k = 1:length(V)
    % Define iterated values
    C_T_k = C_T_glau(k);
    Dpar_k = Dpar(k);
    V_k = V(k);
    adv_k = adv_ratio_m(k);

    % Sovle for lambda_i
    glauert_residual = @(lambda_i) 2*lambda_i * ...
                        sqrt((adv_k * cos(Dpar_k/W))^2 ...
                        + (adv_k * sin(Dpar_k/W) + lambda_i)^2) - C_T_k;

    lambda_i_k = fzero(glauert_residual, sqrt(C_T_k/2));
    
    % Solve matrix system
    A = [1 + 3/2 * adv_k^2, -8/3 * adv_k
        -adv_k,              2/3 + adv_k^2];

    b = [-2 * adv_k^2 * Dpar_k/W - 2 * adv_k * lambda_i_k
         4/sigma_m * C_T_k/cl_alpha + adv_k * Dpar_k/W + lambda_i_k];
    x = A \ b;
    a1_k = x(1);
    theta0_k = x(2);

    lambda_c_k = adv_k * a1_k + adv_k * Dpar_k/W;

    % Store
    lambda_i(k) = lambda_i_k;
    a1(k)       = a1_k * 180/pi; % converted to deg
    theta0(k)   = theta0_k * 180/pi;% converted to deg
    lambda_c(k) = lambda_c_k;
end


%% Results and graphs
% Check result with other equation
C_T_ver = cl_alpha * sigma_m/4 .* ( 2/3 * theta0.*(1 + 3/2 * adv_ratio_m.^2) ...
                               - (lambda_c + lambda_i) );

max_CT_err = max(abs(C_T_ver - C_T_glau));
fprintf('Max blade-element C_T residual: %.2e  (should be ~0)\n', max_CT_err)

% Plot a1 and theta0 versus flight speed V
if plots
    figure;
    plot(V, a1, 'LineWidth', 1.5); hold on;
    plot(V, theta0, 'LineWidth', 1.5);
    ylabel('Control deflection [deg]', 'FontSize', 14);
    xlabel('Velocity V [m/s]', 'FontSize', 14);
    grid on;
    legend('Long. cyclic, a_1','Collective, \theta_0','Location','best');
    set(gca, 'FontSize', 12);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Numerical Simulation: questio 3.2 %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Setup


%% Non-linear EOM

