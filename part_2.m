clear; close all; clc;

run('constants.m');
rho = 1.225;
a_l = 5.7;

%% Inputs
V = 20;
q = deg2rad(20);
p = deg2rad(10);

theta0  = deg2rad(6);
theta1c = deg2rad(2);
theta1s = deg2rad(1);

vi = 5;

%% Rotor parameters
R = R_main;
c = c_main;
Omega = tip_speed_main / R_main;

mu = V / (Omega * R);

%% Azimuth and pitch
nPsi = 721;
psi = linspace(0, 2*pi, nPsi);
psi_deg = rad2deg(psi);

theta = theta0 + theta1c*cos(psi) + theta1s*sin(psi);

%% Flapping model


I_blade = J_main / N_blade_main;

mu = V / (Omega * R);
lambda_i = vi / (Omega * R);
lambda_c = 0;

gamma = rho * c * a_l * R^4 / I_blade;

a0 = (gamma/8) * ( theta0*(1 + mu^2) - (4/3)*(lambda_i + lambda_c) );

a1 = ( (8/3)*mu*theta0 - mu*(lambda_i + lambda_c) - q/Omega ) / (1 - mu^2/2);

b1 = ( -(4/3)*mu*(q/Omega) ) / (1 + mu^2/2);


beta = a0 - a1*cos(psi) - b1*sin(psi);
dbeta_dpsi = a1*sin(psi) - b1*cos(psi);
beta_dot = dbeta_dpsi * Omega;

%% Flapping angle plot
figure;
plot(psi_deg, rad2deg(beta), 'LineWidth', 1.8);
grid on;
xlabel('\psi [deg]');
ylabel('\beta [deg]');
title('Blade flapping angle');

%% Angle of attack
rbar_vec = linspace(0.2, 0.95, 30);
alpha_mat = zeros(length(rbar_vec), length(psi));

for ir = 1:length(rbar_vec)
    r = rbar_vec(ir) * R;
    UT = Omega*r + V*sin(psi);
    UP = vi + r*beta_dot + q*r*cos(psi) + p*r*sin(psi);
    phi = atan2(UP, UT);
    alpha_mat(ir,:) = theta - phi;
end

%% Angle of attack plot
figure; hold on; grid on;
rbar_plot = [0.3 0.5 0.7 0.9];

for k = 1:length(rbar_plot)
    [~, idx] = min(abs(rbar_vec - rbar_plot(k)));
    plot(psi_deg, rad2deg(alpha_mat(idx,:)), 'LineWidth', 1.5);
end

xlabel('\psi [deg]');
ylabel('\alpha [deg]');
title('Angle of attack variation');
legend('r/R = 0.3','r/R = 0.5','r/R = 0.7','r/R = 0.9','Location','best');

%% Angle of attack contours
[PsiGrid, RbarGrid] = meshgrid(psi, rbar_vec);
X = - RbarGrid .* cos(PsiGrid);
Y = RbarGrid .* sin(PsiGrid);

alpha_deg = rad2deg(alpha_mat);
levels = floor(min(alpha_deg(:))) : ceil(max(alpha_deg(:)));

figure;
contourf(X, Y, alpha_deg, 20, 'LineColor', 'none');
colorbar;
axis equal;
xlabel('x/R');
ylabel('y/R');
title('Angle of attack contours');

figure;
[C,h] = contour(X, Y, alpha_deg, levels, 'LineWidth', 1.5);
clabel(C, h, 'FontSize', 10, 'Color', 'k');

colormap(jet);
colorbar;

axis equal;
xlabel('x/R');
ylabel('y/R');
title('Angle of attack contour lines');

%% Fourier fit
A0 = mean(beta);
Ac = 2/length(psi) * sum(beta .* cos(psi));
As = 2/length(psi) * sum(beta .* sin(psi));

a0_fit = A0;
a1_fit = -Ac;
b1_fit = -As;

beta_fit = a0_fit - a1_fit*cos(psi) - b1_fit*sin(psi);

%% Fourier fit plot
figure;
plot(psi_deg, rad2deg(beta), 'LineWidth', 1.8); hold on;
plot(psi_deg, rad2deg(beta_fit), '--', 'LineWidth', 1.5);
grid on;
xlabel('\psi [deg]');
ylabel('\beta [deg]');
title('Fourier fit');
legend('\beta(\psi)', 'Fit', 'Location', 'best');

%% Results
fprintf('mu = %.4f\n', mu);
fprintf('a0 = %.4f deg\n', rad2deg(a0_fit));
fprintf('a1 = %.4f deg\n', rad2deg(a1_fit));
fprintf('b1 = %.4f deg\n', rad2deg(b1_fit));



