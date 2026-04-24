clear; close all; clc;

run('constants.m');
rho = 1.225;
a_l = 5.73;

%% Inputs
V = 20;
q = deg2rad(20);
p = deg2rad(10);

theta0  = deg2rad(6);
theta1c = deg2rad(2);
theta1s = deg2rad(1);

vi = 6.3;

%% Rotor parameters
R = R_main;
c = c_main;
Omega = tip_speed_main / R_main;

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
lambda = lambda_c + lambda_i;

Lock = rho * c * a_l * R^4 / I_blade;

a0 = (Lock/2) * ((1+mu^2)/4 * theta0 + lambda/3 - ...
        (1/3)*mu*theta1s + (mu/6)* (p/Omega));
a1 = 1/(1-mu^2/2) * (4*mu * ((2/3) * theta0 + lambda/2) - ...
        (1 + (3/2)*mu^2)*theta1s - (p/Omega) - (16/Lock)*(q/Omega));
b1 = 1/(1+mu^2/2) * ((4/3)*mu*a0 + (1+mu^2/2)*theta1c + ...
        (q/Omega) - (16/Lock)*(p/Omega));



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
rbar_vec = linspace(hinge_offset, 0.95, 30);
alpha_mat = zeros(length(rbar_vec), length(psi));

for ir = 1:length(rbar_vec)
    r = rbar_vec(ir) * R;
    UT = Omega*r + V*sin(psi);
    UP = vi - r*beta_dot - V*beta.*cos(psi) + q*r*cos(psi) + p*r*sin(psi);
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
X = RbarGrid .* cos(PsiGrid);
Y = RbarGrid .* sin(PsiGrid);

alpha_deg = rad2deg(alpha_mat);

% Clip values below -10 so they are not color coded (set to -10 for plotting)
alpha_plot = max(alpha_deg, -10);

levels = -10 : ceil(max(alpha_plot(:)));

figure;
contourf(X, Y, alpha_plot, 20, 'LineColor', 'none');
c = colorbar;
clim([-10, max(alpha_plot(:))]); % ensure colorbar starts at -10
axis equal;
xlabel('x/R');
ylabel('y/R');
title('Angle of attack contours (values < -10 clipped)');

figure;
[C,h] = contour(X, Y, alpha_plot, levels, 'LineWidth', 1.5);
clabel(C, h, 'FontSize', 10, 'Color', 'k');

colormap(jet);
colorbar;
clim([-10, max(alpha_plot(:))]);

axis equal;
xlabel('x/R');
ylabel('y/R');
title('Angle of attack contour lines (values < -10 clipped)');

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

% %% Results
% fprintf('mu = %.4f\n', mu);
% fprintf('a0 = %.4f deg\n', rad2deg(a0_fit));
% fprintf('a1 = %.4f deg\n', rad2deg(a1_fit));
% fprintf('b1 = %.4f deg\n', rad2deg(b1_fit));
% 
% 
% 
% % Blade pitch angle contour vs radial position
% rbar_vec_full = linspace(0.2, 0.99, 80);      % finer radial grid for plotting
% [PsiGrid_full, RbarGrid_full] = meshgrid(psi, rbar_vec_full);
% 
% % compute theta at each psi and rbar (theta is independent of r here but depends on psi)
% Theta_mat = theta(:).' .* ones(length(rbar_vec_full), length(psi)); % size: nr x nPsi
% 
% % convert to degrees for plotting
% Theta_deg = rad2deg(Theta_mat);
% 
% % map to Cartesian for contour in rotor disk coordinates
% X_theta = RbarGrid_full .* cos(PsiGrid_full);
% Y_theta = RbarGrid_full .* sin(PsiGrid_full);
% 
% figure;
% contourf(X_theta, Y_theta, Theta_deg, 20, 'LineColor', 'none');
% colorbar;
% axis equal;
% xlabel('x/R');
% ylabel('y/R');
% title('Blade pitch angle \theta [deg]');
% 
% % compute inflow angle phi over the same grid used for angles of attack
% % reuse UT and UP expressions but evaluate over rbar_vec_full and psi
% [PsiPhiGrid, RbarPhiGrid] = meshgrid(psi, rbar_vec_full);
% r_mat = RbarPhiGrid * R;
% 
% % Local velocities
% UT_mat = Omega .* r_mat + V .* sin(PsiPhiGrid);
% % beta_dot needs beta derivative at each psi; interpolate dbeta_dpsi (same length as psi)
% dbeta_dpsi_interp = interp1(psi, dbeta_dpsi, PsiPhiGrid(1,:), 'linear', 'extrap');
% % repeat for each radial row
% dbeta_dpsi_mat = repmat(dbeta_dpsi_interp, size(RbarPhiGrid,1), 1);
% beta_dot_mat = dbeta_dpsi_mat * Omega;
% 
% UP_mat = vi - r_mat .* beta_dot_mat + q .* r_mat .* cos(PsiPhiGrid) + p .* r_mat .* sin(PsiPhiGrid);
% 
% % inflow angle phi
% phi_mat = atan2(UP_mat, UT_mat); % radians
% phi_deg_mat = rad2deg(phi_mat);
% 
% % map to Cartesian (using Rbar grid already in RbarPhiGrid)
% X_phi = RbarPhiGrid .* cos(PsiPhiGrid);
% Y_phi = RbarPhiGrid .* sin(PsiPhiGrid);
% 
% % contour plot
% figure;
% contourf(X_phi, Y_phi, phi_deg_mat, 20, 'LineColor', 'none');
% colorbar;
% axis equal;
% xlabel('x/R');
% ylabel('y/R');
% title('Inflow angle \phi [deg]');
