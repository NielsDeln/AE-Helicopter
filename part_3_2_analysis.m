close all; clc;
% Load the simulation data
load('Sim_run.mat');   % loads variable called 'data' (a Simulink.SimulationData.Dataset)

% Extract each signal's time and values
theta_0_sig = data{1};   % returns a timeseries or timetable
theta_c_sig = data{4};
V_sig       = data{22};
Vr_sig      = data{23};
Vz_sig      = data{11};
h_sig       = data{12};
pos_sig     = data{13};


%% Plot control inputs

% Extract signals and squeeze data
theta_0_time = theta_0_sig.Values.Time;
theta_0_data = rad2deg(squeeze(theta_0_sig.Values.Data));

theta_c_time = theta_c_sig.Values.Time;
theta_c_data = rad2deg(squeeze(theta_c_sig.Values.Data));

% Crop to first 120 seconds
mask_0 = theta_0_time <= 100;
mask_c = theta_c_time <= 100;

theta_0_time = theta_0_time(mask_0);
theta_0_data = theta_0_data(mask_0);
theta_c_time = theta_c_time(mask_c);
theta_c_data = theta_c_data(mask_c);

% Compute y-limits equal to saturation limits
t0_ylim = [4, 22.7];
tc_ylim = [-14.75, 19.25];

% Plot
figure;

ax1 = subplot(2,1,1);
plot(theta_0_time, theta_0_data, 'b');
ylabel('$\theta_0$ [deg]', 'Interpreter', 'latex', 'FontSize', 12);
ylim(t0_ylim);
grid on;

ax2 = subplot(2,1,2);
plot(theta_c_time, theta_c_data, 'r');
ylabel('$\theta_c$ [deg]', 'Interpreter', 'latex', 'FontSize', 12);
ylim(tc_ylim);
grid on;
xlabel('Time (s)', 'Interpreter', 'latex', 'FontSize', 14);

% Link x-axes so zoom/pan stays in sync
linkaxes([ax1, ax2], 'x');


%% Plot position, height, and velocity during hover
% Extract position and height signals
pos_time = pos_sig.Values.Time;
pos_data = squeeze(pos_sig.Values.Data);
h_time = h_sig.Values.Time;
h_data = squeeze(h_sig.Values.Data);
V_time = V_sig.Values.Time;
V_data = V_sig.Values.Data;
Vr_time = Vr_sig.Values.Time;
Vr_data = squeeze(Vr_sig.Values.Data);
Vz_time = Vz_sig.Values.Time;
Vz_data = squeeze(Vz_sig.Values.Data);


% Position: plot from 30 s to end, y-limits 500 to 550 m
pos_mask = pos_time >= 30;
pos_time_plot = pos_time(pos_mask);
pos_data_plot = pos_data(pos_mask);

figure;
plot(pos_time_plot, pos_data_plot, 'b');
xlabel('Time [s]', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Position [m]', 'Interpreter', 'latex', 'FontSize', 14);
ylim([500 540]);
xlim([min(pos_time_plot) max(pos_time_plot)]);
grid on;

% Height: first 100 seconds
h_mask = h_time <= 100;

figure;
plot(h_time(h_mask), h_data(h_mask)*1e3, 'g');
xlabel('Time [s]', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Height [mm]', 'Interpreter', 'latex', 'FontSize', 14);
xlim([0 100]);
grid on;

% Velocity: first 100 seconds
vel_mask = V_time <= 100;

figure;
plot(V_time(vel_mask), V_data(vel_mask), 'm');
xlabel('Time [s]', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Velocity [m/s]', 'Interpreter', 'latex', 'FontSize', 14);
xlim([0 100]);
grid on;

% Velocity reference
figure;
plot(Vr_time, Vr_data, 'm', 'Color', 'blue'); hold on;
plot(Vz_time, Vz_data, 'm', 'Color', 'red');
xlabel('Time [s]', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Reference signal [m/s]', 'Interpreter', 'latex', 'FontSize', 14);
legend('Total Velocity', 'Vertical Velocity')
ylim([-1, 25])
grid on;