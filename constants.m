%% Conversion Constants %%
lbs_to_kg = 0.45359237;
ft_to_m = 0.3048;
slugft2_to_kgm2 = 1.35581795;
Hp_to_W = 0.735499 * 1e3;
kts_to_ms = 0.5144444444;
% atmospheric / flight assumptions
rho = 1.225;         % Air density at sea level [kg/m^3]
g   = 9.80665;          % Gravity [m/s^2]

%% Helicopter constants %%
MTOW = 10300 * lbs_to_kg; % kg
Mempty = 5546 * lbs_to_kg; % kg
Mfuel = 1876  * lbs_to_kg; % kg

W = MTOW * g;

N_engine = 2;
P_TO = 1300 * Hp_to_W; % W
P_const = 1300 * Hp_to_W; % W

l_LOA = 16; % m (total length)
A_eq = 12.6 * ft_to_m^2; % Equivalent plate area Found in https://doi.org/10.2514/6.2024-1117

% Rotor Parameters
hinge_offset = 0.038;

N_blade_main = 4;
solidity_main = 0.0747;
R_main = 22 * ft_to_m; % m
A_main = pi * R_main^2; % m^2
c_main = 1.29 * ft_to_m; % m
tip_speed_main = 675 * ft_to_m; % m/s
lin_twist_main = deg2rad(-10); % rad
J_main = 1890 * slugft2_to_kgm2; % kg m^2

N_blade_tail = 4;
solidity_tail = 0.1719;
R_tail = 4 * ft_to_m; % m
A_tail = pi * R_tail^2; % m^2
c_tail = 0.54 * ft_to_m; % m
tip_speed_tail = 674 * ft_to_m; % m/s
lin_twist_tail = deg2rad(-8); % rad
J_tail = 4.04 * slugft2_to_kgm2; % kg m^2