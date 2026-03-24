%% Conversion Constants %%
lbs_to_kg = 0.45359237;
ft_to_m = 0.3048;
slugft2_to_kgm2 = 1.35581795;
Hp_to_KW = 0.735499;

%% Helicopter constants %%
MTOW = 10300 * lbs_to_kg; % kg
Mempty = 5546 * lbs_to_kg; % kg
Mfuel = 1876  * lbs_to_kg; % kg

N_engine = 2;
P_TO = 1300 * Hp_to_KW; % KW
P_const = 1300 * Hp_to_KW; % Kw

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