function I = edt_currentOML(h,A,V_emf)

%% h is the height of EDT above earth's surface in meter
%% A is the m^2 cross-section area of tether
%% V_emf tether voltage 

%% Constants
  k = 1.380649e-23;   % Boltzmann constant [J/K]
  e = 1.602e-19;      % Electron charge [C]
  m_e = 9.11e-31;     % Electron mass [kg]
  T_e = 2000;         % K, electron temperature
  V_P =0;             % Plasma potential
%% get the value of n_e electron density at a given height h

if h>150e3    %height in m
 
  %% Chapman Model Parameters  within range 900km to 150km
  h_max = 300e3;      %height at n_max 
  H = 50e3;           % Scale height [m]
  n_max = 1e11;       % Max electron density [electrons/m^3]

  %% Compute electron density (Chapman)
  z = (h - h_max)/H;
  n_e = n_max * exp(0.5 * (1 - z - exp(-z)));


elseif h<=150e3
%% Exponential D-layer parameters if h<150 km D layer of ionosphere
n0 = 1e10;      % electrons/m^3 at h_ref
h_ref = 100e3;  % Refrence height in m
H = 8e3;        % Scale height in m

%% Compute electron density
n_e = n0 * exp(-(h - h_ref)/H);
end
v_th = sqrt(2*k*T_e/m_e);      % thermal speed of electrons
I = e * n_e * A * v_th * (1 + e*V_emf/(k*T_e));

end
