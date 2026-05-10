function [dx,F_D,tau_LVLH] = Drag_Deorbit(~, x, p,L)
%%%Drag onlyyy
% Extract states
R  = x(1:3);
V  = x(4:6);
th  = x(7);
thd = x(8);

al  = x(9);
ald = x(10);
rnorm = norm(R);
rhat = R/rnorm;
mr = p.m1*p.m2/(p.m1+p.m2);
Im = mr*L^2;          %%moment of inertia
M = p.m1 + p.m2;
%%LVLH Frame
% Build local orbital frame
z_orb = rhat;
y_orb = cross(R,V);
y_orb = y_orb/norm(y_orb);
x_orb = cross(y_orb,z_orb);

% Rotation from angles
e_tether = cos(al)*cos(th)*z_orb + ...
           cos(al)*sin(th)*x_orb + ...
           sin(al)*y_orb;

L_vec = L * e_tether;
% Recover body positions
r1 = R - (p.m2/M)*L_vec;
r2 = R + (p.m1/M)*L_vec;

% Gravity
a_g1 = -p.mu * r1 / norm(r1)^3;
a_g2 = -p.mu * r2 / norm(r2)^3;

% % Magnetic field at COM (approximation)
% %%Get igrf parameters
% theta_e = acos(R(1)/rnorm);
% psie = atan2(R(2),R(1));
% lat = 90- theta_e*(180/pi);
% long = psie*(180/pi);
% alt = (rnorm )/1000;
% 
% [Bx,By,Bz] = igrf('01.Jan.2020',lat,long,alt,'geocentric');
% B = [Bx By -Bz]'*1e-9;
% 
% %% EMF & Current
% omegaE = [0;0;7.2921159e-5];
% Vplasma = cross(omegaE,R);
% V_rel = V- Vplasma ;                   %% To add Earth's rotation  
% E = cross(V_rel,B);
% V_emf = dot(E,e_tether)*L;
% A = 0.01;              %%  m^2, cross-section area of tether
% 
% %Calculating the current with OML Theory
% I = edt_currentOML_eV(rnorm-p.Re,A,V_emf);
% 
% % Lorentz force
% F_L = I * cross(L_vec,B);

%Calculating Drag Force
dpar.Re=6378e3;  %earth radius;
dpar.Cd = 2.2;
dpar.A = 0.05;
format long e
FDX = dragforce(x,dpar);
F_D = FDX(4:6);
%%angular velocity calculation
% omega = sqrt(p.mu / rnorm^3); %%for circular orbit
hh=norm(cross(R,V));
omega = hh/rnorm^2;

% Coupled nonlinear dynamics
%% Rotational equations of motion
%% -2*(ldot/L)*thd  لاكن الطول ثابت فالتغير في الطول =0 
thdd = -3*omega^2*sin(th)*cos(th)...
+ omega^2*sin(th)*(sin(al))^2; 
       
aldd = -omega^2*sin(al)*cos(al) ...
       + omega^2*sin(al)*cos(al)*(sin(th))^2 ;
      
% COM acceleration due to gravity
a_COM = (a_g1*p.m1 + a_g2*p.m2)/M +F_D/M;

dx = [ V; a_COM;thd;thdd;ald;aldd ];
end