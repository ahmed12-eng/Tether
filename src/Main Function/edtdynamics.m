function [dx,F_L,F_D,I,V_emf,tau_LVLH] = edtdynamics(~, x, p, L)

%% ================= STATES =================
R   = x(1:3);
V   = x(4:6);
th  = x(7);     % theta (in-plane)
thd = x(8);

al  = x(9);     % phi (out-of-plane)
ald = x(10);

%% ================= CONSTANTS =================
rnorm = norm(R);
rhat  = R/rnorm;

m1 = p.m1;
m2 = p.m2;
M  = m1 + m2;
mr = (m1*m2)/M;

Im = mr * L^2;

%% ================= LVLH FRAME =================
z_orb = rhat;

y_orb = cross(R,V);
y_orb = y_orb / norm(y_orb);

x_orb = cross(y_orb,z_orb);

%% ================= TETHER DIRECTION (REFERENCE CONSISTENT) =================
e_local = [ ...
    cos(al)*sin(th);   % x_orb
    sin(al);           % y_orb
    cos(al)*cos(th)    % z_orb
];

e_tether = e_local(1)*x_orb + e_local(2)*y_orb + e_local(3)*z_orb;
e_tether = e_tether / norm(e_tether);

L_vec = L * e_tether;

%% ================= BODY POSITIONS =================
r1 = R - (m2/M)*L_vec;
r2 = R + (m1/M)*L_vec;

%% ================= GRAVITY =================
a_g1 = -p.mu * r1 / norm(r1)^3;
a_g2 = -p.mu * r2 / norm(r2)^3;

%% ================= MAGNETIC FIELD =================
theta_e = acos(R(1)/rnorm);
psie    = atan2(R(2),R(1));

lat  = 90 - theta_e*(180/pi);
long = psie*(180/pi);
alt  = rnorm/1000;

[Bx,By,Bz] = igrf('01.Jan.2020',lat,long,alt,'geocentric');
B = [Bx By -Bz]'*1e-9;

%% ================= EMF =================
omegaE = [0;0;7.2921159e-5];
Vplasma = cross(omegaE,R);

V_rel = V - Vplasma;
E     = cross(V_rel,B);

V_emf = dot(E,e_tether)*L;

%% ================= CURRENT (OML) =================
A = 0.01;
I = edt_currentOML_eV(rnorm - p.Re, A, V_emf);

%% ================= LORENTZ FORCE =================
F_L = I * cross(L_vec,B);

%% ================= DRAG =================
dpar.Re = 6378e3;
dpar.Cd = 2.2;
dpar.A  = 0.05;

FDX = dragforce(x,dpar);
F_D = FDX(4:6);

%% ================= ORBITAL RATE =================
h     = norm(cross(R,V));
omega = h / rnorm^2;



%% ================= COM ACCELERATION =================
a_COM = (a_g1*m1 + a_g2*m2)/M + F_L/M + F_D/M;

%% ================= TORQUE FROM LORENTZ =================
tau = cross(r1-R, F_L*m1/M) + cross(r2-R, F_L*m2/M);

Tmat = [x_orb y_orb z_orb];
tau_LVLH = Tmat' * tau;

tau_theta = tau_LVLH(2);   % حول y_orb
tau_alpha = tau_LVLH(1);   % حول x_orb

%% ================= CONTROL =================
xstate = [th; thd; al; ald];
u = LQREDTCurrent(xstate, omega, Im);

tau_theta = tau_theta + u(1);
tau_alpha = tau_alpha + u(2);

%% ================= ROTATIONAL DYNAMICS =================
thdd = -3*omega^2*sin(th)*cos(th) ...
       + 2*ald*thd*tan(al) ...
       + tau_theta/(Im*(cos(al)^2));

aldd = -omega^2*sin(al)*cos(al)*(4 + 3*cos(th)^2) ...
       + tau_alpha/Im;

%% ================= STATE VECTOR =================
dx = [ V;
       a_COM;
       thd;
       thdd;
       ald;
       aldd ];

end