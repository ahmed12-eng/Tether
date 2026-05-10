function dx=dragforce(x,dpar)
%T = dpar.T;
%g0 = dpar.g0  ;
%M = dpar.M ; %Kg/mol
%R = dpar.R ;   %j/(mol*K);
Re = dpar.Re;  %earth radius;
Cd = dpar.Cd ;
A= dpar.A  ;


%mue=3.986e14;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
omega=7.29211e-5;
omegavec = [0;0;omega];
r = x(1:3);
v = x(4:6);
v_rel = v-cross(omegavec,r);

vrel_norm=norm(v_rel);
r_norm = norm(r);
h = r_norm-Re;

if h < 200e3
    rho = 2e-10;
elseif h < 300e3
    rho = 2e-11;
elseif h < 400e3
    rho = 1e-12;
else
    rho = 1e-13;
end

Fd = - (1/(2)) *rho*Cd*A*vrel_norm.*v_rel;
%a = (-mue/r_norm^3).*r - (1/(2*m)) *rho*Cd*A*vrel_norm.*v_rel;

dx=[v;Fd];
end