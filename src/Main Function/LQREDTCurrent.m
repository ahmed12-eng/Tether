function [u] = LQREDTCurrent(xstate,omega,Im)
%% Im is moment of inertia
x=xstate;  %%[theta; thetad; phi ;phid] 
A = [0 1 0 0;
   -3*omega^2 0 0 0;
      0  0  0  1;
      0  0  -7*omega^2 0];
B = [0  0;
    1/Im 0 ;
    0  0;
    0  1/Im];
Q = diag([100 1 100 1]);
R = diag([0.1 0.1]);

k = lqr(A,B,Q,R);  %%gain matrix

u = -k*x;   %%control
end