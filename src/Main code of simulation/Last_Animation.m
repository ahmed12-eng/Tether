clear; clc; close all;

%% ================= PARAMETERS =================
L = 200;   
p.m1 = 5;   
p.m2 = 3;  
p.mu = 3.986e14;
p.Re = 6371e3;
p.R_tether = 30;       
p.R_plasma = 10000;    

%% ================= INITIAL CONDITIONS =================
h0 = 400e3; 

R0 = [p.Re+h0; 0; 0]; 
V0 = [0; sqrt(p.mu/norm(R0)); 0]; 

theta0 = deg2rad(5); 
thetadot0 = deg2rad(0.1); 
alpha0 = deg2rad(2); 
alphadot0 = deg2rad(0); 

X0 = [R0; V0; theta0; thetadot0; alpha0; alphadot0];

tspan = [0 1000];

options = odeset('RelTol',1e-8,'AbsTol',1e-10);
[t,x] = ode45(@(t,x) edtdynamics(t,x,p,L), tspan, X0, options);

R = x(:,1:3);
V = x(:,4:6);
rnorm = vecnorm(R,2,2);

r1 = zeros(length(t),3);
r2 = zeros(length(t),3);

M = p.m1 + p.m2;
%% Liberation angles ploting
title(' \theta , \alpha vs time')
hold on 
plot(t/(90*60),x(:,7)*180/pi,'r')
hold on
plot(t/(90*60),x(:,9)*180/pi,'g')
grid on
legend('\theta','\alpha')
xlabel('Time in orbits')
ylabel('Angle(degree)')
%%
figure
title('\theta dot and \alpha dot vs time')
plot(t/(90*60),x(:,8)*180/pi,'r')
hold on
plot(t/(90*60),x(:,10)*180/pi,'g')
grid on 
legend('\theta dot','\alpha dot')
xlabel('Time in orbits')
ylabel('Anglular velocities (deg/sec)')


% %%TORQUE
% 
% Torque = zeros(length(t),1);
% 
% for i=1:length(t)
%     [~,~,~,~,~,tau] = edtdynamics(t(i),x(i,:)',p,L);
%     Torque(i) = norm(tau);
% end
% 
% figure
% plot(t, Torque)
% title('Torque vs Time')
% ylabel('Torque (N)')
% xlabel('Time (s)')
% grid on



%% ================= RECONSTRUCT TETHER =================
for j = 1:length(t)

    % ===== LVLH =====
    z_orb = R(j,:)/rnorm(j);
    y_orb = cross(R(j,:),V(j,:)); 
    y_orb = y_orb/norm(y_orb);
    x_orb = cross(y_orb,z_orb);

    % ===== angles =====
    th = x(j,7);
    al = x(j,9);

    % ===== CORRECT e_tether (MATCH DYNAMICS) =====
    e_local = [ ...
        cos(al)*sin(th);   % x_orb
        sin(al);           % y_orb
        cos(al)*cos(th)    % z_orb
    ];

    e_tether = e_local(1)*x_orb + ...
               e_local(2)*y_orb + ...
               e_local(3)*z_orb;

    e_tether = e_tether / norm(e_tether); % normalize

    % ===== tether vector =====
    L_vec = L * e_tether;

    % ===== satellites positions =====
    r1(j,:) = R(j,:) - (p.m2/M)*L_vec;
    r2(j,:) = R(j,:) + (p.m1/M)*L_vec;
end

%% ================= LOAD STL =================
TR = stlread('CubeSat 1U Ukr v3 frame.stl');

F = double(TR.ConnectivityList);
Vsat = double(TR.Points);

Vsat = Vsat - mean(Vsat,1);
Vsat = Vsat * 2.5;

%% ================= FIGURE =================
figure
tiledlayout(1,2)

%% ================= ORBIT VIEW =================
ax1 = nexttile;
hold(ax1,'on')
set(ax1,'Color','k')

[X,Y,Z] = sphere(100);

earth = surf(ax1,...
    p.Re*X,p.Re*Y,-p.Re*Z,...
    'FaceColor','texturemap',...
    'CData',imread('earth.jpg'),...
    'EdgeColor','none');

numStars = 2000;
starR = p.Re * 6;

xs = (rand(numStars,1)-0.5)*2*starR;
ys = (rand(numStars,1)-0.5)*2*starR;
zs = (rand(numStars,1)-0.5)*2*starR;

star_phase = 2*pi*rand(numStars,1);
star_freq  = 0.5 + 2*rand(numStars,1);

hStars = scatter3(ax1,xs,ys,zs,2,[1 1 1],'filled');

COM = scatter3(ax1,R(1,1),R(1,2),R(1,3),'r','filled');
traj = plot3(ax1,R(1,1),R(1,2),R(1,3),'w');

axis(ax1,'equal')
view(ax1,3)
camproj(ax1,'perspective')

%% ================= TETHER VIEW =================
ax2 = nexttile;
hold(ax2,'on')
set(ax2,'Color','k')

sat1 = patch(ax2,'Faces',F,'Vertices',Vsat,...
    'FaceColor','red','EdgeColor','none');

sat2 = patch(ax2,'Faces',F,'Vertices',Vsat,...
    'FaceColor','blue','EdgeColor','none');

tether_line = plot3(ax2,...
    [r1(1,1) r2(1,1)],...
    [r1(1,2) r2(1,2)],...
    [r1(1,3) r2(1,3)],...
    'b','LineWidth',1.5);

axis(ax2,'equal')
view(ax2,3)
camproj(ax2,'perspective')

zoom_range = 250;

%% ================= LIGHTING =================
lighting phong

light('Position',[1 0.5 1]*1e7,'Style','infinite');
light('Position',[-1 -0.3 0.5]*1e7,'Style','infinite');

set(sat1,'FaceLighting','phong','SpecularStrength',0.9);
set(sat2,'FaceLighting','phong','SpecularStrength',0.9);

%% ================= LVLH =================
scale = 150;
orb_x = quiver3(ax2,0,0,0,0,0,0,'r','LineWidth',2);
orb_y = quiver3(ax2,0,0,0,0,0,0,'g','LineWidth',2);
orb_z = quiver3(ax2,0,0,0,0,0,0,'b','LineWidth',2);

%% ================= ANIMATION =================
for i = 1:length(t)

    % LVLH
    z_orb = R(i,:)/norm(R(i,:));
    y_orb = cross(R(i,:),V(i,:)); y_orb = y_orb/norm(y_orb);
    x_orb = cross(y_orb,z_orb);

    R_mat = [x_orb(:), y_orb(:), z_orb(:)];

    % satellites
    V1 = (R_mat * Vsat')' + r1(i,:);
    V2 = (R_mat * Vsat')' + r2(i,:);

    set(sat1,'Vertices',V1);
    set(sat2,'Vertices',V2);

    % tether
    set(tether_line,'XData',[r1(i,1) r2(i,1)],...
                    'YData',[r1(i,2) r2(i,2)],...
                    'ZData',[r1(i,3) r2(i,3)]);

    % orbit
    set(COM,'XData',R(i,1),'YData',R(i,2),'ZData',R(i,3));
    set(traj,'XData',R(1:i,1),'YData',R(1:i,2),'ZData',R(1:i,3));

    % LVLH axes
    set(orb_x,'XData',R(i,1),'YData',R(i,2),'ZData',R(i,3),...
        'UData',scale*x_orb(1),'VData',scale*x_orb(2),'WData',scale*x_orb(3));

    set(orb_y,'XData',R(i,1),'YData',R(i,2),'ZData',R(i,3),...
        'UData',scale*y_orb(1),'VData',scale*y_orb(2),'WData',scale*y_orb(3));

    set(orb_z,'XData',R(i,1),'YData',R(i,2),'ZData',R(i,3),...
        'UData',scale*z_orb(1),'VData',scale*z_orb(2),'WData',scale*z_orb(3));

    xlim(ax2,[R(i,1)-zoom_range R(i,1)+zoom_range]);
    ylim(ax2,[R(i,2)-zoom_range R(i,2)+zoom_range]);
    zlim(ax2,[R(i,3)-zoom_range R(i,3)+zoom_range]);

    drawnow
end