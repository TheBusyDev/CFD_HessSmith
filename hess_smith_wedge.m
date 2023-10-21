clc; clear; close all

U = 300;    % undisturbed velocity
U_inf = [U; 0];
N = 100;

% generate wedge
% straight line parameters and equation
m = 0.5;
c = 0.5;
line = @(x) m*x + c;

% find circumference center
xc = m*c;

% circumference radius and equation
R = norm ([xc, c]);
circ = @(x) sqrt(R^2 - (x-xc)^2);

start_x = -c/m;
end_x = xc+R;

i = linspace (0, 1, N);
xd = start_x + (end_x-start_x) * sqrt(i);
xd = xd';
yd = zeros (N, 1);

for i = 1:N
    if xd(i) < 0
        yd(i) = line (xd(i));
    else
        yd(i) = circ (xd(i));
    end
end

xv = xd;
yv = -yd;

xe = [xd(1:end-1); flip(xv)];
ye = [yd(1:end-1); flip(yv)];

Np = length (xe); % number of total points

figure(1)
plot (xd, yd, xv, yv, 'LineWidth', 1.5)
grid on; hold on; axis equal

% initialize variables
xp = zeros (1, Np-1);
yp = zeros (1, Np-1);
l = zeros (1, Np-1);
n = zeros (2, Np-1);
t = zeros (2, Np-1);

for j = 1:Np-1
    xp(j) = (xe(j) + xe(j+1)) / 2;  % x-component of control point
    yp(j) = (ye(j) + ye(j+1)) / 2;  % y-component of control point

    l(j) = sqrt ((ye(j+1) - ye(j))^2 + (xe(j+1) - xe(j))^2);    % length of each panel
    n(:,j) = [-(ye(j+1) - ye(j)); (xe(j+1) - xe(j))] / l(j);    % normal vector
    t(:,j) = [n(2,j); -n(1,j)];                                 % tangent vector
end

% initialize variables
A = zeros (Np, Np);
vel_v = zeros (Np-1, 2);

for i = 1:Np-1
    for j = 1:Np-1
        % source/vortex induced velocity
        [vel_s, vel_v(j,:)] = SourceVortex (xp(i), yp(i), xe(j), ye(j), xe(j+1), ye(j+1), l(j), t(:,j));

        A(i,j) = dot (vel_s, n(:,i));
    end

    A(i,Np) = dot (sum(vel_v), n(:,i));
end

% Kutta condition
% initialize variables
vel_v_1 = zeros (Np-1, 2);
vel_v_Np = zeros (Np-1, 2);

for j = 1:Np-1
    % source/vortex induced velocity
    [vel_s_1, vel_v_1(j,:)] = SourceVortex (xp(1), yp(1), xe(j), ye(j), xe(j+1), ye(j+1), l(j), t(:,j));
    [vel_s_Np, vel_v_Np(j,:)] = SourceVortex (xp(end), yp(end), xe(j), ye(j), xe(j+1), ye(j+1), l(j), t(:,j));
    
    A(Np,j) = dot (vel_s_1, t(:,1)) + dot (vel_s_Np, t(:,end));
end

A(Np,Np) = dot (sum(vel_v_1), t(:,1)) + dot (sum(vel_v_Np), t(:,end));

% right-hand side
% initialize variables
b = zeros (Np, 1);

for i = 1:Np-1
    b(i) = -dot (U_inf, n(:,i));
end

b(Np) = -dot (U_inf, t(:,1) + t(:,end));

% solve linear system: A*sol = b
sol = A \ b;
q = sol(1:end-1);
gamma = sol(end);

% compute velocity and pressure coefficient on each panel
% initialize variables
u = zeros (2, Np-1);
cp = zeros (1, Np-1);
vel_v = [];
U_inf_norm = norm (U_inf);

for i = 1:Np-1
    u(:,i) = U_inf;

    for j = 1:Np-1
        % source/vortex induced velocity
        [vel_s, vel_v] = SourceVortex (xp(i), yp(i), xe(j), ye(j), xe(j+1), ye(j+1), l(j), t(:,j));

        u(:,i) = u(:,i) + vel_s*q(j) + vel_v*gamma;
    end

    % pressure coefficient
    cp(i) = 1 - (dot (u(:,i), t(:,i)) / U_inf_norm)^2;
end

% compute lift coefficient: cl = -sum (cp_i * l_i* cos(beta_i) ), where cos(beta_i) = tangent(1)
cl = -sum (cp.*l.*t(1,:));

% plot pressure coefficient
figure(2)
plot (xp, cp, 'b', 'LineWidth', 1.5)
hold on; grid on
xlabel('x'); ylabel('c_p')

% compute velocity field outside the wedge
Nx = 300; Ny = 300; % number of discretisation points

x = linspace (-4, +4, Nx);
y = linspace (-1.5, +1.5, Ny);

ydplus = interp1 (xd, yd, x);
yvplus = interp1 (xv, yv, x);

% initialize variables
U_field = zeros (Nx, Ny);
V_field = zeros (Nx, Ny);
min_xe = min(xe);
max_xe = max(xe);

for i = 1:Nx
    for j = 1:Ny
        % if the point is located outside the wedge, compute velocity field
        if ( x(i) < min_xe || x(i) > max_xe || y(j) > ydplus(i) || y(j) < yvplus(i) )
            U_field(i,j) = U_inf(1); 
            V_field(i,j) = U_inf(2);

            for k = 1:Np-1
                % source/vortex induced velocity
                [vel_s, vel_v] = SourceVortex (x(i), y(j), xe(k), ye(k), xe(k+1), ye(k+1), l(k), t(:,k));
        
                U_field(i,j) = U_field(i,j) + vel_s(1)*q(k) + vel_v(1)*gamma;
                V_field(i,j) = V_field(i,j) + vel_s(2)*q(k) + vel_v(2)*gamma;
            end
        end
    end
end

% plot velocity field outside the wedge
[X, Y] = meshgrid (x, y);
figure(3); hold on; axis equal; title ('U-component of velocity field')
contourf (X', Y', U_field, 100, 'LineStyle', 'None')
plot (xe, ye, 'Linewidth', 2, 'Color', 'k')

figure(4); hold on; axis equal; title ('V-component of velocity field')
contourf (X', Y', V_field, 100, 'LineStyle', 'None')
plot (xe, ye, 'Linewidth', 2, 'Color', 'k')





function [vel_s, vel_v] = SourceVortex (xc, yc, xe, ye, xe1, ye1, l, t)

    csi = (xc - xe) * (xe1 - xe) + (yc - ye) * (ye1 - ye);
    csi = csi / l;

    eta = -(xc - xe) * (ye1 - ye) + (yc - ye) * (xe1 - xe);
    eta = eta / l;

    r1 = sqrt(csi^2 + eta^2);
    r2 = sqrt((csi - l)^2 + eta^2);

    theta1 = atan2(eta, csi);
    theta2 = atan2(eta, csi - l);

    % solve numerical conflicts that can possibly lead to errors
    if (abs(theta1) < 1e-12 && abs(theta2) > 3)
        theta1 = 0;
        theta2 = pi;
    elseif (abs(theta2) < 1e-12 && abs(theta1) > 3)
        theta1 = pi;
        theta2 = 0;
    end

    u = -log(r2/r1) / (2*pi);
    v = (theta2 - theta1) / (2*pi);

    % return to global reference system
    % Q = [
    %     cos(beta), -sin(beta)
    %     sin(beta), cos(beta)
    % ];
    %
    % Q = [
    %     t(1), -t(2)
    %     t(2), t(1)
    % ];
    %
    % source induced velocity
    % vel_s = Q*[u v]';

    vel_s = [
        t(1)*u - t(2)*v
        t(2)*u + t(1)*v
    ];

    % vortex induced velocity
    vel_v(1,1) = vel_s(2);
    vel_v(2,1) = -vel_s(1);

end
