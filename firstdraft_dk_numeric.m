clear; clc;

%% Constants
deg = pi/180;

alpha0 = pi/2;
alpha1 = pi/2;
alpha2 = pi/2;
alpha3 = pi/2;
alpha4 = pi/2;

beta1  = 75*deg;
beta2  = 90*deg;


%% Input joint angles
theta1 = 40*deg
theta2 = 40*deg

%% Rotation matrices
RotDH = @(th, al) [ cos(th), -sin(th)*cos(al),  sin(th)*sin(al);
                    sin(th),  cos(th)*cos(al), -cos(th)*sin(al);
                    0,        sin(al),          cos(al) ];

Q0 = [ cos(alpha0), 0, sin(alpha0);
       0,           1, 0;
      -sin(alpha0), 0, cos(alpha0)];

Q1 = RotDH(theta1, alpha1);
Q2 = RotDH(theta2, alpha2);

%% Local vectors
e3_3 = [0;0;1];
e4_4 = [0;0;1];
u5   = [0;0;1];

v5 = [sin(beta1)*cos(beta2);
      cos(beta1);
      sin(beta1)*sin(beta2)];

%% Compute e3_1 and e4_1
e3_1 = Q1 * e3_3;
e4_1 = Q0 * Q2 * e4_4;

disp('e3_1 = ');
disp(e3_1);

disp('e4_1 = ');
disp(e4_1);

%% Solve u1 numerically from geometric constraints
syms ux uy uz real

eq1 = dot(e3_1, [ux;uy;uz]) == cos(alpha3);
eq2 = dot(e4_1, [ux;uy;uz]) == cos(alpha4);
eq3 = ux^2 + uy^2 + uz^2 == 1;

sol_u = vpasolve([eq1, eq2, eq3], [ux, uy, uz]);

u1_solutions = [double(sol_u.ux), double(sol_u.uy), double(sol_u.uz)];

disp('u1 solutions:');
disp(u1_solutions);

%% Direct theta4 extraction for alpha0=alpha2=alpha4=pi/2
% u1 = [-cos(theta4);
%        sin(theta2)*sin(theta4);
%       -cos(theta2)*sin(theta4)]

theta4_solutions = [];

for i = 1:size(u1_solutions,1)
    ux_val = u1_solutions(i,1);
    uy_val = u1_solutions(i,2);
    uz_val = u1_solutions(i,3);

    c4 = -ux_val;
    s4 = uy_val / sin(theta2);

    % consistency check with uz
    uz_check = -cos(theta2)*s4;

    if abs(uz_check - uz_val) < 1e-6
        theta4_val = atan2(s4, c4);
        theta4_solutions = [theta4_solutions; theta4_val];
    else
        fprintf('Branch %d failed uz check\n', i);
    end
end

disp('theta4 solutions (deg):');
disp(theta4_solutions/deg);

%% Compute v1 for each theta4
for i = 1:length(theta4_solutions)
    theta4_val = theta4_solutions(i);
    Q4 = RotDH(theta4_val, alpha4);

    v1 = Q0 * Q2 * Q4 * v5;

    vx = v1(1);
    vy = v1(2);
    vz = v1(3);

    % Your actual convention:
    % v1 = [cos(phi1)*cos(phi2);
    %       sin(phi2);
    %      -sin(phi1)*cos(phi2)]
    phi2 = atan2(vy, sqrt(vx^2 + vz^2));
    phi1 = atan2(-vz, vx);

    fprintf('\n--- Solution %d ---\n', i);
    fprintf('theta4 = %.6f deg\n', theta4_val/deg);
    fprintf('v1 = [%.6f %.6f %.6f]\n', vx, vy, vz);
    fprintf('phi1 = %.6f deg\n', phi1/deg);
    fprintf('phi2 = %.6f deg\n', phi2/deg);
end
