clc; clear;

%% =========================================================
%% INPUTS
%% =========================================================
deg = pi/180;

alpha0 = pi/2;
alpha1 = pi/2;
alpha2 = pi/2;
alpha3 = pi/2;
alpha4 = pi/2;

beta1  = 75*deg;
beta2  = 90*deg;

% phi from the dk solved
phi1 = 59.837262 * deg
phi2 = 6.599236 * deg
% =========================================================
%% THETA2 (ACTIVE JOINT)
%% From:
%% e4_5.v5 = cos(gamma)
%% e4_1.v1 = cos(gamma)
%% -> A2*sin(theta2) + B2*cos(theta2) = D2
%% =========================================================

cg5 = sin(alpha4)*cos(beta1) + cos(alpha4)*sin(beta1)*sin(beta2);

A2 = sin(alpha2)*cos(phi2)*(cos(alpha0)*cos(phi1) + sin(alpha0)*sin(phi1));
B2 = -sin(alpha2)*sin(phi2);
C2 = cos(alpha2)*cos(phi2)*(sin(alpha0)*cos(phi1) - cos(alpha0)*sin(phi1));
D2 = cg5 - C2;

R2 = sqrt(A2^2 + B2^2);
delta2 = atan2(B2, A2);

arg2 = D2 / R2;

if abs(arg2) > 1
    error('No real solution for theta2: |D2/R2| > 1');
end

theta2_sol1 = -delta2 + asin(arg2);
theta2_sol2 = -delta2 + (pi - asin(arg2));

% wrap to [-pi, pi]
theta2_sol1 = atan2(sin(theta2_sol1), cos(theta2_sol1));
theta2_sol2 = atan2(sin(theta2_sol2), cos(theta2_sol2));

% choose branch
theta2_chosen = theta2_sol1;



%% THETA4 (PASSIVE JOINT)
%% v1 = Q4.Q2.Q0.v5
%% Q2^T.Q0^T.v1 = Q4.v5

Rot_Y = @(a)[cos(a) 0 sin(a);
             0      1 0;
            -sin(a) 0 cos(a)];

Rot_DH = @(a,t)[cos(t), -sin(t)*cos(a),  sin(t)*sin(a);
                sin(t),  cos(t)*cos(a), -cos(t)*sin(a);
                0,       sin(a),         cos(a)];

Q0 = Rot_Y(alpha0);
Q2 = Rot_DH(alpha2, theta2_chosen);

v1 = [cos(phi1)*cos(phi2);
      sin(phi2);
     -sin(phi1)*cos(phi2)];

w4 = Q2.' * Q0.' * v1;

d = w4(1);
e = w4(2);
f = w4(3);

a = sin(beta1)*cos(beta2);
b = -cos(alpha4)*cos(beta1) + sin(alpha4)*sin(beta1)*sin(beta2);
c =  sin(alpha4)*cos(beta1) + cos(alpha4)*sin(beta1)*sin(beta2);

den4 = a^2 + b^2;
if abs(den4) < 1e-12
    error('Degenerate passive-joint equation: a^2+b^2 is too small');
end

s4 = (a*e + b*d) / den4;
c4 = (a*d - b*e) / den4;

% normalize to avoid small numerical drift
norm4 = sqrt(s4^2 + c4^2);
if norm4 > 1e-12
    s4 = s4 / norm4;
    c4 = c4 / norm4;
end

theta4_sol = atan2(s4, c4);

check_d = d - (a*c4 - b*s4);
check_e = e - (a*s4 + b*c4);
check_f = f - c;

%% =========================================================
%% THETA1 (ACTIVE JOINT OF LEG 2)
%% Constraint:
%% e3_2 . u2 = cos(alpha3)
%% u2 = Q2.Q4.e5_5
%% e3_2 = Q0^T.Q1.e3_3 
%% =========================================================

%Q1 = Rot_DH(alpha1, 0); 
Q2 = Rot_DH(alpha2, theta2_chosen);
Q4 = Rot_DH(alpha4, theta4_sol);

e5_5 = [0;0;1];
u2 = Q2 * Q4 * e5_5;

ux = u2(1);
uy = u2(2);
uz = u2(3);

A1 = sin(alpha1) * (cos(alpha0)*ux + sin(alpha0)*uz);
B1 = -sin(alpha1) * uy;
C1 = cos(alpha1) * (-sin(alpha0)*ux + cos(alpha0)*uz);
D1 = cos(alpha3) - C1;

R1 = sqrt(A1^2 + B1^2);
delta1 = atan2(B1, A1);

arg1 = D1 / R1;

if abs(arg1) > 1
    error('No real solution for theta1: |D1/R1| > 1');
end

theta1_sol1 = -delta1 + asin(arg1);
theta1_sol2 = -delta1 + (pi - asin(arg1));

% wrap to [-pi, pi]
theta1_sol1 = atan2(sin(theta1_sol1), cos(theta1_sol1));
theta1_sol2 = atan2(sin(theta1_sol2), cos(theta1_sol2));

%% =========================================================
%% DISPLAY
%% =========================================================
disp('================ FINAL RESULTS (DEGREES) ================')

disp('theta2 branch 1 (deg) = ')
disp(theta2_sol1 * 180/pi)

disp('theta2 branch 2 (deg) = ')
disp(theta2_sol2 * 180/pi)

disp('chosen theta2 (deg) = ')
disp(theta2_chosen * 180/pi)

disp('theta4 (deg) = ')
disp(theta4_sol * 180/pi)

disp('theta1 branch 1 (deg) = ')
disp(theta1_sol1 * 180/pi)

disp('theta1 branch 2 (deg) = ')
disp(theta1_sol2 * 180/pi)

disp('================ FINAL RESULTS (RADIANS) ================')

disp('theta2 branch 1 (rad) = ')
disp(theta2_sol1)

disp('theta2 branch 2 (rad) = ')
disp(theta2_sol2)

disp('chosen theta2 (rad) = ')
disp(theta2_chosen)

disp('theta4 (rad) = ')
disp(theta4_sol)

disp('theta1 branch 1 (rad) = ')
disp(theta1_sol1)

disp('theta1 branch 2 (rad) = ')
disp(theta1_sol2)

disp('================ CHECKS ================')
disp('check_d = ')
disp(check_d)
disp('check_e = ')
disp(check_e)
disp('check_f = ')
disp(check_f)