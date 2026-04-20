clc;
clear;
deg = pi/180;
theta1 = 30*deg;
theta2 = 20*deg;
theta4 = 120.64*deg;
phi1   = 60.48*deg;
phi2   = 2.35*deg;

out = mechanism_jacobian(theta1, theta2, theta4, phi1, phi2);

disp('e3_1 ='); disp(out.e3_1)
disp('e4_1 ='); disp(out.e4_1)
disp('e5_1 ='); disp(out.e5_1)
disp('v_1  ='); disp(out.v_1)

disp('K =');  disp(out.K)
disp('J0 ='); disp(out.J0)
disp('S =');  disp(out.S)
disp('J =');  disp(out.J)

