clc; clear; close all;

deg = pi/180;

%deg_to_rad

theta1 = 0 * deg;
theta2 = 30 * deg;
alpha0 = pi/2;
alpha1 = pi/2;
alpha2 = pi/2;
alpha3 = pi/2;
alpha4 = pi/2;

beta1 = 0;
beta2 = 75 * deg;

theta3_ref = 15 * deg;

tolCross  = 1e-10;
tolTheta  = 1e-10;
tolTheta1 = 1e-8;

previousTheta3 = [];



dk = DK_solver(theta1, theta2, ...
               theta3_ref, ...
               alpha0, alpha1, alpha2, alpha3, alpha4, ...
               beta1, beta2, ...
               tolCross, tolTheta, tolTheta1, previousTheta3);


fprintf('\n================ DIRECT KINEMATICS RESULT ================\n');

fprintf('\nInput motor angles:\n');
fprintf('theta1 = %.6f deg\n', rad2deg(theta1));
fprintf('theta2 = %.6f deg\n', rad2deg(theta2));

fprintf('\nSolved passive angle:\n');
fprintf('theta3 = %.6f deg\n', rad2deg(dk.theta3));

fprintf('\nBranch information:\n');
fprintf('branch = %d\n', dk.branch);
fprintf('branch name = %s\n', dk.branchName);

fprintf('\nConstraint axes:\n');
fprintf('e31 = \n');
disp(dk.e31)

fprintf('e41 = \n');
disp(dk.e41)

fprintf('crossNorm = %.6e\n', dk.crossNorm);

fprintf('\nPassive direction u = Q1*Q3*e55:\n');
disp(dk.u)

fprintf('\nInternal output vector v1 = Q1*Q3*v5:\n');
disp(dk.v1)

fprintf('\nFinal-frame direction d:\n');
disp(dk.d)

fprintf('\nOutput orientation:\n');
fprintf('phi1 = %.6f deg\n', rad2deg(dk.phi1));
fprintf('phi2 = %.6f deg\n', rad2deg(dk.phi2));

fprintf('\nConstraint checks:\n');
fprintf('e31.u = %.6e\n', dk.constraint1);
fprintf('e41.u - cos(alpha4) = %.6e\n', dk.constraint2);

fprintf('\nReconstruction check:\n');
fprintf('norm(d - direction_from_orientation(phi1,phi2)) = %.6e\n', ...
    norm(dk.d - direction_from_orientation(dk.phi1, dk.phi2)));

if dk.branch == 0
    fprintf('\nSpecial CAD/theta1=0 branch check:\n');
    fprintf('expected theta3 = theta3_ref + theta2 = %.6f deg\n', ...
        rad2deg(theta3_ref + theta2));

    fprintf('expected v1 = [cos(theta2); 0; sin(theta2)] = \n');
    disp([cos(theta2); 0; sin(theta2)])

    fprintf('expected d = [0; cos(theta2); sin(theta2)] = \n');
    disp([0; cos(theta2); sin(theta2)])

    fprintf('error norm(v1 - expected_v1) = %.6e\n', ...
        norm(dk.v1 - [cos(theta2); 0; sin(theta2)]));

    fprintf('error norm(d - expected_d) = %.6e\n', ...
        norm(dk.d - [0; cos(theta2); sin(theta2)]));
end



function dk = DK_solver(theta1, theta2,theta3_ref,alpha0, alpha1, alpha2, alpha3, alpha4, beta1, beta2,tolCross, tolTheta, tolTheta1, previousTheta3)

    Rot_Y = @(alpha) [ cos(alpha), 0, sin(alpha);
                       0,          1, 0;
                      -sin(alpha), 0, cos(alpha)];

    Rot_DH = @(alpha,theta) [ cos(theta), -sin(theta)*cos(alpha),  sin(theta)*sin(alpha);
                              sin(theta),  cos(theta)*cos(alpha), -cos(theta)*sin(alpha);
                              0,           sin(alpha),             cos(alpha)];

    e33 = [0;0;1];
    e44 = [0;0;1];
    e55 = [0;0;1];

    v5 = [ sin(beta2);
           sin(beta1)*cos(beta2);
           cos(beta1)*cos(beta2)];

    v5 = v5 / norm(v5);

    Q0 = Rot_Y(alpha0);
    Q1 = Rot_DH(alpha1, theta1);
    Q2 = Rot_DH(alpha2, theta2);

    e31 = Q1 * e33;
    e41 = Q0 * Q2 * e44;

    cross_e = cross(e31, e41);
    crossNorm = norm(cross_e);

    if crossNorm < tolCross || abs(theta1) < tolTheta1

        theta3 = theta3_ref + theta2;
        branch = 0;
        branchName = 'special CAD/theta1=0 branch: theta3 = theta3_ref + theta2';

    else

        u_plus  =  cross_e / crossNorm;
        u_minus = -cross_e / crossNorm;

        theta3_plus  = theta3_from_u(u_plus,  theta1, tolTheta);
        theta3_minus = theta3_from_u(u_minus, theta1, tolTheta);

        if isempty(previousTheta3)
            err_plus  = abs(wrapToPiLocal(theta3_plus  - theta3_ref));
            err_minus = abs(wrapToPiLocal(theta3_minus - theta3_ref));
        else
            err_plus  = abs(wrapToPiLocal(theta3_plus  - previousTheta3));
            err_minus = abs(wrapToPiLocal(theta3_minus - previousTheta3));
        end

        if err_plus <= err_minus
            theta3 = theta3_plus;
            branch = 1;
            branchName = 'general constraint DK plus branch';
        else
            theta3 = theta3_minus;
            branch = 2;
            branchName = 'general constraint DK minus branch';
        end
    end

    Q3 = Rot_DH(alpha3, theta3);

    u = Q1 * Q3 * e55;
    u = u / norm(u);

    v1 = Q1 * Q3 * v5;
    v1 = v1 / norm(v1);

    d = final_direction_from_v1(v1);

    [phi1, phi2] = orientation_from_direction(d);

    constraint1 = dot(e31, u);
    constraint2 = dot(e41, u) - cos(alpha4);

    dk.theta1 = theta1;
    dk.theta2 = theta2;
    dk.theta3 = theta3;

    dk.branch = branch;
    dk.branchName = branchName;

    dk.e31 = e31;
    dk.e41 = e41;

    dk.crossNorm = crossNorm;

    dk.u = u;
    dk.v1 = v1;
    dk.d = d;

    dk.phi1 = phi1;
    dk.phi2 = phi2;

    dk.constraint1 = constraint1;
    dk.constraint2 = constraint2;
end

% SOLVE THETA3 FROM u


function theta3 = theta3_from_u(u, theta1, tol)

    ux = u(1);
    uy = u(2);
    uz = u(3);

    if abs(cos(theta1)) > abs(sin(theta1))
        sin_theta3 = ux / cos(theta1);
    else
        if abs(sin(theta1)) < tol
            error('Cannot solve theta3: theta1 too close to singular.');
        end
        sin_theta3 = uy / sin(theta1);
    end

    cos_theta3 = -uz;

    sin_theta3 = clamp(sin_theta3, -1, 1);
    cos_theta3 = clamp(cos_theta3, -1, 1);

    theta3 = atan2(sin_theta3, cos_theta3);
end


% FINAL-FRAME DIRECTION


function d = final_direction_from_v1(v1)

    d = [-v1(2);
          v1(1);
          v1(3)];

    d = d / norm(d);
end


% DIRECTION FROM ORIENTATION


function d = direction_from_orientation(phi1, phi2)

    d = [
        -cos(phi2)*sin(phi1);
         cos(phi2)*cos(phi1);
         sin(phi2)
    ];

    d = d / norm(d);
end


% ORIENTATION FROM DIRECTION


function [phi1, phi2] = orientation_from_direction(d)

    d = d / norm(d);

    dx = d(1);
    dy = d(2);
    dz = d(3);

    dz = clamp(dz, -1, 1);

    phi1 = atan2(-dx, dy);
    phi2 = asin(dz);
end



function x = clamp(x, lo, hi)

    x = min(max(x, lo), hi);
end

function a = wrapToPiLocal(a)

    a = atan2(sin(a), cos(a));
end