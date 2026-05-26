clear; clc; close all;

deg = pi/180;
%deg_to_rad
phi1 = 0 * deg;
phi2 = 0 * deg;

alpha0 = pi/2;
alpha1 = pi/2;
alpha2 = pi/2;
alpha3 = pi/2;
alpha4 = pi/2;

beta1 = 0;
beta2 = 75 * deg;

theta3_ref = 15 * deg;

tolPhi1 = 1e-8;



d_des = direction_from_orientation(phi1, phi2);

solutions = IK_solns_updated(d_des, phi1, phi2, beta1, beta2, alpha0, alpha1, alpha2, alpha3, alpha4, theta3_ref, tolPhi1);

fprintf('\nDesired orientation:\n');
fprintf('phi1 = %.6f deg\n', rad2deg(phi1));
fprintf('phi2 = %.6f deg\n', rad2deg(phi2));

fprintf('\nDesired d vector in final frame:\n');
disp(d_des)

fprintf('\nNumber of IK solutions found: %d\n\n', length(solutions));

for k = 1:length(solutions)

    fprintf('===== Solution %d =====\n', k);

    fprintf('solution type = %s\n', solutions(k).type);

    fprintf('theta1 = %.6f rad = %.6f deg\n', ...
        solutions(k).theta1, rad2deg(solutions(k).theta1));

    fprintf('theta2 = %.6f rad = %.6f deg\n', ...
        solutions(k).theta2, rad2deg(solutions(k).theta2));

    fprintf('theta3 = %.6f rad = %.6f deg\n', ...
        solutions(k).theta3, rad2deg(solutions(k).theta3));

    fprintf('Forward d_check in final frame = \n');
    disp(solutions(k).d_check);

    fprintf('Recovered orientation from d_check:\n');
    fprintf('phi1_check = %.6f deg\n', rad2deg(solutions(k).phi1_check));
    fprintf('phi2_check = %.6f deg\n', rad2deg(solutions(k).phi2_check));

    fprintf('Error norm = %.3e\n\n', solutions(k).error_norm);
end



chosen = choose_preferred_solution(solutions);

fprintf('\n================ CHOSEN IK SOLUTION ================\n');
fprintf('type = %s\n', chosen.type);
fprintf('theta1 = %.6f deg\n', rad2deg(chosen.theta1));
fprintf('theta2 = %.6f deg\n', rad2deg(chosen.theta2));
fprintf('theta3 = %.6f deg\n', rad2deg(chosen.theta3));



function solutions = IK_solns_updated(d_des, phi1_des, phi2_des, beta1, beta2, alpha0, alpha1, alpha2, alpha3, alpha4,theta3_ref, tolPhi1)

    solutions = struct([]);
    sol_index = 0;

    %% ========================================================
    % SPECIAL CAD / DK-CONSISTENT BRANCH
    %% ========================================================
    %
    % Updated DK says:
    %
    % theta1 = 0  ->  theta2 = phi2
    %                  theta3 = theta3_ref + theta2
    %
    % This branch is used when phi1 = 0.

    if abs(phi1_des) < tolPhi1

        theta1 = 0;
        theta2 = phi2_des;
        theta3 = theta3_ref + theta2;

        [d_check, phi1_check, phi2_check] = FK_check(theta1, theta3, ...
                                                     beta1, beta2, ...
                                                     alpha1, alpha3);

        err = norm(d_check - d_des);

        sol_index = sol_index + 1;

        solutions(sol_index).theta1 = wrapToPiLocal(theta1);
        solutions(sol_index).theta2 = wrapToPiLocal(theta2);
        solutions(sol_index).theta3 = wrapToPiLocal(theta3);

        solutions(sol_index).d_check = d_check;
        solutions(sol_index).phi1_check = phi1_check;
        solutions(sol_index).phi2_check = phi2_check;

        solutions(sol_index).error_norm = err;
        solutions(sol_index).constraint_error = NaN;
        solutions(sol_index).type = '';

        % If this branch exactly matches the desired direction, return it first.
        % Normal IK solutions are still added below for comparison.
    end


    normal_solutions = IK_solns_normal(d_des, beta1, beta2, alpha0, alpha1, alpha2, alpha3, alpha4);

    for k = 1:length(normal_solutions)

        sol_index = sol_index + 1;

        solutions(sol_index).theta1 = normal_solutions(k).theta1;
        solutions(sol_index).theta2 = normal_solutions(k).theta2;
        solutions(sol_index).theta3 = normal_solutions(k).theta3;

        solutions(sol_index).d_check = normal_solutions(k).d_check;
        solutions(sol_index).phi1_check = normal_solutions(k).phi1_check;
        solutions(sol_index).phi2_check = normal_solutions(k).phi2_check;

        solutions(sol_index).error_norm = normal_solutions(k).error_norm;
        solutions(sol_index).constraint_error = normal_solutions(k).constraint_error;
        solutions(sol_index).type = 'normal constraint IK branch';
    end

    if isempty(solutions)
        error('No IK solution found.');
    end
end


function solutions = IK_solns_normal(d_des, beta1, beta2, alpha0, alpha1, alpha2, alpha3, alpha4)

    tol = 1e-9;
    solutions = struct([]);

    Rot_Y = @(alpha) [ cos(alpha), 0, sin(alpha);
                       0,          1, 0;
                      -sin(alpha), 0, cos(alpha)];

    Rot_DH = @(alpha,theta) [ cos(theta), -sin(theta)*cos(alpha),  sin(theta)*sin(alpha);
                              sin(theta),  cos(theta)*cos(alpha), -cos(theta)*sin(alpha);
                              0,           sin(alpha),             cos(alpha)];

    Q0 = Rot_Y(alpha0);

    e11 = [0;0;1];
    e44 = [0;0;1];

    v5 = [ sin(beta2);
           sin(beta1)*cos(beta2);
           cos(beta1)*cos(beta2)];

    v5 = v5 / norm(v5);

    d_des = d_des / norm(d_des);

    v = [ d_des(2);
         -d_des(1);
          d_des(3)];

    v = v / norm(v);

    vx = v(1);
    vy = v(2);

    %% Solve theta1

    b = sin(beta1)*cos(beta2);

    R = sqrt(vx^2 + vy^2);

    if R < tol
        error('Cannot solve theta1: vx and vy are both zero.');
    end

    if abs(b) > R + tol
        error('No real theta1 solution.');
    end

    ratio = clamp(b/R, -1, 1);
    psi = atan2(vy, vx);

    theta1_list = [
        psi + asin(ratio);
        psi + pi - asin(ratio)
    ];

    theta1_list = uniqueAngles(theta1_list);

    a = sin(beta2);
    c = cos(beta1)*cos(beta2);

    sol_index = 0;

    for i = 1:length(theta1_list)

        theta1 = wrapToPiLocal(theta1_list(i));
        Q1 = Rot_DH(alpha1, theta1);

        %% Solve theta3

        A = Q1.' * v;

        denom = a^2 + c^2;

        if denom < tol
            continue;
        end

        cos_theta3 = (a*A(1) - c*A(2)) / denom;
        sin_theta3 = (c*A(1) + a*A(2)) / denom;

        cos_theta3 = clamp(cos_theta3, -1, 1);
        sin_theta3 = clamp(sin_theta3, -1, 1);

        theta3_list = atan2(sin_theta3, cos_theta3);
        theta3_list = uniqueAngles(theta3_list);

        for j = 1:length(theta3_list)

            theta3 = wrapToPiLocal(theta3_list(j));
            Q3 = Rot_DH(alpha3, theta3);

            if norm(Q1*Q3*v5 - v) > 1e-7
                continue;
            end

            %% Solve theta2 from constraint

            u = Q1 * Q3 * e11;
            u = u / norm(u);

            C = cos(alpha4);

            A2 = -u(3);
            B2 = -u(2);

            R2 = sqrt(A2^2 + B2^2);

            if R2 < tol
                continue;
            end

            if abs(C) > R2 + tol
                continue;
            end

            ratio2 = clamp(C/R2, -1, 1);
            delta = atan2(B2, A2);

            theta2_list = [
                asin(ratio2) - delta;
                pi - asin(ratio2) - delta
            ];

            theta2_list = uniqueAngles(theta2_list);

            for m = 1:length(theta2_list)

                theta2 = wrapToPiLocal(theta2_list(m));
                Q2 = Rot_DH(alpha2, theta2);

                e41 = Q0 * Q2 * e44;

                constraint_error = dot(u,e41) - cos(alpha4);

                if abs(constraint_error) > 1e-7
                    continue;
                end

                v_check = Q1 * Q3 * v5;
                v_check = v_check / norm(v_check);

                d_check = final_direction_from_v1(v_check);

                err = norm(d_check - d_des);

                [phi1_check, phi2_check] = orientation_from_direction(d_check);

                sol_index = sol_index + 1;

                solutions(sol_index).theta1 = theta1;
                solutions(sol_index).theta2 = theta2;
                solutions(sol_index).theta3 = theta3;

                solutions(sol_index).d_check = d_check;
                solutions(sol_index).phi1_check = phi1_check;
                solutions(sol_index).phi2_check = phi2_check;

                solutions(sol_index).error_norm = err;
                solutions(sol_index).constraint_error = constraint_error;
            end
        end
    end

    if isempty(solutions)
        error('No normal IK solution found.');
    end
end

%% ============================================================
% FK CHECK USED BY SPECIAL BRANCH
%% ============================================================

function [d_check, phi1_check, phi2_check] = FK_check(theta1, theta3, beta1, beta2, alpha1, alpha3)

    Rot_DH = @(alpha,theta) [ cos(theta), -sin(theta)*cos(alpha),  sin(theta)*sin(alpha);
                              sin(theta),  cos(theta)*cos(alpha), -cos(theta)*sin(alpha);
                              0,           sin(alpha),             cos(alpha)];

    Q1 = Rot_DH(alpha1, theta1);
    Q3 = Rot_DH(alpha3, theta3);

    v5 = [ sin(beta2);
           sin(beta1)*cos(beta2);
           cos(beta1)*cos(beta2)];

    v5 = v5 / norm(v5);

    v_check = Q1 * Q3 * v5;
    v_check = v_check / norm(v_check);

    d_check = final_direction_from_v1(v_check);

    [phi1_check, phi2_check] = orientation_from_direction(d_check);
end

% SOLUTION SELECTION


function chosen = choose_preferred_solution(solutions)

    % Prefer the special DK-consistent branch if it exists and error is small.
    for k = 1:length(solutions)
        if contains(solutions(k).type, 'special') && solutions(k).error_norm < 1e-8
            chosen = solutions(k);
            return;
        end
    end

    % Otherwise choose solution closest to home.
    deg = pi/180;

    theta_home = [0; 0; 15*deg];

    cost = zeros(length(solutions),1);

    for k = 1:length(solutions)

        theta_k = [
            solutions(k).theta1;
            solutions(k).theta2;
            solutions(k).theta3
        ];

        cost(k) = norm(wrapToPiLocal(theta_k - theta_home));
    end

    [~,idx] = min(cost);

    chosen = solutions(idx);
end


% DIRECTION FUNCTIONS


function d = direction_from_orientation(phi1, phi2)

    d = [
        -cos(phi2)*sin(phi1);
         cos(phi2)*cos(phi1);
         sin(phi2)
    ];

    d = d / norm(d);
end

function d = final_direction_from_v1(v1)

    d = [-v1(2);
          v1(1);
          v1(3)];

    d = d / norm(d);
end

function [phi1, phi2] = orientation_from_direction(d)

    d = d / norm(d);

    dx = d(1);
    dy = d(2);
    dz = d(3);

    dz = clamp(dz, -1, 1);

    phi1 = atan2(-dx, dy);
    phi2 = asin(dz);
end

% HELPERS


function x = clamp(x, lo, hi)

    x = min(max(x, lo), hi);
end

function a = wrapToPiLocal(a)

    a = atan2(sin(a), cos(a));
end

function list_unique = uniqueAngles(list)

    tol = 1e-8;
    list_unique = [];

    for i = 1:length(list)

        a = wrapToPiLocal(list(i));

        if isempty(list_unique)
            list_unique = a;
        else
            diffs = abs(wrapToPiLocal(a - list_unique));
            if all(diffs > tol)
                list_unique(end+1,1) = a; %#ok<AGROW>
            end
        end
    end
end