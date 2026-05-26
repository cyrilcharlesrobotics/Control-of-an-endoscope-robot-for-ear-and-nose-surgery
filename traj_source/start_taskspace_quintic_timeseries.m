clc; clear; close all;

deg = pi/180;

%% ============================================================
% USER SETTINGS
%% ============================================================

% Simulation time
tf = 2.0;          % total motion time [s]
dt = 0.001;        % sample time [s]
t = (0:dt:tf)';

% Initial task-space orientation
phi0 = [-40;
        -40] * deg;

% Final task-space orientation
phif = [-20;
        25] * deg;

% Initial and final task-space velocity
phidot0 = [0;
           0] * deg;

phidotf = [0;
           0] * deg;

% Initial and final task-space acceleration
phiddot0 = [0;
            0] * deg;

phiddotf = [0;
            0] * deg;

%% ============================================================
% ROBOT PARAMETERS
%% ============================================================

alpha0 = pi/2;
alpha1 = pi/2;
alpha2 = pi/2;
alpha3 = pi/2;
alpha4 = pi/2;

beta1 = 0;
beta2 = 75 * deg;

theta3_home = 15 * deg;

% Gear ratio
N = 70;

% Joint-side safety limit
theta_lim = 45 * deg;

%% ============================================================
% GENERATE QUINTIC TASK-SPACE TRAJECTORY
%% ============================================================

[phi, phidot, phiddot] = quintic_taskspace_traj( ...
    t, tf, phi0, phif, phidot0, phidotf, phiddot0, phiddotf);

%% ============================================================
% INITIAL IK FOR PLANT INITIAL CONDITION
%% ============================================================

d0 = direction_from_orientation(phi0(1), phi0(2));

solutions0 = IK_solns(d0, beta1, beta2, ...
                      alpha0, alpha1, alpha2, alpha3, alpha4);

fprintf('\n================ INITIAL IK SOLUTIONS ================\n');

for k = 1:length(solutions0)
    fprintf('Solution %d:\n', k);
    fprintf('theta1 = %.6f deg\n', rad2deg(solutions0(k).theta1));
    fprintf('theta2 = %.6f deg\n', rad2deg(solutions0(k).theta2));
    fprintf('theta3 = %.6f deg\n\n', rad2deg(solutions0(k).theta3));
end

% Choose safe initial IK solution closest to zero
theta_prev0 = [0; 0];

[theta0_from_IK, theta3_0_from_IK, ok0] = chooseSafeIKSolution( ...
    solutions0, theta_prev0, theta_lim);

if ~ok0
    error('Initial pose is outside the safe joint range [-45 deg, +45 deg].');
end

theta_dot0_from_IK = [0;
                      0];

theta3_dot0_from_IK = 0;

% Motor-side initial conditions
theta_m0_from_IK    = N * theta0_from_IK;
theta_mdot0_from_IK = N * theta_dot0_from_IK;

fprintf('\nChosen initial IK solution inside safe range.\n');

fprintf('\nInitial joint position for Simulink plant [deg]:\n');
disp(rad2deg(theta0_from_IK))

fprintf('Initial motor position for Simulink plant [deg]:\n');
disp(rad2deg(theta_m0_from_IK))

fprintf('Initial theta3 for Simulink [deg]:\n');
disp(rad2deg(theta3_0_from_IK))

%% ============================================================
% COMPUTE DESIRED JOINT TRAJECTORY FROM IK
%% ============================================================

Nt = length(t);

theta_d_all  = nan(2, Nt);
theta3_d_all = nan(1, Nt);

theta_prev = theta0_from_IK;

for ii = 1:Nt

    d_i = direction_from_orientation(phi(ii,1), phi(ii,2));

    sol_i = IK_solns(d_i, beta1, beta2, ...
                     alpha0, alpha1, alpha2, alpha3, alpha4);

    [theta_i, theta3_i, ok_i] = chooseSafeIKSolution( ...
        sol_i, theta_prev, theta_lim);

    if ~ok_i
        fprintf('\n====================================================\n');
        fprintf('TRAJECTORY REJECTED: NO SAFE IK BRANCH\n');
        fprintf('====================================================\n');
        fprintf('time = %.4f s\n', t(ii));
        fprintf('phi1 = %.3f deg\n', phi(ii,1)/deg);
        fprintf('phi2 = %.3f deg\n', phi(ii,2)/deg);
        error('No IK solution found inside [-45 deg, +45 deg]. Change phi0 or phif.');
    end

    theta_d_all(:,ii) = theta_i;
    theta3_d_all(ii)  = theta3_i;

    theta_prev = theta_i;
end

% Motor-side desired joint trajectory for checking or plotting
theta_m_d_all = N * theta_d_all;

%% ============================================================
% OFFLINE JOINT WORKSPACE SAFETY CHECK
%% ============================================================

theta1_d_all = theta_d_all(1,:);
theta2_d_all = theta_d_all(2,:);

unsafe_idx = abs(theta1_d_all) > theta_lim | abs(theta2_d_all) > theta_lim;

if any(unsafe_idx)

    fprintf('\n====================================================\n');
    fprintf('TRAJECTORY REJECTED: JOINT LIMIT EXCEEDED\n');
    fprintf('====================================================\n');

    fprintf('Allowed joint range:\n');
    fprintf('theta1, theta2 must stay within [-45 deg, +45 deg]\n\n');

    fprintf('theta1 min = %.3f deg\n', min(theta1_d_all)/deg);
    fprintf('theta1 max = %.3f deg\n', max(theta1_d_all)/deg);

    fprintf('theta2 min = %.3f deg\n', min(theta2_d_all)/deg);
    fprintf('theta2 max = %.3f deg\n', max(theta2_d_all)/deg);

    first_bad = find(unsafe_idx,1,'first');

    fprintf('\nFirst unsafe point:\n');
    fprintf('time = %.4f s\n', t(first_bad));
    fprintf('theta1 = %.3f deg\n', theta1_d_all(first_bad)/deg);
    fprintf('theta2 = %.3f deg\n', theta2_d_all(first_bad)/deg);

    error('Unsafe trajectory. Change phi0, phif, or reduce motion range.');

else

    fprintf('\n====================================================\n');
    fprintf('TRAJECTORY ACCEPTED\n');
    fprintf('====================================================\n');

    fprintf('theta1 min = %.3f deg\n', min(theta1_d_all)/deg);
    fprintf('theta1 max = %.3f deg\n', max(theta1_d_all)/deg);

    fprintf('theta2 min = %.3f deg\n', min(theta2_d_all)/deg);
    fprintf('theta2 max = %.3f deg\n', max(theta2_d_all)/deg);

    fprintf('\nBoth joints remain within [-45 deg, +45 deg].\n');

end

%% ============================================================
% CREATE SIMULINK TIMESERIES
%% ============================================================

% Each task-space signal is Nt x 2:
% column 1 = phi1
% column 2 = phi2

phi_d_ts     = timeseries(phi, t);
phidot_d_ts  = timeseries(phidot, t);
phiddot_d_ts = timeseries(phiddot, t);

phi_d_ts.Name     = 'phi_d';
phidot_d_ts.Name  = 'phidot_d';
phiddot_d_ts.Name = 'phiddot_d';

phi_d_ts     = setinterpmethod(phi_d_ts, 'linear');
phidot_d_ts  = setinterpmethod(phidot_d_ts, 'linear');
phiddot_d_ts = setinterpmethod(phiddot_d_ts, 'linear');

%% ============================================================
% OPTIONAL SEPARATE TIMESERIES
%% ============================================================

phi1_d_ts      = timeseries(phi(:,1), t);
phi2_d_ts      = timeseries(phi(:,2), t);

phi1dot_d_ts   = timeseries(phidot(:,1), t);
phi2dot_d_ts   = timeseries(phidot(:,2), t);

phi1ddot_d_ts  = timeseries(phiddot(:,1), t);
phi2ddot_d_ts  = timeseries(phiddot(:,2), t);

phi1_d_ts.Name      = 'phi1_d';
phi2_d_ts.Name      = 'phi2_d';

phi1dot_d_ts.Name   = 'phi1dot_d';
phi2dot_d_ts.Name   = 'phi2dot_d';

phi1ddot_d_ts.Name  = 'phi1ddot_d';
phi2ddot_d_ts.Name  = 'phi2ddot_d';

%% ============================================================
% SAVE FOR SIMULINK
%% ============================================================

save('taskspace_quintic_timeseries.mat', ...
     't', 'tf', 'dt', ...
     'deg', ...
     'N', ...
     'theta_lim', ...
     'alpha0', 'alpha1', 'alpha2', 'alpha3', 'alpha4', ...
     'beta1', 'beta2', 'theta3_home', ...
     'phi0', 'phif', ...
     'phi', 'phidot', 'phiddot', ...
     'theta_d_all', 'theta3_d_all', 'theta_m_d_all', ...
     'phi_d_ts', 'phidot_d_ts', 'phiddot_d_ts', ...
     'phi1_d_ts', 'phi2_d_ts', ...
     'phi1dot_d_ts', 'phi2dot_d_ts', ...
     'phi1ddot_d_ts', 'phi2ddot_d_ts', ...
     'theta0_from_IK', 'theta_dot0_from_IK', ...
     'theta_m0_from_IK', 'theta_mdot0_from_IK', ...
     'theta3_0_from_IK', 'theta3_dot0_from_IK');

disp(' ')
disp('Saved task-space trajectory, IK safety check, and initial values to:')
disp('taskspace_quintic_timeseries.mat')

%% ============================================================
% DISPLAY FINAL VALUES
%% ============================================================

fprintf('\n================ TRAJECTORY CHECK ================\n');

fprintf('Initial phi [deg] = [%.6f %.6f]\n', ...
    phi(1,1)/deg, phi(1,2)/deg);

fprintf('Final phi [deg]   = [%.6f %.6f]\n', ...
    phi(end,1)/deg, phi(end,2)/deg);

fprintf('Initial phidot [deg/s] = [%.6f %.6f]\n', ...
    phidot(1,1)/deg, phidot(1,2)/deg);

fprintf('Final phidot [deg/s]   = [%.6f %.6f]\n', ...
    phidot(end,1)/deg, phidot(end,2)/deg);

fprintf('Initial phiddot [deg/s^2] = [%.6f %.6f]\n', ...
    phiddot(1,1)/deg, phiddot(1,2)/deg);

fprintf('Final phiddot [deg/s^2]   = [%.6f %.6f]\n', ...
    phiddot(end,1)/deg, phiddot(end,2)/deg);

fprintf('\nInitial theta_d [deg] = [%.6f %.6f]\n', ...
    theta_d_all(1,1)/deg, theta_d_all(2,1)/deg);

fprintf('Final theta_d [deg]   = [%.6f %.6f]\n', ...
    theta_d_all(1,end)/deg, theta_d_all(2,end)/deg);

fprintf('\nInitial theta_m_d [deg] = [%.6f %.6f]\n', ...
    theta_m_d_all(1,1)/deg, theta_m_d_all(2,1)/deg);

fprintf('Final theta_m_d [deg]   = [%.6f %.6f]\n', ...
    theta_m_d_all(1,end)/deg, theta_m_d_all(2,end)/deg);

%% ============================================================
% PLOTS
%% ============================================================

figure;
plot(t, phi(:,1)/deg, 'LineWidth', 1.5); hold on;
plot(t, phi(:,2)/deg, 'LineWidth', 1.5);
xlabel('Time [s]');
ylabel('\phi [deg]');
title('Task-space quintic position');
legend('\phi_1','\phi_2');
grid on;

figure;
plot(t, phidot(:,1)/deg, 'LineWidth', 1.5); hold on;
plot(t, phidot(:,2)/deg, 'LineWidth', 1.5);
xlabel('Time [s]');
ylabel('\dot{\phi} [deg/s]');
title('Task-space quintic velocity');
legend('$\dot{\phi}_1$', '$\dot{\phi}_2$', 'Interpreter', 'latex');
grid on;

figure;
plot(t, phiddot(:,1)/deg, 'LineWidth', 1.5); hold on;
plot(t, phiddot(:,2)/deg, 'LineWidth', 1.5);
xlabel('Time [s]');
ylabel('\ddot{\phi} [deg/s^2]');
title('Task-space quintic acceleration');
legend('$\ddot{\phi}_1$', '$\ddot{\phi}_2$', 'Interpreter', 'latex');
grid on;

figure;
plot(t, theta_d_all(1,:)/deg, 'LineWidth', 1.5); hold on;
plot(t, theta_d_all(2,:)/deg, 'LineWidth', 1.5);
yline(45,'r--','+45 deg');
yline(-45,'r--','-45 deg');
xlabel('Time [s]');
ylabel('\theta_d [deg]');
title('Desired joint trajectory with safety limits');
legend('\theta_{1d}','\theta_{2d}','Location','best');
grid on;

figure;
plot(t, theta_m_d_all(1,:)/deg, 'LineWidth', 1.5); hold on;
plot(t, theta_m_d_all(2,:)/deg, 'LineWidth', 1.5);
yline(N*45,'r--','+3150 deg');
yline(-N*45,'r--','-3150 deg');
xlabel('Time [s]');
ylabel('\theta_{m,d} [deg]');
title('Desired motor-side trajectory');
legend('\theta_{m1,d}','\theta_{m2,d}','Location','best');
grid on;

%% ============================================================
% IMPORTANT SIMULINK NOTES
%% ============================================================

disp(' ')
disp('================ SIMULINK INITIAL CONDITIONS ================')
disp('For MOTOR-SIDE plant:')
disp('Set theta_m position integrator initial condition to:')
disp('theta_m0_from_IK')
disp(' ')
disp('Set theta_m_dot velocity integrator initial condition to:')
disp('theta_mdot0_from_IK')
disp(' ')
disp('After motor-side integrators:')
disp('theta_j     = theta_m / N')
disp('theta_j_dot = theta_m_dot / N')
disp(' ')
disp('If you use integral error, set its initial condition to:')
disp('[0;0]')
disp(' ')
disp('Set Simulink Stop Time to:')
disp(num2str(tf))
disp(' ')
disp('Offline safety check:')
disp('theta1_d and theta2_d must stay within [-45 deg, +45 deg].')

%% ============================================================
% LOCAL FUNCTIONS
%% ============================================================

function [q, qdot, qddot] = quintic_taskspace_traj( ...
    t, tf, q0, qf, qdot0, qdotf, qddot0, qddotf)

    n = length(q0);
    Nt = length(t);

    q     = zeros(Nt,n);
    qdot  = zeros(Nt,n);
    qddot = zeros(Nt,n);

    for i = 1:n

        A = [1, 0,      0,        0,          0,           0;
             0, 1,      0,        0,          0,           0;
             0, 0,      2,        0,          0,           0;
             1, tf,     tf^2,     tf^3,       tf^4,        tf^5;
             0, 1,      2*tf,     3*tf^2,     4*tf^3,      5*tf^4;
             0, 0,      2,        6*tf,       12*tf^2,     20*tf^3];

        b = [q0(i);
             qdot0(i);
             qddot0(i);
             qf(i);
             qdotf(i);
             qddotf(i)];

        a = A\b;

        q(:,i) = ...
            a(1) + ...
            a(2)*t + ...
            a(3)*t.^2 + ...
            a(4)*t.^3 + ...
            a(5)*t.^4 + ...
            a(6)*t.^5;

        qdot(:,i) = ...
            a(2) + ...
            2*a(3)*t + ...
            3*a(4)*t.^2 + ...
            4*a(5)*t.^3 + ...
            5*a(6)*t.^4;

        qddot(:,i) = ...
            2*a(3) + ...
            6*a(4)*t + ...
            12*a(5)*t.^2 + ...
            20*a(6)*t.^3;
    end
end

function [theta_sel, theta3_sel, ok] = chooseSafeIKSolution(solutions, theta_prev, theta_lim)

    ok = false;
    theta_sel = [NaN; NaN];
    theta3_sel = NaN;

    best_cost = inf;

    for k = 1:length(solutions)

        theta_k = [solutions(k).theta1;
                   solutions(k).theta2];

        theta3_k = solutions(k).theta3;

        inside_limit = abs(theta_k(1)) <= theta_lim && ...
                       abs(theta_k(2)) <= theta_lim;

        if inside_limit

            cost = norm(wrapToPiVec(theta_k - theta_prev));

            if cost < best_cost
                best_cost = cost;
                theta_sel = theta_k;
                theta3_sel = theta3_k;
                ok = true;
            end
        end
    end
end

function solutions = IK_solns(d_des, beta1, beta2, alpha0, alpha1, alpha2, alpha3, alpha4)

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

    % Final-frame inverse conversion:
    % d = [-v(2); v(1); v(3)]
    % therefore:
    % v = [d(2); -d(1); d(3)]
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
        error('No real theta1 solution: abs(b) > sqrt(vx^2 + vy^2).');
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
            if norm(A - [0;0;b]) > 1e-7
                continue;
            end
            theta3_list = 0;
        else
            cos_theta3 = (a*A(1) - c*A(2)) / denom;
            sin_theta3 = (c*A(1) + a*A(2)) / denom;

            cos_theta3 = clamp(cos_theta3, -1, 1);
            sin_theta3 = clamp(sin_theta3, -1, 1);

            theta3_list = atan2(sin_theta3, cos_theta3);
        end

        theta3_list = uniqueAngles(theta3_list);

        for j = 1:length(theta3_list)

            theta3 = wrapToPiLocal(theta3_list(j));
            Q3 = Rot_DH(alpha3, theta3);

            if norm(Q1*Q3*v5 - v) > 1e-7
                continue;
            end

            %% Solve theta2

            u = Q1 * Q3 * e11;

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

                %% Forward check

                v_check = Q1 * Q3 * v5;
                v_check = v_check / norm(v_check);

                d_check = [-v_check(2);
                            v_check(1);
                            v_check(3)];

                d_check = d_check / norm(d_check);

                err = norm(d_check - d_des);

                [phi1_check, phi2_check] = orientation_from_direction(d_check);

                sol_index = sol_index + 1;

                solutions(sol_index).theta1 = theta1;
                solutions(sol_index).theta3 = theta3;
                solutions(sol_index).theta2 = theta2;
                solutions(sol_index).d_check = d_check;
                solutions(sol_index).phi1_check = phi1_check;
                solutions(sol_index).phi2_check = phi2_check;
                solutions(sol_index).error_norm = err;
                solutions(sol_index).constraint_error = constraint_error;
            end
        end
    end

    if isempty(solutions)
        error('No IK solution found.');
    end
end

function d = direction_from_orientation(phi1, phi2)

    d = [
        -cos(phi2)*sin(phi1);
         cos(phi2)*cos(phi1);
         sin(phi2)
    ];

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

function x = clamp(x, lo, hi)

    x = min(max(x, lo), hi);
end

function a = wrapToPiLocal(a)

    a = atan2(sin(a), cos(a));
end

function a = wrapToPiVec(a)

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
                list_unique(end+1,1) = a; 
            end
        end
    end
end