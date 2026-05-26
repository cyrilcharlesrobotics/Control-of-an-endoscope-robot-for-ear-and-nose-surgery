clc; clear; close all;

deg = pi/180;

%% ============================================================
% USER SETTINGS
%% ============================================================

beta1 = 0*deg;
beta2 = 75*deg;

theta1_vals = -50:0.2:50;
theta2_vals = -50:0.2:50;

theta3_home = 15*deg;

tolTheta  = 1e-10;
tolTheta1 = 1e-8;      
tolCross  = 1e-10;
tolK      = 1e-12;
tolSigma  = 1e-3;
tolDet    = 1e-6;

sigma_safe = 1e-2;
cond_safe  = 1000;

theta1_lim = 50;
theta2_lim = 50;

%% ============================================================
% INITIALIZATION
%% ============================================================

n1 = numel(theta1_vals);
n2 = numel(theta2_vals);
maxRows = n1*n2;

data = nan(maxRows,19);
row = 0;

alpha0 = pi/2;
alpha1 = pi/2;
alpha2 = pi/2;
alpha3 = pi/2;
alpha4 = pi/2;

e33 = [0;0;1];
e44 = [0;0;1];
e55 = [0;0;1];

e1_1 = [0;0;1];
e2_1 = [1;0;0];

Rot_Y = @(alpha) [ cos(alpha), 0, sin(alpha);
                   0,          1, 0;
                  -sin(alpha), 0, cos(alpha)];

Rot_DH = @(alpha,theta) [ cos(theta), -sin(theta)*cos(alpha),  sin(theta)*sin(alpha);
                          sin(theta),  cos(theta)*cos(alpha), -cos(theta)*sin(alpha);
                          0,           sin(alpha),             cos(alpha)];

Q0 = Rot_Y(alpha0);

v5 = [ sin(beta2);
       sin(beta1)*cos(beta2);
       cos(beta1)*cos(beta2)];

v5 = v5 / norm(v5);

%% ============================================================
% MAIN DK + SINGULARITY SCAN
%% ============================================================

for t1d = theta1_vals

    theta1 = t1d*deg;
    Q1 = Rot_DH(alpha1, theta1);

    e3_1 = Q1 * e33;

    for t2d = theta2_vals

        theta2 = t2d*deg;
        Q2 = Rot_DH(alpha2, theta2);

        e4_1 = Q0 * Q2 * e44;

        cross_e = cross(e3_1, e4_1);
        cross_norm = norm(cross_e);

        %% ====================================================
        % DK BRANCH LOGIC
        %% ====================================================
        %
        % branch:
        %  0 = theta1=0 physical branch
        % -1 = true cross-product singular fallback
        %  1 = general plus branch
        %  2 = general minus branch
        %
        % branchFlag:
        % 0 = normal cross-product branch
        % 1 = theta1=0 physical branch
        % 2 = true DK singular, cross product zero

        if abs(theta1) < tolTheta1

            % Physical branch:
            % theta3 = theta3_home + theta2
            %
            % At theta1 = 0:
            % v1 = [cos(theta2); 0; sin(theta2)]
            % d  = [0; cos(theta2); sin(theta2)]
            % phi1 = 0
            % phi2 = theta2

            theta3 = theta3_home + theta2;
            Q3 = Rot_DH(alpha3, theta3);

            u = Q1 * Q3 * e55;
            u = u / norm(u);

            branch = 0;
            branchFlag = 1;

        elseif cross_norm < tolCross

            % True DK singular case: u cannot be computed from cross product

            theta3 = theta3_home;
            Q3 = Rot_DH(alpha3, theta3);

            u = Q1 * Q3 * e55;
            u = u / norm(u);

            branch = -1;
            branchFlag = 2;

        else

            % General constraint DK case

            u_plus  =  cross_e / cross_norm;
            u_minus = -cross_e / cross_norm;

            theta3_plus  = theta3_from_u(u_plus,  theta1, tolTheta);
            theta3_minus = theta3_from_u(u_minus, theta1, tolTheta);

            err_plus  = abs(wrapToPiLocal(theta3_plus  - theta3_home));
            err_minus = abs(wrapToPiLocal(theta3_minus - theta3_home));

            if err_plus <= err_minus
                u = u_plus;
                theta3 = theta3_plus;
                branch = 1;
            else
                u = u_minus;
                theta3 = theta3_minus;
                branch = 2;
            end

            Q3 = Rot_DH(alpha3, theta3);
            branchFlag = 0;
        end

        %% ====================================================
        % OUTPUT VECTOR
        %% ====================================================

        e5_1 = Q1 * Q3 * e55;
        e5_1 = e5_1 / norm(e5_1);

        v_internal = Q1 * Q3 * v5;
        v_internal = v_internal / norm(v_internal);

        d = [-v_internal(2);
              v_internal(1);
              v_internal(3)];

        d = d / norm(d);

        vx = d(1);
        vy = d(2);
        vz = d(3);

        phi2 = asin(max(min(vz,1),-1));
        phi1 = atan2(-vx, vy);

        v1 = v_internal;

        %% ====================================================
        % JACOBIAN / SINGULARITY METRICS
        %% ====================================================

        S_omega_1 = [0,  sin(phi1);
                     0, -cos(phi1);
                     1,  0];

        k11 = dot(v1,   cross(e1_1, e3_1));
        k22 = dot(e5_1, cross(e2_1, e4_1));

        K = [k11, 0;
             0,   k22];

        J0 = [cross(e3_1, v1).';
              cross(e4_1, e5_1).'];

        J = J0 * S_omega_1;

        %% ====================================================
        % SINGULARITY CLASSIFICATION
        %% ====================================================
        %
        % singType:
        % 0 = regular
        % 1 = k11 = 0
        % 2 = k22 = 0
        % 3 = both K terms zero
        % 4 = near-singular M
        % 5 = true DK singular, cross product zero
        % 7 = theta1 = 0 branch, but Jacobian regular
        %
        % Important:
        % singType = 7 is treated as regular in the plots.

        k11_zero = abs(k11) < tolK;
        k22_zero = abs(k22) < tolK;

        singType = 0;

        if branchFlag == 2

            singType = 5;
            detM = 0;
            sigmaMin = 0;
            condM = inf;

        elseif k11_zero || k22_zero

            detM = 0;
            sigmaMin = 0;
            condM = inf;

            if k11_zero && ~k22_zero
                singType = 1;
            elseif ~k11_zero && k22_zero
                singType = 2;
            else
                singType = 3;
            end

        else

            M = K \ J;

            detM = det(M);
            sv = svd(M);
            sigmaMin = min(sv);
            condM = cond(M);

            if abs(detM) < tolDet || sigmaMin < tolSigma
                singType = 4;
            else
                if branchFlag == 1
                    singType = 7;
                else
                    singType = 0;
                end
            end
        end

        row = row + 1;

        data(row,:) = [ ...
            t1d, t2d, branch, branchFlag, theta3/deg, ...
            phi1/deg, phi2/deg, ...
            detM, sigmaMin, condM, ...
            k11, k22, ...
            vx, vy, vz, ...
            e5_1(1), e5_1(2), e5_1(3), ...
            singType];
    end
end

data = data(1:row,:);

%% ============================================================
% TABLE
%% ============================================================

T = array2table(data, ...
    'VariableNames', {'theta1_deg','theta2_deg','branch','branchFlag','theta3_deg', ...
                      'phi1_deg','phi2_deg', ...
                      'detM','sigmaMin','condM', ...
                      'k11','k22', ...
                      'vx','vy','vz', ...
                      'e5x','e5y','e5z', ...
                      'singType'});

disp('Number of reachable configurations:')
disp(height(T))

disp('Number of reachable task-space points:')
P = round([T.phi1_deg, T.phi2_deg], 3);
P_unique = unique(P, 'rows');
disp(size(P_unique,1))

%% ============================================================
% SINGULARITY INDEXES
%% ============================================================

idx_regular = T.singType == 0 | T.singType == 7;
idx_k11     = T.singType == 1;
idx_k22     = T.singType == 2;
idx_both    = T.singType == 3;
idx_nearM   = T.singType == 4;
idx_trueDK  = T.singType == 5;

% Do not count singType = 7 as singular
Ts = T(T.sigmaMin < tolSigma | ismember(T.singType,[1 2 3 4 5]), :);
Ts_sorted = sortrows(Ts,'sigmaMin','ascend');

Ts_clean = Ts(abs(Ts.theta1_deg) > 1.0 | abs(Ts.theta2_deg) > 1.0, :);

disp('Number of singular / near-singular points found:')
disp(height(Ts))

disp('Most singular points first:')
if height(Ts_sorted) > 0
    disp(Ts_sorted(1:min(20,height(Ts_sorted)),:))
else
    disp('No near-singular points found.')
end

%% ============================================================
% SAFE WORKSPACE
%% ============================================================

T_safe = T( ...
    abs(T.theta1_deg) <= theta1_lim & ...
    abs(T.theta2_deg) <= theta2_lim & ...
    T.sigmaMin >= sigma_safe & ...
    T.condM <= cond_safe & ...
    (T.singType == 0 | T.singType == 7), :);

disp('Total reachable configurations:')
disp(height(T))

disp('Safe/free configurations:')
disp(height(T_safe))

%% ============================================================
% EXPORT
%% ============================================================

writetable(T,'DK_full_workspace_scan.csv');
writetable(Ts_sorted,'DK_near_singular_points.csv');
writetable(Ts_clean,'DK_near_singular_points_clean.csv');
writetable(T_safe,'DK_safe_workspace.csv');

disp(' ')
disp('Exported:')
disp('DK_full_workspace_scan.csv')
disp('DK_near_singular_points.csv')
disp('DK_near_singular_points_clean.csv')
disp('DK_safe_workspace.csv')

%% ============================================================
% PLOTS
%% ============================================================

figure;
scatter(T.theta1_deg,T.theta2_deg,10,log10(max(T.sigmaMin,1e-16)),'filled');
colorbar;
xlabel('\theta_1 (deg)');
ylabel('\theta_2 (deg)');
title('DK branch: log_{10}(\sigma_{min})');
grid on;
xlim([min(theta1_vals) max(theta1_vals)]);
ylim([min(theta2_vals) max(theta2_vals)]);

figure; hold on;

scatter(T.theta1_deg,T.theta2_deg,4,[0.85 0.85 0.85],'filled');

scatter(T.theta1_deg(idx_k11),    T.theta2_deg(idx_k11),    20, 'r', 'filled');
scatter(T.theta1_deg(idx_k22),    T.theta2_deg(idx_k22),    20, 'b', 'filled');
scatter(T.theta1_deg(idx_both),   T.theta2_deg(idx_both),   30, 'm', 'filled');
scatter(T.theta1_deg(idx_nearM),  T.theta2_deg(idx_nearM),  15, 'k');
scatter(T.theta1_deg(idx_trueDK), T.theta2_deg(idx_trueDK), 35, 'y', 'filled');

xlabel('\theta_1 (deg)');
ylabel('\theta_2 (deg)');
title('DK branch: singularity types in joint space');

legend('Reachable', ...
       'k_{11}=0','k_{22}=0','both K terms zero', ...
       'near-singular M','true DK singular');

grid on;
axis equal;
xlim([min(theta1_vals) max(theta1_vals)]);
ylim([min(theta2_vals) max(theta2_vals)]);

figure; hold on;

scatter(T.phi1_deg,T.phi2_deg,5,[0.8 0.8 0.8],'filled');
scatter(Ts.phi1_deg,Ts.phi2_deg,12,'r','filled');

xlabel('\phi_1 (deg)');
ylabel('\phi_2 (deg)');
title('DK branch: task-space singular boundary');
legend('Workspace','Singular / near-singular');
grid on;

figure; hold on;

scatter(T.phi1_deg,T.phi2_deg,5,[0.85 0.85 0.85],'filled');

scatter(T.phi1_deg(idx_k11),    T.phi2_deg(idx_k11),    18, 'r', 'filled');
scatter(T.phi1_deg(idx_k22),    T.phi2_deg(idx_k22),    18, 'b', 'filled');
scatter(T.phi1_deg(idx_both),   T.phi2_deg(idx_both),   25, 'm', 'filled');
scatter(T.phi1_deg(idx_nearM),  T.phi2_deg(idx_nearM),  10, 'k');
scatter(T.phi1_deg(idx_trueDK), T.phi2_deg(idx_trueDK), 30, 'y', 'filled');

xlabel('\phi_1 (deg)');
ylabel('\phi_2 (deg)');
title('DK branch: task-space singularity types');

legend('Workspace','k_{11}=0','k_{22}=0','both K terms zero', ...
       'near-singular M','true DK singular');

grid on;

figure; hold on;

scatter3(T.vx,T.vy,T.vz,4,[0.85 0.85 0.85],'filled');
scatter3(Ts.vx,Ts.vy,Ts.vz,15,'r','filled');

axis equal;
xlabel('X_f');
ylabel('Y_f');
zlabel('Z_f');
title('DK branch: orientation workspace');
legend('Workspace','Singular / near-singular');
grid on;
view(3);

figure; hold on;

scatter3(T.vx,T.vy,T.vz,4,[0.85 0.85 0.85],'filled');

scatter3(T.vx(idx_k11),    T.vy(idx_k11),    T.vz(idx_k11),    18, 'r', 'filled');
scatter3(T.vx(idx_k22),    T.vy(idx_k22),    T.vz(idx_k22),    18, 'b', 'filled');
scatter3(T.vx(idx_both),   T.vy(idx_both),   T.vz(idx_both),   25, 'm', 'filled');
scatter3(T.vx(idx_nearM),  T.vy(idx_nearM),  T.vz(idx_nearM),  10, 'k');
scatter3(T.vx(idx_trueDK), T.vy(idx_trueDK), T.vz(idx_trueDK), 30, 'y', 'filled');

axis equal;
xlabel('X_f');
ylabel('Y_f');
zlabel('Z_f');
title('DK branch: separated singularities on sphere');

legend('Workspace','k_{11}=0','k_{22}=0','both K terms zero', ...
       'near-singular M','true DK singular');

grid on;
view(3);

figure; hold on;

scatter(T.phi1_deg, T.phi2_deg, 4, [0.85 0.85 0.85], 'filled');
scatter(T_safe.phi1_deg, T_safe.phi2_deg, 8, 'g', 'filled');

xlabel('\phi_1 (deg)');
ylabel('\phi_2 (deg)');
title('DK reachable workspace and safe workspace');
legend('Reachable','Safe/free');
grid on;

figure; hold on;

scatter3(T.vx, T.vy, T.vz, 4, [0.85 0.85 0.85], 'filled');
scatter3(T_safe.vx, T_safe.vy, T_safe.vz, 8, 'g', 'filled');

axis equal;
xlabel('X_f');
ylabel('Y_f');
zlabel('Z_f');
title('DK safe free orientation workspace');
legend('Reachable','Safe/free');
grid on;
view(3);

%% ============================================================
% HELPER FUNCTIONS
%% ============================================================

function theta3 = theta3_from_u(u, theta1, tol)

    ux = u(1);
    uy = u(2);
    uz = u(3);

    % Solve theta3 from:
    %
    % Q1 Q3 e55 =
    % [cos(theta1)*sin(theta3);
    %  sin(theta1)*sin(theta3);
    % -cos(theta3)]

    if abs(cos(theta1)) > abs(sin(theta1))
        sin_theta3 = ux / cos(theta1);
    else
        if abs(sin(theta1)) < tol
            error('Cannot compute sin(theta3): theta1 is singular.');
        end
        sin_theta3 = uy / sin(theta1);
    end

    cos_theta3 = -uz;

    sin_theta3 = max(min(sin_theta3,1),-1);
    cos_theta3 = max(min(cos_theta3,1),-1);

    theta3 = atan2(sin_theta3, cos_theta3);
end

function a = wrapToPiLocal(a)

    a = atan2(sin(a), cos(a));
end