function out = mechanism_jacobian(theta1, theta2, theta4, phi1, phi2)
% mechanism_jacobian
% Inputs:
%   theta1, theta2, theta4, phi1, phi2  [rad]
%
% Output:
%   out.e3_1, out.e4_1, out.e5_1, out.v_1
%   out.K, out.J0, out.S, out.J

    alpha0 = pi/2;
    alpha1 = pi/2;
    alpha2 = pi/2;
    alpha4 = pi/2;

    z = [0;0;1];

    Q0 = RotY(alpha0);
    Q1 = RotDH(alpha1, theta1);
    Q2 = RotDH(alpha2, theta2);
    Q4 = RotDH(alpha4, theta4);

    e3_1 = Q1 * z;
    e4_1 = Q0 * Q2 * z;
    e5_1 = Q0 * Q2 * Q4 * z;

    e1_1 = z;
    e2_1 = Q0 * z;

    v_1 = [ cos(phi1)*cos(phi2);
            sin(phi2);
           -sin(phi1)*cos(phi2) ];

    K = [ dot(e5_1, cross(e1_1, e3_1)), 0;
          0, dot(v_1, cross(e2_1, e4_1)) ];

    J0 = [ cross(e3_1, e5_1).';
           cross(e4_1, v_1).' ];

    S = [ 0, -sin(phi1);
          0,  cos(phi1);
          1,  0 ];

    J = J0 * S;

    out.e3_1 = e3_1;
    out.e1_1 = e1_1;
    out.e4_1 = e4_1;
    out.e5_1 = e5_1;
    out.v_1  = v_1;
    out.K    = K;
    out.J0   = J0;
    out.S    = S;
    out.J    = J;
    out.e2_1 = e2_1;
end

function R = RotY(a)
    R = [ cos(a), 0, sin(a);
          0,      1, 0;
         -sin(a), 0, cos(a) ];
end

function R = RotDH(a, t)
    R = [ cos(t), -sin(t)*cos(a),  sin(t)*sin(a);
          sin(t),  cos(t)*cos(a), -cos(t)*sin(a);
          0,       sin(a),         cos(a) ];
end