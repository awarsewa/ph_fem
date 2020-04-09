function [M, K, B, Lout] = get_element_matrix_ext(dim, elem_type, props)
%GET_ELEMENT_MATRIX Summary of this function goes here
%   Detailed explanation goes here

E = props(1);
rho = props(2);
A = props(3);
Iy = props(4);
Iz = props(5);
Ip = (Iz + Iy)/A;
G = props(6);
It = props(7);
L = props(8);

% Element type:
et = elem_type(1);

% Actuation type:
at = elem_type(2);

% Load type defines:
load elem_types type_*

switch at
    case {0, 1}
        % no actuator: passive element, B = 0
        % Parallel actuator: passive element, B = [1, 0, ... , 0, -1, 0, ... , 0]
        switch et
            % Link Element (supports normal forces only):
            case 1
                % Mass:
                M = rho*A*L/420*[140, 0, 0,  70, 0, 0; ...
                                   0, 0, 0,   0, 0, 0; ...
                                   0, 0, 0,   0, 0, 0; ...
                                  70, 0, 0, 140, 0, 0; ...
                                   0, 0, 0,   0, 0, 0; ...
                                   0, 0, 0,   0, 0, 0];

                % Stiffness
                K = E*A/L*[ 1, 0, 0, -1, 0, 0; ...
                            0, 0, 0,  0, 0, 0; ...
                            0, 0, 0,  0, 0, 0; ...
                           -1, 0, 0,  1, 0, 0; ...
                            0, 0, 0,  0, 0, 0; ...
                            0, 0, 0,  0, 0, 0];

                % Input matrix
                if at == type_actparallel
                    B = [1; 0; 0; -1; 0; 0];
                else
                    B = double.empty(6,0);
                end

                % Output of element length:
                Lout = [-1, 0, 0, 1, 0, 0];

                % Delete z-direction for 2-dimension matrix:
                if dim == 2
                    M([3,6],:) = [];
                    M(:,[3,6]) = [];

                    K([3,6],:) = [];
                    K(:,[3,6]) = [];

                    B([3,6],:) = [];
                    Lout([3,6]) = [];
                end

            % Euler-Bernoulli beam element (supports normal forces, bending and torsion) für kubische Ansatzfunktionen:
            case 2
                % Mass:
                M = rho*A*L/420 * [ 140,     0,     0,      0,      0,      0,  70,     0,     0,      0,      0,      0; ...
                                      0,   156,     0,      0,      0,  -22*L,   0,    54,     0,      0,      0,   13*L; ...
                                      0,     0,   156,      0,   22*L,      0,   0,     0,    54,      0,  -13*L,      0; ...
                                      0,     0,     0, 140*Ip,      0,      0,   0,     0,     0,  70*Ip,      0,      0; ...
                                      0,     0,  22*L,      0,  4*L^2,      0,   0,     0,  13*L,      0, -3*L^2,      0; ...
                                      0, -22*L,     0,      0,      0,  4*L^2,   0, -13*L,     0,      0,      0, -3*L^2; ...
                                     70,     0,     0,      0,      0,      0, 140,     0,     0,      0,      0,      0; ...
                                      0,    54,     0,      0,      0,  -13*L,   0,   156,     0,      0,      0,   22*L; ...
                                      0,     0,    54,      0,   13*L,      0,   0,     0,   156,      0,  -22*L,      0; ...
                                      0,     0,     0,  70*Ip,      0,      0,   0,     0,     0, 140*Ip,      0,      0; ...
                                      0,     0, -13*L,      0, -3*L^2,      0,   0,     0, -22*L,      0,  4*L^2,      0; ...
                                      0,  13*L,     0,      0,      0, -3*L^2,   0,  22*L,     0,      0,      0,  4*L^2];

                % Stiffness:
                K = E * [  A/L,            0,            0,           0,           0,           0, -A/L,            0,            0,           0,           0,          0; ...
                             0,  (12*Iy)/L^3,            0,           0,           0, -(6*Iy)/L^2,    0, -(12*Iy)/L^3,            0,           0,           0, -(6*Iy)/L^2; ...
                             0,            0,  (12*Iz)/L^3,           0,  (6*Iz)/L^2,           0,    0,            0, -(12*Iz)/L^3,           0,  (6*Iz)/L^2,          0; ...
                             0,            0,            0,  (G/E)*It/L,           0,           0,    0,            0,            0, -(G/E)*It/L,           0,          0; ...
                             0,            0,   (6*Iz)/L^2,           0,    (4*Iz)/L,           0,    0,            0,  -(6*Iz)/L^2,           0,    (2*Iz)/L,          0; ...
                             0,  -(6*Iy)/L^2,            0,           0,           0,    (4*Iy)/L,    0,   (6*Iy)/L^2,            0,           0,           0,    (2*Iy)/L; ...
                          -A/L,            0,            0,           0,           0,           0,  A/L,            0,            0,           0,           0,          0; ...
                             0, -(12*Iy)/L^3,            0,           0,           0,  (6*Iy)/L^2,    0,  (12*Iy)/L^3,            0,           0,           0,  (6*Iy)/L^2; ...
                             0,            0, -(12*Iz)/L^3,           0, -(6*Iz)/L^2,           0,    0,            0,  (12*Iz)/L^3,           0, -(6*Iz)/L^2,          0; ...
                             0,            0,            0, -(G/E)*It/L,           0,           0,    0,            0,            0,  (G/E)*It/L,           0,          0; ...
                             0,            0,   (6*Iz)/L^2,           0,    (2*Iz)/L,           0,    0,            0,  -(6*Iz)/L^2,           0,    (4*Iz)/L,          0; ...
                             0,  -(6*Iy)/L^2,            0,           0,           0,    (2*Iy)/L,    0,   (6*Iy)/L^2,            0,           0,           0,   (4*Iy)/L];

                % Input matrix
                if at == 1
                    B = [1; 0; 0; 0; 0; 0; -1; 0; 0; 0; 0; 0];
                else
                    B = double.empty(12,0);
                end
                
                % Output of element length:
                Lout = [-1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0];

                % Delete z-, phi_x- and phi_y-directions for 2-dimension matrix:
                if dim == 2
                    M([3,4,5,9,10,11],:) = [];
                    M(:,[3,4,5,9,10,11]) = [];

                    K([3,4,5,9,10,11],:) = [];
                    K(:,[3,4,5,9,10,11]) = [];

                    B([3,4,5,9,10,11],:) = [];
                    Lout([3,4,5,9,10,11]) = [];
                end

        end
        
    case 2
        % Serial actuator element: Two nodes connected via a force employing actuator
        % Mass:
        M = zeros(6, 6);

        % Stiffness:
        K = zeros(6, 6);

        % Input matrix:
        B = [1; 0; 0; -1; 0; 0];

        % Output of element length:
        Lout = [-1, 0, 0, 1, 0, 0];

        % Delete z-direction for 2-dimension matrix:
        if dim == 2
            M([3,6],:) = [];
            M(:,[3,6]) = [];

            K([3,6],:) = [];
            K(:,[3,6]) = [];

            B([3,6]) = [];
            Lout([3,6]) = [];
        end
end



end

