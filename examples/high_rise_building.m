clear all;
close all;

E = 210e9;                          % Young's modulus (for steel)
rho = 7850;                         % Density
rho_dh = 7900;                      % Density horizontal diagonals
nu = 0.3;                           % Poisson ratio
G = E/(2*(1+nu));                 % Shear modulus

w_v     = 0.3;      % Vertical column cross sectional width
w_h     = 0.2;      % Horizontal bar cross sectional width
h_h     = 0.12;     % Horizontal bar cross sectional height
w_dh    = 0.01;
h_dh    = 0.06;
w_dv    = 0.012;
h_dv    = 0.15;
 
t_v1    = 0.01;     % Lower column wall thickness
t_v2    = 0.008;    % Upper column wall thickness
t_h     = 0.008;    % Horizontal bar wall thickness

h_floor = 3;
w_floor = 4.75;

L_v  = h_floor;           % Length of vertical columns
L_h  = w_floor;    
L_dv = sqrt((3*h_floor)^2 + w_floor^2);
L_dh = sqrt(2*w_floor^2);

% Cross sectional area
A_v = w_v^2 - (w_v - 2*t_v1)^2;
A_h  = 2*((w_h * h_h) - (w_h - 2*t_h)*(h_h - 2*t_h));
A_dh = w_dh * h_dh;
A_dv = w_dv * h_dv;

myu_v = A_v * rho;
myu_h  = A_h * rho; 
myu_dh = A_dh * rho_dh;
myu_dv = A_dv * rho;

% 2nd moment of area
I_v = 1/12 * (w_v^4 - (w_v - 2*t_v1)^4);

% Torsion constant
It_v = w_v*(w_v-t_v1)^3;

% Node table
n_stories = 12;
nodes = zeros(n_stories+1*4, 3);
for i = 0:n_stories
    nodes(i*4+1:(i+1)*4,:) = [0 0 i*h_floor; w_floor 0 i*h_floor; w_floor w_floor i*h_floor; 0 w_floor i*h_floor];
end
    
% Elements of a single module
elems_mod = [1 5 1; 2 6 1; 3 7 1; 4 8 1; ...    % Vertical columns
         5 9 1; 6 10 1; 7 11 1; 8 12 1; ...
         9 13 1; 10 14 1; 11 15 1; 12 16 1; ...
         5 6 2; 6 7 2; 7 8 2; 8 5 2; ...        % Horizontal bars
         9 10 2; 10 11 2; 11 12 2; 12 9 2; ...
         13 14 2; 14 15 2; 15 16 2; 16 13 2; ...
         5 7 3; 6 8 3; ...                      % Horizontal diagonals
         9 11 3; 10 12 3; ...
         13 15 3; 14 16 3; ... 
         1 14 4; 2 13 4; 2 15 4; 3 14 4; ...    % Vertical diagonals
         3 16 4; 4 15 4; 4 13 4; 1 16 4];       
% Full element table
elems = [elems_mod; ...
        [elems_mod(:, 1:2)+12 elems_mod(:, 3)]; ...
        [elems_mod(:, 1:2)+24 elems_mod(:, 3)]; ...
        [elems_mod(:, 1:2)+36 elems_mod(:, 3)]];
                     

% Element types
% Euler-Bernouli beam
column_v    = PH_FEM_Beam(4, 2, 2, myu_v, E, A_v, G, It_v, I_v, I_v, L_v);
% Rod elements
bar_h       = PH_FEM_Link(2, myu_h, E, A_h, L_h); 
bar_dh      = PH_FEM_Link(2, myu_dh, E, A_dh, L_dh);
diagonal    = PH_FEM_Link(2, myu_dv, E, A_dv, L_dv); 

% FEM mesh
mesh = PH_FEM_mesh(nodes, elems, {column_v, bar_h, bar_dh, diagonal});
% Lock DOFs of lowermost nodes
mesh.fixNodeDOFs(nodes(1,:), [1 1 1 0 0 1]);
mesh.fixNodeDOFs(nodes(2,:), [1 1 1 0 0 1]);
mesh.fixNodeDOFs(nodes(3,:), [1 1 1 0 0 1]);
mesh.fixNodeDOFs(nodes(4,:), [1 1 1 0 0 1]);

% Add external forces acting on global DOFs
mesh.addExternalInputsAtNodes();
% Generate constraints
mesh.generateConstraints();
% Assemble constraint matrix B
mesh.assembleDAESystem();
% Eliminate algebraic constraints and linearly dependent states
mesh.eliminateAlgebraicConstraints();
% Transform coordinate space to global DOFs
mesh.transformToGlobalDOFs();

% Get mass and stiffness matrix
M_ph = mesh.Q(1:mesh.n/2, 1:mesh.n/2)^-1;
K_ph = mesh.Q(mesh.n/2+1:end, mesh.n/2+1:end);

%% Conventional FE
load elem_types type_*
nodes_fe = [linspace(1, length(nodes), length(nodes))', nodes];
elems_fe = [linspace(1, length(elems_mod), length(elems_mod))', elems_mod, ones(length(elems_mod), 1)*type_beam, ones(length(elems_mod),1)*type_actparallel]; 
elems_fe(13:end, 5) = type_link;
elems_fe(13:end, 6) = type_actnone; 

elems_fe = [elems_fe; ...
            elems_fe(:,1:3)+12, elems_fe(:, 4:end); ... 
            elems_fe(:,1:3)+24, elems_fe(:, 4:end); ... 
            elems_fe(:,1:3)+36, elems_fe(:, 4:end)];

global_dof_order = nodes_fe(:,1)';

Re  = 355e6;          % yield strength of S355 steel
props = [1, E, rho     , A_v , I_v, I_v, G, It_v, Re; ... % vertical columns
         2, E, rho     , A_h , 0, 0, 0, 0, Re; ...  % horizontal bars
         3, E, rho_dh  , A_dh , 0, 0, 0, 0, Re; ... % horizonta diagonals
         4, E, rho     , A_dv , 0, 0, 0, 0, Re];
[M_fe, K_fe, ~, ~, ~] = truss3ext(nodes_fe, elems_fe, props, global_dof_order, 'srel');
M_fe = full(M_fe);
K_fe = full(K_fe);

fixed_dof = [1:3, 6, 7:9, 12, 13:15, 18, 19:21, 24];

M_fe(:, fixed_dof) = [];
K_fe(:, fixed_dof) = [];
M_fe(fixed_dof, :) = [];
K_fe(fixed_dof, :) = [];

%% Modal Analysis

% Enforce symmetry (not symmetric due to numerical errors)
M_ph = M_ph - (M_ph - M_ph')/2; 
K_ph = K_ph - (K_ph - K_ph')/2;

highestMode = 10;

opts.v0 = ones(size(M_ph, 1),1);
opts.isreal = 1;
[~,Lambda] = eigs(K_ph, M_ph, highestMode, 'smallestabs',opts); 
[omega,~] = sort(sqrt(diag(Lambda)));
omega_Hz_ph = omega/2/pi;

[~,Lambda] = eigs(K_fe, M_fe, highestMode, 'smallestabs',opts); 
[omega,~] = sort(sqrt(diag(Lambda)));
omega_Hz_fe = omega/2/pi;



%% Simulation

% Add Rayleigh damping
D_ph = M_ph * 0.05 + K_ph * 0.005;
J = mesh.J;
Q = mesh.Q; 
R = [D_ph, zeros(mesh.n/2);
    zeros(mesh.n/2), zeros(mesh.n/2)];

A = (J-R)*Q;

n_dofs = mesh.n/2;
D_fe = M_fe * 0.05 + K_fe * 0.005;
A_fe = [-M_fe\D_fe,  -M_fe\K_fe; ...
        eye(n_dofs), zeros(n_dofs)];


% Initial displacement due to wind in x-direction
wind_dofs_x = sort([(5:4:49).*6-5 (8:4:52).*6-5]-24+8);
E_wind = zeros(n_dofs, 1);
for i=1:12
    E_wind(wind_dofs_x(2*i-1), 1) = i/12;
    E_wind(wind_dofs_x(2*i), 1) = i/12;
end
E_wind_fe = [M_fe\E_wind; zeros(mesh.n/2, 1)];
E_wind = [E_wind; zeros(mesh.n/2, 1)];
x0_ph = -A\E_wind*0.2e5;
x0_fe = -A_fe\E_wind_fe*0.2e5;

time = 1:0.01:5; 
odefun_ph = @(t, x) A*x;
odefun_fe = @(t, x) A_fe*x;
[~, y_ph] = ode15s(odefun_ph, time, x0_ph);
%[~, y_fe] = ode15s(odefun_fe, time, x0_fe);

%% Plots
c_x = 1:6:6*13*4;
c_x(1:4) = [];
c_x = c_x - 16;

x_ph = y_ph(:, n_dofs+1:end);
x_ph = x_ph(:, c_x);

figure();
plot(time, 1e3.*x_ph(:, end));
xlabel('time in s');
ylabel('x-displacement in mm');
title('Response to an initial displacement');

figure();
H = zeros(length(time), 1);
for k=1:length(time)
    H(k) = 0.5*y_ph(k,:)*mesh.Q*y_ph(k,:)';
end
plot(time, 1e-3.*H);
xlabel('time in s');
ylabel('energy in kJ');
title('System Hamiltonian');