clear all;
close all;

E = 210e9;                          % Young's modulus (for steel)
rho = 7850;                         % Density
rho_dh = 7900;                      % Density horizontal diagonals
nu = 0.3;                           % Poisson ratio
G = 1/(2*(1+nu))*E;                 % Shear modulus

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

% Cross sectional areas
A_v1 = w_v^2 - (w_v - 2*t_v1)^2;
A_v2 = w_v^2 - (w_v - 2*t_v2)^2;
A_h  = 2*((w_h * h_h) - (w_h - 2*t_h)*(h_h - 2*t_h));
A_dh = w_dh * h_dh;
A_dv = w_dv * h_dv;

myu_v1 = A_v1 * rho;
myu_v2 = A_v2 * rho;
myu_h  = A_h * rho; 
myu_dh = A_dh * rho_dh;
myu_dv = A_dv * rho;

% 2nd moment of area
I_v1 = 1/12 * (w_v^4 - (w_v - 2*t_v1)^4);
I_v2 = 1/12 * (w_v^4 - (w_v - 2*t_v2)^4);

It_v1 = w_v*(w_v-t_v1)^3;
It_v2 = w_v*(w_v-t_v2)^3;

kappa = 5/6;

% node and element table
nodes = [0 0 0; w_floor 0 0; w_floor w_floor 0; 0 w_floor 0; ...
         0 0 h_floor; w_floor 0 h_floor; w_floor w_floor h_floor; 0 w_floor h_floor; ...
         0 0 2*h_floor; w_floor 0 2*h_floor; w_floor w_floor 2*h_floor; 0 w_floor 2*h_floor; ...
         0 0 3*h_floor; w_floor 0 3*h_floor; w_floor w_floor 3*h_floor; 0 w_floor 3*h_floor];
elems = [1 5 1; 2 6 1; 3 7 1; 4 8 1; ...        % Vertical columns
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

     
% Element types
% Timoshenko beam
%column_v    = PH_FEM_Beam(3, 2, 2, myu_v1, E, A_v1, G, It_v1, I_v1, I_v1, L_v, kappa);
% Bernouli beam
column_v    = PH_FEM_Beam(4, 2, 2, myu_v1, E, A_v1, G, It_v1, I_v1, I_v1, L_v);

bar_h       = PH_FEM_Link(2, myu_h, E, A_h, L_h, 'free-free'); 
bar_dh      = PH_FEM_Link(2, myu_dh, E, A_dh, L_dh, 'free-free');
diagonal    = PH_FEM_Link(2, myu_dv, E, A_dv, L_dv, 'free-free'); 

% FEM-Mesh
mesh = PH_FEM_mesh(nodes, elems, {column_v, bar_h, bar_dh, diagonal});
% Lock ground node DOFs
mesh.fixNodeDOFs(1, [1 1 1 1 1 1]);
mesh.fixNodeDOFs(2, [1 1 1 1 1 1]);
mesh.fixNodeDOFs(3, [1 1 1 1 1 1]);
mesh.fixNodeDOFs(4, [1 1 1 1 1 1]);

% Assemble systems
mesh.generateConstraints()
mesh.assemble();
mesh.eliminateAlgebraicConstraints();


%% Simulation
% Single input force in x-direction on top node
B = zeros(mesh.n_d, 1);
B(49) = 1; 

% Time
t = 0:0.01:3;

% Setup odefun 0for the pH-system
A = sparse(mesh.A);
ph_odefun = @(t, x) A*[x; B*(50000/3*t)];

% Simulate the excited system
[time, x_ph] = ode15s(ph_odefun, t, zeros(size(mesh.A,1),1)); 

%% Output
c_x = 1:6:67;
c_y = 2:6:68;
c_z = 3:6:69;
z_ph = zeros(length(t), mesh.n_d);
for k=1:length(t)
    z_ph(k, :) = mesh.getInteractionOutput(x_ph(k,:)');
end
% Integrate to get positions from velocities
z_ph = cumtrapz(z_ph).*0.01;
figure();
plot(time, z_ph(:, c_x));

% legend entries will be velocities instead of positions
legend(mesh.externalOutputNames{c_x});