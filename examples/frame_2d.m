close all

E = 210e9;                          % Young's modulus (for steel)
rho = 7850;                         % Density
nu = 0.3;                           % Poisson ratio
G = 1/(2*(1+nu))*E;                 % Shear modulus
kappa = 0.099;

A_v = 0.3^2 - (0.3 - 0.02)^2; 
A_h = 0.2 * 0.12 - (0.2 - 0.016) * (0.12 - 0.016);
A_h = A_h * 2;
A_d = 0.15*0.015;
mu_v = A_v * rho;
mu_h = A_h * rho;
mu_d = A_d * rho;

I_v = 1/12 * (0.3^4-(0.3-0.02)^4);    % 2nd moment of area
I_h = 1/12 * (0.12*0.2^3 - (0.12 - 0.016)*(0.2 -0.016)^3);
I_d = 1/12 * A_d * 0.15^2;
L_v = 3;                              % Beam height
L_h = 4.75;
L_d = sqrt(9^2 + 4.75^2);

% node and element table
nodes = [0 0 0; 4.75 0 0; 0 0 3; 4.75 0 3; ...
         0 0 6; 4.75 0 6; 0 0 9; 4.75 0 9];
elems = [1 3 1; 3 5 1; 5 7 1; 2 4 1; 4 6 1; 6 8 1; ... 
         1 3 4; 3 5 4; 5 7 4; 2 4 4; 4 6 4; 6 8 4; ...
         3 4 2; 5 6 2; 7 8 2; ...
         1 8 3; 2 7 3]; 
attribs = cell(size(elems, 1), 1);
for e=1:size(elems, 1)
    if elems(e, 3) == 1 || elems(e, 3) == 6
        % bending axis is the x-axis
        attribs{e} = [0 1 0];
    end
    if elems(e, 3) == 5
        attribs{e} = [0 1 0];
    end
end   

% Element types
column      = PH_FEM_Bernoulli(4, mu_v, E, I_v, L_v);
column_rod  = PH_FEM_Link(2, mu_v, E, A_v, L_v); 
bar         = PH_FEM_Link(2, mu_h, E, A_h, L_h);
diagonal    = PH_FEM_Link(2, mu_d, E, A_d, L_d);

% FEM-Mesh
mesh = PH_FEM_mesh(nodes, elems, {column, bar, diagonal, column_rod}, attribs);
mesh.fixNodeDOFs([0 0 0], [1 0 1 0 1 0]);
mesh.fixNodeDOFs([4.75 0 0], [1 0 1 0 1 0]);

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

%% Simulation
t = 0:0.01:3;

B = zeros(mesh.n_u, 1);
B(16) = 1;
% Setup state derivative function for the pH-system
ph_odefun = @(t, x) mesh.A*[x; B*(50000/3*t)];

% Simulate the excited system
[time, x_ph] = ode45(ph_odefun, t, zeros(size(mesh.A,1),1)); 

%% Output
y_ph = zeros(length(t), mesh.n_u);
for k=1:length(t)
    % Get global nodal velocities
    y_ph(k, :) = mesh.getSystemOutput(x_ph(k,:)', B*(50000/3*t(k)));
end
% Integrate to get nodal displacements
y_ph = cumtrapz(y_ph).*0.01;
figure;
% Plot upper 2 nodes only
plot(time, y_ph(:, 13:18));
xlabel('time in s');
ylabel('displacement'); 
title('displacements by output integration');

% Displacements can also be obtained from the state vector 
figure;
plot(time, x_ph(:, 13+mesh.n/2:18+mesh.n/2));
xlabel('time in s');
ylabel('displacement'); 
title('displacements from state vector');