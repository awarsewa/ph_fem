clear all;

% Parameters
E = 210e9;                          % Young's modulus (for steel)
rho = 7850;                         % Density
nu = 0.3;                           % Poisson ratio
G = 1/(2*(1+nu))*E;                 % Shear modulus
kappa = 5/6;

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

% element types
column      = PH_FEM_Timoshenko(3, rho, E, G, A_v, I_v, kappa, L_v);
%column      = PH_FEM_Bernoulli(4, mu_v, E, I_v, L_v, 'free-free');
column_rod  = PH_FEM_Link(2, mu_v, E, A_v, L_v, 'free-free'); 
bar         = PH_FEM_Link(2, mu_h, E, A_h, L_h, 'free-free');
diagonal    = PH_FEM_Link(2, mu_d, E, A_d, L_d, 'free-free'); 

% FEM-Mesh
mesh = PH_FEM_mesh(nodes, elems, {column, bar, diagonal, column_rod});
mesh.fixNodeDOFs(1, [1 0 1 0 1 0]);
mesh.fixNodeDOFs(2, [1 0 1 0 1 0]);

% Assembly of concatenated system
mesh.generateConstraints();
mesh.assemble();
mesh.eliminateAlgebraicConstraints();


%% Simulation
t = 0:0.01:3;

% Setup state derivative function for the pH-system
ph_odefun = @(t, x) mesh.A*[x; [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0]'*(50000/3*t)];

% Simulate the excited system
[time, x_ph] = ode45(ph_odefun, t, zeros(size(mesh.A,1),1)); 

%% Output
y_ph = zeros(length(t), mesh.n_d);
for k=1:length(t)
    y_ph(k, :) = mesh.getInteractionOutput(x_ph(k,:)');
end
y_ph = cumtrapz(y_ph).*0.01;
figure();
plot(time, y_ph(:, 13:18));
legend(mesh.externalOutputNames{13:18});
