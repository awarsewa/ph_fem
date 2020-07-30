clear
close all

%% Passive frame
E = 210e9;                          % Young's modulus (for steel)
rho = 7850;                         % Density
nu = 0.3;                           % Poisson ratio

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
% Eliminate algebraic constraints 
mesh.eliminateAlgebraicConstraints();
% Eliminate linearly dependent states/efforts
mesh.eliminateLinearDependencies();
% Transform coordinate space to global DOFs
mesh.transformToGlobalDOFs();
% And add damping
mesh.addRayleighDamping(0.05, 0.005);


%% Hydraulic actuator
beta = 1.3e9;      % Fluid bulk modulus in Pa
A1 =  pi*(130e-3/2)^2;                     % Cross-sectional area chamber
A2 =  pi*(130e-3/2)^2 - pi*(70e-3/2)^2;    % Cross-sectional area ring
L = 350e-3;        % Cylinder length
m_p = 10;          % Piston mass
m_c = 20;          % Cylinder mass
p_s = 200e5;       % Supply pressure 
p_t = 0;           % Tank pressure
kv = 2/60000;      % Valve coefficient (m^3/s)
C_int = 0;         % Internal leakage flow coefficient
C_ext = 0;         % Externa leakage flow coefficient
b_p = 0;            % Linear damping coefficient piston
b_c = b_p;         % Linear damping coefficient cylinder

% PH_HydraulicActuator(beta, m, L, A1, A2, kv, p_s, p_t)
cylinder = PH_HydraulicActuator(beta, m_p, m_c, b_p, b_c, L, A1, A2, kv, p_s, p_t, C_int, C_ext);
cylinder.setPosition([nodes(1,:); nodes(7,:)]);


%% System interconnection
frame = PH_NonlinearSystem.toNonlinearPHSystem(mesh);
% frame with cylinder
proto = frame.copy(); 
proto.add(cylinder);
proto.fixNodeDOFs([0 0 0], [1 0 1 0 1 0]);
proto.fixNodeDOFs([4.75 0 0], [1 0 1 0 1 0]);
proto.addExternalInputsAtNodes();
proto.generateConstraints();
proto.assembleDAESystem();


%% Compilation
% Precompile functions to achieve ~ x2 speedup
% If you change any of the parameters, you need to recompile
if ~exist('prototype.mexw64', 'file')
    [f, df] = proto.getSystemODEfun();
    
    opts = struct('mex', true, 'verbose', false);
    C = casadi.CodeGenerator('prototype.c', opts);
    C.add(f);
    C.add(df);
    C.generate();
    % Takes several minutes...
    mex prototype.c -LargeArrayDims;
end

%% Simulation
x0 = [zeros(mesh.n, 1); 0.1; 0; 0; 0; 1e5; 1e5];
dt = 0.01;
time = 0:dt:2;
u = zeros(proto.n_u, length(time));
u(2,:) = 0.5e-3*sin(2*pi*time);
u(3, u(2,:) < 0) = u(2, u(2,:) < 0);
u(2, u(2,:) < 0) = 0;
  
x_ph = zeros(proto.n, length(time));
h_ph = zeros(1, length(time));
warning('off', 'MATLAB:nearlySingularMatrix');
for k=2:length(time)
    if mod(k, 10) == 0
        disp([num2str(k/100) 's elapsed...']);
    end
    x_ph(:, k-1:k) = gls(@(t, x) full(prototype('f',x, u(:, k-1))), @(t, x) full(prototype('dfdx', x, u(:,k-1))), time(k-1:k), x0, 1)'; 
    x0 = x_ph(:, k);
end
warning('on', 'MATLAB:nearlySingularMatrix');

%% Plots
close all

% Integrate to get nodal displacements
y = proto.getSystemOutputFun();
y_ph = zeros(proto.n_u, length(time)); 
for k=1:length(time)
    y_ph(:, k) = full(y(x_ph(:,k))); 
end

d = cumtrapz(y_ph([2 3 5 6 8 9 11 12 14 15 17 18]+3, :), 2)*dt;
figure;
ch = plot(time, d([9 10 11 12], :)*1e3);
xlabel('time in s');
ylabel('displacement in mm'); 
legend('n7_x', 'n7_z', 'n8_x', 'n8_z');
ax = gca;
ax.Box = 'off';
ax.FontName = 'Arial';
ax.FontSize = 12;
ax.LineWidth = 1;
ch(1).LineWidth = 1;
ch(2).LineWidth = 1;
ch(3).LineWidth = 1;
ch(3).LineStyle = '--';
ch(4).LineWidth = 1;
ch(4).LineStyle = '--';

figure
ax = plotyy(time, 1e3*x_ph(37,:)-100, time, 1e-5*x_ph(41:42, :));
ch1 = ax(1).Children;
ch2 = ax(2).Children;
ylabel('piston displacement in mm');
ax(1).Box = 'off';
ax(1).FontName = 'Arial';
ax(1).FontSize = 12;
ax(1).LineWidth = 1;
xlabel(ax(1), 'time in s');
ylabel(ax(1), 'piston displacement in mm');
ylabel(ax(2),'pressure in bar');
ax(2).YAxis.Color = 'b';
ax(1).YAxis.Color = 'k';
ax(1).YTickMode = 'auto';
ax(2).YTickMode = 'auto';
legend('\Delta s', 'p_1', 'p_2');
ch1(1).LineWidth = 1;
ch2(1).LineWidth = 1;
ch2(2).LineWidth = 1;
ch1(1).LineStyle = '--';