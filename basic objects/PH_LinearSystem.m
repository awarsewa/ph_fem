classdef PH_LinearSystem < PH_System
    %LINEARPHSYSTEM Linear port-Hamiltonian System in input-state-output
    %form
    %   Detailed explanation goes here
    
    properties(SetAccess = protected, GetAccess = public)
        J       % Structure matrix J
        R       % Resistance matrix R (dissipation)
        Q       % Energy matrix Q
        G       % Input matrix symmetric part G
        P       % Input matrix skew symmetric part P
        K       % Interaction matrix K
        M       % Feedthrough matrix skew symmetric part M
        S       % Feedthrough matrix symmetric part S (dissipation)

        C_u     % Matrix for I/O-coupling 
        C_y     % Matrix for I/O-coupling
        C_d     % Matrix for coupling with external forces
        
        B       % Algebraic constraints on efforts
        A       % System matrix for simulation
        E       % Jacobian matrix for simulation
    end
    
    properties(SetAccess = public, GetAccess = public)
        inputSigns
        outputSigns
    end 
    
    methods(Access = public)
        %              PH_LinearSystem(n, E, J, Q, G, R, K, P, S, M, B)
        function obj = PH_LinearSystem(n, E, J, Q, varargin)
            % Sanity checks
            obj = obj@PH_System('Linear port-Hamiltonian system');
           
            % Strict type checking
            if ~isscalar(n) || n < 0
                error('n must be >= 0');
            end
            obj.n = n;
            
            
            if isempty(E)
                obj.E = eye(n);
            elseif ~ismatrix(E) || any(size(E) ~= n) 
                error('E must be a nxn matrix');
            else
                obj.E = E;
            end

            
            if ~ismatrix(J) || any(size(J) ~= n) || ~issymmetric(J, 'skew')
                error('J must be a skew-symmetric nxn matrix');
            end
            obj.J = J;
            
            if ~ismatrix(Q) || any(size(Q) ~= n)
                error('Q must be a nxn matrix');
            end
            obj.Q = Q;
            
            % Inputs and Outputs
            obj.n_u = 0; 
            obj.G = zeros(obj.n, 0);
            if nargin > 4 && ~isempty(varargin{1})
                G = varargin{1};
                if ~ismatrix(G) || size(G, 1) ~= n
                    error('G must be a matrix with n rows');
                end            
                obj.n_u = size(G, 2);
                obj.G = G; 
            end
            
            % Dissipation
            obj.R = zeros(n);
            if nargin > 5
                R = varargin{2};
                if ~isempty(R)
                    if ~ismatrix(R) || any(size(R) ~= n) || ~issymmetric(R) || any(round(eig(R),17) < 0)
                        error('R must be a symmetric positive semidefinite nxn matrix');
                    end
                    obj.R = R;
                end
            end
                       
            % Interaction
            obj.n_d = 0;
            obj.K = zeros(obj.n, 0);
            if nargin > 6
                K = varargin{3};
                if ~isempty(K)
                    if ~ismatrix(K) || size(K, 1) ~= n
                        error('K must be a matrix with n rows');
                    end
                    obj.K = K;
                    obj.n_d = size(K, 2);
                end
            end
            
            % System with direct input feedthrough
            obj.P = zeros(obj.n, obj.n_u); 
            obj.M = zeros(obj.n_u);
            obj.S = zeros(obj.n_u);
            if nargin > 9
                P = varargin{4};
                S = varargin{5};
                M = varargin{6};
                
                if ~isempty(P) && ~isempty(S) && ~isempty(M)
                    if ~ismatrix(P) || size(P, 1) ~= obj.n || size(P, 2) ~= obj.n_u
                        error('P must be a n x n_u matrix');
                    end
                    if ~ismatrix(S) || size(S, 1) ~= obj.n_u || size(S, 2) ~= obj.n_u || ~issymmetric(S)
                        error('S must be a symmetric n_u x n_u matrix');
                    end

                    if ~ismatrix(M) || size(M, 1) ~= obj.n_u || size(M, 2) ~= obj.n_u || ~issymmetric(M, 'skew')
                        error('M must be a skew-symmetric n_u x n_u matrix');
                    end

                    Z = [obj.R P; P' S];
                    if any(round(eig(Z),17) < 0)
                        error('P,S and R must satisfy Z = [R P; P'' S] is positive semidefinite');
                    end

                    obj.P = P;
                    obj.S = S;
                    obj.M = M;
                end
            end
            
            obj.C_u = zeros(0, obj.n_u);
            obj.C_y = zeros(0, obj.n_u);
            obj.C_d = zeros(0, obj.n_d);
            
            obj.B = zeros(obj.n, 0);
            if nargin > 10 
                B = varargin{7};
                if size(B, 1) ~= obj.n 
                    error('B must be a n x n_c matrix');
                end
                
                obj.B = B;
            end
            obj.n_c = size(obj.B, 2);
            
            obj.A = [(obj.J - obj.R)*obj.Q, obj.G - obj.P, obj.K;
                    obj.B', zeros(obj.n_c, obj.n_u + obj.n_d)];
            obj.E = [obj.E, zeros(obj.n, obj.n_u+obj.n_d);
                    zeros(obj.n_c, obj.n+obj.n_u+obj.n_d)];
                
            obj.n_elements = 1;
            obj.elements{1} = PH_Element('Linear PH system', []);
            
            obj.inputSigns = ones(obj.n_u, 1);
            obj.outputSigns = obj.inputSigns;
            for i=1:obj.n_u
                obj.inputNames{i} = ['u' num2str(i)];
                obj.outputNames{i} = ['y' num2str(i)];
            end
            
            for i=1:obj.n_d
                obj.externalInputNames{i} = ['d' num2str(i)];
                obj.externalOutputNames{i} = ['z' num2str(i)];
            end
        end
        
        function addCouplingConstraints(obj, C_u, C_y, varargin)   
            if ~isempty(C_u)
                if ~ismatrix(C_u) || size(C_u, 2) ~= obj.n_u 
                    error(['C_u must be a matrix with ', num2str(obj.n_u) ' columns']);
                end
                if ~ismatrix(C_y) || size(C_y, 2) ~= obj.n_u 
                    error(['C_y must be a matrix with ', num2str(obj.n_u) ' columns']);
                end
                if ~(size(C_u, 1) == size(C_y, 1))
                    error('C_u and C_y must be matrices of equal size!'); 
                end
                
               
                obj.C_u = [obj.C_u; C_u];
                obj.C_y = [obj.C_y; C_y];
            end
            
            if nargin > 3 
                Cd = varargin{1};
                if ~isempty(Cd)
                    if ~ismatrix(Cd) || size(Cd, 1) ~= size(C_u, 1) || size(Cd, 2) ~= obj.n_d
                            error('C_d must be n_c x n_d!'); 
                    end
                    obj.C_d = [obj.C_d; Cd];
                else
                    obj.C_d = [obj.C_d; zeros(size(C_u, 1), obj.n_d)];
                end
            else 
                obj.C_d = [obj.C_d; zeros(size(C_u, 1), obj.n_d)];
            end 
            obj.n_c = size(obj.C_u,1) + size(obj.B, 2);
        end
        
        % Add another linear PH system 
        function add(obj, system)
            if ~isa(system, 'PH_LinearSystem') || ~system.isValidPHSystem()
                error('system must be a valid ''PH_LinearSystem''');
            end
            
            % System does not have two equally sized energy domains?
            % This needs to be checked more thoroughly..
            if mod(system.n, 2)
                obj.add_single(system);
                return;
            end
            
            N_b2 = system.n/2; 
            N_a = obj.n;            
                                  
            Ea = obj.E(1:obj.n, 1:obj.n);
            Eb = system.E(1:system.n, 1:system.n);
            
            obj.E = [Eb(1:N_b2,1:N_b2) zeros(N_b2, N_a) Eb(1:N_b2, N_b2+1:end);
                  zeros(N_a, N_b2) Ea zeros(N_a, N_b2);
                  Eb(N_b2+1:end, 1:N_b2) zeros(N_b2, N_a) Eb(N_b2+1:end, N_b2+1:end)];
            obj.J = [system.J(1:N_b2,1:N_b2) zeros(N_b2, N_a) system.J(1:N_b2, N_b2+1:end);
                  zeros(N_a, N_b2) obj.J zeros(N_a, N_b2);
                  system.J(N_b2+1:end, 1:N_b2) zeros(N_b2, N_a) system.J(N_b2+1:end, N_b2+1:end)];
            obj.Q = [system.Q(1:N_b2,1:N_b2) zeros(N_b2, N_a) system.Q(1:N_b2, N_b2+1:end);
                  zeros(N_a, N_b2) obj.Q zeros(N_a, N_b2);
                  system.Q(N_b2+1:end, 1:N_b2) zeros(N_b2, N_a) system.Q(N_b2+1:end, N_b2+1:end)];
            
            obj.G = [zeros(N_b2, obj.n_u), system.G(1:N_b2,:);
                  obj.G, zeros(N_a, system.n_u);
                  zeros(N_b2, obj.n_u)  system.G(N_b2+1:end, :)];
            obj.P = [zeros(N_b2, obj.n_u), system.P(1:N_b2,:);
                  obj.P, zeros(N_a, system.n_u);
                  zeros(N_b2, obj.n_u)  system.P(N_b2+1:end, :)];
            obj.B = [zeros(N_b2, size(obj.B, 2)), system.B(1:N_b2,:);
                  obj.B, zeros(N_a, size(system.B,2));
                  zeros(N_b2, size(obj.B, 2))  system.B(N_b2+1:end, :)];
            obj.K = [zeros(N_b2, obj.n_d), system.K(1:N_b2,:);
                  obj.K, zeros(N_a, system.n_d);
                  zeros(N_b2, obj.n_d)  system.K(N_b2+1:end, :)];
              
            obj.M = blkdiag(obj.M, system.M);
            obj.S = blkdiag(obj.S, system.S); 
            
            if ~isempty(system.R)
                obj.R = [system.R(1:N_b2,1:N_b2) zeros(N_b2, N_a) system.R(1:N_b2, N_b2+1:end);
                      zeros(N_a, N_b2) obj.R zeros(N_a, N_b2);
                      system.R(N_b2+1:end, 1:N_b2) zeros(N_b2, N_a) system.R(N_b2+1:end, N_b2+1:end)];
            else
                obj.R = [zeros(N_b2, 2*N_b2+N_a);
                      zeros(N_a, N_b2) obj.R zeros(N_a, N_b2);
                      zeros(N_b2, 2*N_b2+N_a)];
            end
            
            obj.C_u = [obj.C_u, zeros(size(obj.C_u, 1), system.n_u); ...
                       zeros(size(system.C_u, 1), obj.n_u) system.C_u];
            obj.C_y = [obj.C_y, zeros(size(obj.C_y, 1), system.n_u); ...
                       zeros(size(system.C_y, 1), obj.n_u) system.C_y];
            obj.C_d = [obj.C_d, zeros(size(obj.C_d, 1), system.n_d); ...
                       zeros(size(system.C_d, 1), obj.n_d) system.C_d];
            obj.n_c = size(obj.C_u, 1) + size(obj.B, 2);   
            
            for p = 1:system.n_ports
                obj.ports{end+1} = copy(system.ports{p});
                obj.ports{end}.node = obj.ports{end}.node + obj.n_nodes; 
                obj.ports{end}.IOPair = obj.ports{end}.IOPair + obj.n_u;
            end
            
            for n = 1:system.n_nodes
                obj.nodes{end+1} = copy(system.nodes{n});
                obj.nodes{end}.ports = obj.nodes{end}.ports + obj.n_ports;
            end
            
            for e = 1:system.n_elements
                obj.elements{end+1} = copy(system.elements{e});
                obj.elements{end}.ports = obj.elements{end}.ports + obj.n_ports;
            end
            
            obj.n_ports = obj.n_ports + system.n_ports;
            obj.n_nodes = obj.n_nodes + system.n_nodes;
            obj.n_elements = obj.n_elements + system.n_elements; 
            obj.n = obj.n + system.n;
            obj.n_u = obj.n_u + system.n_u;
            obj.n_d = obj.n_d + system.n_d;
            
            obj.A = [(obj.J - obj.R)*obj.Q, obj.G - obj.P, obj.K;
                     obj.B', zeros(obj.n_c, obj.n_u + obj.n_d)];
            obj.E = [obj.E, zeros(obj.n, obj.n_u+obj.n_d);
                    zeros(obj.n_c, obj.n+obj.n_u+obj.n_d)];
            
            
            for i = 1:system.n_u
                obj.inputNames{end+1} = system.inputNames{i};
                obj.outputNames{end+1} = system.outputNames{i};
                obj.inputSigns(end+1) = system.inputSigns(i);
                obj.outputSigns(end+1) = system.outputSigns(i);
            end
            
            for i = 1:system.n_d
                obj.externalInputNames{end+1} = system.externalInputNames{i};
                obj.externalOutputNames{end+1} = system.externalOutputNames{i};
            end
        end
                
        % Generate constraints of the form Cy * y + Cu * u + Cd * d = 0.
        % Constraints are generated automatically for each node. Only the
        % mechanical domain is supported so far.
        function generateConstraints(obj)
            % Mechanical constraint generation
            [Cu_v, Cy_v, Cd_v] = obj.generateVelocityConstraints();           
            [Cu_f, Cy_f, Cd_f] = obj.generateForceConstraints();
            
            % Ensure matrices have the same size
            if size(Cd_v, 2) < obj.n_d
                Cd_v(1, obj.n_d) = 0; 
            end
            if size(Cd_f, 2) < obj.n_d
                Cd_f(1, obj.n_d) = 0; 
            end
            
            
            % Add constraint matrices
            obj.C_u = [Cu_v; Cu_f];
            obj.C_y = [Cy_v; Cy_f];
            obj.C_d = [Cd_v; Cd_f];
            
            % Adjust the number of constraints
            obj.n_c = size(obj.C_u, 1) + size(obj.B, 2); 
            
            % Constraints were generated successfully. To obtain an ODE system,
            % call assemble() followed by eliminateAlgebraicConstraints()
        end
        
        % Eliminates the Cu, Cy, Cd constraints by solving for u.
        % In the process, a set of algebraic constraints of the form B'*e = 0
        % are usually generated. They can be eliminated by a call to
        % eliminateAlgebraicConstraints
        function assemble(obj)
            if obj.n_c == 0
                % No constraints --> nothing to do
                return
            end
            
            % Ignore external inputs mapping to outputs..
            % Find single y and d only constraints
            idx_yd = any(obj.C_y, 2) & any(obj.C_d, 2) & ~any(obj.C_u, 2);
            idx_d = any(obj.C_d(idx_yd, :), 1);
            idx_yd = find(sum(double(obj.C_y(idx_yd, :) ~= 0), 2) == 1);

            % Those can be deleted...
            obj.C_u(idx_yd, :) = [];
            obj.C_y(idx_yd, :) = [];
            obj.C_d(idx_yd, :) = [];
            obj.C_d(:, idx_d) = [];
            obj.n_d = size(obj.C_d, 2);
            
            % Get algebraic constraints
            idx_ac = ~any(obj.C_u') & ~any(obj.C_d');
            B_t = (obj.outputSigns.*obj.C_y(idx_ac, :)) * (obj.G+obj.P)';         
            if ~isempty(obj.B)
                    obj.B = [obj.B; rref(B_t)']; 
            else
                obj.B = rref(B_t)';
            end
            obj.B(:, ~any(obj.B)) = [];    
            obj.C_y(idx_ac, :) = [];
            obj.C_d(idx_ac, :) = [];
            obj.C_u(idx_ac, :) = [];
                        
            C = rref([obj.C_u, obj.C_y, obj.C_d]);
            obj.C_u = C(:, 1:obj.n_u);
            obj.C_y = C(:, obj.n_u+1:2*obj.n_u);
            obj.C_d = C(:, 2*obj.n_u+1:end);
            
            % Solve remaining equations for u
            u_f = ones(1, obj.n_u);
            A_uu = eye(obj.n_u);
            A_ue = zeros(obj.n_u, obj.n);
            A_ud = zeros(obj.n_u, obj.n_d);
            for i = 1:size(obj.C_u, 1)
                idx_u = find(obj.C_u(i, :), 1);
                if ~isempty(idx_u)
                    % Constrain this input
                    % TODO: make sure to ignore contributions in orthogonal
                    % directions... (this shouldn't be done here)
                    u_f(idx_u) = 0;
                    a_uu = obj.C_u(i, :)*A_uu + obj.C_y(i,:)*(obj.M+obj.S)*A_uu;
                    b = a_uu(idx_u);
                    a_uu(idx_u) = 0;
                    A_uu(idx_u, :) = -b .* a_uu;
                    A_ue(idx_u, :) = -b .* obj.C_y(i,:)*(obj.G+obj.P)';
                    A_ud(idx_u, :) = -b .* obj.C_d(i,:);

                    % Update dependent inputs
                    idx_uc = find(A_uu(:, idx_u));
                    for k=1:length(idx_uc)
                        b = A_uu(idx_uc(k), idx_u);
                        A_uu(idx_uc(k), idx_u) = 0;
                        A_uu(idx_uc(k), :) = A_uu(idx_uc(k), :) + b*A_uu(idx_u, :);
                        A_ue(idx_uc(k), :) = A_ue(idx_uc(k), :) + b*A_ue(idx_u, :);
                        A_ud(idx_uc(k), :) = A_ud(idx_uc(k), :) + b*A_ud(idx_u, :);
                    end

                    % Update constraints
                    b = obj.C_u(i, idx_u);
                    c_uu = -b*obj.C_u(i, :);
                    c_uu(idx_u) = 0;
                    c_uy = -b*obj.C_y(i, :);
                    c_ud = -b*obj.C_d(i, :);
                    idx_uc = find(obj.C_u(i+1:end, idx_u)) + i;
                    for k=1:length(idx_uc)
                        b = obj.C_u(idx_uc(k), idx_u);
                        obj.C_u(idx_uc(k), idx_u) = 0;
                        obj.C_u(idx_uc(k), :) = obj.C_u(idx_uc(k), :) + b * c_uu;
                        obj.C_y(idx_uc(k), :) = obj.C_y(idx_uc(k), :) + b * c_uy;
                        obj.C_d(idx_uc(k), :) = obj.C_d(idx_uc(k), :) + b * c_ud;
                    end
                else
                    warning(['ignoring constraint ' num2str(i)]);
                end
            end
            A_uu(:, ~any(A_uu)) = [];       
            u_new = find(u_f); 
            
            J_ = obj.J + (obj.G - obj.P) * A_ue;
            R_ = obj.R; 
            G_ = (obj.G - obj.P) * A_uu; 
            K_ = (obj.G - obj.P) * A_ud;
            
            P_ = zeros(obj.n, size(G_, 2));
            D_ = A_uu' * (obj.M + obj.S) * A_uu;
            if ~issymmetric(D_, 'skew')
                error('Splitting D into S and M needs to be implemented');
            else
                M_ = D_; 
                S_ = zeros(size(M_, 1), size(M_, 2));
            end
            
            % Alter the decimal precision in case of numerical problems
            decimalPrecision = 14;
            
            % Assign new system matrices
            obj.J = round(J_, decimalPrecision);
            obj.R = round(R_, decimalPrecision);
            obj.G = round(G_, decimalPrecision);
            obj.P = round(P_, decimalPrecision); 
            obj.M = round(M_, decimalPrecision); 
            obj.S = round(S_, decimalPrecision);
            obj.K = round(K_, decimalPrecision); 
            
            % Check for zero colums in G (inputs that are zero) and remove
            % them
            idx_u0 = ~any(obj.G-obj.P);
            obj.G(:, idx_u0) = [];
            obj.P(:, idx_u0) = [];
            obj.M(:, idx_u0) = [];
            obj.M(idx_u0, :) = [];
            obj.S(:, idx_u0) = [];
            obj.S(idx_u0, :) = [];
            
            u_new(idx_u0) = [];
            u_new_idx = u_new;
            
            obj.n_u = size(obj.G, 2);           
            obj.n_d = size(obj.K, 2);
            
            % Update input and output names
            inputNames_ = cell(length(u_new), 1); 
            outputNames_ = cell(length(u_new), 1); 
            for i=1:obj.n_u
                inputNames_{i} = obj.inputNames{u_new_idx(i)};
                outputNames_{i} = obj.outputNames{u_new_idx(i)};
            end         
            obj.inputNames = inputNames_;
            obj.outputNames = outputNames_;
            obj.inputSigns = obj.inputSigns(u_new_idx);
            obj.outputSigns = obj.outputSigns(u_new_idx);
            
            % System check... 
            if ~obj.isValidPHSystem()
                error('Something is messed up... check your constraints');
            end
            
            % All I/O coupling constraints are now eliminated
            obj.C_u = [];
            obj.C_y = [];
            obj.C_d = [];
            
            % Some algebraic constraints might remain
            % Call eliminateAlgebraicConstraints() to eliminate them
            % This also eliminates constraint forces (inputs)
            obj.n_c = size(obj.B, 2);
            
            % Apply changes to A and E
            obj.A = [(obj.J - obj.R)*obj.Q, zeros(obj.n, obj.n_c), obj.G - obj.P, obj.K;
                    obj.B', zeros(obj.n_c, obj.n_c + obj.n_u + obj.n_d)];
            obj.E = [obj.E(1:obj.n, 1:obj.n), zeros(obj.n, obj.n_c+obj.n_u+obj.n_d);
                    zeros(obj.n_c, obj.n+obj.n_c+obj.n_u+obj.n_d)];
        end
        
        % Check whether our System is a valid PH system
        function b = isValidPHSystem(obj)
            b = 1; 
            
            if ~issymmetric(obj.J, 'skew') 
                disp('J is not skew-symmetric');
                b = 0;
            end
            
            if ~issymmetric(obj.R) || any(round(eig(obj.R),17) < 0)
                disp('R is not symmetric and positive semidefinite');
                b = 0;
            end
            
            if ~issymmetric(obj.S)
                disp('S is not a symmetric  matrix');
                b = 0;
            end
                
            if ~issymmetric(obj.M, 'skew')
                disp('M is not a skew-symmetric matrix');
                b = 0;
            end
                
            Z = [obj.R obj.P; obj.P' obj.S];
            if any(round(eig(Z),17) < 0)
                disp('Z = [R P; P'' S] is not positive semidefinite');
                b = 0;
            end
        end
        
        % Eliminates all constraints of the form B'*e = 0. Such constraints
        % are usually generated by a call to assemble()
        function [V_xz] = eliminateAlgebraicConstraints(obj)
            if size(obj.B, 2) == 0
                % No algebraic constraints --> nothing to do!
                V_xz = eye(obj.n);
                return
            end
            % Get number of algebraic constraints
            n_ac = size(obj.B, 2);
            % Construct a left annihilator such that G_a*B = 0
            G_a = null(obj.B')';
            
            % Compute transformation matrix
            V = [G_a; (obj.B'*obj.B)\obj.B'];
            
            % Transform system matrices
            J_ = V*obj.J*V';
            R_ = V*obj.R*V';
            Q_ = inv(V)'*obj.Q*inv(V);
            % Q_ must be positive definite 
            if any(eig(Q_) <= 0)
                error('Check your system matrices');
            end
            
            G_ = V*obj.G;
            P_ = V*obj.P;
            K_ = V*obj.K;
            
            % Partition Q (Q2x will be eliminated) 
            Q_11 = Q_(1:obj.n-n_ac, 1:obj.n-n_ac);
            Q_12 = Q_(1:obj.n-n_ac, obj.n-n_ac+1:end);
            Q_22 = Q_(obj.n-n_ac+1:end, obj.n-n_ac+1:end);
            Q_21 = Q_(obj.n-n_ac+1:end, 1:obj.n-n_ac); 
            
            % System order is reduced by the number of constraints
            obj.n = obj.n - n_ac;
            
            % Map from new to old state variables
            V_xz = obj.Q\V'*Q_*[eye(obj.n); -Q_22\Q_21];
            
            % Alter the decimal precision in case of numerical problems
            decimalPrecision = 14;

            obj.J = round(J_(1:obj.n, 1:obj.n), decimalPrecision);
            obj.R = round(R_(1:obj.n, 1:obj.n), decimalPrecision);
            obj.Q = round(Q_11 - Q_12*(Q_22\Q_21), decimalPrecision);
            obj.G = round(G_(1:obj.n, :), decimalPrecision);     
            obj.P = round(P_(1:obj.n, :), decimalPrecision); 
            obj.K = round(K_(1:obj.n, :), decimalPrecision);

            % Did we eliminate constraint forces?
            idx_u0 = find(~any(round(obj.G-obj.P, 10)));
            % If yes, remove them
            obj.n_u = obj.n_u - length(idx_u0);
            obj.G(:, idx_u0) = [];
            obj.P(:, idx_u0) = [];
            obj.M(idx_u0, :) = [];
            obj.M(:, idx_u0) = [];
            obj.S(idx_u0, :) = [];
            obj.S(:, idx_u0) = [];
            obj.inputNames(idx_u0) = [];
            obj.outputNames(idx_u0) = [];
            obj.inputSigns(idx_u0) = [];
            obj.outputSigns(idx_u0) = [];
            
            % G_c is now empty
            obj.B = zeros(obj.n, 0);
            % The number of constraints is reduced by a_c
            obj.n_c = obj.n_c - n_ac; 
            
            % Update simulation matrices
            obj.A = [(obj.J - obj.R)*obj.Q, obj.G - obj.P, obj.K];
            obj.E = [eye(obj.n), zeros(obj.n, obj.n_u+obj.n_d)];
            
            % Save the relation between old and new state variables
            %{
            for i=1:obj.n_ss
                obj.xMap{i} = obj.xMap{i} * V_xz; 
            end
            %}
            
            % Finally, check the validity of the system
            if ~obj.isValidPHSystem()
                error('Whoops, something went wrong');
            end
        end 
        
        function dxdt = getStateDerivative(obj, x, u, d)
            if nargin < 3
                error('x and u must be given to calculate dxdt');
            end
            if ~iscolumn(x) || length(x) ~= obj.n
                error(['x must be a column vector of length ' num2str(obj.n)]);
            end
            if ~iscolumn(u) || length(u) ~= obj.n_u
                error(['u must be a column vector of length ' num2str(obj.n_u)]);
            end
            
            interaction_term = 0;
            if obj.hasInteractionPort() 
                if nargin < 4
                    error('d must be given for systems with interaction port');
                end

                if ~iscolumn(d) || length(d) ~= obj.n_d
                    error(['d must be a column vector of length ' num2str(obj.n_z)]);
                end
                interaction_term = obj.K*d;
            end
            
            dxdt = obj.A * [x; u.*obj.inputSigns; interaction_term];
        end
        
        function dxdt = getStateDerivativeNoChecks(obj, x, u, d)           
            interaction_term = 0;
            if obj.hasInteractionPort() 
                interaction_term = obj.K*d;
            end
            
            dxdt = obj.A * [x; u.*obj.inputSigns; interaction_term];
        end
        
        function y = getSystemOutput(obj, x, u)
            if nargin < 2
                error('x must be given to calculate y');
            end
            if ~iscolumn(x) || length(x) ~= obj.n
                error(['x must be a column vector of length ' num2str(obj.n)]);
            end
                        
            y = obj.outputSigns .* ((obj.G + obj.P)' * obj.Q * x + (obj.M + obj.S) * u.*obj.inputSigns);
        end
        
        function z = getInteractionOutput(obj, x)
            if nargin < 2
                error('x must be given to calculate z');
            end
            if ~iscolumn(x) || length(x) ~= obj.n
                error(['x must be a column vector of length ' num2str(obj.n)]);
            end
            z = obj.K'*obj.Q*x;
        end
        
        function H = getHamiltonian(obj, x)
            if nargin < 2
                error('x must be given to calculate z');
            end
            if ~iscolumn(x) || length(x) ~= obj.n
                error(['x must be a column vector of length ' num2str(obj.n)]);
            end
            H = (x' * obj.Q * x) / 2;
        end
        
        % Call this function to get the ports numbers of mechanical ports
        % connected to a certain node. You can specify a port type and an
        % orientation
        function portNumbers = getExternalPorts(obj)
            nPortsFound = 0;
            portNumbers = zeros(obj.n_ports, 1); 
            for p = 1:obj.n_ports
                port = obj.ports{p};
                if strcmp(port.scope, 'external')
                    nPortsFound = nPortsFound + 1;
                    portNumbers(nPortsFound) = p;
                end
                
            end
            portNumbers = portNumbers(1:nPortsFound);
        end
        
        function portNumbers = getMechanicalPorts(obj, node, type, orientation)
            checkType = 0;
            checkOrientation = 0;
            if nargin > 2
                checkType = 1;
            end
            if nargin > 3
                checkOrientation = 1;
            end
            
            nPortsFound = 0;
            portNumbers = zeros(obj.n_ports, 1); 
            for p = 1:length(obj.nodes{node}.ports)
                port = obj.ports{obj.nodes{node}.ports(p)};
                matchFlag = 1;
                if ~isa(port, 'PH_MechanicalPort')
                    matchFlag = 0;
                end
                if checkType && ~strcmp(port.type, type)
                    matchFlag = 0;
                end
                if checkOrientation && dot(port.orientation, orientation) == 0
                    matchFlag = 0;
                end
                if matchFlag
                    nPortsFound = nPortsFound + 1;
                    portNumbers(nPortsFound) = obj.nodes{node}.ports(p);
                end
            end
            portNumbers = portNumbers(1:nPortsFound);
        end
        
        function port = getPortForIOPair(obj, pair)
            port = 0;
            for p = 1:obj.n_ports
                if obj.ports{p}.IOPair == pair
                    port = p;
                    break;
                end
            end
        end
        
    end
    
    methods (Access = protected)
        function addExternalPort(obj, domain, type, inputType, outputType, node, inputName, outputName, orientation)
            if strcmp(domain, 'mechanical')
                obj.ports{obj.n_ports+1} = PH_MechanicalPort('external', type, obj.n_d+1, inputType, outputType, node, orientation);
            else
                obj.ports{obj.n_ports+1} = PH_Port(domain, 'external', type, obj.n_d+1, inputType, outputType, node); 
            end
            obj.externalInputNames{obj.n_d+1} = inputName;
            obj.externalOutputNames{obj.n_d+1} = outputName;
            obj.n_ports = obj.n_ports + 1;
            obj.n_d = obj.n_d + 1;
        end
        
        % This function is invoked by subclasses to copy properties
        function copyData(obj, cp)
            cp.C_u = obj.C_u;
            cp.C_y = obj.C_y;
            cp.C_d = obj.C_d;
            cp.n_c = size(cp.C_u, 1) + size(cp.B, 2);

            cp.n_elements = obj.n_elements;
            cp.elements = cell(obj.n_elements, 1);
            for e=1:obj.n_elements
                cp.elements{e} = copy(obj.elements{e});
            end

            cp.n_nodes = obj.n_nodes;
            cp.nodes = cell(obj.n_nodes, 1);
            for n=1:obj.n_nodes
                cp.nodes{n} = copy(obj.nodes{n});
            end

            cp.n_ports = obj.n_ports;
            cp.ports = cell(obj.n_ports, 1);
            for p=1:obj.n_ports
                cp.ports{p} = copy(obj.ports{p});
            end

            for i=1:cp.n_u
                cp.inputNames{i} = obj.inputNames{i};
                cp.outputNames{i} = obj.outputNames{i};
                cp.inputSigns(i) = obj.inputSigns(i);
                cp.outputSigns(i) = obj.outputSigns(i);
            end
            for i=1:cp.n_d
                cp.externalInputNames{i} = obj.externalInputNames{i};
                cp.externalOutputNames{i} = obj.externalOutputNames{i};
            end
        end
        
        % Returns a deep copy of a PH_LinearSystem
        function cp = copyElement(obj)
            % PH_LinearSystem(n, E, J, Q, G, R, K, P, S, M, B)
            cp = PH_LinearSystem(obj.n, obj.E(1:obj.n, 1:obj.n), obj.J, obj.Q, obj.G, ...
                                 obj.R, obj.K, obj.P, obj.S, obj.M, obj.B);
            % copyData can also be called by subclasses
            obj.copyData(cp);
        end
        
        % Generates constraints in the form v_a = v_b and w_a = w_b at each
        % node, where v refers to a velocity in x/y/z direction and w to an 
        % angular velocity in x/y/z direction.
        function [Cu, Cy, Cd] = generateVelocityConstraints(obj)
            nConstraints = 0;
            V = zeros(nConstraints, obj.n_nodes*6);
            Cu = zeros(nConstraints, obj.n_u);
            Cy = zeros(nConstraints, obj.n_u); 
            Cd = zeros(nConstraints, obj.n_d);
            
            % Loop through the system's nodes
            for n = 1:obj.n_nodes
                node = obj.nodes{n};
                if ~isa(node, 'PH_MechanicalNode')
                    continue;
                end
                ports_mech = obj.getMechanicalPorts(n);
                for p = 1:length(ports_mech)
                    port = obj.ports{ports_mech(p)};
                    nConstraints = nConstraints + 1;
                    if strcmp(port.type, 'force/velocity')
                        V(nConstraints, n*6-5:n*6-3) = port.orientation .* double(~node.lockedDOFs(1:3))';
                        if strcmp(port.inputType, 'velocity')
                            Cu(nConstraints, port.IOPair) = obj.inputSigns(port.IOPair); 
                        else
                            Cy(nConstraints, port.IOPair) = obj.outputSigns(port.IOPair); 
                        end
                    end
                    if strcmp(port.type, 'torque/angular velocity')
                        V(nConstraints, n*6-2:n*6) = port.orientation .* double(~node.lockedDOFs(4:6))';
                        if strcmp(port.inputType, 'angular velocity')
                            Cu(nConstraints, port.IOPair) = obj.inputSigns(port.IOPair); 
                        else
                            Cy(nConstraints, port.IOPair) = obj.outputSigns(port.IOPair); 
                        end
                    end
                end
            end  
                       
            if size(Cu, 1) < nConstraints && obj.n_u > 0
                Cu(nConstraints,1) = 0;
            end
            if size(Cy, 1) < nConstraints && obj.n_u > 0 
                Cy(nConstraints,1) = 0;
            end
            if size(Cd, 1) && obj.n_d > 0
                Cd(nConstraints, 1) = 0;
            end
            
            idx = ~any([Cu Cy Cd]'); 
            Cu(idx,:) = [];
            Cy(idx,:) = [];
            Cd(idx,:) = [];
            V(idx, :) = [];
            V(:, ~any(V)) = [];
            G_v = null(V')';
            C = round(rref(G_v*[Cu, Cy, Cd]), 14);
            Cu = C(:, 1:obj.n_u);
            Cy = C(:, obj.n_u+1:2*obj.n_u);
            Cd = C(:, 2*obj.n_u+1:end);
        end
        
        % Generates constraints of the form sum(M) = 0 / sum(F) = 0 at each
        % node
        function [Cu, Cy, Cd] = generateForceConstraints(obj)
            nConstraints = 0; 
            Cu = zeros(nConstraints, obj.n_u);
            Cy = zeros(nConstraints, obj.n_u); 
            Cd = zeros(nConstraints, obj.n_d);
            dir_str = 'xyz';
            for n = 1:obj.n_nodes
                node = obj.nodes{n};
                if ~isa(node, 'PH_MechanicalNode')
                    continue;
                end
                for d = 1:3
                    dir = zeros(3, 1);
                    dir(d) = 1;
                    ports_fv = obj.getMechanicalPorts(n, 'force/velocity', dir);
                    % All forces sum up to zero
                    if ~node.lockedDOFs(d) && ~isempty(ports_fv)
                        nConstraints = nConstraints+1;
                        for p=1:length(ports_fv)
                            port = obj.ports{ports_fv(p)};
                            if strcmp(port.inputType, 'force')
                                Cu(nConstraints, port.IOPair) = obj.inputSigns(port.IOPair) * port.orientation(d);
                            else
                                Cy(nConstraints, port.IOPair) = obj.outputSigns(port.IOPair) * port.orientation(d);
                            end
                        end
                        obj.addExternalPort('mechanical', 'force/velocity', 'force', 'velocity', n, ...
                                            ['n' num2str(n) '.F_' dir_str(d)], ['n' num2str(n) '.v_' dir_str(d)], dir);
                        Cd(nConstraints, obj.n_d) = 1;       
                    end
                end
                for d = 1:3
                    dir = zeros(3, 1);
                    dir(d) = 1;
                    ports_tr = obj.getMechanicalPorts(n, 'torque/angular velocity', dir);
                    % All torques sum up to zero
                    if ~node.lockedDOFs(d+3) && ~isempty(ports_tr)
                        nConstraints = nConstraints+1;
                        for p=1:length(ports_tr)
                            port = obj.ports{ports_tr(p)};
                            if strcmp(port.inputType, 'torque')
                                Cu(nConstraints, port.IOPair) = obj.inputSigns(port.IOPair) * port.orientation(d);
                            else
                                Cy(nConstraints, port.IOPair) = obj.outputSigns(port.IOPair) * port.orientation(d);
                            end
                        end
                        obj.addExternalPort('mechanical', 'torque/angular velocity', 'torque', 'angular velocity', n, ...
                                            ['n' num2str(n) '.M_' dir_str(d)], ['n' num2str(n) '.w_' dir_str(d)], dir);
                        Cd(nConstraints, obj.n_d) = 1;       
                    end
                end
            end

            if size(Cu, 1) < nConstraints && obj.n_u > 0
                Cu(nConstraints,1) = 0;
            end
            if size(Cy, 1) < nConstraints && obj.n_u > 0
                Cy(nConstraints, 1) = 0;
            end
            if size(Cd, 1) < nConstraints && obj.n_d > 0 
                Cd(nConstraints, 1) = 0;
            end
            idx = ~any([Cu Cy]'); 
            Cu(idx,:) = [];
            Cy(idx,:) = [];
            Cd(idx,:) = [];
        end
    end 
end