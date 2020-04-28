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
        M       % Feedthrough matrix skew symmetric part M
        S       % Feedthrough matrix symmetric part S (dissipation)
        
        C_e     % Matrix for effort coupling
        C_f     % Matrix for flow coupling
        C_u     % Matrix for I/O coupling 
        C_y     % Matrix for I/O coupling
        
        B       % Algebraic constraints on efforts
    end
    
    methods(Access = public)
        %              PH_LinearSystem(n, J, Q, G, R, P, S, M, B)
        function obj = PH_LinearSystem(n, J, Q, varargin)
            % Sanity checks
            obj = obj@PH_System('Linear port-Hamiltonian system');
           
            % Strict type checking
            if ~isscalar(n) || n < 0
                error('n must be >= 0');
            end
            obj.n = n;
                       
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
            if nargin > 3 && ~isempty(varargin{1})
                G = varargin{1};
                if ~ismatrix(G) || size(G, 1) ~= n
                    error('G must be a matrix with n rows');
                end            
                obj.n_u = size(G, 2);
                obj.G = G; 
            end
            
            % Dissipation
            obj.R = zeros(n);
            if nargin > 4
                R = varargin{2};
                if ~isempty(R)
                    if ~ismatrix(R) || any(size(R) ~= n) || ~issymmetric(R) || any(round(eig(R),17) < 0)
                        error('R must be a symmetric positive semidefinite nxn matrix');
                    end
                    obj.R = R;
                end
            end
                                  
            % System with direct input feedthrough
            obj.P = zeros(obj.n, obj.n_u); 
            obj.M = zeros(obj.n_u);
            obj.S = zeros(obj.n_u);
            if nargin > 7
                P = varargin{3};
                S = varargin{4};
                M = varargin{5};
                
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
            
            obj.C_e = zeros(0, obj.n);
            obj.C_f = zeros(0, obj.n);
            obj.C_u = zeros(0, obj.n_u);
            obj.C_y = zeros(0, obj.n_u);
            
            obj.B = zeros(obj.n, 0);
            if nargin > 8 
                B = varargin{6};
                if size(B, 1) ~= obj.n 
                    error('B must be a n x n_c matrix');
                end
                
                obj.B = B;
            end
            obj.n_c = size(obj.B, 2);
                
            obj.n_elements = 1;
            obj.elements{1} = PH_Element('Linear PH system', [], []);
        end
        
        % Coupling constraints for flows, effors, inputs & outputs
        function addCouplingConstraints(obj, C_e, C_f, C_u, C_y) 
            if ~isempty(C_e)
                if ~ismatrix(C_e) || size(C_e, 2) ~= obj.n 
                    error(['C_e must be a matrix with ', num2str(obj.n) ' columns']);
                end
                if ~ismatrix(C_f) || size(C_f, 2) ~= obj.n 
                    error(['C_f must be a matrix with ', num2str(obj.n) ' columns']);
                end
                if ~(size(C_e, 1) == size(C_f, 1))
                    error('C_e and C_f must be matrices of equal size!'); 
                end
                
                % remove 'empty' constraints
                idx_0 = ~any(C_e') & ~any(C_f');
                C_e(idx_0,:) = [];
                C_f(idx_0,:) = [];
                
                obj.C_e = [obj.C_e; C_e];
                obj.C_f = [obj.C_f; C_f];
                obj.C_u = [obj.C_u; zeros(size(obj.C_e, 1), obj.n_u)];
                obj.C_y = [obj.C_y; zeros(size(obj.C_e, 1), obj.n_u)];
            end
            
            if nargin > 4
                if ~ismatrix(C_u) || size(C_u, 2) ~= obj.n_u 
                    error(['C_u must be a matrix with ', num2str(obj.n_u) ' columns']);
                end
                if ~ismatrix(C_y) || size(C_y, 2) ~= obj.n_u 
                    error(['C_y must be a matrix with ', num2str(obj.n_u) ' columns']);
                end
                if ~(size(C_u, 1) == size(C_y, 1))
                    error('C_u and C_y must be matrices of equal size!'); 
                end
                
                % remove 'empty' constraints
                idx_0 = ~any(C_u') & ~any(C_y');
                C_u(idx_0,:) = [];
                C_y(idx_0,:) = [];
                
                obj.C_u = [obj.C_u; C_u];
                obj.C_y = [obj.C_y; C_y];
                obj.C_e = [obj.C_e; zeros(size(obj.C_u, 1), obj.n)];
                obj.C_f = [obj.C_f; zeros(size(obj.C_u, 1), obj.n)];
            end
            
            obj.n_c = size(obj.C_e,1) + size(obj.B, 2);
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
            
            obj.C_e = [system.C_e(:,1:N_b2) zeros(size(system.C_e, 1), N_a) system.C_e(:, N_b2+1:end);
                        zeros(size(obj.C_e, 1), N_b2) obj.C_e zeros(size(obj.C_e, 1), N_b2)];
            obj.C_f = [system.C_f(:,1:N_b2) zeros(size(system.C_f, 1), N_a) system.C_f(:, N_b2+1:end);
                        zeros(size(obj.C_f, 1), N_b2) obj.C_f zeros(size(obj.C_f, 1), N_b2)];
            obj.C_u = [obj.C_u, zeros(size(obj.C_u, 1), system.n_u); ...
                        zeros(size(system.C_u, 1), obj.n_u) system.C_u];
            obj.C_y = [obj.C_y, zeros(size(obj.C_y, 1), system.n_u); ...
                        zeros(size(system.C_y, 1), obj.n_u) system.C_y];
            obj.n_c = size(obj.C_e, 1) + size(obj.B, 2);   
            
            for p = 1:system.n_ports
                obj.ports{end+1} = copy(system.ports{p});
                obj.ports{end}.nodes = obj.ports{end}.nodes + obj.n_nodes;
                if isa(obj.ports{end}, 'PH_Port_external') 
                    obj.ports{end}.IOPair = obj.ports{end}.IOPair + obj.n_u;
                end
                if isa(obj.ports{end}, 'PH_Port_boundary')
                    obj.ports{end}.IOPair = obj.ports{end}.IOPair + obj.n_u;
                    obj.ports{end}.element = obj.ports{end}.element + obj.n_elements;
                end
            end
            
            for n = 1:system.n_nodes
                obj.nodes{end+1} = copy(system.nodes{n});
                obj.nodes{end}.ports = obj.nodes{end}.ports + obj.n_ports;
                obj.nodes{end}.elements = obj.nodes{end}.elements + obj.n_elements;
            end
            
            for e = 1:system.n_elements
                obj.elements{end+1} = copy(system.elements{e});
                obj.elements{end}.ports = obj.elements{end}.ports + obj.n_ports;
                obj.elements{end}.nodes = obj.elements{end}.nodes + obj.n_nodes;
            end
            
            obj.n_ports = obj.n_ports + system.n_ports;
            obj.n_nodes = obj.n_nodes + system.n_nodes;
            obj.n_elements = obj.n_elements + system.n_elements; 
            obj.n = obj.n + system.n;
            obj.n_u = obj.n_u + system.n_u;
        end
        
        function add_single(obj, system)
            if ~isa(system, 'PH_LinearSystem') || ~system.isValidPHSystem()
                error('system must be a valid ''PH_LinearSystem''');
            end

            obj.J = blkdiag(obj.J, system.J);
            obj.R = blkdiag(obj.R, system.R);
            obj.Q = blkdiag(obj.Q, system.Q);
            obj.G = blkdiag(obj.G, system.G);
            obj.P = blkdiag(obj.P, system.P);
            obj.M = blkdiag(obj.M, system.M);
            obj.S = blkdiag(obj.S, system.S);
            obj.B = blkdiag(obj.B, system.B);
            
            obj.C_e = [obj.C_e, zeros(size(obj.C_e, 1), system.n); ...
                       zeros(size(system.C_e, 1), obj.n_u) system.C_e];
            obj.C_f = [obj.C_f, zeros(size(obj.C_f, 1), system.n); ...
                       zeros(size(system.C_f, 1), obj.n_u) system.C_f];
            obj.C_u = [obj.C_u, zeros(size(obj.C_u, 1), system.n_u); ...
                       zeros(size(system.C_u, 1), obj.n_u) system.C_u];
            obj.C_y = [obj.C_y, zeros(size(obj.C_y, 1), system.n_u); ...
                       zeros(size(system.C_y, 1), obj.n_u) system.C_y];
            obj.n_c = size(obj.C_e, 1) + size(obj.B, 2);   
            
            for p = 1:system.n_ports
                obj.ports{end+1} = copy(system.ports{p});
                obj.ports{end}.node = obj.ports{end}.node + obj.n_nodes; 
                if isa(obj.ports{end}, 'PH_Port_external') 
                    obj.ports{end}.IOPair = obj.ports{end}.IOPair + obj.n_u;
                end
                if isa(obj.ports{end}, 'PH_Port_boundary')
                    obj.ports{end}.IOPair = obj.ports{end}.IOPair + obj.n_u;
                    obj.ports{end}.element = obj.ports{end}.element + obj.n_elements;
                end
            end
            
            for n = 1:system.n_nodes
                obj.nodes{end+1} = copy(system.nodes{n});
                obj.nodes{end}.ports = obj.nodes{end}.ports + obj.n_ports;
                obj.nodes{end}.elements = obj.nodes{end}.elements + obj.n_elements;
            end
            
            for e = 1:system.n_elements
                obj.elements{end+1} = copy(system.elements{e});
                obj.elements{end}.ports = obj.elements{end}.ports + obj.n_ports;
                obj.elements{end}.nodes = obj.elements{end}.nodes + obj.n_nodes;
            end
            
            obj.n_ports = obj.n_ports + system.n_ports;
            obj.n_nodes = obj.n_nodes + system.n_nodes;
            obj.n_elements = obj.n_elements + system.n_elements; 
            obj.n = obj.n + system.n;
            obj.n_u = obj.n_u + system.n_u;
        end
        
        % Generate constraints of the form Cy * y + Cu * u = 0.
        % Constraints are generated automatically for each node. Only the
        % mechanical domain is supported so far.
        function generateConstraints(obj)
            % Mechanical constraint generation
            [Ce_v, Cf_v, Cu_v, Cy_v] = obj.generateVelocityConstraints();           
            [Ce_f, Cf_f, Cu_f, Cy_f] = obj.generateForceConstraints();

            % Add constraint matrices
            obj.C_e = [obj.C_e; Ce_v; Ce_f];
            obj.C_f = [obj.C_f; Cf_v; Cf_f];
            obj.C_u = [obj.C_u; Cu_v; Cu_f];
            obj.C_y = [obj.C_y; Cy_v; Cy_f];
            
            % Adjust the number of constraints
            obj.n_c = size(obj.C_e, 1) + size(obj.B, 2); 
            
            % Constraints were generated successfully. To obtain an ODE system,
            % call assemble() followed by eliminateAlgebraicConstraints()
        end

        % Eliminates the Cu, Cy constraints by solving for u 
        % In the process, a set of algebraic constraints of the form B'*e = 0
        % are usually generated. They can be eliminated by a call to
        % eliminateAlgebraicConstraints
        function assembleDAESystem(obj)
            if obj.n_c == 0
                % No constraints --> nothing to do
                return
            end
                        
            if any(any(obj.C_e)) || any(any(obj.C_f))
                warning('Constraints on efforts and flows are not handled... setting to zero');
                obj.C_e = [];
                obj.C_f = [];
            end
            
            
            % Get algebraic constraints
            idx_ac = ~any(obj.C_u');
            B_y = obj.C_y(idx_ac, :) * (obj.G+obj.P)';         
            if ~isempty(obj.B)
                obj.B = [obj.B, B_y'];
            else
                obj.B = B_y';
            end
            obj.B(:, ~any(obj.B)) = [];    
            obj.C_y(idx_ac, :) = [];
            obj.C_u(idx_ac, :) = [];
                        
            
            % Solve remaining equations for u
            u_f = ones(1, obj.n_u);
            A_uu = eye(obj.n_u);
            A_ue = zeros(obj.n_u, obj.n);
            for i = 1:size(obj.C_u, 1)
                idx_u = find(obj.C_u(i, :), 1);
                if ~isempty(idx_u)
                    % Constrain this input
                    u_f(idx_u) = 0;
                    a_uu = obj.C_u(i, :)*A_uu + obj.C_y(i,:)*(obj.M+obj.S)*A_uu;
                    b = a_uu(idx_u);
                    a_uu(idx_u) = 0;
                    A_uu(idx_u, :) = -b .* a_uu;
                    A_ue(idx_u, :) = -b .* obj.C_y(i,:)*(obj.G+obj.P)';

                    % Update dependent inputs
                    idx_uc = find(A_uu(:, idx_u));
                    for k=1:length(idx_uc)
                        b = A_uu(idx_uc(k), idx_u);
                        A_uu(idx_uc(k), idx_u) = 0;
                        A_uu(idx_uc(k), :) = A_uu(idx_uc(k), :) + b*A_uu(idx_u, :);
                        A_ue(idx_uc(k), :) = A_ue(idx_uc(k), :) + b*A_ue(idx_u, :);
                    end

                    % Update constraints
                    b = obj.C_u(i, idx_u);
                    c_uu = -b*obj.C_u(i, :);
                    c_uu(idx_u) = 0;
                    c_uy = -b*obj.C_y(i, :);
                    idx_uc = find(obj.C_u(i+1:end, idx_u)) + i;
                    for k=1:length(idx_uc)
                        b = obj.C_u(idx_uc(k), idx_u);
                        obj.C_u(idx_uc(k), idx_u) = 0;
                        obj.C_u(idx_uc(k), :) = obj.C_u(idx_uc(k), :) + b * c_uu;
                        obj.C_y(idx_uc(k), :) = obj.C_y(idx_uc(k), :) + b * c_uy;
                    end
                else
                    warning(['ignoring constraint ' num2str(i)]);
                end
            end
            A_uu(:, ~any(A_uu)) = []; 
            u_new = find(u_f); 
            u_dep = find(~u_f);
            
            J_ = obj.J + (obj.G - obj.P) * A_ue;
            R_ = obj.R; 
            G_ = (obj.G - obj.P) * A_uu; 
            P_ = zeros(obj.n, size(G_, 2));
            D_ = A_uu' * (obj.M + obj.S) * A_uu;
            if ~issymmetric(D_, 'skew')
                error('Splitting D into S and M needs to be implemented');
            else
                M_ = D_; 
                S_ = zeros(size(M_, 1), size(M_, 2));
            end
            
            % remove dependent inputs
            u_dep = sort(u_dep, 'descend');
            ports = obj.getPortsForIOPairs(u_dep);
            obj.deletePorts(ports);
            u_new = 1:length(u_new);
            
            % Assign new system matrices
            obj.J = J_;
            obj.R = R_;
            obj.G = G_;
            obj.P = P_; 
            obj.M = M_; 
            obj.S = S_; 
            
            % Check for zero colums in G (inputs that are zero) and remove
            % them
            idx_u0 = ~any(obj.G-obj.P);
          
            u_dep = u_new(idx_u0);
            u_dep = sort(u_dep, 'descend');
            ports = obj.getPortsForIOPairs(u_dep);
            obj.deletePorts(ports);
            
            % System check... 
            if ~obj.isValidPHSystem()
                error('Something is messed up... check your constraints');
            end
            
            % All I/O coupling constraints are now eliminated
            obj.C_u = zeros(0, obj.n_u);
            obj.C_y = zeros(0, obj.n_u);
            obj.C_e = zeros(0, obj.n);
            obj.C_f = zeros(0, obj.n);
            
            % Some algebraic constraints might remain
            % Call eliminateAlgebraicConstraints() to eliminate them
            % This also eliminates constraint forces (inputs)
            obj.n_c = size(obj.B, 2);
        end
        
        % Check whether our System is a valid PH system
        function b = isValidPHSystem(obj)
            b = 1; 
            
            if ~issymmetric(round(obj.J, 10), 'skew') 
                disp('J is not skew-symmetric');
                b = 0;
            end
            
            if ~issymmetric(round(obj.R, 10)) || any(round(eig(obj.R),17) < 0)
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
        % are usually generated by a call to assembleDAESystem()
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

            % Partition Q (Q2x will be eliminated) 
            Q_11 = Q_(1:obj.n-n_ac, 1:obj.n-n_ac);
            Q_12 = Q_(1:obj.n-n_ac, obj.n-n_ac+1:end);
            Q_22 = Q_(obj.n-n_ac+1:end, obj.n-n_ac+1:end);
            Q_21 = Q_(obj.n-n_ac+1:end, 1:obj.n-n_ac); 
            
            % System order is reduced by the number of constraints
            obj.n = obj.n - n_ac;
            
            % Map from new to old state variables
            V_xz = obj.Q\V'*Q_*[eye(obj.n); -Q_22\Q_21];
                
            obj.J = J_(1:obj.n, 1:obj.n);
            obj.R = R_(1:obj.n, 1:obj.n);
            obj.Q = Q_11 - Q_12*(Q_22\Q_21);
            obj.G = G_(1:obj.n, :);     
            obj.P = P_(1:obj.n, :); 

            % Did we eliminate constraint forces?
            idx_u0 = find(~any(round(obj.G-obj.P, 10)));
            % If yes, remove them
            u_dep = sort(idx_u0, 'descend');
            ports = obj.getPortsForIOPairs(u_dep);
            obj.deletePorts(ports);

            % G_c is now empty
            obj.B = zeros(obj.n, 0);
            % The number of constraints is reduced by a_c
            obj.n_c = obj.n_c - n_ac; 
            
            % Check whether there is a block of linearly dependent flows
            row_idx = ~any(obj.G'-obj.P');
            col_idx = any(obj.J(row_idx, :)-obj.R(row_idx,:));
            
            J_sorted = [obj.J(:, col_idx == 1), obj.J(:, col_idx ~=1)]; 
            J_sorted = [J_sorted(row_idx ~= 1, :); J_sorted(row_idx == 1, :)];
            R_sorted = [obj.R(:, col_idx == 1), obj.R(:, col_idx ~=1)]; 
            R_sorted = [R_sorted(row_idx ~= 1, :); R_sorted(row_idx == 1, :)];
            G_sorted = [obj.G(row_idx ~= 1, :); obj.G(row_idx == 1, :)];
            P_sorted = [obj.P(row_idx ~= 1, :); obj.P(row_idx == 1, :)];
            Q_sorted = [obj.Q(:, col_idx == 1), obj.Q(:, col_idx ~= 1)];
            Q_sorted = [Q_sorted(row_idx ~= 1, :); Q_sorted(row_idx == 1, :)]; 
            
            row_idx = ~any(G_sorted'-P_sorted');
            col_idx = any(J_sorted(row_idx, :) - R_sorted(row_idx, :));
            
            if sum(row_idx) > sum(col_idx)
                % If such a block was found, construct an orthonormal basis 
                % to eliminate dependent states
                A_b = J_sorted(row_idx == 1, :) - R_sorted(row_idx == 1, :);
                T_b = orth(A_b);
                T = blkdiag(eye(obj.n-size(T_b, 1)), T_b);
                
                % Construct algebraic constraints B_ for linearly dependent
                % states
                B_ = null((obj.Q*T)');
                G_a = null(B_')';
                
                % Proceed as above 
                V = [G_a; (B_'*B_)\B_'];
            
                % Transform system matrices
                J_ = V*J_sorted*V';
                R_ = V*R_sorted*V';
                Q_ = (V')\(Q_sorted/V);
                % Q_ must be positive definite 
                if any(eig(Q_) <= 0)
                    error('Check your system matrices');
                end
                G_ = V*G_sorted;
                P_ = V*P_sorted;

                % Partition Q (Q2x will be eliminated) 
                Q_11 = Q_(1:obj.n-n_ac, 1:obj.n-n_ac);
                Q_12 = Q_(1:obj.n-n_ac, obj.n-n_ac+1:end);
                Q_22 = Q_(obj.n-n_ac+1:end, obj.n-n_ac+1:end);
                Q_21 = Q_(obj.n-n_ac+1:end, 1:obj.n-n_ac); 

                % System order is reduced by the number of constraints
                obj.n = obj.n - n_ac;
                % Ensure correct size of B
                obj.B = zeros(obj.n, 0);

                % Map from new to old state variables
                % TODO: fix this mapping... 
                V_xz = obj.Q\V'*Q_*[eye(obj.n); -Q_22\Q_21];


                obj.J = J_(1:obj.n, 1:obj.n);
                obj.R = R_(1:obj.n, 1:obj.n);
                obj.Q = Q_11 - Q_12*(Q_22\Q_21);
                obj.G = G_(1:obj.n, :);     
                obj.P = P_(1:obj.n, :); 
            end
            
            % Ensure correct size of constraint matrices 
            obj.C_u = zeros(0, obj.n_u);
            obj.C_y = zeros(0, obj.n_u);
            obj.C_e = zeros(0, obj.n);
            obj.C_f = zeros(0, obj.n);
            
            % Finally, check the validity of the system
            if ~obj.isValidPHSystem()
                error('Whoops, something went wrong. Try decreasing arithmetic precision');
            end
        end 
        
        function roundSystemMatrices(obj, decimalPlaces)
            obj.J = round(obj.J, decimalPlaces);
            obj.R = round(obj.R, decimalPlaces);
            obj.Q = round(obj.Q, decimalPlaces);
            obj.G = round(obj.G, decimalPlaces);
            obj.P = round(obj.P, decimalPlaces);
            obj.M = round(obj.M, decimalPlaces);
            obj.S = round(obj.S, decimalPlaces);
        end
        
        function transformToGlobalDOFs(obj)
            row_idx = any(obj.G'-obj.P');
            G_a = null(obj.G(row_idx,:)');
            T = blkdiag([(obj.G(row_idx, :)'*obj.G(row_idx,:))\obj.G(row_idx,:)'; G_a'], eye(obj.n-sum(row_idx)));
            
            % Left mulitply the system with T. J and R need to be right
            % multiplied with T' to preserve symmetry properties
            obj.J = T*obj.J*T';
            obj.R = T*obj.R*T';
            % This transformation affects e_p only --> turns to M^-1
            % Q needs to be left multiplied with the inverse of T' and
            % right multiplied with the inverse of T
            obj.Q = T'\(obj.Q/T);
            obj.G = T*obj.G;
            obj.P = T*obj.P;
            
            T = blkdiag(eye(obj.n/2), pinv(obj.J(obj.n/2+1:end, 1:obj.n/2)));

            % Again, J and R needs to be right multiplied by T'
            obj.J = T*obj.J*T';
            obj.R = T*obj.R*T';
            % This transformation affects e_q only --> turns to K
            obj.Q = T'\(obj.Q/T);
            obj.G = T*obj.G;
            obj.P = T*obj.P;
        end
        
        function names = getInputNames(obj, IOPairs)
            names = cell(1, length(IOPairs));
            ports = obj.getPortsForIOPairs(IOPairs);
            for p = 1:length(ports)
                names{p} = obj.ports{ports(p)}.inputName;
            end
        end
        
        function dxdt = getStateDerivative(obj, x, u)
            if nargin < 3
                error('x and u must be given to calculate dxdt');
            end
            if ~iscolumn(x) || length(x) ~= obj.n
                error(['x must be a column vector of length ' num2str(obj.n)]);
            end
            if ~iscolumn(u) || length(u) ~= obj.n_u
                error(['u must be a column vector of length ' num2str(obj.n_u)]);
            end            
            dxdt = [(obj.J - obj.R)*obj.Q, obj.G - obj.P] * [x; u];
        end
        
        function dxdt = getStateDerivativeNoChecks(obj, x, u)           
            dxdt = [(obj.J - obj.R)*obj.Q, obj.G - obj.P] * [x; u];
        end
        
        function y = getSystemOutput(obj, x, u)
            if nargin < 2
                error('x must be given to calculate y');
            end
            if ~iscolumn(x) || length(x) ~= obj.n
                error(['x must be a column vector of length ' num2str(obj.n)]);
            end
                        
            y = ((obj.G + obj.P)' * obj.Q * x + (obj.M + obj.S) * u);
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
        
        % Returns all external ports of the system
        function portNumbers = getExternalPorts(obj)
            nPortsFound = 0;
            portNumbers = zeros(obj.n_ports, 1); 
            for p = 1:obj.n_ports
                port = obj.ports{p};
                if isa(port, 'PH_Port_external')
                    nPortsFound = nPortsFound + 1;
                    portNumbers(nPortsFound) = p;
                end
                
            end
            portNumbers = portNumbers(1:nPortsFound);
        end
        
        % Call this function to get the ports numbers of mechanical ports
        % connected to a certain node. You can specify a port type and an
        % orientation
        function portNumbers = getMechanicalPorts(obj, node, types, orientation)
            checkTypes = 0;
            checkOrientation = 0;
            if nargin > 2 && ~isempty(types) && iscell(types)
                checkTypes = 1;
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
                    continue
                end
                if checkOrientation
                    if dot(port.orientation(:, find(port.nodes == node, 1)), orientation) == 0
                        continue
                    end
                end
                if checkTypes
                    matchFlag = 0;
                    for t = 1:length(types)
                        if strcmp(port.type, types{t}) 
                            matchFlag = 1;
                        end
                    end 
                end

                if matchFlag
                    nPortsFound = nPortsFound + 1;
                    portNumbers(nPortsFound) = obj.nodes{node}.ports(p);
                end
            end
            portNumbers = portNumbers(1:nPortsFound);
        end
        
        function ports = getPortsForIOPairs(obj, pairs)
            ports = zeros(1, length(pairs));
            for p = 1:obj.n_ports
                if (isa(obj.ports{p}, 'PH_Port_external') || isa(obj.ports{p}, 'PH_Port_boundary'))
                    for i = 1:length(pairs)
                        if obj.ports{p}.IOPair == pairs(i)
                            ports(i) = p;
                        end
                    end
                end
            end
        end
        
        function deletePorts(obj, ports)
            if ~isvector(ports) || any(ports > obj.n_ports) 
                error('ports must be given as an array of scalars <= n_ports');
            end
            
            ports = sort(ports, 'descend');

            IOPairs = zeros(1, length(ports));
            for p = 1:length(ports)
                port = obj.ports{ports(p)};
                if isa(port, 'PH_Port_external') || isa(port, 'PH_Port_boundary')
                    IOPairs(p) = port.IOPair;
                end
                port.delete();
                obj.ports(ports(p)) = [];
            end
            obj.n_ports = obj.n_ports - length(ports);
            IOPairs = sort(IOPairs, 'descend');
            IOPairs(IOPairs == 0) = [];
            
            for p=1:obj.n_ports
                port = obj.ports{p};
                if isa(port, 'PH_Port_external') || isa(port, 'PH_Port_boundary')
                    if(any(IOPairs < port.IOPair))
                        port.IOPair = port.IOPair - sum(double(IOPairs < port.IOPair));
                    end
                end
            end
            
            obj.G(:, IOPairs) = [];
            obj.P(:, IOPairs) = [];
            obj.M(IOPairs, :) = [];
            obj.M(:, IOPairs) = [];
            obj.S(IOPairs, :) = [];
            obj.S(:, IOPairs) = [];
            obj.n_u = obj.n_u - length(ports); 
            
            for n=1:obj.n_nodes
                node = obj.nodes{n};
                for p = 1:length(ports)
                    node.ports(node.ports == ports(p)) = [];
                    node.ports(node.ports > ports(p)) = node.ports(node.ports > ports(p)) - 1; 
                end
            end
            for e=1:obj.n_elements
                element = obj.elements{e};
                for p = 1:length(ports)
                    element.ports(element.ports == ports(p)) = [];
                    element.ports(element.ports > ports(p)) = element.ports(element.ports > ports(p)) - 1;
                end
            end
        end
        
    end
    
    methods (Access = protected)
        function addExternalPort(obj, domain, type, nodes, IOPair, inputName, outputName, orientation)
            if strcmp(domain, 'mechanical')
                obj.ports{obj.n_ports+1} = PH_MechanicalPort_external(type, nodes, orientation, IOPair, inputName, outputName);
            else
                obj.ports{obj.n_ports+1} = PH_Port_external(domain, type, nodes, IOPair, inputName, outputName); 
            end
            for n = 1:length(nodes)
                obj.nodes{nodes(n)}.ports = [obj.nodes{nodes(n)}.ports obj.n_ports+1];
            end
            obj.n_ports = obj.n_ports + 1;
        end
        
        % This function is invoked by subclasses to copy properties
        function copyData(obj, cp)
            cp.C_u = obj.C_u;
            cp.C_y = obj.C_y;
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
        end
        
        % Returns a deep copy of a PH_LinearSystem
        function cp = copyElement(obj)
            % PH_LinearSystem(n, J, Q, G, R, P, S, M, B)
            cp = PH_LinearSystem(obj.n, obj.J, obj.Q, obj.G, ...
                                 obj.R, obj.P, obj.S, obj.M, obj.B);
            % copyData can also be called by subclasses
            obj.copyData(cp);
        end
        
        % Generates constraints in the form v_a = v_b and w_a = w_b at each
        % node, where v refers to a velocity in x/y/z direction and w to an 
        % angular velocity in x/y/z direction.
        function [Ce, Cf, Cu, Cy] = generateVelocityConstraints(obj)
            nConstraints = 0;
            V = zeros(nConstraints, obj.n_nodes*6);
            % TODO: include effort and flow constraints
            Ce = zeros(nConstraints, obj.n);
            Cf = zeros(nConstraints, obj.n);
            Cu = zeros(nConstraints, obj.n_u);
            Cy = zeros(nConstraints, obj.n_u); 
                        
            % Loop through the system's nodes
            for n = 1:obj.n_nodes
                node = obj.nodes{n};
                if ~isa(node, 'PH_MechanicalNode') 
                    continue;
                end
                
                % Get mechanical ports at current node
                ports_mech = obj.getMechanicalPorts(n);
                duplicates = [];
                % Loop through ports at current node
                for p = 1:length(ports_mech)
                    % Get next port and its orientation at current node
                    port = obj.ports{ports_mech(p)};
                    orientation = port.orientation(:,find(port.nodes==n,1));
                    
                    % External ports are not constrained as long as there
                    % are no duplicates (one constraint per duplicate)
                    if isa(port, 'PH_Port_external')
                        duplicate = 0;
                        for pp = 1:length(ports_mech)
                            if pp ~= p  && isa(obj.ports{ports_mech(pp)}, 'PH_Port_external') ...
                                    && ~any(obj.ports{ports_mech(pp)}.orientation ~= orientation)...
                                    && ~any(duplicates == ports_mech(pp))
                                duplicate = 1;
                                duplicates(end+1) = ports_mech(p);
                            end
                        end
                        if ~duplicate
                            continue
                        end
                    end
        
                    % Add one velocity constraint for each port
                    nConstraints = nConstraints + 1;
                    % Define how port is related to global DOFs
                    if strcmp(port.type, 'force') || strcmp(port.type, 'velocity')
                        V(nConstraints, n*6-5:n*6-3) = orientation .* double(~node.lockedDOFs(1:3))';   
                    else
                        V(nConstraints, n*6-2:n*6) = orientation .* double(~node.lockedDOFs(4:6))';
                    end
                    % Add entry in respective constraint matrix
                    if isa(port, 'PH_Port_storage')
                        % For storage ports, add flow/effort constraints
                        if strcmp(port.type, 'force') || strcmp(port.type, 'torque')
                            % effort is force/torque -> flow is velocity
                            Cf(nConstraints, port.flow(find(port.nodes == n,1))) = 1;
                        else
                            % effort is velocity
                            Ce(nConstraints, port.effort(find(port.nodes == n, 1))) = 1; 
                        end
                    elseif isa(port, 'PH_Port_external') || isa(port, 'PH_Port_boundary')
                        % For external/boundary ports, add input/output constraints
                        if strcmp(port.type, 'force') || strcmp(port.type, 'torque')
                            % input is force/torque -> output is velocity
                            Cy(nConstraints, port.IOPair) = 1;
                        else
                            % input is velocity
                            Cu(nConstraints, port.IOPair) = 1; 
                        end
                    end
                end
            end  
            
            % Ensure all arrays are of the same size 
            if size(Ce, 1) < nConstraints
                Ce(nConstraints, 1) = 0;
            end
            if size(Cf, 1) < nConstraints
                Cf(nConstraints, 1) = 0;
            end
            if size(Cu, 1) < nConstraints && obj.n_u > 0
                Cu(nConstraints, 1) = 0;
            end
            if size(Cy, 1) < nConstraints && obj.n_u > 0 
                Cy(nConstraints, 1) = 0;
            end
            
            % Remove constraints that are zero
            idx = ~any([Ce Cf Cu Cy]');
            Ce(idx,:) = [];
            Cf(idx,:) = [];
            Cu(idx,:) = [];
            Cy(idx,:) = [];
            V(idx, :) = [];
            V(:, ~any(V)) = [];
            G_v = null(V')';
            %C = round(G_v*[Ce, Cf, Cu, Cy], 14);
            C = G_v*[Ce, Cf, Cu, Cy];
            Ce = C(:, 1:obj.n);
            Cf = C(:, obj.n+1:2*obj.n);
            Cu = C(:, 2*obj.n+1:2*obj.n+obj.n_u);
            Cy = C(:, 2*obj.n+obj.n_u+1:2*obj.n+2*obj.n_u);
        end
        
        % Generates constraints of the form sum(M) = 0 / sum(F) = 0 at each
        % node
        function [Ce, Cf, Cu, Cy] = generateForceConstraints(obj)
            nConstraints = 0; 
            Ce = zeros(nConstraints, obj.n);
            Cf = zeros(nConstraints, obj.n); 
            Cu = zeros(nConstraints, obj.n_u);
            Cy = zeros(nConstraints, obj.n_u); 

            % Loop through nodes
            for n = 1:obj.n_nodes
                % Get next node
                node = obj.nodes{n};
                % Only continue in case of a mechanical node
                if ~isa(node, 'PH_MechanicalNode') 
                    continue;
                end
                
                %% Translatory part
                % Separate constraints in each direction (x/y/z)
                for d = 1:3
                    dir = zeros(3, 1);
                    dir(d) = 1;
                    % Get all translational ports
                    ports_fv = obj.getMechanicalPorts(n, {'force', 'velocity'}, dir);
                    % Do not continue if the current DOF is locked or
                    % when no ports are connected to the current node
                    if  ~isempty(ports_fv) %&& ~node.lockedDOFs(d)
                        % If the only port connected at the node is an
                        % external one, no constraint is added
                        if length(ports_fv) == 1 && isa(obj.ports{ports_fv(1)}, 'PH_Port_external')
                            continue
                        end
                        
                        % All forces sum up to zero -> one constraint                      
                        nConstraints = nConstraints+1;
                        for p=1:length(ports_fv)
                            port = obj.ports{ports_fv(p)};
                            orientation = port.orientation(:,find(port.nodes==n,1));
                            
                            if isa(port, 'PH_Port_storage')
                                % For storage ports, add flow/effort constraints
                                if strcmp(port.type, 'force')
                                    % effort is force 
                                    Ce(nConstraints, port.effort(find(port.nodes == n,1))) = orientation(d);
                                else
                                    % flow is force
                                    Cf(nConstraints, port.flow(find(port.nodes == n, 1))) = orientation(d); 
                                end
                            elseif isa(port, 'PH_Port_external') || isa(port, 'PH_Port_boundary')
                                % For external/boundary ports, add input/output constraints
                                if strcmp(port.type, 'force')
                                    % input is force
                                    Cu(nConstraints, port.IOPair) = orientation(d);
                                else
                                    % ouput is velocity
                                    Cy(nConstraints, port.IOPair) = orientation(d); 
                                end
                            end
                        end      
                    end
                end
                %% Rotational part
                for d = 1:3
                    dir = zeros(3, 1);
                    dir(d) = 1;
                    % Get all rotational ports
                    ports_tr = obj.getMechanicalPorts(n, {'torque', 'angular velocity'}, dir);
                    % Do not continue if the current DOF is locked or
                    % when no ports are connected to the current node
                    if ~isempty(ports_tr) %&& ~node.lockedDOFs(d+3) 
                        % If the only port connected at the node is an
                        % external one, no constraint is added
                        if length(ports_tr) == 1 && isa(obj.ports{ports_tr(1)}, 'PH_Port_external')
                            continue
                        end
                        % All torques sum up to zero
                        nConstraints = nConstraints+1;
                        for p=1:length(ports_tr)
                            port = obj.ports{ports_tr(p)};
                            orientation = port.orientation(:,find(port.nodes==n,1));
                            if isa(port, 'PH_Port_storage')
                                % For storage ports, add flow/effort constraints
                                if strcmp(port.type, 'torque')
                                    % effort is torque 
                                    Ce(nConstraints, port.effort(find(port.nodes == n,1))) = orientation(d);
                                else
                                    % flow is torque
                                    Cf(nConstraints, port.flow(find(port.nodes == n, 1))) = orientation(d); 
                                end
                            elseif isa(port, 'PH_Port_external') || isa(port, 'PH_Port_boundary')
                                % For externa/boundary ports, add input/output constraints
                                if strcmp(port.type, 'torque')
                                    % input is torque
                                    Cu(nConstraints, port.IOPair) = orientation(d);
                                else
                                    % ouput is torque
                                    Cy(nConstraints, port.IOPair) = orientation(d); 
                                end
                            end
                        end
                    end
                end
            end

            % Ensure all arrays are of the same size 
            if size(Ce, 1) < nConstraints
                Ce(nConstraints, 1) = 0;
            end
            if size(Cf, 1) < nConstraints
                Cf(nConstraints, 1) = 0;
            end
            if size(Cu, 1) < nConstraints && obj.n_u > 0
                Cu(nConstraints, 1) = 0;
            end
            if size(Cy, 1) < nConstraints && obj.n_u > 0 
                Cy(nConstraints, 1) = 0;
            end
            % Remove constraints that are zero
            idx = ~any([Ce Cf Cu Cy]');
            Ce(idx,:) = [];
            Cf(idx,:) = [];
            Cu(idx,:) = [];
            Cy(idx,:) = [];
        end
    end 
end