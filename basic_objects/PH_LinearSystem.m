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
        B       % Algebraic constraints on efforts / Lagrange multiplier input mapping 
        
        C       % Non-collocated output mapping
        
        C_e     % Matrix for effort coupling
        C_f     % Matrix for flow coupling
        C_u     % Matrix for I/O coupling 
        C_y     % Matrix for I/O coupling
        
        x_p     % State variables associated with the generalized momenta
        x_q     % State variables associated with the generalized displacements
    end
    
    methods(Access = public)
        %              PH_LinearSystem(n, J, Q, G, R, P, S, M, B, x_p, x_q)
        function obj = PH_LinearSystem(n, J, Q, varargin)
            % Sanity checks
            obj = obj@PH_System('Linear port-Hamiltonian system');
           
            % Strict type checking
            if ~isscalar(n) || n < 0
                error('n must be >= 0');
            end
            obj.n = n;
            
            if ~ismatrix(J) || any(size(J) ~= n) || ~issymmetric(round(J, 16), 'skew')
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
                    if ~ismatrix(R) || any(size(R) ~= n) || ~issymmetric(round(R, 16)) || any(round(eig(R),17) < 0)
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
            
            % Constraints
            obj.C_e = zeros(0, obj.n);
            obj.C_f = zeros(0, obj.n);
            obj.C_u = zeros(0, obj.n_u);
            obj.C_y = zeros(0, obj.n_u);
            
            % Lagrange multiplier input matrix
            obj.B = zeros(obj.n, 0);
            if nargin > 8 
                B = varargin{6};
                if ~isempty(B)
                    if size(B, 1) ~= obj.n 
                        error('B must be a n x n_c matrix');
                    end
                    obj.B = B;
                end
            end
            obj.n_c = size(obj.B, 2);
            
            % Which states are generalized momenta and which
            % generalized displacements
            if nargin > 10
                x_p = varargin{7};
                x_q = varargin{8};
                if ~isvector(x_p) || ~isvector(x_q) || length(x_p) + length(x_q) ~= n
                    error('x_p and x_q must have a combined length of n')
                end
                obj.x_p = x_p;
                obj.x_q = x_q;
            else
                % We assume that the number of generalized momenta equals the
                % number of generalized displacements and that the generalized
                % momenta are given first
                obj.x_p = [1:floor(obj.n/2)]'; 
                obj.x_q = [floor(obj.n/2)+1:obj.n]'; 
            end

            obj.n_elements = 1;
            obj.elements{1} = PH_Element('Linear PH system', [], []);
        end
        
        % Check whether our System is a valid PH system
        function b = isValidPHSystem(obj)
            b = 1; 
            
            if ~issymmetric(round(obj.J, 10), 'skew') 
                disp('J is not skew-symmetric');
                b = 0;
            end
            
            E = eig(obj.R);
            if ~issymmetric(round(obj.R, 10)) || any(E < -1e-12 * max(abs(E)))
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
            E = eig(Z);
            if any(E < -1e-12*max(abs(E)))
                disp('Z = [R P; P'' S] is not positive semidefinite');
                b = 0;
            end
        end
        
        % Add another linear PH system 
        function add(obj, system)
            if ~isa(system, 'PH_LinearSystem') || ~system.isValidPHSystem()
                error('system must be a valid ''PH_LinearSystem''');
            end
            
            states = ones(system.n, 1);
            if (sum(states(system.x_p)) + sum(states(system.x_q))) ~= system.n
                error('system must be separable into two energy domains');
            end
            
            N_p = length(system.x_p);
            N = obj.n;            
            N_q = length(system.x_q);
            
            obj.x_p = [system.x_p; obj.x_p + N_p];
            obj.x_q = [obj.x_q + N_p; system.x_q + N];
            
            obj.J = [system.J(system.x_p,system.x_p) zeros(N_p, N) system.J(system.x_p, system.x_q);
                  zeros(N, N_p) obj.J zeros(N, N_q);
                  system.J(system.x_q, system.x_p) zeros(N_q, N) system.J(system.x_q, system.x_q)];
            obj.Q = [system.Q(system.x_p,system.x_p) zeros(N_p, N) system.Q(system.x_p, system.x_q);
                  zeros(N, N_p) obj.Q zeros(N, N_q);
                  system.Q(system.x_q, system.x_p) zeros(N_q, N) system.Q(system.x_q, system.x_q)];
            
            obj.G = [zeros(N_p, obj.n_u), system.G(system.x_p,:);
                  obj.G, zeros(N, system.n_u);
                  zeros(N_q, obj.n_u)  system.G(system.x_q, :)];
            obj.P = [zeros(N_p, obj.n_u), system.P(system.x_p,:);
                  obj.P, zeros(N, system.n_u);
                  zeros(N_q, obj.n_u)  system.P(system.x_q, :)];
            obj.B = [zeros(N_p, size(obj.B,2)), system.B(system.x_p,:);
                  obj.B, zeros(N, size(system.B, 2));
                  zeros(N_q, size(obj.B, 2))  system.B(system.x_q, :)];
              
            obj.M = blkdiag(obj.M, system.M);
            obj.S = blkdiag(obj.S, system.S); 
            
            if ~isempty(system.R)
                obj.R = [system.R(system.x_p,system.x_p) zeros(N_p, N) system.R(system.x_p, system.x_q);
                  zeros(N, N_p) obj.R zeros(N, N_q);
                  system.R(system.x_q, system.x_p) zeros(N_q, N) system.R(system.x_q, system.x_q)];
            else
                obj.R = [zeros(N_p, N_p+N_q+N);
                      	zeros(N, N_p) obj.R zeros(N, N_q);
                        zeros(N_q, N_p+N_q+N)];
            end
            
            obj.C_e = [system.C_e(:,system.x_p) zeros(size(system.C_e, 1), N) system.C_e(:, system.x_q);
                        zeros(size(obj.C_e, 1), N_p) obj.C_e zeros(size(obj.C_e, 1), N_q)];
            obj.C_f = [system.C_f(:,system.x_p) zeros(size(system.C_f, 1), N) system.C_f(:, system.x_q);
                        zeros(size(obj.C_f, 1), N_p) obj.C_f zeros(size(obj.C_f, 1), N_q)];
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
            
            % Look for identical nodes and delete duplicates
            duplicateNodes = zeros(0,2);
            for n = 1:obj.n_nodes
                current = obj.nodes{n};
                if ~isa(current, 'PH_MechanicalNode') || any(any(duplicateNodes == n)) 
                    continue;
                else
                    for i=1:obj.n_nodes
                        if ~(i == n) && isa(obj.nodes{i}, 'PH_MechanicalNode')
                            if sum(abs(obj.nodes{i}.location - current.location)) < 1e-14
                                duplicateNodes(end+1, 1) = n;
                                duplicateNodes(end, 2) = i;
                                break;
                            end
                        end
                    end
                end
            end
            
            if ~isempty(duplicateNodes)
                for n = 1:size(duplicateNodes, 1)
                    current = obj.nodes{duplicateNodes(n, 1)};
                    current.ports = [current.ports obj.nodes{duplicateNodes(n, 2)}.ports];
                    current.elements = [current.elements obj.nodes{duplicateNodes(n, 2)}.elements];
                    current.lockedDOFs = double(current.lockedDOFs | obj.nodes{duplicateNodes(n, 2)}.lockedDOFs);
                    for p = 1:length(current.ports)
                        for pn = 1:length(obj.ports{current.ports(p)}.nodes)
                            if(obj.ports{current.ports(p)}.nodes(pn) == duplicateNodes(n, 2))
                                obj.ports{current.ports(p)}.nodes(pn) = duplicateNodes(n, 1);
                            end
                        end
                    end
                    for e = 1:length(current.elements)
                        for en = 1:length(obj.elements{current.elements(e)}.nodes)
                            if(obj.elements{current.elements(e)}.nodes(en) == duplicateNodes(n, 2))
                                obj.elements{current.elements(e)}.nodes(en) = duplicateNodes(n, 1);
                            end
                        end
                    end
                    obj.nodes{duplicateNodes(n, 2)}.delete();
                end
                
                for n = 1:size(duplicateNodes, 1)            
                    for p = 1:obj.n_ports
                        port = obj.ports{p};
                        for pn = 1:length(port.nodes)
                            if port.nodes(pn) >= duplicateNodes(n, 2)
                                port.nodes(pn) = port.nodes(pn) -1; 
                            end
                        end
                    end
                    for e = 1:obj.n_elements
                        element = obj.elements{e};
                        for en = 1:length(element.nodes)
                            if element.nodes(en) > duplicateNodes(n, 2)
                                element.nodes(en) = element.nodes(en) -1; 
                            end
                        end
                    end
                    obj.nodes(duplicateNodes(n, 2)) = [];
                    obj.n_nodes = length(obj.nodes);
                    for i=n+1:size(duplicateNodes,1)
                        if duplicateNodes(i, 2) > duplicateNodes(n, 2) 
                            duplicateNodes(i, 2) = duplicateNodes(i, 2) -1;
                        end
                    end
                end 
            end
        end
               
        % Add Rayleigh damping to the system
        % Method can be found in the PhD thesis of Flavio Luiz Cardoso-Ribeiro
        % NOTE: works only if there is no direct coupling b/w energy variables
        function addRayleighDamping(obj, alpha, beta)
            if any(obj.x_p == obj.x_q)
                error('not applicable when there is direct coupling between the energy variables');
            end
            
            Q_1 = obj.Q(obj.x_p, obj.x_p);
            Q_2 = obj.Q(obj.x_q, obj.x_q);
            J_12 = obj.J(obj.x_p, obj.x_q);
            J_21 = obj.J(obj.x_q, obj.x_p);
            
            D = alpha*(Q_1)^-1 - beta * J_12*Q_2*J_21;
            D = D - 0.5*(D - D');
            obj.R(obj.x_p, obj.x_p) = D;       
        end
        
        % Get the system's first n smallest magnitude eigenvalues
        % NOTE: works only if there is no direct coupling b/w energy variables
        function [lambda, K] = getSmallestMagnitudeEigenvalues(obj, n)
            % TODO: check the system matrices instead of the p/q indices...
            if any(ismember(obj.x_p, obj.x_q))
                error('not applicable when there is direct coupling between the energy variables');
            end
 
            Q_1 = obj.Q(obj.x_p, obj.x_p);
            Q_2 = obj.Q(obj.x_q, obj.x_q);
            J_12 = obj.J(obj.x_p, obj.x_q);
            J_21 = obj.J(obj.x_q, obj.x_p);
            
            K = -J_21*Q_1*J_12*Q_2;
            
            K = K - 0.5*(K-K');
            
            if rank(K) < size(K, 1)
                error('system is singular - eigenvalues cannot be computed');
            end
            
            %[L,U,p] = lu(K,'vector');
            %Afun = @(x) U\(L\(x(p)));
            opts.v0 = ones(size(K, 1),1);
            opts.isreal = 1;
            %opts.issym = 1;
            %opts.tol = 1e-20;
            %opts.p = floor(size(K, 1)/2);
            %opts.disp = 1;
            [~, lambda] = eigs(K, n, 'smallestabs',opts); 
            lambda = diag(lambda);
        end
        
        % Get a subsystem with the given states x_p and x_q.
        % Used for decentralized approaches (substructuring)
        function sys = getSubsystem(obj, x_p, x_q)
            sys = obj.copy();
            states = [x_p; x_q];
            sys.n = length(states);
            sys.J = sys.J(states, states);
            sys.R = sys.R(states, states);
            sys.Q = sys.Q(states, states);
            sys.G = sys.G(states, :);
            sys.P = sys.P(states, :);            
            % Remove inputs that do not influence the subsystem
            idx_u0 = find(~any(round(sys.G-sys.P, 12)));
            u_dep = sort(idx_u0, 'descend');
            ports = sys.getPortsForIOPairs(u_dep);
            sys.deletePorts(ports);
            
            [~, idx_a] = unique(sys.G' - sys.P', 'rows');
            idx_u0 = find(~ismember(1:size(sys.G,2), idx_a));
            u_dep = sort(idx_u0, 'descend');
            ports = sys.getPortsForIOPairs(u_dep);
            sys.deletePorts(ports);
            
            sys.B = sys.B(states, :);
            sys.C_e = sys.C_e(:, states);
            sys.C_f = sys.C_f(:, states); 
            sys.n_c = size(sys.C_u, 1) + size(sys.B, 2);
            
            sys.x_p = [1:length(x_p)]';
            sys.x_q = [length(x_p)+1:length(x_p)+length(x_q)]';
            
            % Delete all nodes and elements that do not belong to the
            % subsystem
            for n = sys.n_nodes:-1:1
                ports = sys.getMechanicalPorts(n);
                if isempty(ports)
                    % Search elements that belong to this node
                    for e = sys.n_elements:-1:1
                        element = sys.elements{e};
                        if any(element.nodes == n)
                            sys.elements{e}.delete();
                            sys.elements(e) = [];
                            sys.n_elements = sys.n_elements-1;
                        end
                    end
                    sys.nodes{n}.delete();
                    sys.nodes(n) = [];
                    sys.n_nodes = sys.n_nodes - 1;
                end
            end
        end
        
        % Round all dynamic system matrices (for visualization purposes)
        function roundSystemMatrices(obj, decimalPlaces)
            obj.J = round(obj.J, decimalPlaces);
            obj.R = round(obj.R, decimalPlaces);
            obj.Q = round(obj.Q, decimalPlaces);
            obj.G = round(obj.G, decimalPlaces);
            obj.P = round(obj.P, decimalPlaces);
            obj.M = round(obj.M, decimalPlaces);
            obj.S = round(obj.S, decimalPlaces);
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
        
        % Generate constraints of the form Cy * y + Cu * u = 0.
        % Constraints are generated automatically for each node. Only the
        % mechanical domain is supported so far. Hydraulic actuators are
        % coupled via their mechanical ports.
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
       
        % Eliminates all constraints of the form B'*e = 0. Such constraints
        % are usually generated by a call to assembleDAESystem()
        % Setting useRationalBasis to 1 uses a rational basis for the null
        % space of B. This preserves separation of energy domains in some
        % cases but is less numerically stable... see for yourself.
        function [V_xz] = eliminateAlgebraicConstraints(obj, useRationalBasis)
            if size(obj.B, 2) == 0
                % No algebraic constraints --> nothing to do!
                V_xz = eye(obj.n);
                return
            end
            
            % Get number of algebraic constraints
            n_ac = size(obj.B, 2);
            if rank(obj.B) ~= n_ac
                obj.B = rref(obj.B')';
                obj.B(:,~any(obj.B)) = [];
                n_ac = size(obj.B, 2);
            end
            % Construct a left annihilator such that G_a*B = 0
            if nargin > 1 && useRationalBasis
                G_a = null(obj.B', 'r')';
            else
                G_a = null(obj.B')';
            end
            
            % Compute transformation matrix
            V = [G_a; (obj.B'*obj.B)\obj.B'];
            
            % Transform system matrices
            J_ = V*obj.J*V';
            R_ = V*obj.R*V';
            Q_ = inv(V)'*obj.Q*inv(V);
            G_ = V*obj.G;
            P_ = V*obj.P;
            
            % Q_ must be positive definite 
            if any(eig(Q_) <= 0)
                error('Positive definiteness of Q lost after constraint elimination... aborting');
            end

            % Partition Q (Q2x will be eliminated) 
            Q_11 = Q_(1:obj.n-n_ac, 1:obj.n-n_ac);
            Q_12 = Q_(1:obj.n-n_ac, obj.n-n_ac+1:end);
            Q_22 = Q_(obj.n-n_ac+1:end, obj.n-n_ac+1:end);
            Q_21 = Q_(obj.n-n_ac+1:end, 1:obj.n-n_ac); 
            
            x_p_ = zeros(obj.n, 1);
            x_p_(obj.x_p) = 1;
            x_q_ = zeros(obj.n, 1);
            x_q_(obj.x_q) = 1;           
            
            % System order is reduced by the number of constraints
            obj.n = obj.n - n_ac;
            
            % Map from new to old state variables
            V_xz = V\[eye(obj.n); -Q_22\Q_21];
            obj.x_p = find(abs(V_xz\x_p_) > 1e-10);
            obj.x_q = find(abs(V_xz\x_q_) > 1e-10); 
            
            obj.J = J_(1:obj.n, 1:obj.n);
            obj.R = R_(1:obj.n, 1:obj.n);
            obj.Q = Q_11 - Q_12*(Q_22\Q_21);
            obj.G = G_(1:obj.n, :);     
            obj.P = P_(1:obj.n, :); 
            % We know that this method preserves structural properties...
            % make sure they are not lost due to the numerics involved
            obj.J = obj.J - 0.5*(obj.J+obj.J');
            obj.R = obj.R - 0.5*(obj.R-obj.R');
            obj.Q = obj.Q - 0.5*(obj.Q - obj.Q');
            
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
        
        % Check whether there exist linear dependencies between sets of
        % effort variables. If yes, derive an orthonormal basis for and
        % eliminate them by calling eliminateAlgebraicConstraints
        function [V_xz] = eliminateLinearDependencies(obj, useRationalBasis)
            V_xz = eye(obj.n);
            if size(obj.B, 2) > 0
                error('There are still algebraic constraints on the system ... eliminate those first');
            end
            
            % Check whether there is a block of linearly dependent flows
            row_idx = ~any(round(obj.G'-obj.P', 12),1);
            col_idx = any(obj.J(row_idx, :)-obj.R(row_idx,:));
            
            J_sorted = [obj.J(:, col_idx == 1), obj.J(:, col_idx ~=1)]; 
            J_sorted = [J_sorted(row_idx ~= 1, :); J_sorted(row_idx == 1, :)];
            R_sorted = [obj.R(:, col_idx == 1), obj.R(:, col_idx ~=1)]; 
            R_sorted = [R_sorted(row_idx ~= 1, :); R_sorted(row_idx == 1, :)];
            G_sorted = [obj.G(row_idx ~= 1, :); obj.G(row_idx == 1, :)];
            P_sorted = [obj.P(row_idx ~= 1, :); obj.P(row_idx == 1, :)];
            
            row_idx = ~any(round(G_sorted'-P_sorted', 12), 1);
            col_idx = any(round(J_sorted(row_idx, :) - R_sorted(row_idx, :), 12));
            
            if sum(row_idx) > sum(col_idx)
                % If such a block was found, construct an orthonormal basis 
                % to eliminate dependent states
                A_b = J_sorted(row_idx == 1, :) - R_sorted(row_idx == 1, :);
                [U,~] = svd(A_b,'econ');
                T_b = U(:, 1:sum(row_idx == 0));
                
                %T_b = orth(A_b);
                T = blkdiag(eye(obj.n-size(T_b, 1)), T_b);
                
                % Construct algebraic constraints B_ for linearly dependent
                % states
                [~,lambda,V] = svd((obj.Q*T)', 0);
                lambda = diag(lambda);
                lambda = lambda / max(lambda);
                if length(lambda) < size(V, 1)
                    lambda(size(V,1)) = 0;
                end
                obj.B = V(:, lambda < 1e-18);
                obj.n_c = obj.n_c + size(obj.B, 2);
                if nargin > 1 && useRationalBasis
                    V_xz = obj.eliminateAlgebraicConstraints(1);
                else 
                    V_xz = obj.eliminateAlgebraicConstraints();
                end
            end
        end
        
        % Perform A modal analysis on the system and truncate high
        % frequency eigenmodes
        function [T, V] = transformGeneralizedCoordinates(obj, U)
            M_ = inv(obj.Q(obj.x_p, obj.x_p));
            M_ = M_ - 0.5*(M_-M_');
            %obj.Q(obj.x_p, obj.x_q) = 0;
            %obj.Q(obj.x_q, obj.x_p) = 0;
            
            T = pinv(blkdiag(M_*U, U));
            V = blkdiag(M_*U, U);
            obj.Q = pinv(blkdiag(U, M_*U))*obj.Q*V;
            obj.Q = obj.Q - 0.5*(obj.Q-obj.Q');
            obj.J = T*obj.J*blkdiag(U, M_*U);
            obj.J = obj.J - 0.5*(obj.J+obj.J');
            obj.R = T*obj.R*blkdiag(U, M_*U);
            obj.R = obj.R - 0.5*(obj.R-obj.R');
            obj.G = T*obj.G;
            obj.P = T*obj.P;
            obj.B = T*obj.B;
            obj.n = 2*size(U, 2);
            obj.x_p = [1:obj.n/2]';
            obj.x_q = [obj.n/2+1:obj.n]';
        end
        
        % Transform the state space by formulation of algebraic
        % constraints. T is the transformation matrix between state spaces.
        function V_xz = applyStateTransformation(obj, T)
            % Construct algebraic constraints B_ for linearly dependent
            % states
            [~,lambda,V] = svd(T*obj.Q, 0);
            lambda = diag(lambda);
            lambda = lambda / max(lambda);
            if length(lambda) < size(V, 1)
                lambda(size(V,1)) = 0;
            end
            obj.B = V(:, lambda < 1e-18);
            obj.n_c = obj.n_c + size(obj.B, 2);
            V_xz = obj.eliminateAlgebraicConstraints();
        end
        
        % Designed for mechanical systems where the global dofs are given
        % as x = [M dq,  q] for a system M ddq + D dq + K q = B u
        % NOTE: doens't work for hypostatic systems
        function transformToGlobalDOFs(obj, skip2ndStep)      
            row_idx = any(obj.G'-obj.P', 1);
            G_a = null(obj.G(row_idx,:)');
            T = blkdiag([(obj.G(row_idx, :)'*obj.G(row_idx,:))\obj.G(row_idx,:)'; G_a'], eye(obj.n-sum(row_idx)));
            %T = blkdiag([pinv(obj.G); G_a'], eye(obj.n-sum(row_idx)));
            if rank(T) < obj.n
                error('system is singular (maybe hypostatic?)... aborting')
            end

            % Left mulitply the system with T. J and R need to be right
            % multiplied with T' to preserve symmetry properties
            obj.J = T*obj.J*T';
            obj.R = T*obj.R*T';
            % This transformation affects e_p only --> turns to M^-1
            % Q needs to be left multiplied with the inverse of T' and
            % right multiplied with the inverse of T
            obj.Q = inv(T)'*obj.Q*inv(T);
            obj.G = T*obj.G;
            obj.P = T*obj.P;
            obj.B = T*obj.B;
            
            obj.x_p = [1:obj.n/2]';
            obj.x_q = [obj.n/2+1:obj.n]';
            
            if nargin > 1 && skip2ndStep
                return
            end
            
            T = blkdiag(eye(obj.n/2), pinv(obj.J(obj.n/2+1:end, 1:obj.n/2)));

            if rank(T) < obj.n
                error('system is singular (maybe hypostatic?)... aborting')
            end

            % Again, J and R needs to be right multiplied by T'
            obj.J = T*obj.J*T';
            obj.R = T*obj.R*T';
            % This transformation affects e_q only --> turns to K
            obj.Q = inv(T)'*obj.Q*inv(T);
            obj.G = T*obj.G;
            obj.P = T*obj.P;
            obj.B = T*obj.B;

            % We know that this method preserves structural properties...
            % make sure they are not lost due to the numerics involved
            obj.J = obj.J - 0.5*(obj.J+obj.J');
            obj.R = obj.R - 0.5*(obj.R-obj.R');
            obj.Q = obj.Q - 0.5*(obj.Q - obj.Q');
        end
        
        % Assuming the system is in some kind of balanced form, use the
        % effort-constraint method to eliminate the efforcts e_c
        % NOTE: implemented according to Polyuga & v. d. Schaft 2011
        % NOTE: MOR is short for model order reduction
        function effortConstraintMOR(obj, e_2)
            if islogical(e_2)
                e_1 = ~e_2;
            else
                e_1 = 1:obj.n;
                e_1(ismember(e_1, e_2)) = [];
            end
            
            obj.J = obj.J(e_1, e_1);
            obj.R = obj.R(e_1, e_1);
            obj.Q = obj.Q(e_1, e_1)-obj.Q(e_1, e_2)*(obj.Q(e_2, e_2)\obj.Q(e_2, e_1));
            obj.G = obj.G(e_1, :);
            obj.P = obj.P(e_1, :);
            obj.M = obj.M(e_1, e_1);
            obj.S = obj.S(e_1, e_1);
        end
        
        % Assuming the system is in some kind of balanced form, use the
        % flow-constraint method to elminate the flows f_c
        % NOTE: implemented according to Polyuga & v. d. Schaft 2011
        % NOTE: MOR is short for model order reduction
        function flowConstraintMOR(obj, f_2)
            if nnz(obj.P) || nnz(obj.M) || nnz(obj.S)
                error('not implemented for systems with input feedthrough');
            end
            
            if islogical(f_2)
                f_1 = ~f_2;
            else
                f_1 = 1:obj.n;
                f_1(ismember(f_1, f_2)) = [];
            end
            
            J_11 = obj.J(f_1, f_1);
            J_12 = obj.J(f_1, f_2);
            J_21 = obj.J(f_2, f_1);
            J_22 = obj.J(f_2, f_2); 
            G_1 = obj.G(f_1, :);
            G_2 = obj.G(f_2, :);
            
            alpha = G_2'*J_22\J_21 - G_1';
            beta = J_22\J_21 - eye(size(J_11, 1));
            gamma = G_2'*J_22^-1;
            delta = J_22^-1;
            eta = G_2'*J_22\G_2;
            
            Z = obj.R*(eye(obj.n) - delta*obj.R)^-1;
            Z_sym = 0.5*(Z+Z');
            Z_sk = 0.5*(Z-Z');
            J_s = J_11 - J_12*J_22\J_21;
            
            obj.J = J_s - beta'*Z_sk*beta;
            obj.R = beta'*Z_sym*beta;
            obj.G = -alpha'+beta'*Z_sk*gamma';
            
            % Flow-constraint method produces feedthrough for G_2 ~= 0
            obj.P = -beta'*Z_sym*gamma';
            obj.M = -eta+gamma*Z_sk*gamma';
            obj.S = gamma*Z_sym*gamma';
        end
        
        % Provides the possibility to replace the input matrix G by a
        % user-specified one. Deletes all existing inputs.
        function setInputMatrix(obj, G_)
            if size(G_, 1) ~= obj.n
                error('input matrix must have n rows');
            end
            % Delete all ports associated with current input matrix 
            obj.deletePorts(obj.getPortsForIOPairs(1:obj.n_u));
            obj.G = G_; 
            obj.n_u = size(obj.G, 2);
            
            obj.P = zeros(size(obj.G));
            obj.M = zeros(obj.n_u);
            obj.S = zeros(obj.n_u);
        end
        
        % Returns the input names of the inputs specified by IOPairs in the
        % form of a cell array of strings
        function names = getInputNames(obj, IOPairs)
            names = cell(1, length(IOPairs));
            ports = obj.getPortsForIOPairs(IOPairs);
            for p = 1:length(ports)
                names{p} = obj.ports{ports(p)}.inputName;
            end
        end
        
        % Returns the system ODEs as an anonymous function. In the presence
        % of algebraic constraints, the lagrange multipliers are calculated
        % to obtain an explicit representation.
        function [f, dfdx] = getSystemODEfun(obj)
            if obj.isConstrained()
                lambda = @(x, u) -(obj.B'*obj.Q*obj.B)\obj.B'*obj.Q*((obj.J-obj.R)*obj.Q*x + (obj.G - obj.P)*u);
                f = @(x, u) (obj.J-obj.R)*obj.Q*x + (obj.G - obj.P)*u + obj.B*lambda(x, u);
                dfdx = (obj.J-obj.R)*obj.Q - obj.B*((obj.B'*obj.Q*obj.B)\obj.B')*obj.Q*(obj.J-obj.R)*obj.Q;
            else
                f = @(x, u) (obj.J-obj.R)*obj.Q*x + (obj.G - obj.P)*u;
                dfdx = (obj.J-obj.R)*obj.Q;
            end
        end
        
        % Implicit version of the system dynamics
        function [f, dfdt, x0, dx0] = getSystemDAEfun(obj, x0, u0)
            if obj.isConstrained()
                f = @(x, u) [(obj.J-obj.R)*obj.Q*x(1:obj.n) + (obj.G - obj.P)*u + obj.B*x(obj.n+1:obj.n+size(obj.B,2)); obj.B'*obj.Q*obj.x(1:obj.n)];
                dfdt = blkdiag(eye(obj.n), zeros(size(obj.B,2)));
                % Consistent initial condition
                lambda0 = -(obj.B'*obj.Q*obj.B)\obj.B'*obj.Q*((obj.J-obj.R)*obj.Q*x0 + (obj.G - obj.P)*u0);
                x0 = [x0; lambda0];
                dx0 = f(x0, u0);
            else
                f = @(x, u) (obj.J-obj.R)*obj.Q*x + (obj.G - obj.P)*u;
                dfdt = eye(obj.n);
                dx0 = f(x0, u0);
            end
        end
        
        % For symplectic integrators (symplectic Euler)
        function [f_p, f_q] = getSymplecticDEs(obj)
            % TODO: check the system matrices instead of the p/q indices...
            if any(ismember(obj.x_p, obj.x_q))
                error('Not applicable when there is direct coupling between the energy variables');
            end
            if obj.isConstrained()
                error('Separation into energy domains not implemented for constrained systems')
            end
            
            f_p = @(q, u) (obj.J(obj.x_p, obj.x_q)-obj.R(obj.x_p, obj.x_q))*obj.Q(obj.x_q, obj.x_q)*q + (obj.G(obj.x_p, :) - obj.P(obj.x_p, :))*u;
            f_q = @(p, u) (obj.J(obj.x_q, obj.x_p)-obj.R(obj.x_q, obj.x_p))*obj.Q(obj.x_p, obj.x_p)*p + (obj.G(obj.x_q, :) - obj.P(obj.x_q, :))*u;
        end
        
        % Returns the system output for a given state and input vector
        function y = getSystemOutput(obj, x, u)
            if nargin < 2
                error('x must be given to calculate y');
            end
            if ~iscolumn(x) || length(x) ~= obj.n
                error(['x must be a column vector of length ' num2str(obj.n)]);
            end
            if obj.isConstrained()
                % TODO: is it really correct to include the Lagrange
                % multipliers here?
                y = (obj.G' + obj.P' - ((obj.B'*obj.Q*obj.B)\obj.B'*obj.Q*obj.G)'*obj.B')*obj.Q*x + (obj.M + obj.S)*u;
            else
                y = ((obj.G + obj.P)' * obj.Q * x + (obj.M + obj.S) * u);
            end
        end
        
        % Computs the Hamiltonian for a given state vector
        function H = getHamiltonian(obj, x)
            if nargin < 2
                error('x must be given to calculate z');
            end
            if ~iscolumn(x) || length(x) ~= obj.n
                error(['x must be a column vector of length ' num2str(obj.n)]);
            end
            H = (x' * obj.Q * x) / 2;
        end
        
        % Searches the system's mechanical nodes for the one at the
        % specified location
        function node = getNodeAtLocation(obj, location)
            node = 0;
            for n=1:obj.n_nodes
                if isa(obj.nodes{n}, 'PH_MechanicalNode') && sum(abs(obj.nodes{n}.location - location)) < 1e-14
                    node = n;
                    break;
                end
            end
        end
        
        % Fix/lock the specified DOFs of the node indicated by location
        function fixNodeDOFs(obj, location, dofs)
            if ~ismatrix(location) || size(location, 1) ~= 1 || size(location, 2) ~= 3
                error('location must be supplied as a 1x3 vector');
            end
            if size(dofs, 1) ~= 1 || size(dofs, 2) ~= 6
                error('dofs must be an 1x6 vector'); 
            end
            
            node = obj.getNodeAtLocation(location);
            obj.nodes{node}.lockedDOFs = dofs;
        end 
        
        % Adds an external input (force) at each mechanical node
        function addExternalInputsAtNodes(obj)
            % Original number of inputs
            n_u0 = obj.n_u;
            for n=1:obj.n_nodes
                node = obj.nodes{n};
                if ~isa(node, 'PH_MechanicalNode') || node.internal
                    continue
                end
                
                % Get all ports at this node
                dir_str = 'xyz';
                for d = 1:3
                    dir = zeros(1, 3);
                    dir(d) = 1;
                    % ports are connected at this port in this direction?
                    ports_f = obj.getMechanicalPorts(n, {'force'}, dir);
                    if ~isempty(ports_f)
                        obj.addExternalPort('mechanical', 'force', n, obj.n_u+1, ...
                                            ['n' num2str(n) '_F' dir_str(d)], ['n' num2str(n) '_v' dir_str(d)], dir');
                        obj.n_u = obj.n_u +1;
                    end
                end
                for d = 1:3
                    dir = zeros(1, 3);
                    dir(d) = -1;

                    ports_t = obj.getMechanicalPorts(n, {'torque'}, dir);
                    if ~isempty(ports_t)
                        obj.addExternalPort('mechanical', 'torque', n, obj.n_u+1, ...
                                            ['n' num2str(n) '_M' dir_str(d)], ['n' num2str(n) '_w' dir_str(d)], dir');
                        obj.n_u = obj.n_u +1;
                    end
                end
            end
            
            obj.G = [obj.G, zeros(obj.n, obj.n_u - n_u0)];
            obj.P = [obj.P, zeros(obj.n, obj.n_u - n_u0)];
            obj.M = blkdiag(obj.M, zeros(obj.n_u - n_u0));
            obj.S = blkdiag(obj.S, zeros(obj.n_u - n_u0));
            obj.C_u = [obj.C_u, zeros(size(obj.C_u, 1), obj.n_u - n_u0)];
            obj.C_y = [obj.C_y, zeros(size(obj.C_y, 1), obj.n_u - n_u0)];
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
        
        % Returns the ports corresponding to the inputs specified by pairs
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
        
        % You can add multiply ports of the same type at once
        function addPorts(obj, ports)
            if ~isa(ports, 'PH_Port')
                error('port has to be a member of class PH_Port');
            end
            for p = 1:length(ports)
                port = ports(p);
                obj.ports{obj.n_ports+1} = port; 

                % Append to node's ports
                for n = 1:length(port.nodes)
                    obj.nodes{port.nodes(n)}.ports = [obj.nodes{port.nodes(n)}.ports obj.n_ports+1];
                end

                obj.n_ports = obj.n_ports+1;
            end
        end
        
        % Delete ports by numbers
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
            obj.C_u(:, IOPairs) = [];
            obj.C_y(:, IOPairs) = [];
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
        % Internal function used for adding external ports
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
            cp.C_e = obj.C_e;
            cp.C_f = obj.C_f;
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
                    % skip ports located at mutliple nodes...
                    if(length(port.nodes) > 1)
                        continue
                    end
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
            C_ = G_v*[Ce, Cf, Cu, Cy];
            Ce = C_(:, 1:obj.n);
            Cf = C_(:, obj.n+1:2*obj.n);
            Cu = C_(:, 2*obj.n+1:2*obj.n+obj.n_u);
            Cy = C_(:, 2*obj.n+obj.n_u+1:2*obj.n+2*obj.n_u);
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