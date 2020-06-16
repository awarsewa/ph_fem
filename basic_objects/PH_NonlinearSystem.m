classdef PH_NonlinearSystem < PH_System
    %PH_NONLINEARSYSTEM Nonlinear port-Hamiltonian System in input-state-output
    %form
    
    properties(SetAccess = protected, GetAccess = public)
        x       % System state vector
        J       % Structure matrix J(x)
        R       % Resistance matrix R(x)
        H       % System Hamiltonian H(x)
        G       % Input matrix symmetric part G(x)
        B       % Algebraic constraints on efforts / Lagrange multiplier input mapping 
        
        C_u     % Matrix for I/O coupling 
        C_y     % Matrix for I/O coupling
    end
    
    methods(Access = public)
        %              PH_NonlinearSystem(n, x, J, H, G, R, B)
        function obj = PH_NonlinearSystem(n, x, J, H, varargin)
            % Sanity checks
            obj = obj@PH_System('Nonlinear port-Hamiltonian system');
           
            % Strict type checking
            if ~isscalar(n) || n < 0
                error('n must be >= 0');
            end
            obj.n = n;       
            
            if ~isa(x, 'casadi.MX') || length(x) ~= n
                error('x must be a casadi.MX object with length n');
            end
            obj.x = x;
            
            % TODO: Skew-symmetry check...
            if ~isa(J, 'casadi.MX') || any(size(J) ~= n)  
                error('J must be a skew-symmetric nxn casadi matrix function');
            end
            obj.J = J;
            
            if ~isa(H, 'casadi.MX') || any(size(H) ~= 1)
                error('H must be a scalar casadi.MX object');
            end
            obj.H = H;
            
            % Inputs and Outputs
            obj.G = casadi.MX.zeros(obj.n, 0);
            if nargin > 4 && ~isempty(varargin{1})
                G = varargin{1};
                if ~isa(G, 'casadi.MX') || size(G, 1) ~= n 
                    error('G must be a casadi.MX object of size n x n_u');
                end            
                obj.G = G; 
            end
            obj.n_u = size(obj.G,2);
            
            % Dissipation
            obj.R = casadi.MX.zeros(n, n);
            if nargin > 5 && ~isempty(varargin{2})
                R = varargin{2};
                if ~isa(R, 'casadi.MX') || any(size(R) ~= n)
                    error('R must be a symmetric positive semidefinite nxn casadi.MX object');
                end
                obj.R = R;
            end
            
            % Algebraic constraints
            obj.B = casadi.MX.zeros(n, 0);
            if nargin > 6 && ~isempty(varargin{3})
                B = varargin{3};
                if ~isa(B, 'casadi.MX') || size(B, 1) ~= n
                    error('B must be a casadi.MX matrix with n rows');
                end
                obj.B = B;
            end
            obj.n_c =  size(obj.B, 2);
                                  
            obj.C_u = zeros(0, obj.n_u);
            obj.C_y = zeros(0, obj.n_u);

            obj.n_elements = 1;
            obj.elements{1} = PH_Element('Nonlinear PH system', [], []);
        end
                
        % Add another linear PH system 
        function add(obj, system)
            if ~isa(system, 'PH_NonlinearSystem')
                error('system must be a ''PH_NonlinearSystem''');
            end
            
            % New x
            obj.x = [obj.x; system.x];
            % Concatenate matrix functions
            obj.J = blkdiag(obj.J, system.J);
            obj.R = blkdiag(obj.R, system.R); 
            obj.H = obj.H + system.H;
            obj.G = blkdiag(obj.G, system.G);
            obj.B = blkdiag(obj.B, system.B);
            
            obj.C_u = blkdiag(obj.C_u, system.C_u);
            obj.C_y = blkdiag(obj.C_y, system.C_y);
            obj.n_c = size(obj.C_u, 1) + size(obj.B, 2);  
            
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
                if ~isa(current, 'PH_MechanicalNode') || any(any(duplicateNodes == n)) || current.internal
                    continue;
                else
                    for i=1:obj.n_nodes
                        if ~(i == n) && isa(obj.nodes{i}, 'PH_MechanicalNode')
                            if ~any(obj.nodes{i}.location ~= current.location)
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
        
        function generateConstraints(obj)
            % Mechanical constraint generation
            [Cu_v, Cy_v] = obj.generateVelocityConstraints();           
            [Cu_f, Cy_f] = obj.generateForceConstraints();

            % Add constraint matrices
            obj.C_u = [obj.C_u; Cu_v; Cu_f];
            obj.C_y = [obj.C_y; Cy_v; Cy_f];
            
            % Adjust the number of constraints
            obj.n_c = size(obj.C_u, 1) + size(obj.B, 2); 
            
            % Constraints were generated successfully. To obtain an ODE system,
            % call assemble() followed by eliminateAlgebraicConstraints()
        end
        
        function assembleDAESystem(obj)
            if obj.n_c == 0
                % No constraints --> nothing to do
                return
            end

            % Get algebraic constraints
            idx_ac = ~any(obj.C_u');
            if ~isempty(obj.B)
                obj.B = [obj.B'; obj.C_y(idx_ac, :) * obj.G']';
            else
                obj.B = (obj.C_y(idx_ac, :) * obj.G')';
            end
            obj.C_y(idx_ac, :) = [];
            obj.C_u(idx_ac, :) = [];
            
            if nnz(obj.C_y)
                error('mixed input/output constraints not yet supported');
            end
            
            % Solve remaining equations for u
            u_f = ones(1, obj.n_u);
            A_uu = eye(obj.n_u);
            for i = 1:size(obj.C_u, 1)
                idx_u = find(obj.C_u(i, :), 1);
                if ~isempty(idx_u)
                    % Constrain this input
                    u_f(idx_u) = 0;
                    a_uu = obj.C_u(i, :)*A_uu;
                    b = a_uu(idx_u);
                    a_uu(idx_u) = 0;
                    A_uu(idx_u, :) = -b .* a_uu;

                    % Update dependent inputs
                    idx_uc = find(A_uu(:, idx_u));
                    for k=1:length(idx_uc)
                        b = A_uu(idx_uc(k), idx_u);
                        A_uu(idx_uc(k), idx_u) = 0;
                        A_uu(idx_uc(k), :) = A_uu(idx_uc(k), :) + b*A_uu(idx_u, :);
                    end

                    % Update constraints
                    b = obj.C_u(i, idx_u);
                    c_uu = -b*obj.C_u(i, :);
                    c_uu(idx_u) = 0;
                    idx_uc = find(obj.C_u(i+1:end, idx_u)) + i;
                    for k=1:length(idx_uc)
                        b = obj.C_u(idx_uc(k), idx_u);
                        obj.C_u(idx_uc(k), idx_u) = 0;
                        obj.C_u(idx_uc(k), :) = obj.C_u(idx_uc(k), :) + b * c_uu;
                    end
                else
                    warning(['ignoring constraint ' num2str(i)]);
                end
            end
            A_uu(:, ~any(A_uu)) = []; 
            u_dep = find(~u_f);
             
            obj.G = obj.G * A_uu; 

            % remove dependent inputs
            u_dep = sort(u_dep, 'descend');
            ports = obj.getPortsForIOPairs(u_dep);
            obj.deletePorts(ports);
            
            % All I/O coupling constraints are now eliminated
            obj.C_u = zeros(0, obj.n_u);
            obj.C_y = zeros(0, obj.n_u);
            
            % Some algebraic constraints might remain
            % Call eliminateAlgebraicConstraints() to eliminate them
            % This also eliminates constraint forces (inputs)
            obj.n_c = size(obj.B, 2);
        end
        
        function names = getInputNames(obj, IOPairs)
            names = cell(1, length(IOPairs));
            ports = obj.getPortsForIOPairs(IOPairs);
            for p = 1:length(ports)
                names{p} = obj.ports{ports(p)}.inputName;
            end
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
        
        function node = getNodeAtLocation(obj, location)
            node = 0;
            for n=1:obj.n_nodes
                if isa(obj.nodes{n}, 'PH_MechanicalNode') && ~any(obj.nodes{n}.location ~= location)
                    node = n;
                    break;
                end
            end
        end
        
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
        
        function addExternalInputsAtNodes(obj)
            % Original number of inputs
            n_u0 = obj.n_u;
            for n=1:obj.n_nodes
                node = obj.nodes{n};
                if ~isa(node, 'PH_MechanicalNode') %|| node.internal
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
            
            obj.G = [obj.G, casadi.MX.zeros(obj.n, obj.n_u - n_u0)];
            obj.C_u = [obj.C_u, zeros(size(obj.C_u, 1), obj.n_u - n_u0)];
            obj.C_y = [obj.C_y, zeros(size(obj.C_y, 1), obj.n_u - n_u0)];
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
        
        % Explicit version of the state derivative and its jacobian
        function [f, dfdx] = getSystemODEfun(obj)
            u = casadi.MX.sym('u', obj.n_u, 1);
            if obj.isConstrained()
                lambda = casadi.Function('lambda', {obj.x, u}, {-(obj.B'*hessian(obj.H, obj.x)*obj.B)^-1* ... 
                                                obj.B'*hessian(obj.H, obj.x)*((obj.J-obj.R)*jacobian(obj.H, obj.x)' + obj.G*u)});
                f_sym = (obj.J-obj.R)*jacobian(obj.H, obj.x)' + obj.G*u + obj.B*lambda(obj.x, u);
            else
                f_sym = (obj.J-obj.R)*jacobian(obj.H, obj.x)' + obj.G*u;
            end
            f = casadi.Function('f', {obj.x, u}, {f_sym});
            if nargout > 1
                dfdx_sym = jacobian(f_sym, obj.x);
                dfdx = casadi.Function('dfdx', {obj.x, u}, {dfdx_sym});
            end
        end
        
        % Implicit version of the system dynamics
        function [f, dfdt, x0, dx0] = getSystemDAEfun(obj, x0, u0)
            u = casadi.MX.sym('u', obj.n_u, 1);
            if obj.isConstrained()
                lambda = casadi.MX.sym('lambda', size(obj.B, 2), 1);
                f_sym = [(obj.J-obj.R)*jacobian(obj.H, obj.x)' + obj.G*u + obj.B*lambda; obj.B'*jacobian(obj.H, obj.x)'];
                f = casadi.Function('f', {[obj.x; lambda], u}, {f_sym});
                dfdt = blkdiag(eye(obj.n), zeros(size(obj.B,2)));
                % Consistent initial condition
                lambda0 = casadi.Function('lambda', {obj.x, u}, {-(obj.B'*hessian(obj.H, obj.x)*obj.B)^-1* ... 
                                                obj.B'*hessian(obj.H, obj.x)*((obj.J-obj.R)*jacobian(obj.H, obj.x)' + obj.G*u)});
                x0 = [x0; full(lambda0(x0, u0))];
                dx0 = full(f(x0, u0));
            else
                f_sym = (obj.J-obj.R)*jacobian(obj.H, obj.x)' + obj.G*u;
                f = casadi.Function('f', {obj.x, u}, {f_sym});
                dfdt = eye(obj.n);
                dx0 = f(x0, u0);
            end
        end
        
        function y = getSystemOutputFun(obj)
            if obj.isConstrained()
                y = casadi.Function('y', {obj.x}, {(obj.G' - ((obj.B'*hessian(obj.H, obj.x)*obj.B)^-1* ...
                                                    obj.B'*hessian(obj.H, obj.x)*obj.G)'*obj.B')*jacobian(obj.H, obj.x)'});
            else
                y = casadi.Function('y', {obj.x}, {obj.G'*jacobian(obj.H, obj.x)'});
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
            cp = PH_NonlinearSystem(obj.n, obj.x, obj.J, obj.H, obj.G, obj.R, obj.B);
            % copyData can also be called by subclasses
            obj.copyData(cp);
        end
        
        % Generates constraints in the form v_a = v_b and w_a = w_b at each
        % node, where v refers to a velocity in x/y/z direction and w to an 
        % angular velocity in x/y/z direction.
        function [Cu, Cy] = generateVelocityConstraints(obj)
            nConstraints = 0;
            V = zeros(nConstraints, obj.n_nodes*6);
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
                    if strcmp(port.type, 'force') || strcmp(port.type, 'torque')
                        % input is force/torque -> output is velocity
                        Cy(nConstraints, port.IOPair) = 1;
                    else
                        % input is velocity
                        Cu(nConstraints, port.IOPair) = 1; 
                    end
                end
            end  
            
            if size(Cu, 1) < nConstraints && obj.n_u > 0
                Cu(nConstraints, 1) = 0;
            end
            if size(Cy, 1) < nConstraints && obj.n_u > 0 
                Cy(nConstraints, 1) = 0;
            end
            
            % Remove constraints that are zero
            idx = ~any([Cu Cy]');
            Cu(idx,:) = [];
            Cy(idx,:) = [];
            V(idx, :) = [];
            V(:, ~any(V)) = [];
            G_v = null(V')';
            C = G_v*[Cu, Cy];
            Cu = C(:, 1:obj.n_u);
            Cy = C(:, obj.n_u+1:2*obj.n_u);
        end
        
        % Generates constraints of the form sum(M) = 0 / sum(F) = 0 at each
        % node
        function [Cu, Cy] = generateForceConstraints(obj)
            nConstraints = 0; 
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
            
            if size(Cu, 1) < nConstraints && obj.n_u > 0
                Cu(nConstraints, 1) = 0;
            end
            if size(Cy, 1) < nConstraints && obj.n_u > 0 
                Cy(nConstraints, 1) = 0;
            end
            % Remove constraints that are zero
            idx = ~any([Cu Cy]');
            Cu(idx,:) = [];
            Cy(idx,:) = [];
        end
    end 
    
    methods(Static, Access = public)
        function obj = toNonlinearPHSystem(system)
            if ~isa(system, 'PH_LinearSystem')
                error('Conversion only happening for PH_LinearSystemObjects...');
            end
            if nnz(system.P) || nnz(system.M)
                error('Output feedthrough is not yet supported for PH_NonlinearSystem');
            end
            
            x = casadi.MX.sym('x', system.n, 1); 
            J = casadi.MX.zeros(system.n, system.n);
            J(:,:) = system.J;
            R = casadi.MX.zeros(system.n, system.n);
            R(:,:) = system.R;
            H = x'*system.Q*x;
            G = casadi.MX.zeros(system.n, system.n_u);
            G(:,:) = system.G;
            B = casadi.MX.zeros(system.n, size(system.B,2));
            B(:,:) = system.B;
            
            obj = PH_NonlinearSystem(system.n, x, J, H, G, R, B);
            obj.C_u = system.C_u;
            obj.C_y = system.C_y; 
            obj.n_c = obj.n_c + size(system.C_u, 1); 
            
            obj.elements{1}.delete();
            obj.elements(1) = [];
            for e=1:system.n_elements
                obj.elements{e} = system.elements{e};
            end
            obj.n_elements = system.n_elements;
            
            for n=1:system.n_nodes
                obj.nodes{n} = system.nodes{n};
            end
            obj.n_nodes = system.n_nodes;
            
            for p=1:system.n_ports
                obj.ports{p} = system.ports{p};
            end
            obj.n_ports = system.n_ports;
            
            obj.name = system.name;
        end
    end
end