classdef PH_HydraulicActuator < PH_NonlinearSystem
    properties(SetAccess = protected, GetAccess = public)
        beta    % fluid bulk modulus
        m_p     % piston mass
        m_c     % cylinder mass
        b_p     % linear damping coefficient piston
        b_c     % linear damping coefficient cylinder
        L       % cylinder length
        A1      % piston-side cross-sectional area
        A2      % full cross-sectional area
        kv      % valve coefficient
        p_s     % supply pressure
        p_t     % tank pressure
        C_int   % internal leakage flow coefficient
        C_ext   % external leakage flow coefficient
        action  % 'push' or 'pull' ...  returns different input matrices
    end
    
    methods(Access = public)
        %              PH_NonlinearSystem(n, J, H, G, R)
        function obj = PH_HydraulicActuator(beta, m_p, m_c, b_p, b_c, L, A1, A2, kv, p_s, p_t, C_int, C_ext, action)
            x_ = casadi.MX.sym('x', 6, 1);
            J_ = [0     0   1                       0                       0                   0; ...
                  0     0   0                       1                       0                   0; ...
                  -1    0   0                       0                       beta/(x_(1)-x_(2)) -beta/(L-x_(1)+x_(2)); ...
                  0     -1  0                       0                       -beta/(x_(1)-x_(2)) beta/(L-x_(1)+x_(2)); ...
                  0     0   -beta/(x_(1)-x_(2))     beta/(x_(1)-x_(2))      0                   0; ... 
                  0     0   beta/(L-x_(1)+x_(2))    beta/(L-x_(1)+x_(2))   0                   0];
            R_ = [casadi.MX.zeros(2, 6); ...
                  0 0 b_p 0 0 0; ...
                  0 0 0 b_c 0 0; ... 
                  0 0 0 0 (C_ext+C_int)*beta/(A1*(x_(1)-x_(2)))-C_int*x_(6) 0; ...
                  0 0 0 0 0 (C_ext+C_int)*(beta/(A2*(L-x_(1)+x_(2))))-C_int*x_(5)];
            
            H_ = 0.5/m_p*x_(3)^2 + 0.5/m_c*x_(4)^2 + ...                    % kinetic energy
                 A1*(x_(1)-x_(2))*(beta*(exp(x_(5)/beta)-1)-x_(5)) + ...    % chamber 1 energy
                 A2*(L - x_(1) + x_(2))*(beta*(exp(x_(6)/beta)-1)-x_(6));   % chamber 2 energy
            
            if strcmp(action, 'pull')
                G_ = [0 0 0; ...
                      0 0 0; ... 
                      -1 0 0; ...
                      0 -1 0; ...
                      0 0 beta/(A1*x_(1)-x_(2))*kv*sqrt(x_(5)-p_t); ...
                      0 0 -beta/(A2*(L-x_(1)+x_(2)))*kv*sqrt(p_s-x_(6))];
            else
                G_ = [0 0 0; ...
                      0 0 0; ... 
                      -1 0 0; ...
                      0 -1 0; ...
                      0 0 beta/(A1*x_(1)-x_(2))*kv*sqrt(p_s-x_(5)); ...
                      0 0 -beta/(A2*(L-x_(1)+x_(2)))*kv*sqrt(x_(6)-p_t)];
            end
            % PH_NonlinearSystem(n, x, J, H, u, G, R)
            obj = obj@PH_NonlinearSystem(6, x_, J_, H_, G_, R_);

            obj.n_ports = 2;
            obj.name = 'double acting hydraulic piston actuator';
            
            % Physical parameters
            obj.beta = beta;
            obj.m_p = m_p;
            obj.m_c = m_c;
            obj.b_p = b_p;
            obj.b_c = b_c; 
            obj.L = L;
            obj.A1 = A1;
            obj.A2 = A2;
            obj.kv = kv;
            obj.p_s = p_s;
            obj.p_t = p_t; 
            obj.C_int = C_int;
            obj.C_ext = C_ext;
            obj.action = action;
            
            % Nodes
            obj.n_nodes = 2;   
            % PH_MechanicalNode(ports, elements, location, lockedDOFs)
            obj.nodes{1} = PH_MechanicalNode(2, 1, [0 0 0]);
            obj.nodes{2} = PH_MechanicalNode(1, 1, [L 0 0]);
            
            obj.elements{1}.delete();
                             %PH_Element(type, nodes, ports, varargin)
            obj.elements{1} = PH_Element('PH_HydraulicActuator', 1:obj.n_nodes, 1:obj.n_nodes); 
            
            % Ports
            obj.n_ports = 2;

            %PH_MechanicalPort_external(type, nodes, orientation, IOPair, inputName, outputName) 
            obj.ports{1} = PH_MechanicalPort_external('force', 2, [1; 0; 0], 1, 'Fp_x', 'vp_x');
            obj.ports{2} = PH_MechanicalPort_external('force', 1, [1; 0; 0], 2, 'Fc_x', 'vc_x');
        end
        
        function setPosition(obj, nodePositions)
            if ~ismatrix(nodePositions) || ...
                size(nodePositions, 1) ~= 2 || ...
                size(nodePositions, 2) ~= 3 
                error('NodePositions must be a 2x3 matrix');
            end
            dir = obj.nodes{2}.location - obj.nodes{1}.location;
            dir = dir./norm(dir);
            nDir = nodePositions(2,:) - nodePositions(1,:);
            nDir = nDir./norm(nDir);
            
            Tr = round(vrrotvec2mat(vrrotvec(dir,nDir)), 12);
            
            for p=1:obj.n_ports
                obj.ports{p}.orientation = Tr*obj.ports{p}.orientation;
                obj.ports{p}.orientation(abs(obj.ports{p}.orientation) < 1e-9) = 0;
                obj.ports{p}.orientation = obj.ports{p}.orientation/norm(obj.ports{p}.orientation);
            end    

            for n=1:obj.n_nodes
                obj.nodes{n}.location = nodePositions(n,:);
            end
        end
    end
    
    methods(Access = protected)   
        function copyData(obj, cp)
            copyData@PH_NonlinearSystem(obj, cp);
        end
        function cp = copyElement(obj)
               % PH_HydraulicActuator(beta, m_p, m_c, b_p, b_c, L, A1, A2, kv, p_s, p_t, C_int, C_ext)
            cp = PH_HydraulicActuator(obj.beta, obj.m_p, obj.m_c, obj.b_p, obj.b_c, obj.L, obj.A1, obj.A2, ...
                                obj.kv, obj.p_s, obj.p_t, obj.C_int, obj.C_ext, obj.action);
            obj.copyData(cp);
        end
    end
end