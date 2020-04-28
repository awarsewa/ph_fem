classdef PH_FEM_Bernoulli < PH_FEM_Beam_element
    properties(SetAccess = protected, GetAccess = public)
        N
        massPerUnitLength
        youngsModulus
        areaMomentOfInertia
        length
    end
    
    methods(Access = public)
        function obj = PH_FEM_Bernoulli(N, myu, E, I, L)
            % Get system matrices
            [n_, J_, Q_, G_] = pH_BernoulliBeam_PFEM_p(N, N, myu, E, I, L);
            % LinearPHSystem(n, J, Q, G, R, P, S, M, C_u, C_y) 
            obj = obj@PH_FEM_Beam_element(n_, J_, Q_, G_);
            
            obj.name = 'PH-FEM Euler-Bernoulli beam';
            
            % Physical parameters
            obj.N = N;
            obj.youngsModulus = E;
            obj.areaMomentOfInertia = I;
            obj.length = L;
            obj.massPerUnitLength = myu;
            
            % Nodes
            obj.n_nodes = size(obj.G, 2)/2;
                         
            node_pos = linspace(0, L, obj.n_nodes);
            for n = 1:obj.n_nodes
                obj.nodes{n} = PH_MechanicalNode([n, obj.n_nodes+n], 1, [node_pos(n) 0 0]);
                obj.nodes{n}.internal = 1;
            end
            obj.nodes{1}.internal = 0;
            obj.nodes{end}.internal = 0;
            %{
            obj.n_nodes = 2; 
            obj.nodes{1} = PH_MechanicalNode([1, 3], 1, [0 0 0]);
            obj.nodes{2} = PH_MechanicalNode([2, 4], 1, [L 0 0]);
            %}
            obj.elements{1}.delete();
                             %PH_Element(type, nodes, ports, varargin)
            obj.elements{1} = PH_Element('PH_FEM_Bernoulli', 1:obj.n_nodes, 1:obj.n_nodes*2, ...
                                        'Young''s modulus', E, ...
                                        'Area moment of inertia', I, ...
                                        'Length', L, ...
                                        'Mass per unit length', myu); 
            obj.n_ports = size(obj.G, 2); 
            for n = 1:obj.n_nodes
                obj.ports{n} = PH_MechanicalPort_external('torque', n, [0; 1; 0], n, ['M' num2str(n) '_y'], ['w' num2str(n) '_y']);
                obj.ports{obj.n_nodes+n} = PH_MechanicalPort_external('force', n, [0; 0; 1], obj.n_nodes+n, ['F' num2str(n) '_z'], ['v' num2str(n) '_z']);
            end
            %{
            % Ports
            obj.n_ports = 4;

            % PH_MechanicalPort_external(type, nodes, orientation, IOPair, inputName, outputName)     
            obj.ports{1} = PH_MechanicalPort_external('torque', 1, [0; 1; 0], 1, 'Ma_y', 'wa_y');
            obj.ports{2} = PH_MechanicalPort_external('torque', 2, [0; 1; 0], 2, 'Mb_y', 'wb_y');
            obj.ports{3} = PH_MechanicalPort_external('force', 1, [0; 0; 1], 3, 'Fa_z', 'va_z');
            obj.ports{4} = PH_MechanicalPort_external('force', 2, [0; 0; 1], 4, 'Fb_z', 'vb_z'); 
            %}
        end
    end
    
    methods(Access = protected)  
        function copyData(obj, cp)
            copyData@PH_LinearSystem(obj, cp);
        end
                
        function cp = copyElement(obj)
            % PH_FEM_Bernoulli(N, myu, E, I, L, boundary)
            cp = PH_FEM_Bernoulli(obj.N, obj.massPerUnitLength, obj.youngsModulus, ...
                                  obj.areaMomentOfInertia, obj.length);
            obj.copyData(cp);
        end
    end
end