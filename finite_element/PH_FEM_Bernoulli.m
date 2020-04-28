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
            end

            obj.elements{1}.delete();
                             %PH_Element(type, nodes, ports, varargin)
            obj.elements{1} = PH_Element('PH_FEM_Bernoulli', 1:obj.n_nodes, 1:obj.n_nodes*2, ...
                                        'Young''s modulus', E, ...
                                        'Area moment of inertia', I, ...
                                        'Length', L, ...
                                        'Mass per unit length', myu); 
            obj.n_ports = size(obj.G, 2); 
            for n = 1:obj.n_nodes
                %PH_MechanicalPort_boundary(type, nodes, orientation, element, IOPair, inputName, outputName)  
                obj.ports{n} = PH_MechanicalPort_boundary('torque', n, [0; 1; 0], 1, n, ['M' num2str(n) '_y'], ['w' num2str(n) '_y']);
                obj.ports{obj.n_nodes+n} = PH_MechanicalPort_boundary('force', n, [0; 0; 1], 1, obj.n_nodes+n, ['F' num2str(n) '_z'], ['v' num2str(n) '_z']);
            end
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