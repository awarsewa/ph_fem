classdef PH_FEM_Link < PH_FEM_element
    properties(Access = protected)
        N
        massPerUnitLength
        youngsModulus
        crossSectionalArea
        length
        inputType
    end
    
    methods(Access = public)
        function obj = PH_FEM_Link(N, myu, E, A, L, inputType)   
            type = 'force'; 
            if nargin > 5
                if ~strcmp(inputType, 'velocity') && ~strcmp(inputType, 'force')
                    error('inputType must be either ''force'' (default) or ''velocity''');
                end
                type = inputType;
            end

            if strcmp(type, 'force')
                [n_, J_, Q_, G_, C_] = pH_Link_PFEM_p(N, N, myu, E, A, L);
            else
                [n_, J_, Q_, G_, C_] = pH_Link_PFEM_q(N, N, myu, E, A, L);
            end

            % PH_LinearSystem(n, J, Q, G, R, P, S, M, B, x_p, x_q, C)
            obj = obj@PH_FEM_element(n_, J_, Q_, G_, [], [], [], [], [], [], [], C_);
            
            obj.name = 'PH-FEM link element';
            obj.inputType = type;
            
            % Physical parameters
            obj.N = N;
            obj.youngsModulus = E;
            obj.crossSectionalArea = A;
            obj.length = L;
            obj.massPerUnitLength = myu;
            
            % Nodes
            obj.n_nodes = size(obj.G, 2); 
            
            node_pos = linspace(0, L, obj.n_nodes);
            for n = 1:obj.n_nodes
                obj.nodes{n} = PH_MechanicalNode(n, 1, [node_pos(n) 0 0]);
            end
            
            obj.elements{1}.delete();
            obj.elements{1} = PH_Element('PH_FEM_Link', 1:obj.n_nodes, 1:obj.n_nodes, ...
                                        'Young''s modulus', E, ...
                                        'Cross sectional area', A, ...
                                        'Length', L, ...
                                        'Mass per unit length', myu); 
            
            % Ports
            obj.n_ports = size(obj.G, 2);
            for n = 1:obj.n_nodes
                if strcmp(type, 'force')
                    %PH_MechanicalPort_boundary(type, nodes, orientation, element, IOPair, inputName, outputName) 
                    obj.ports{n} = PH_MechanicalPort_boundary('force', n, [1; 0; 0], 1, n, ['F' num2str(n) '_x'], ['v' num2str(n) '_x']);
                else
                    obj.ports{n} = PH_MechanicalPort_boundary('velocity', n, [1; 0; 0], 1, n, ['v' num2str(n) '_x'], ['F' num2str(n) '_x']);
                end
            end
        end
    end
    
    methods(Access = protected)
        function copyData(obj, cp)
            copyData@PH_LinearSystem(obj, cp);
        end
         
        function cp = copyElement(obj)
            % PH_FEM_Link(N, myu, E, A, L)
            cp = PH_FEM_Link(obj.N, obj.massPerUnitLength, obj.youngsModulus, ...
                                    obj.crossSectionalArea, obj.length, obj.inputType);
            obj.copyData(cp);
        end
    end
end