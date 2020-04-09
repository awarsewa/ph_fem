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
                [n_, J_, Q_, G_] = pH_Link_PFEM_p(N, N, myu, E, A, L);
            else
                [n_, J_, Q_, G_] = pH_Link_PFEM_q(N, N, myu, E, A, L);
            end

            % LinearPHSystem(n, E, J, Q, G, R, P, S, M, C_u, C_y) 
            obj = obj@PH_FEM_element(n_, [], J_, Q_, G_);
            
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
                obj.nodes{n}.internal = 1;
            end
            obj.nodes{1}.internal = 0;
            obj.nodes{end}.internal = 0;
            
            %obj.nodes{1} = PH_MechanicalNode(1, 1, [0 0 0]);
            %obj.nodes{2} = PH_MechanicalNode(2, 1, [L 0 0]);
            
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
                    obj.ports{n} = PH_MechanicalPort_external('force', n, [1; 0; 0], n, ['F' num2str(n) '_x'], ['v' num2str(n) '_x']);
                else
                    obj.ports{n} = PH_MechanicalPort_external('velocity', n, [1; 0; 0], n, ['v' num2str(n) '_x'], ['F' num2str(n) '_x']);
                end
            end
            %{
            if strcmp(type, 'force')
                % PH_MechanicalPort_external(type, nodes, orientation, IOPair, inputName, outputName)     
                obj.ports{1} = PH_MechanicalPort_external('force', 1, [1; 0; 0], 1, 'Fa_x', 'va_x');
                obj.ports{2} = PH_MechanicalPort_external('force', 2, [1; 0; 0], 2, 'Fb_x', 'vb_x');
            else
                obj.ports{1} = PH_MechanicalPort_external('velocity', 1, [1; 0; 0], 1, 'va_x', 'Fa_x');
                obj.ports{2} = PH_MechanicalPort_external('velocity', 2, [1; 0; 0], 2, 'vb_x', 'Fb_x');
            end
            %}
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