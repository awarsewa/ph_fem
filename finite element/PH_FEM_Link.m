classdef PH_FEM_Link < PH_FEM_element
    properties(Access = protected)
        N
        massPerUnitLength
        youngsModulus
        crossSectionalArea
        length
        boundaryType
        inputType
    end
    
    methods(Access = public)
        function obj = PH_FEM_Link(N, myu, E, A, L, boundary, inputType)
            % Check whether valid boundary conditions were selected
            if ~ischar(boundary) || (~strcmp(boundary, 'clamped-free') && ~strcmp(boundary, 'free-free'))
                error('Given boundary conditions are not valid. Please enter either ''clamped-free'' or ''free-free''.');
            end
            
            type = 'force'; 
            if nargin > 6
                if ~strcmp(inputType, 'velocity') && ~strcmp(inputType, 'force')
                    error('inputType must be either ''force'' (default) or ''velocity''');
                end
                type = inputType;
            end
            
            if strcmp(boundary, 'clamped-free')
                [n_, J_, Q_, G_, P_, S_, M_] = pH_Link_PSM(N, myu, E, A, L);
            else
                if strcmp(type, 'force')
                    [n_, J_, Q_, G_] = pH_Link_PFEM_p(N, N, myu, E, A, L);
                else
                    [n_, J_, Q_, G_] = pH_Link_PFEM_q(N, N, myu, E, A, L);
                end
                P_ =  [];
                S_ =  [];
                M_ =  [];
            end
            % LinearPHSystem(n, E, J, Q, G, R, K, P, S, M, C_u, C_y) 
            obj = obj@PH_FEM_element(n_, [], J_, Q_, G_, [], [], P_, S_, M_);
            
            obj.name = [boundary ' link element'];
            obj.boundaryType = boundary;
            obj.inputType = type;
            
            % Physical parameters
            obj.N = N;
            obj.youngsModulus = E;
            obj.crossSectionalArea = A;
            obj.length = L;
            obj.massPerUnitLength = myu;
            
            % Nodes
            obj.n_nodes = 2; 
            obj.nodes{1} = PH_MechanicalNode(1, [0 0 0]);
            obj.nodes{2} = PH_MechanicalNode(2, [L 0 0]);
            
            obj.elements{1}.delete();
            obj.elements{1} = PH_Element('PH_FEM_Link', [1 2], ...
                                        'Young''s modulus', E, ...
                                        'Cross sectional area', A, ...
                                        'Length', L, ...
                                        'Mass per unit length', myu); 
            
            % Ports
            obj.n_ports = 2;
            if strcmp(boundary, 'clamped-free')   
                % Create ports
                             % PH_MechanicalPort(scope, type, IOPair, inputType, outputType, node, orientation)
                obj.ports{1} = PH_MechanicalPort('internal', 'force/velocity', 1, 'velocity', 'force', 1, [1; 0; 0]);
                obj.ports{2} = PH_MechanicalPort('internal', 'force/velocity', 2, 'force', 'velocity', 2, [1; 0; 0]);
                
                % I/O names
                obj.inputNames{1} = 'va_x';
                obj.inputNames{2} = 'Fb_x';
                obj.outputNames{1} =  'Fa_x';
                obj.outputNames{2} = 'vb_x';
                obj.inputSigns =  [1; 1];
                obj.outputSigns = [1; 1];
                
                % First node's x-axis is locked
                obj.nodes{1}.lockedDOFs = [1 0 0 0 0 0];
            else
                if strcmp(type, 'force')
                    obj.ports{1} = PH_MechanicalPort('internal', 'force/velocity', 1, 'force', 'velocity', 1, [1; 0; 0]);
                    obj.ports{2} = PH_MechanicalPort('internal', 'force/velocity', 2, 'force', 'velocity', 2, [1; 0; 0]);

                    % I/O names
                    obj.inputNames{1} = 'Fa_x';
                    obj.inputNames{2} = 'Fb_x';
                    obj.outputNames{1} = 'va_x';
                    obj.outputNames{2} = 'vb_x';
                else
                    obj.ports{1} = PH_MechanicalPort('internal', 'force/velocity', 1, 'velocity', 'force', 1, [1; 0; 0]);
                    obj.ports{2} = PH_MechanicalPort('internal', 'force/velocity', 2, 'velocity', 'force', 2, [1; 0; 0]);

                    % I/O names
                    obj.inputNames{1} = 'va_x';
                    obj.inputNames{2} = 'vb_x';
                    obj.outputNames{1} = 'Fa_x';
                    obj.outputNames{2} = 'Fb_x';
                end
                
                obj.inputSigns =  [1; 1];
                obj.outputSigns = [1; 1];
            end           
            
        end
    end
    
    methods(Access = protected)
        function copyData(obj, cp)
            copyData@PH_LinearSystem(obj, cp);
        end
         
        function cp = copyElement(obj)
            % PH_FEM_Link(N, myu, E, A, L, boundary)
            cp = PH_FEM_Link(obj.N, obj.massPerUnitLength, obj.youngsModulus, ...
                                    obj.crossSectionalArea, obj.length, obj.boundaryType, obj.inputType);
            obj.copyData(cp);
        end
    end
end