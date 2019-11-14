classdef PH_FEM_Bernoulli < PH_FEM_element
    properties(SetAccess = protected, GetAccess = public)
        N
        massPerUnitLength
        youngsModulus
        areaMomentOfInertia
        length
        boundaryType
    end
    
    methods(Access = public)
        function obj = PH_FEM_Bernoulli(N, myu, E, I, L, boundary)
            % Check whether valid boundary conditions were selected
            if ~ischar(boundary) || (~strcmp(boundary, 'clamped-free') && ~strcmp(boundary, 'free-free'))
                error('Given boundary conditions are not valid. Please enter either ''clamped-free'' or ''free-free''.');
            end
            
            if strcmp(boundary, 'clamped-free')
                [n_, J_, Q_, G_, P_, S_, M_] = pH_BernoulliBeam_PSM(N, myu, E, I, L);
            else
                [n_, J_, Q_, G_] = pH_BernoulliBeam_PFEM_p(N, N, myu, E, I, L);
                P_ =  [];
                S_ =  [];
                M_ =  [];
            end
            % LinearPHSystem(n, E, J, Q, G, R, K, P, S, M, C_u, C_y) 
            obj = obj@PH_FEM_element(n_, [], J_, Q_, G_, [], [], P_, S_, M_);
            
            obj.name = [boundary ' Euler-Bernoulli beam element'];
            obj.boundaryType = boundary;
            
            % Physical parameters
            obj.N = N;
            obj.youngsModulus = E;
            obj.areaMomentOfInertia = I;
            obj.length = L;
            obj.massPerUnitLength = myu;
            
            % Nodes
            obj.n_nodes = 2; 
            obj.nodes{1} = PH_MechanicalNode([1, 3], [0 0 0]);
            obj.nodes{2} = PH_MechanicalNode([2, 4], [L 0 0]);
            
            obj.elements{1}.delete();
            obj.elements{1} = PH_Element('PH_FEM_Bernoulli', [1 2 3 4], ...
                                        'Young''s modulus', E, ...
                                        'Area moment of inertia', I, ...
                                        'Length', L, ...
                                        'Mass per unit length', myu); 
            
            % Ports
            obj.n_ports = 4;
            if strcmp(boundary, 'clamped-free')          
                % PH_MechanicalPort(scope, type, IOPair, inputType, outputType, node, orientation)
                obj.ports{1}  = PH_MechanicalPort('internal', 'force/velocity', 2, 'velocity', 'force', 1, [0; 0; 1]);
                obj.ports{2} = PH_MechanicalPort('internal', 'force/velocity', 3, 'force', 'velocity', 2, [0; 0; 1]);
                obj.ports{3} = PH_MechanicalPort('internal', 'torque/angular velocity', 1, 'angular velocity', 'torque', 1, [0; 1; 0]);
                obj.ports{4} = PH_MechanicalPort('internal', 'torque/angular velocity', 4, 'torque', 'angular velocity', 2, [0; 1; 0]);
                
                % I/O names
                obj.inputNames{1} = 'wa_z';
                obj.inputNames{2} = 'va_y';
                obj.inputNames{3} = 'Fb_z';
                obj.inputNames{4} = 'Mb_y';
                obj.outputNames{1} = 'Ma_z';
                obj.outputNames{2} = 'Fa_y';
                obj.outputNames{3} = 'vb_y';
                obj.outputNames{4} = 'wb_z';
                obj.inputSigns =  [1; 1; 1; 1];
                obj.outputSigns = [1; 1; 1; 1];
                
                obj.nodes{1}.lockedDOFs = [0 1 0 0 0 1];
            else
                % PH_MechanicalPort(scope, type, IOPair, inputType, outputType, node, orientation)
                obj.ports{1} = PH_MechanicalPort('internal', 'torque/angular velocity', 1, 'torque', 'angular velocity', 1, [0; 1; 0]);
                obj.ports{2} = PH_MechanicalPort('internal', 'torque/angular velocity', 2, 'torque', 'angular velocity', 2, [0; 1; 0]);
                obj.ports{3} = PH_MechanicalPort('internal', 'force/velocity', 3, 'force', 'velocity', 1, [0; 0; 1]);
                obj.ports{4} = PH_MechanicalPort('internal', 'force/velocity', 4, 'force', 'velocity', 2, [0; 0; 1]);
                
                % I/O names
                obj.inputNames{1} = 'Ma_y';
                obj.inputNames{2} = 'Mb_y';
                obj.inputNames{3} = 'Fa_z';
                obj.inputNames{4} = 'Fb_z';
                obj.outputNames{1} = 'wa_y';
                obj.outputNames{2} = 'wb_y';
                obj.outputNames{3} = 'va_z';
                obj.outputNames{4} = 'vb_z';
                obj.inputSigns =  [1; 1; 1; 1];
                obj.outputSigns = [1; 1; 1; 1];
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
                                  obj.areaMomentOfInertia, obj.length, obj.boundaryType);
            obj.copyData(cp);
        end
    end
end