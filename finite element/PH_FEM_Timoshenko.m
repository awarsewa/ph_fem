classdef PH_FEM_Timoshenko < PH_FEM_element
    properties(SetAccess = protected, GetAccess = public)
        N
        density
        youngsModulus
        shearModulus
        areaMomentOfInertia
        crossSectionalArea
        timoshenkoShearCoefficient
        length
    end
    
    methods(Access = public)
        function obj = PH_FEM_Timoshenko(N, rho, E, G, A, I, kappa, L)
            % pH_TimoshenkoBeam_PFEM_p(N_p, N_q, rho, E, G, A, I, kappa, L) 
            [n_, J_, Q_, G_] = pH_TimoshenkoBeam_PFEM_p(N, N, rho, E, G, A, I, kappa, L);
            % LinearPHSystem(n, E, J, Q, G, R, K, P, S, M, C_u, C_y) 
            obj = obj@PH_FEM_element(n_, [], J_, Q_, G_);
            
            obj.name = 'Timoshenko beam element';
            
            % Physical parameters
            obj.N = N;
            obj.youngsModulus = E;
            obj.shearModulus = G;
            obj.areaMomentOfInertia = I;
            obj.crossSectionalArea = A;
            obj.length = L;
            obj.density = rho;
            obj.timoshenkoShearCoefficient = kappa;
            
            % Nodes
            obj.n_nodes = 2; 
            obj.nodes{1} = PH_MechanicalNode([1, 3], [0 0 0]);
            obj.nodes{2} = PH_MechanicalNode([2, 4], [L 0 0]);
            
            obj.elements{1}.delete();
            obj.elements{1} = PH_Element('PH_FEM_Timoshenko', [1 2 3 4], ...
                                        'Young''s modulus', E, ...
                                        'Shear modulus', G, ...
                                        'Area moment of inertia', I, ...
                                        'Cross sectional area', A, ...
                                        'Length', L, ...
                                        'Density', rho, ...
                                        'Timoshenko shear coefficient', kappa); 
            
            % Ports
            obj.n_ports = 4;

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
    
    methods(Access = protected)  
        function copyData(obj, cp)
            copyData@PH_LinearSystem(obj, cp);
        end
         
        function cp = copyElement(obj)
            % PH_FEM_Timoshenko(N, rho, E, G, A, I, kappa, L)
            cp = PH_FEM_Timoshenko(obj.N, obj.density, obj.youngsModulus, ...
                                  obj.shearModulus, obj.crossSectionalArea, ...
                                  obj.areaMomentOfInertia, obj.timoshenkoShearCoefficient,...
                                  obj.length);
            obj.copyData(cp);
        end
    end
end