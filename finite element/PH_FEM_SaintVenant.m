classdef PH_FEM_SaintVenant < PH_FEM_element
    properties(Access = protected)
        N
        polarMomentOfInertia
        shearModulus
        torsionConstant
        length
    end
    
    methods(Access = public)
        function obj = PH_FEM_SaintVenant(N, I_p, G_t, I_t, L)
            % Check whether valid boundary conditions were selected
            [n_, J_, Q_, G_] = pH_SaintVenant_PFEM_p(N, N, I_p, G_t, I_t, L);
            P_ =  [];
            S_ =  [];
            M_ =  [];
            % LinearPHSystem(n, E, J, Q, G, R, K, P, S, M, C_u, C_y) 
            obj = obj@PH_FEM_element(n_, [], J_, Q_, G_, [], [], P_, S_, M_);
            
            obj.name = 'Saint-Venant torsion element';
            
            % Physical parameters
            obj.N = N;
            obj.shearModulus = G_t;
            obj.torsionConstant = I_t;
            obj.length = L;
            obj.polarMomentOfInertia = I_p;
            
            % Nodes
            obj.n_nodes = 2; 
            obj.nodes{1} = PH_MechanicalNode(1, [0 0 0]);
            obj.nodes{2} = PH_MechanicalNode(2, [L 0 0]);
            
            obj.elements{1}.delete();
            obj.elements{1} = PH_Element('PH_FEM_SaintVenant', [1 2], ...
                                        'Shear modulus', G_t, ...
                                        'Torsion constant', I_t, ...
                                        'Length', L, ...
                                        'Polar moment of inertia', I_p); 
            
            % Ports
            obj.n_ports = 2;
            % PH_MechanicalPort(scope, type, IOPair, inputType, outputType, node, orientation)
            obj.ports{1} = PH_MechanicalPort('internal', 'torque/angular velocity', 1, 'torque', 'angular velocity', 1, [1; 0; 0]);
            obj.ports{2} = PH_MechanicalPort('internal', 'torque/angular velocity', 2, 'torque', 'angular velocity', 2, [1; 0; 0]);

            % I/O names
            obj.inputNames{1} = 'Ma_x';
            obj.inputNames{2} = 'Mb_x';
            obj.outputNames{1} = 'wa_x';
            obj.outputNames{2} = 'wb_x';
            obj.inputSigns =  [1; 1];
            obj.outputSigns = [1; 1];    
        end
    end
    
    methods(Access = protected)  
        function copyData(obj, cp)
            copyData@PH_LinearSystem(obj, cp);
        end
         
        function cp = copyElement(obj)
            % PH_FEM_SaintVenant(N, I_p, G_t, I_t, L)
            cp = PH_FEM_SaintVenant(obj.N, obj.polarMomentOfInertia, obj.shearModulus, ...
                                  obj.torsionConstant, obj.length);
            obj.copyData(cp);
        end
    end
end