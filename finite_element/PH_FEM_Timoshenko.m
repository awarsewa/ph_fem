classdef PH_FEM_Timoshenko < PH_FEM_Beam_element
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
            % Get system matrices
            % pH_TimoshenkoBeam_PFEM_p(N_p, N_q, rho, E, G, A, I, kappa, L) 
            [n_, J_, Q_, G_] = pH_TimoshenkoBeam_PFEM_p(N, N, rho, E, G, A, I, kappa, L);
            % LinearPHSystem(n, J, Q, G, R, K, P, S, M, C_u, C_y) 
            obj = obj@PH_FEM_Beam_element(n_, J_, Q_, G_);
            
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
            obj.n_nodes = size(obj.G, 2)/2;
                         
            node_pos = linspace(0, L, obj.n_nodes);
            for n = 1:obj.n_nodes
                obj.nodes{n} = PH_MechanicalNode([n, obj.n_nodes+n], 1, [node_pos(n) 0 0]);
            end
            
            obj.elements{1}.delete();
            obj.elements{1} = PH_Element('PH_FEM_Timoshenko', 1:obj.n_nodes, 1:obj.n_nodes*2, ...
                                        'Young''s modulus', E, ...
                                        'Shear modulus', G, ...
                                        'Area moment of inertia', I, ...
                                        'Cross sectional area', A, ...
                                        'Length', L, ...
                                        'Density', rho, ...
                                        'Timoshenko shear coefficient', kappa); 
            
            obj.n_ports = size(obj.G, 2); 
            for n = 1:obj.n_nodes
                %PH_MechanicalPort_boundary(type, nodes, orientation, element, IOPair, inputName, outputName) 
                obj.ports{n} = PH_MechanicalPort_boundary('torque', n, [0; 1; 0], 1, n, ['M' num2str(n) '_y'], ['w' num2str(n) '_y']);
                obj.ports{obj.n_nodes+n} = PH_MechanicalPort_boundary('force', 1, n, [0; 0; 1], obj.n_nodes+n, ['F' num2str(n) '_z'], ['v' num2str(n) '_z']);
            end
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