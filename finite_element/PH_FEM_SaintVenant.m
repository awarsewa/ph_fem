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
            % Get system matrices
            [n_, J_, Q_, G_] = pH_SaintVenant_PFEM_p(N, N, I_p, G_t, I_t, L);
            % LinearPHSystem(n, J, Q, G, R, P, S, M, C_u, C_y) 
            obj = obj@PH_FEM_element(n_, J_, Q_, G_);
            
            obj.name = 'Saint-Venant torsion element';
            
            % Physical parameters
            obj.N = N;
            obj.shearModulus = G_t;
            obj.torsionConstant = I_t;
            obj.length = L;
            obj.polarMomentOfInertia = I_p;
            
            % Nodes
            obj.n_nodes = size(obj.G, 2);   
            node_pos = linspace(0, L, obj.n_nodes);
            for n = 1:obj.n_nodes
                obj.nodes{n} = PH_MechanicalNode(n, 1, [node_pos(n) 0 0]);
                obj.nodes{n}.internal = 1;
            end
            obj.nodes{1}.internal = 0;
            obj.nodes{end}.internal = 0;
            %{
            obj.n_nodes = 2; 
            obj.nodes{1} = PH_MechanicalNode(1, 1, [0 0 0]);
            obj.nodes{2} = PH_MechanicalNode(2, 1, [L 0 0]);
            %}
            obj.elements{1}.delete();
            obj.elements{1} = PH_Element('PH_FEM_SaintVenant', 1:obj.n_nodes, 1:obj.n_nodes, ...
                                        'Shear modulus', G_t, ...
                                        'Torsion constant', I_t, ...
                                        'Length', L, ...
                                        'Polar moment of inertia', I_p); 
            
            % Ports
            obj.n_ports = size(obj.G, 2);
            for n = 1:obj.n_nodes
                obj.ports{n} = PH_MechanicalPort_external('torque', n, [1; 0; 0], n, ['M' num2str(n) '_x'], ['w' num2str(n) '_x']);
            end
            %{
            obj.n_ports = 2;
            % PH_MechanicalPort_external(type, nodes, orientation, IOPair, inputName, outputName)     
            obj.ports{1} = PH_MechanicalPort_external('torque', 1, [1; 0; 0], 1, 'Ma_x', 'wa_x');
            obj.ports{2} = PH_MechanicalPort_external('torque', 2, [1; 0; 0], 2, 'Mb_x', 'wb_x');
            %}
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