% Isoparametric disk element
classdef PH_FEM_Disk < PH_FEM_element
    properties(Access = protected)
        order
        corner_nodes
        youngsModulus
        poissonRatio
        density
        thickness
    end
    
    methods(Access = public)
        function obj = PH_FEM_Disk(order, corner_nodes, rho, E, nu, h)   
            
            [n_, J_, Q_, G_, nodes_xy] = pH_Disk_PFEM_p(order, corner_nodes, rho, E, nu, h);

            N = (order+1)^2;
            % PH_LinearSystem(n, J, Q, G, R, P, S, M, B, x_p, x_q)
            x_p_ = [1:2*N]';
            x_q_ = [2*N+1:5*N]';
            obj = obj@PH_FEM_element(n_, J_, Q_, G_, [], [], [], [], [], x_p_, x_q_);
            
            obj.name = 'PH-FEM disk element';

            % Physical parameters
            obj.order = order;
            obj.corner_nodes = corner_nodes;
            obj.youngsModulus = E;
            obj.poissonRatio = nu;
            obj.density = rho;
            obj.thickness = h;
            
            % Nodes
            obj.n_nodes = N; 
            
            for n=1:obj.n_nodes
                %PH_MechanicalNode(ports, elements, location, lockedDOFs)
                obj.nodes{n} = PH_MechanicalNode([n obj.n_nodes+n], 1, [nodes_xy(n,:) 0]);
            end
            
            obj.elements{1}.delete();
            obj.elements{1} = PH_Element('PH_FEM_Disk', 1:obj.n_nodes, 1:obj.n_nodes, ...
                                        'order', order, ...
                                        'corner nodes', corner_nodes, ...
                                        'Young''s modulus', E, ...
                                        'Poisson ratio', nu, ...
                                        'density', rho, ...
                                        'thickness', h); 
            % Ports
            obj.n_ports = 2*obj.n_nodes;
            
            for n = 1:obj.n_nodes
                %PH_MechanicalPort_boundary(type, nodes, orientation, element, IOPair, inputName, outputName)
                obj.ports{n} = PH_MechanicalPort_boundary('force', n, [1; 0; 0], 1, n, ['F' num2str(n) '_x'], ['v' num2str(n) '_x']);
                obj.ports{obj.n_nodes+n} = PH_MechanicalPort_boundary('force', n, [0; 1; 0], 1, obj.n_nodes + n, ['F' num2str(n) '_y'], ['v' num2str(n) '_y']);
            end
        end
        
        function setPosition(obj, nodePositions)
            if ~ismatrix(nodePositions) || ...
                size(nodePositions, 1) ~= 4 || ...
                size(nodePositions, 2) ~= 3 
                error('NodePositions must be a 4 x 3 matrix');
            end
            
            for n = 1:4
                L0 = norm(obj.nodes{n}.location - obj.nodes{max(mod(n+1, 4),1)}.location);
                L = norm(nodePositions(n, :) - nodePositions(max(mod(n+1, 4),1), :));
                if L - L0 > 1e-12
                    error('When setting the position of a plate, do not change its geometry!');
                end
            end
            
            surfNorm0 = cross(obj.nodes{2}.location-obj.nodes{1}.location, ...
                              obj.nodes{4}.location-obj.nodes{1}.location);
            surfNorm0 = surfNorm0./norm(surfNorm0);
            
            surfNorm = cross(nodePositions(2, :) - nodePositions(1,:), ...
                             nodePositions(4, :) - nodePositions(1,:));
            surfNorm = surfNorm./norm(surfNorm);     
            Tr = round(vrrotvec2mat(vrrotvec(surfNorm0,surfNorm)), 12);
            
            for p=1:obj.n_ports
                obj.ports{p}.orientation = Tr*obj.ports{p}.orientation;
                obj.ports{p}.orientation(abs(obj.ports{p}.orientation) < 1e-9) = 0;
                obj.ports{p}.orientation = obj.ports{p}.orientation/norm(obj.ports{p}.orientation);
            end   
            
            nodes = nodePositions;
            if obj.order > 1
                for i=1:obj.order-1
                    nodes(4*i+1, :) = nodes(1, :) + i/obj.order*(nodes(2, :) - nodes(1, :));
                    nodes(4*i+2, :) = nodes(2, :) + i/obj.order*(nodes(3, :) - nodes(2, :));
                    nodes(4*i+3, :) = nodes(3, :) - i/obj.order*(nodes(3, :) - nodes(4, :));
                    nodes(4*i+4, :) = nodes(4, :) - i/obj.order*(nodes(4, :) - nodes(1, :));
                end
                % number of inner rectangles
                N_in = floor((obj.order+1)/2-1);
                for i = 1:N_in
                    inner(1, :) = nodes(1, :) + i/obj.order*(nodes(2, :) - nodes(1, :)) + i/obj.order*(nodes(4, :) - nodes(1, :));
                    inner(2, :) = nodes(2, :) - i/obj.order*(nodes(2, :) - nodes(1, :)) + i/obj.order*(nodes(3, :) - nodes(2, :));
                    inner(3, :) = nodes(3, :) - i/obj.order*(nodes(3, :) - nodes(4, :)) - i/obj.order*(nodes(3, :) - nodes(2, :));
                    inner(4, :) = nodes(4, :) - i/obj.order*(nodes(4, :) - nodes(1, :)) + i/obj.order*(nodes(3, :) - nodes(4, :));
                    for j=1:obj.order-2*i-1
                        inner(4*j+1, :) = inner(1, :) + j/(obj.order-2*i)*(inner(2, :) - inner(1, :));
                        inner(4*j+2, :) = inner(2, :) + j/(obj.order-2*i)*(inner(3, :) - inner(2, :));
                        inner(4*j+3, :) = inner(3, :) - j/(obj.order-2*i)*(inner(3, :) - inner(4, :));
                        inner(4*j+4, :) = inner(4, :) - j/(obj.order-2*i)*(inner(4, :) - inner(1, :));
                    end
                    nodes = [nodes; inner]; 
                end
                % single node at the center
                if mod(obj.order+1, 2)
                    nodes(end+1, :) = nodes(1,:) + 0.5*(nodes(4,:) + 0.5*(nodes(3,:) - nodes(4,:)) - nodes(1,:) + 0.5*(nodes(2,:) - nodes(1,:)));
                end
            end
            
            for n = 1:obj.n_nodes
                obj.nodes{n}.location = nodes(n, :);
            end
        end
    end
    
    methods(Access = protected)
        function copyData(obj, cp)
            copyData@PH_LinearSystem(obj, cp);
        end
         
        function cp = copyElement(obj)
            % PH_FEM_Link(N, myu, E, A, L)
            cp = PH_FEM_Disk(obj.order, obj.corner_nodes, obj.density, obj.youngsModulus, obj.poissonRatio, ...
                                    obj.thickness);
            obj.copyData(cp);
        end
    end
end