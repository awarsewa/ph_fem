classdef PH_FEM_element < PH_LinearSystem 
    methods(Access = public)
        function obj = PH_FEM_element(n, E, J, Q, varargin)
            obj = obj@PH_LinearSystem(n, E, J, Q, varargin{:});
        end
        
        function setPosition(obj, nodePositions)
            if obj.n_nodes > 2
                error('Setting the position of objects with more than 2 nodes not yet supported');
            end
            if ~ismatrix(nodePositions) || ...
                size(nodePositions, 1) ~= obj.n_nodes || ...
                size(nodePositions, 2) ~= 3 
                error('NodePositions must be a %dx3 matrix', length(obj.n_nodes));
            end
                       
            % Change orientation
            % For a beam we need to calculate a rotation matrix
            nDir = nodePositions(2,:) - nodePositions(1,:);
            nDir = nDir./norm(nDir);
            
            Tr = round(vrrotvec2mat(vrrotvec([1 0 0],nDir)), 12);
            for p=1:obj.n_ports
                obj.ports{p}.orientation = Tr*obj.ports{p}.orientation;
            end    
            
            for i=1:obj.n_nodes
                obj.nodes{i}.location = nodePositions(i, :);
                obj.nodes{1}.lockedDOFs = [(Tr*obj.nodes{1}.lockedDOFs(1:3)')' ...
                                           (Tr*obj.nodes{1}.lockedDOFs(4:6)')'];
            end
        end
        
        % Modified add function 
        function add(obj, system)
            if ~isa(system, 'PH_FEM_element')
                error('Only works if system is another ''PH_FEM_element''.');
            end
                        
            obj.add@PH_LinearSystem(system);
            
            % Collect attributes in a single element
            if obj.n_elements > 0
                obj.elements{1}.attributes = [obj.elements{1}.attributes; obj.elements{2}.attributes];
                [~, idx] = unique(obj.elements{1}.attributes(:,1));
                idx_del = 1:size(obj.elements{1}.attributes, 1);
                idx_del(idx) = [];
                obj.elements{1}.attributes(idx_del, :) = []; 
                obj.elements{1}.ports = [obj.elements{1}.ports obj.elements{2}.ports]; 
                obj.elements{2}.delete();
                obj.elements(2) = []; 
                
                obj.n_elements = 1;
            end  
            
            % Look for identical nodes and delete duplicates
            duplicateNodes = zeros(0,2);
            for n = 1:obj.n_nodes
                current = obj.nodes{n};
                if ~isa(current, 'PH_MechanicalNode') || any(any(duplicateNodes == n))
                    continue;
                else
                    for i=1:obj.n_nodes
                        if ~(i == n) && isa(obj.nodes{i}, 'PH_MechanicalNode')
                            if ~any(obj.nodes{i}.location ~= current.location)
                                duplicateNodes(end+1, 1) = n;
                                duplicateNodes(end, 2) = i;
                                break;
                            end
                        end
                    end
                end
            end
            
            if ~isempty(duplicateNodes)
                for n = 1:size(duplicateNodes, 1)
                    current = obj.nodes{duplicateNodes(n, 1)};
                    current.ports = [current.ports obj.nodes{duplicateNodes(n, 2)}.ports];
                    current.lockedDOFs = double(current.lockedDOFs | obj.nodes{duplicateNodes(n, 2)}.lockedDOFs);
                    for p = 1:length(current.ports)
                        obj.ports{current.ports(p)}.node = duplicateNodes(n, 1);
                    end
                end
                obj.nodes(duplicateNodes(:, 2)) = [];
                obj.n_nodes = length(obj.nodes);
            end
            
            
        end 
    end
    
    methods(Access = protected)   
        function copyData(obj, cp)
            copyData@PH_LinearSystem(obj, cp);
        end
        function cp = copyElement(obj)
               % PH_LinearSystem(n, E, J, Q, G, R, K, P, S, M, B)
            cp = PH_FEM_element(obj.n, obj.E(1:obj.n, 1:obj.n), obj.J, obj.Q, obj.G, ...
                                obj.R, obj.K, obj.P, obj.S, obj.M, obj.B);
            obj.copyData(cp);
        end
    end
end