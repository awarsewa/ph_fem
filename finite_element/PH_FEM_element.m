classdef PH_FEM_element < PH_LinearSystem 
    methods(Access = public)
        function obj = PH_FEM_element(n, E, J, Q, varargin)
            obj = obj@PH_LinearSystem(n, E, J, Q, varargin{:});
        end
        
        function setPosition(obj, nodePositions)
            if ~ismatrix(nodePositions) || ...
                size(nodePositions, 1) ~= 2 || ...
                size(nodePositions, 2) ~= 3 
                error('NodePositions must be a 2x3 matrix');
            end
            dir = obj.nodes{end}.location - obj.nodes{1}.location;
            dir = dir./norm(dir);
            nDir = nodePositions(2,:) - nodePositions(1,:);
            nDir = nDir./norm(nDir);
            
            Tr = round(vrrotvec2mat(vrrotvec([1 0 0],nDir)), 12);
            
            
            % new direction is in some plane
            for p=1:obj.n_ports
                obj.ports{p}.orientation = Tr*obj.ports{p}.orientation;
                obj.ports{p}.orientation(abs(obj.ports{p}.orientation) < 1e-9) = 0;
                obj.ports{p}.orientation = obj.ports{p}.orientation/norm(obj.ports{p}.orientation);
            end    

            pos_x = linspace(nodePositions(1,1), nodePositions(2,1), obj.n_nodes);
            pos_y = linspace(nodePositions(1,2), nodePositions(2,2), obj.n_nodes);
            pos_z = linspace(nodePositions(1,3), nodePositions(2,3), obj.n_nodes);
            for n=1:obj.n_nodes
                obj.nodes{n}.location = [pos_x(n), pos_y(n), pos_z(n)];
                obj.nodes{n}.lockedDOFs = [round(Tr*obj.nodes{n}.lockedDOFs(1:3)', 10)' ...
                                           round(Tr*obj.nodes{n}.lockedDOFs(4:6)', 10)'];
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
                obj.elements{1}.nodes = [obj.elements{1}.nodes obj.elements{2}.nodes];
                obj.elements{2}.delete();
                obj.elements(2) = []; 
                
                for n = 1:obj.n_nodes
                    obj.nodes{n}.elements(obj.nodes{n}.elements == 2) = 1;
                end
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
                        for pn = 1:length(obj.ports{current.ports(p)}.nodes)
                            if(obj.ports{current.ports(p)}.nodes(pn) == duplicateNodes(n, 2))
                                obj.ports{current.ports(p)}.nodes(pn) = duplicateNodes(n, 1);
                            end
                        end
                    end
                    for e = 1:length(current.elements)
                        for en = 1:length(obj.elements{current.elements(e)}.nodes)
                            if(obj.elements{current.elements(e)}.nodes(en) == duplicateNodes(n, 2))
                                obj.elements{current.elements(e)}.nodes(en) = duplicateNodes(n, 1);
                            end
                        end
                    end
                    obj.nodes{duplicateNodes(n, 2)}.delete();
                end
                obj.nodes(duplicateNodes(:, 2)) = [];
                obj.n_nodes = length(obj.nodes);
            end           
        end 
        
        function setAttributes(obj, attribs)
            if ~isempty(attribs{1})
                warning('The element has no attributes to set');
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
                                obj.R, obj.P, obj.S, obj.M, obj.B);
            obj.copyData(cp);
        end
    end
end