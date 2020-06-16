classdef PH_FEM_element < PH_LinearSystem 
    methods(Access = public)
        function obj = PH_FEM_element(n, J, Q, varargin)
            obj = obj@PH_LinearSystem(n, J, Q, varargin{:});
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
            
            Tr = round(vrrotvec2mat(vrrotvec(dir,nDir)), 12);
            
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
               % PH_LinearSystem(n, J, Q, G, R, K, P, S, M, B)
            cp = PH_FEM_element(obj.n, obj.J, obj.Q, obj.G, ...
                                obj.R, obj.P, obj.S, obj.M, obj.B);
            obj.copyData(cp);
        end
    end
end