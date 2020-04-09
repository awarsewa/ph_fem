classdef PH_FEM_mesh < PH_LinearSystem
    properties(SetAccess = protected, GetAccess = public) 
        nodeTable
        elementTable
        elementTypes
    end
    
    methods(Access = public)
        function obj = PH_FEM_mesh(nodeArray, elementArray, elements, attributes)
            % Error checking 
            n_n = size(nodeArray, 1);
            n_e = size(elementArray, 1);
                        
            if n_n < 1 || n_e < 1
                error('At least one node/element has to be given');
            end
            
            if size(nodeArray, 2) ~= 3
                error('Node array must be an nx3 matrix');
            end
            
            if size(elementArray, 2) < 2
                error('Element array must have at least 2 entries per element. One node and the element type');
            end
            
            if ~isa(elements, 'cell')
                error('elements must be given as cell array of elements');
            end
            
            if nargin > 3 && ~isempty(attributes)
                if ~isa(attributes, 'cell') || size(attributes, 1) ~= n_e
                    error('attributes must be given as a cell array');
                end
            else
                attributes = cell(n_e, 1);
            end
            
            obj = obj@PH_LinearSystem(0, [], zeros(0), zeros(0), zeros(0));
            
            
            for n=1:n_n
                obj.nodes{n} = PH_MechanicalNode([], [], nodeArray(n, :));
            end
            obj.n_nodes = n_n;
            
            
            for i=1:n_e
                current = copy(elements{elementArray(i, end)});
                %current.setPosition(nodeArray(elementArray(i, 1:current.n_nodes), :));
                current.setPosition(nodeArray(elementArray(i, 1:2), :));
                current.setAttributes(attributes(i, :));
                
                for k=1:length(current.ports)
                    current.ports{k}.inputName = ['e' num2str(i) '.' current.ports{k}.inputName];
                    current.ports{k}.outputName = ['e' num2str(i) '.' current.ports{k}.outputName];
                end
                % Concatenate the current element to the mesh
                obj.add(current);
                current.delete();
            end
  
            obj.elements{1}.delete();
            obj.elements(1) = [];
            obj.n_elements = obj.n_elements - 1;
            for n = 1:obj.n_nodes
                if ~isempty(obj.nodes{n}.elements)
                    obj.nodes{n}.elements = obj.nodes{n}.elements - 1; 
                end
            end
            
            % Set class properties
            obj.nodeTable = nodeArray;
            obj.elementTable = elementArray;
            obj.elementTypes = elements;
        end
        
        % Modified add function that removes duplicate nodes
        function add(obj, system)
            obj.add@PH_LinearSystem(system);
            
            % Look for identical nodes and delete duplicates
            duplicateNodes = zeros(0,2);
            for n = 1:obj.n_nodes
                current = obj.nodes{n};
                if ~isa(current, 'PH_MechanicalNode') || any(any(duplicateNodes == n)) || current.internal
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
                    current.elements = [current.elements obj.nodes{duplicateNodes(n, 2)}.elements];
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
                
                
                for n = 1:size(duplicateNodes, 1)            
                    for p = 1:obj.n_ports
                        port = obj.ports{p};
                        for pn = 1:length(port.nodes)
                            if port.nodes(pn) > duplicateNodes(n, 2)
                                port.nodes(pn) = port.nodes(pn) -1; 
                            end
                        end
                    end
                    for e = 1:obj.n_elements
                        element = obj.elements{e};
                        for en = 1:length(element.nodes)
                            if element.nodes(en) > duplicateNodes(n, 2)
                                element.nodes(en) = element.nodes(en) -1; 
                            end
                        end
                    end
                    obj.nodes(duplicateNodes(n, 2)) = [];
                    obj.n_nodes = length(obj.nodes);
                    for i=n+1:size(duplicateNodes,1)
                        if duplicateNodes(i, 2) > duplicateNodes(n, 2) 
                            duplicateNodes(i, 2) = duplicateNodes(i, 2) -1;
                        end
                    end
                end

                
            end
        end 
                
        % Check whether a selected DOF of a given node is fixed?
        function ret = isFixed(obj, node, dof)
            n = obj.getNodeAtLocation(obj.nodeTable(node, :));
            ret = obj.nodes{n}.lockedDOFS(dof);
        end
        
        function node = getNodeAtLocation(obj, location)
            node = 0;
            for n=1:obj.n_nodes
                if isa(obj.nodes{n}, 'PH_MechanicalNode') && ~any(obj.nodes{n}.location ~= location)
                    node = n;
                    break;
                end
            end
        end
        
        function fixNodeDOFs(obj, location, dofs)
            if ~ismatrix(location) || size(location, 1) ~= 1 || size(location, 2) ~= 3
                error('location must be supplied as a 1x3 vector');
            end
            if size(dofs, 1) ~= 1 || size(dofs, 2) ~= 6
                error('dofs must be an 1x6 vector'); 
            end
            
            node = obj.getNodeAtLocation(location);
            obj.nodes{node}.lockedDOFs = dofs;
        end 
        
        function addExternalInputsAtNodes(obj)
            for n=1:obj.n_nodes
                node = obj.nodes{n};
                if ~isa(node, 'PH_MechanicalNode') %|| node.internal
                    continue
                end
                
                % Get all ports at this node
                dir_str = 'xyz';
                for d = 1:3
                    dir = zeros(1, 3);
                    dir(d) = 1;
                    % ports are connected at this port in this direction?
                    ports_f = obj.getMechanicalPorts(n, {'force'}, dir);
                    if ~isempty(ports_f)
                        obj.addExternalPort('mechanical', 'force', n, obj.n_u+1, ...
                                            ['n' num2str(n) '_F' dir_str(d)], ['n' num2str(n) '_v' dir_str(d)], dir');
                        obj.n_u = obj.n_u +1;
                        
                        g_new = zeros(obj.n, 1);
                        p_new = zeros(obj.n, 1);
                        
                        obj.G = [obj.G, g_new];
                        obj.P = [obj.P, p_new];
                        obj.M = blkdiag(obj.M, 0);
                        obj.S = blkdiag(obj.S, 0);
                        obj.C_u = [obj.C_u, zeros(size(obj.C_u, 1), 1)];
                        obj.C_y = [obj.C_y, zeros(size(obj.C_y, 1), 1)];
                    end
                end
                for d = 1:3
                    dir = zeros(1, 3);
                    dir(d) = -1;

                    ports_t = obj.getMechanicalPorts(n, {'torque'}, dir);
                    if ~isempty(ports_t)
                        obj.addExternalPort('mechanical', 'torque', n, obj.n_u+1, ...
                                            ['n' num2str(n) '_M' dir_str(d)], ['n' num2str(n) '_w' dir_str(d)], dir');
                        obj.n_u = obj.n_u +1;
                        
                        g_new = zeros(obj.n, 1);
                        p_new = zeros(obj.n, 1);
                        
                        obj.G = [obj.G, g_new];
                        obj.P = [obj.P, p_new];
                        obj.M = blkdiag(obj.M, 0);
                        obj.S = blkdiag(obj.S, 0);
                        obj.C_u = [obj.C_u, zeros(size(obj.C_u, 1), 1)];
                        obj.C_y = [obj.C_y, zeros(size(obj.C_y, 1), 1)];
                    end
                end
            end
        end
    end
end