classdef PH_FEM_mesh < PH_LinearSystem
    properties(SetAccess = protected, GetAccess = public) 
        nodeTable
        elementTable
        elementTypes
    end
    
    methods(Access = public)
        function obj = PH_FEM_mesh(nodeArray, elementArray, elements)
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
            
            obj = obj@PH_LinearSystem(0, [], zeros(0), zeros(0), zeros(0));
            
            for n=1:n_n
                obj.nodes{n} = PH_MechanicalNode([], nodeArray(n, :));
            end
            obj.n_nodes = n_n;
            
            for i=1:n_e
                current = copy(elements{elementArray(i, end)});
                current.setPosition(nodeArray(elementArray(i, 1:current.n_nodes), :));
                
                for k=1:length(current.inputNames)
                    current.inputNames{k} = ['e' num2str(i) '.' current.inputNames{k}];
                    current.outputNames{k} = ['e' num2str(i) '.' current.outputNames{k}];
                end
                % Concatenate the current element to the mesh
                obj.add(current);
                current.delete();
            end
  
            obj.elements{1}.delete();
            obj.elements(1) = [];
            obj.n_elements = obj.n_elements - 1;
            
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
        
        function fixNodeDOFs(obj, node, dofs)
            if ~ismatrix(node)
                error('nodes must be supplied as a column vector');
            end
            if size(dofs, 1) ~= length(node) || size(dofs, 2) ~= 6
                error('dofs must be an nx6 matrix, where n is the number of nodes'); 
            end
            
            for i=1:length(node)
                obj.nodes{node(i)}.lockedDOFs = dofs(i,:);
            end
        end 
    end
end